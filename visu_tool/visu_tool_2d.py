import os
import struct
import sys
import tkinter as tk
from tkinter import filedialog
import numpy as np
import pandas as pd

# HEADLESS BACKEND (must be before pyplot import)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

output_dir = "frames"
Nx, Ny = 256, 256
L = 100.0
fps = 30
save_frames = True

# --- Tkinter File Dialog ---
_root = tk.Tk()
_root.withdraw()
_root.attributes("-topmost", True)
kwf_file = filedialog.askopenfilename(
    title="Open KetCat wavefunction file",
    filetypes=[("KetCat Wavefunction", "*.kwf"), ("All files", "*.*")],
)
_root.destroy()

if not kwf_file:
    print("No file selected. Exiting.")
    sys.exit(0)

if save_frames and not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- Extended .kwf File Parser (Perfectly mirrors RealImag C++ stream) ---
def load_kwf(path):
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic == b"KWF odds":
            raise ValueError("Deprecated format version")
        elif magic == b"KWF\x02":
            mode = struct.unpack("<B", f.read(1))[0]
            num_qubits = struct.unpack("<B", f.read(1))[0]
        else:
            raise ValueError(f"Not a valid .kwf file (bad magic: {magic!r})")

        size = struct.unpack("<Q", f.read(8))[0]

        global Nx, Ny, L
        N = struct.unpack("<Q", f.read(8))[0]
        Nx = Ny = N
        L = struct.unpack("<d", f.read(8))[0]

        state_dim = 2 ** num_qubits
        floats_per_state = 2
        frame_bytes = floats_per_state * size * 8 * num_qubits

        captions, times, frames, state_vectors, bloch_vectors, purities, stirap = [], [], [], [], [], [], []

        while True:
            t_raw = f.read(8)
            if not t_raw:
                break
            t = struct.unpack("<d", t_raw)[0]

            sv_bytes = state_dim * 8
            sv_raw = f.read(sv_bytes)
            if len(sv_raw) < sv_bytes: break
            state_vector = np.array(struct.unpack(f"<{state_dim}d", sv_raw))

            frame_qubit_titles = []
            for _ in range(num_qubits):
                raw_len = f.read(2)
                if not raw_len: break
                cap_len = struct.unpack("<H", raw_len)[0]
                caption = f.read(cap_len).decode("utf-8")
                frame_qubit_titles.append(caption)

            b_bytes = num_qubits * 4 * 8
            b_raw = f.read(b_bytes)
            if len(b_raw) < b_bytes: break
            b_vals =  np.array(struct.unpack(f"<{num_qubits * 4}d", b_raw)).reshape(num_qubits, 4)

            p_bytes = num_qubits * 8
            p_raw = f.read(p_bytes)
            if len(p_raw) < p_bytes: break
            qubit_purities = np.array(struct.unpack(f"<{num_qubits}d", p_raw))

            s_bytes = num_qubits * 4 * 8
            s_raw = f.read(s_bytes)
            if len(s_raw) < s_bytes: break
            s_vals = np.array(struct.unpack(f"<{num_qubits * 4}d", s_raw)).reshape(num_qubits, 4)

            payload = f.read(frame_bytes)
            if len(payload) < frame_bytes: break
            arr = np.frombuffer(payload, dtype="<f8")

            captions.append(frame_qubit_titles)
            times.append(t)
            frames.append(arr)
            state_vectors.append(state_vector)
            bloch_vectors.append(b_vals)    
            purities.append(qubit_purities)
            stirap.append(s_vals)

    return num_qubits, captions, np.array(times), np.array(frames), np.array(state_vectors), np.array(bloch_vectors), np.array(purities), np.array(stirap)

# Load data from file
num_qubits, captions, times, raw, state_vectors, bloch_data, purities, stirap_data = load_kwf(kwf_file)
n_timesteps = raw.shape[0]
state_dim = 2 ** num_qubits

# --- CSV Export ---
csv_headers = ['Time']
for i in range(state_dim):
    bitstr = bin(i)[2:].zfill(num_qubits)
    csv_headers.append(f'Prob_{bitstr}')
for i in range(num_qubits):
    csv_headers.append(f'Purity_Q{i}')

csv_rows = []
for i in range(n_timesteps):
    row = [times[i]]
    row.extend(state_vectors[i])
    row.extend(purities[i])
    csv_rows.append(row)

df_export = pd.DataFrame(csv_rows, columns=csv_headers)
df_export.to_csv("quantum_data.csv", index=False)
print("Data exported to quantum_data.csv")

# Reshape RealImag arrays
reshaped_raw = raw.reshape(n_timesteps, num_qubits, Ny, Nx, 2)
psi_all_atoms = (reshaped_raw[..., 0] + 1j * reshaped_raw[..., 1]).astype(complex)

x_coord = np.linspace(-L, L, Nx)
y_coord = np.linspace(-L, L, Ny)

# --- Phase to RGB Conversion Map ---
def phase_to_rgb(psi):
    phase = (np.angle(psi) + np.pi) / (2 * np.pi)
    phase = phase - np.floor(phase)
    amplitude = np.abs(psi) ** 2
    if amplitude.max() > 0:
        amplitude /= amplitude.max()

    r, g, b = np.zeros_like(phase), np.zeros_like(phase), np.zeros_like(phase)
    m1 = (phase < 0.25)
    f1 = phase[m1] / 0.25
    r[m1], g[m1], b[m1] = 1.0, 0.0, f1

    m2 = (phase >= 0.25) & (phase < 0.50)
    f2 = (phase[m2] - 0.25) / 0.25
    r[m2], g[m2], b[m2] = 1.0 - f2, 0.0, 1.0

    m3 = (phase >= 0.50) & (phase < 0.75)
    f3 = (phase[m3] - 0.50) / 0.25
    r[m3], g[m3], b[m3] = 0.0, f3, 1.0

    m4 = (phase >= 0.75)
    f4 = (phase[m4] - 0.75) / 0.25
    r[m4], g[m4], b[m4] = f4, 1.0 - f4, 1.0 - f4

    return np.clip(np.stack((r * amplitude, g * amplitude, b * amplitude), axis=-1), 0, 1)

# --- Main Canvas & Base Multi-panel Grid Layout ---
ncols = num_qubits
nrows = 1
# Increased height to 6.5 and tightly squeezed padding to maximize atom box sizes
fig, axs = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6.5), squeeze=False)
fig.subplots_adjust(left=0.04, right=0.90, top=0.92, bottom=0.24, wspace=0.15)

spatial_images = []
caption_boxes = []

# --- Custom High-Fidelity Bloch Sphere Insets per Atom Box ---
ax_blochs = []
bloch_vectors = []
bloch_warnings = []

def setup_bloch(axis):
    axis.set_axis_off()

    # High-density wireframe grid mesh matching legacy version
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    axis.plot_wireframe(x, y, z, color="gray", alpha=0.4, linewidth=0.5)

    # Subdued inner projection axes
    axis.plot([-1.1, 1.1], [0, 0], [0, 0], color="white", alpha=0.05)
    axis.plot([0, 0], [-1.1, 1.1], [0, 0], color="white", alpha=0.05)
    axis.plot([0, 0], [0, 0], [-1.1, 1.1], color="white", alpha=0.1)

    # Distinct RGB basis vector lines (X: Red, Y: Green, Z: Blue)
    axis.plot([0, 1], [0, 0], [0, 0], color="red", alpha=0.6)
    axis.plot([0, 0], [0, 1], [0, 0], color="green", alpha=0.6)
    axis.plot([0, 0], [0, 0], [0, 1], color="blue", alpha=0.6)

    # Standard state label decorations surrounding the sphere boundary
    o = 1.25
    axis.text(0, 0, o, r"$|0\rangle$", color="white", ha="center", va="bottom", fontsize=8)
    axis.text(0, 0, -o, r"$|1\rangle$", color="white", ha="center", va="top", fontsize=8)
    axis.text(o, 0, 0, r"$|+\rangle$", color="white", ha="left", va="center", fontsize=7)
    axis.text(-o, 0, 0, r"$|-\rangle$", color="white", ha="right", va="center", fontsize=7)
    axis.text(0, o, 0, r"$|+i\rangle$", color="white", ha="left", va="center", fontsize=7)
    axis.text(0, -o, 0, r"$|-i\rangle$", color="white", ha="right", va="center", fontsize=7)

    axis.view_init(elev=20, azim=45)

for q in range(num_qubits):
    ax_q = axs[0, q]
    img = ax_q.imshow(
        phase_to_rgb(psi_all_atoms[0, q]),
        origin="lower",
        extent=[-L, L, -L, L],
        interpolation="bicubic"
    )
    spatial_images.append(img)
    ax_q.set_xlabel("x (a.u.)", fontsize=9)
    ax_q.set_ylabel("y (a.u.)", fontsize=9)
    ax_q.set_title(f"Atom {q} Subspace", fontsize=11, weight='bold')

    # HUD title box (Top-left corner)
    box = ax_q.text(
        0.02, 0.98, "", transform=ax_q.transAxes, fontsize=8, color="white",
        verticalalignment="top", multialignment="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.6, edgecolor="none")
    )
    caption_boxes.append(box)

    # Embed translucent inset 3D panel into the bottom-right corner
    pos = ax_q.get_position()
    ax_b = fig.add_axes([pos.x1 - 0.2, pos.y0 + 0.01, 0.2, 0.3], projection='3d')
    ax_b.set_facecolor((0, 0, 0, 0))
    
    setup_bloch(ax_b)
    
    # State tracking pointer vector
    b_vec, = ax_b.plot([0, 0], [0, 0], [0, 0], color="cyan", linewidth=2.5, marker='o', markersize=4, zorder=10)
    
    # Entanglement alert box placed at center of sphere when pure state collapses
    b_warn = ax_b.text2D(
        0.5, 0.5, "MIXED\nSTATE", color="coral", ha="center", va="center", fontsize=8, weight="bold",
        transform=ax_b.transAxes, bbox=dict(facecolor='black', alpha=0.7, edgecolor='none', pad=2)
    )
    b_warn.set_visible(False)
    
    ax_blochs.append(ax_b)
    bloch_vectors.append(b_vec)
    bloch_warnings.append(b_warn)

ref_pos = axs[0, 0].get_position()
last_pos = axs[0, -1].get_position()

# --- Qiskit-style Probability Histogram Panel (Right Margin) ---
ax_qiskit = fig.add_axes([last_pos.x1 + 0.02, ref_pos.y0 + 0.22, 0.07, 0.35])
ax_qiskit.set_facecolor("#111111")
ax_qiskit.set_title("Basis Probs", color="gray", fontsize=8, pad=4)
ax_qiskit.set_ylim(0, 1.05)
bitstrings = [f"|{bin(i)[2:].zfill(num_qubits)}⟩" for i in range(state_dim)]
bars = ax_qiskit.bar(bitstrings, state_vectors[0], color="#6f42c1", edgecolor="#563d7c")
ax_qiskit.tick_params(colors='gray', labelsize=7, axis='x', rotation=90 if num_qubits > 2 else 45)
for spine in ax_qiskit.spines.values():
    spine.set_color('gray')

# --- Shared Global Phase Colorbar ---
colorbar_rgb = phase_to_rgb(np.exp(1j * np.linspace(-np.pi, np.pi, 256)).reshape(1, -1))
custom_cmap = ListedColormap(colorbar_rgb[0])
cax = fig.add_axes([last_pos.x1 + 0.02, ref_pos.y0, 0.012, 0.18])
norm = Normalize(vmin=-np.pi, vmax=np.pi)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=custom_cmap), cax=cax)
cb.set_ticks([-np.pi, 0, np.pi])
cb.set_ticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
cb.ax.tick_params(labelsize=8)

# --- Shared Laser Diagnostics Panel with Safe Logarithmic Scaling ---
# Adjusted height and bottom alignment to push the main atom plots larger
ax_stirap = fig.add_axes([ref_pos.x0, 0.05, (last_pos.x1 - ref_pos.x0), 0.12])
ax_stirap.set_facecolor((0, 0, 0, 0.4))
ax_stirap.set_yscale('log')
ax_stirap.tick_params(colors='gray', labelsize=8)
for spine in ax_stirap.spines.values():
    spine.set_color('gray')
ax_stirap.set_xlabel("Time (sec)", color="gray", fontsize=8)
ax_stirap.set_ylabel("Intensity (W/cm²)", color="gray", fontsize=8)

laser1_line, = ax_stirap.plot([], [], color="lime", linewidth=1.5, label="Pump")
laser2_line, = ax_stirap.plot([], [], color="deepskyblue", linewidth=1.5, label="Stokes")
time_marker = ax_stirap.axvline(x=times[0], color="white", linestyle="--", alpha=0.5)

ax_stirap.set_xlim(times[0], times[-1])
max_laser_val = max(np.max(stirap_data[:, :, [1, 3]]), 200.0)
ax_stirap.set_ylim(0.1, max_laser_val * 1.5) 
ax_stirap.legend(loc="upper right", framealpha=0.3, facecolor="black", edgecolor="none", fontsize=7, labelcolor="gray")

current_limits = [L] * num_qubits

# --- Frame Renderer Pipeline ---
def update(frame):
    global current_limits
    current_probs = state_vectors[frame]
    
    for q in range(num_qubits):
        # 1. Update Spatial Wavefunction Cloud
        psi = psi_all_atoms[frame, q]
        spatial_images[q].set_data(phase_to_rgb(psi))
        
        formatted_caption = captions[frame][q].replace("|", "\n")
        caption_boxes[q].set_text(formatted_caption)

        # Per-Atom Independent Adaptive Auto-Zoom Engine
        amplitude_sq = np.abs(psi) ** 2
        max_amp = amplitude_sq.max()
        target_limit = L
        if max_amp > 1e-8:
            threshold = 0.005 * max_amp
            inside_indices = np.argwhere(amplitude_sq >= threshold)
            if len(inside_indices) > 0:
                y_indices, x_indices = inside_indices[:, 0], inside_indices[:, 1]
                max_bound = max(abs(x_coord[x_indices.min()]), abs(x_coord[x_indices.max()]),
                                abs(y_coord[y_indices.min()]), abs(y_coord[y_indices.max()]))
                target_limit = np.clip(max_bound * 1.25, 8.0, L)
        
        current_limits[q] = 0.9 * current_limits[q] + 0.1 * target_limit
        axs[0, q].set_xlim(-current_limits[q], current_limits[q])
        axs[0, q].set_ylim(-current_limits[q], current_limits[q])

        # 2. Update Bloch Sphere Insets
        purity_q = purities[frame][q]
        if num_qubits > 1 and purity_q < 0.95:
            bloch_vectors[q].set_visible(False)
            bloch_warnings[q].set_visible(True)
            ax_blochs[q].set_title(f"Q{q} Mix({purity_q:.2f})", color="coral", fontsize=7, y=0.98, weight="bold")
        else:
            bloch_vectors[q].set_visible(True)
            bloch_warnings[q].set_visible(False)
            ax_blochs[q].set_title(f"Q{q} Pure({purity_q:.2f})", color="darkgray", fontsize=7, y=0.98)
            
            bloch = bloch_data[frame][q]
            a = complex(bloch[0], bloch[1])
            b = complex(bloch[2], bloch[3])
            bx = 2 * (a.real * b.real + a.imag * b.imag)
            by = 2 * (a.real * b.imag - a.imag * b.real)
            bz = (np.abs(a)**2) - (np.abs(b)**2)
            bloch_vectors[q].set_data([0, bx], [0, by])
            bloch_vectors[q].set_3d_properties([0, bz])
            print(a, b, bx, by, bz)


    # 3. Update Qiskit-style Histogram Columns
    for bar, h in zip(bars, current_probs):
        bar.set_height(h)

    # 4. Update Laser Diagnostics Timeline Tracking (With 0.1 floor safety limits)
    y1_vals = np.clip(stirap_data[:frame+1, 0, 1], 0.1, None)
    y2_vals = np.clip(stirap_data[:frame+1, 0, 3], 0.1, None)
    
    laser1_line.set_data(times[:frame+1], y1_vals)
    laser2_line.set_data(times[:frame+1], y2_vals)
    time_marker.set_xdata([times[frame], times[frame]])

# --- Runtime Export Execution ---
print(f"Exporting {n_timesteps} multi-panel frames for a {num_qubits}-atom system...")
for frame in range(n_timesteps):
    update(frame)
    if save_frames:
        filename = os.path.join(output_dir, f"{frame:04d}.png")
        fig.savefig(filename, dpi=130)
        if frame % 20 == 0:
            print(f"Rendered step: {frame}/{n_timesteps}")

print("Done processing all frames successfully.")
plt.close('all')