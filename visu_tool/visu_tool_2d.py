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

# --- Extended .kwf File Parser (Perfectly mirrors final C++ stream) ---
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
        floats_per_state = 2 if mode == 0 else 1
        frame_bytes = floats_per_state * size * 8 * num_qubits

        captions, times, frames, state_vectors, purities, stirap = [], [], [], [], [], []

        while True:
            # 1. Read simulation time (double)
            t_raw = f.read(8)
            if not t_raw:
                break
            t = struct.unpack("<d", t_raw)[0]

            # 2. Read global logical state vector amplitudes
            sv_bytes = state_dim * 8
            sv_raw = f.read(sv_bytes)
            if len(sv_raw) < sv_bytes: break
            state_vector = np.array(struct.unpack(f"<{state_dim}d", sv_raw))

            # 3. Read specific independent title string per qubit
            frame_qubit_titles = []
            for _ in range(num_qubits):
                raw_len = f.read(2)
                if not raw_len: break
                cap_len = struct.unpack("<H", raw_len)[0]
                caption = f.read(cap_len).decode("utf-8")
                frame_qubit_titles.append(caption)

            # 4. Read physical purity scalars per atom
            p_bytes = num_qubits * 8
            p_raw = f.read(p_bytes)
            if len(p_raw) < p_bytes: break
            qubit_purities = np.array(struct.unpack(f"<{num_qubits}d", p_raw))

            # 5. Read dynamic laser wavelength & intensity parameters per atom
            s_bytes = num_qubits * 4 * 8
            s_raw = f.read(s_bytes)
            if len(s_raw) < s_bytes: break
            s_vals = np.array(struct.unpack(f"<{num_qubits * 4}d", s_raw)).reshape(num_qubits, 4)

            # 6. Read consecutive spatial wavefunction grids for ALL atoms
            payload = f.read(frame_bytes)
            if len(payload) < frame_bytes: break
            arr = np.frombuffer(payload, dtype="<f8")

            captions.append(frame_qubit_titles)
            times.append(t)
            frames.append(arr)
            state_vectors.append(state_vector)
            purities.append(qubit_purities)
            stirap.append(s_vals)

    return num_qubits, captions, np.array(times), np.array(frames), np.array(state_vectors), np.array(purities), mode, np.array(stirap)

# Load layout arrays
num_qubits, captions, times, raw, state_vectors, purities, mode, stirap_data = load_kwf(kwf_file)
n_timesteps = raw.shape[0]
state_dim = 2 ** num_qubits

# --- CSV Export Logic ---
csv_headers = ['Time']
for i in range(state_dim):
    bitstr = bin(i)[2:].zfill(num_qubits)
    csv_headers.extend([f'Amp_{bitstr}_Real', f'Amp_{bitstr}_Imag'])
for i in range(num_qubits):
    csv_headers.append(f'Purity_Q{i}')

csv_rows = []
for i in range(n_timesteps):
    row = [times[i]]
    for amp in state_vectors[i]:
        row.extend([amp.real, amp.imag])
    row.extend(purities[i])
    csv_rows.append(row)

df_export = pd.DataFrame(csv_rows, columns=csv_headers)
df_export.to_csv("quantum_data.csv", index=False)
print("Data exported to quantum_data.csv")

# Reshape full spatial configurations into (timesteps, atoms, Ny, Nx)
if mode == 0:
    psi_all_atoms = raw.reshape(n_timesteps, num_qubits, Ny, Nx).astype(complex)
else:
    psi_all_atoms = np.sqrt(raw).reshape(n_timesteps, num_qubits, Ny, Nx).astype(complex)

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

# --- Grid Resolution Engine for Per-Atom Visual Boxes ---
# Automatically computes side-by-side grids depending on atom counts
ncols = num_qubits
nrows = 1
fig, axs = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6), squeeze=False)
fig.subplots_adjust(left=0.06, right=0.88, top=0.85, bottom=0.28)

spatial_images = []
caption_boxes = []

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

    # Assign a unique visual text box inside EACH individual atom panel container
    box = ax_q.text(
        0.02, 0.98, "", transform=ax_q.transAxes, fontsize=10, color="white",
        verticalalignment="top", bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.6, edgecolor="none")
    )
    caption_boxes.append(box)

# Reference positions from the first panel to scale peripheral indicators
ref_pos = axs[0, 0].get_position()
last_pos = axs[0, -1].get_position()

# --- Shared Inset Panel Setup ---
ax_qiskit = fig.add_axes([last_pos.x1 + 0.02, ref_pos.y1 - 0.22, 0.08, 0.22])
ax_bloch = fig.add_axes([last_pos.x1 + 0.02, ref_pos.y0, 0.09, 0.22], projection='3d')
ax_bloch.set_facecolor((0, 0, 0, 0))

bloch_warning = ax_bloch.text2D(
    0.5, 0.5, "ENTANGLED STATE\n\nBloch Sphere\nDisabled.",
    color="coral", ha="center", va="center", fontsize=8, weight="bold",
    transform=ax_bloch.transAxes, bbox=dict(facecolor='black', alpha=0.7, edgecolor='none')
)
bloch_warning.set_visible(False)

# Shared Laser diagnostics graph pinned across bottom boundary
ax_stirap = fig.add_axes([ref_pos.x0, 0.06, (last_pos.x1 - ref_pos.x0), 0.14])
ax_stirap.set_facecolor((0, 0, 0, 0.4))
ax_stirap.tick_params(colors='gray', labelsize=8)
for spine in ax_stirap.spines.values():
    spine.set_color('gray')
ax_stirap.set_xlabel("Time (sec)", color="gray", fontsize=8)
ax_stirap.set_ylabel("Intensity (W/cm²)", color="gray", fontsize=8)

laser1_line, = ax_stirap.plot([], [], color="lime", linewidth=1.5, label="Pump")
laser2_line, = ax_stirap.plot([], [], color="deepskyblue", linewidth=1.5, label="Stokes")
time_marker = ax_stirap.axvline(x=times[0], color="white", linestyle="--", alpha=0.5)
ax_stirap.set_xlim(times[0], times[-1])
ax_stirap.set_ylim(0, np.max(stirap_data[:, :, [1, 3]]) * 1.1)

# 3D Bloch Axis Init
def draw_bloch_skeleton(axis):
    u, v = np.mgrid[0:2*np.pi:14j, 0:np.pi:7j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)
    axis.plot_wireframe(x, y, z, color="gray", alpha=0.2, linewidth=0.5)
    axis.set_axis_off()
    axis.view_init(elev=20, azim=45)

draw_bloch_skeleton(ax_bloch)
bloch_vector, = ax_bloch.plot([0, 0], [0, 0], [0, 0], color="cyan", linewidth=2, marker='o', markersize=3)

# Qiskit Histogram layout configuration
ax_qiskit.set_facecolor("#111111")
ax_qiskit.set_title("Basis Probs", color="gray", fontsize=8, pad=4)
ax_qiskit.set_ylim(0, 1.05)
bitstrings = [f"|{bin(i)[2:].zfill(num_qubits)}⟩" for i in range(state_dim)]
bars = ax_qiskit.bar(bitstrings, np.abs(state_vectors[0]) ** 2, color="#6f42c1", edgecolor="#563d7c")
ax_qiskit.tick_params(colors='gray', labelsize=7, axis='x', rotation=90 if num_qubits > 2 else 45)
for spine in ax_qiskit.spines.values():
    spine.set_color('gray')

# Main Colorbar mapping
colorbar_rgb = phase_to_rgb(np.exp(1j * np.linspace(-np.pi, np.pi, 256)).reshape(1, -1))
custom_cmap = ListedColormap(colorbar_rgb[0])
cax = fig.add_axes([last_pos.x1 + 0.02, ref_pos.y0 + 0.26, 0.015, ref_pos.y1 - ref_pos.y0 - 0.28])
norm = Normalize(vmin=-np.pi, vmax=np.pi)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=custom_cmap), cax=cax)
cb.set_ticks([-np.pi, 0, np.pi])
cb.set_ticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
cb.ax.tick_params(labelsize=8)

current_limits = [L] * num_qubits

# --- Frame Renderer Pipeline ---
def update(frame):
    global current_limits
    
    # 1. Update each individual atom visualization box dynamically
    for q in range(num_qubits):
        psi = psi_all_atoms[frame, q]
        spatial_images[q].set_data(phase_to_rgb(psi))
        
        # Pull distinct dynamic caption saved specifically for this atom
        caption_boxes[q].set_text(captions[frame][q])

        # Individual Box Adaptive Zoom 
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

    # 2. Update Qiskit probability profiles
    current_probs = np.abs(state_vectors[frame]) ** 2
    for bar, h in zip(bars, current_probs):
        bar.set_height(h)

    # 3. Update Subsystem 0 Coherence (Bloch Vector tracking)
    purity_0 = purities[frame][0]
    if num_qubits > 1 and purity_0 < 0.95:
        bloch_vector.set_visible(False)
        bloch_warning.set_visible(True)
        ax_bloch.set_title(f"Q0 Mixed ({purity_0:.2f})", color="coral", fontsize=8, weight="bold")
    else:
        bloch_vector.set_visible(True)
        bloch_warning.set_visible(False)
        ax_bloch.set_title(f"Q0 Pure ({purity_0:.2f})", color="gray", fontsize=8)
        
        sv = state_vectors[frame]
        rho_00 = np.sum(np.abs(sv[0::2])**2)
        rho_11 = np.sum(np.abs(sv[1::2])**2)
        rho_01 = np.sum(sv[0::2] * np.conj(sv[1::2]))
        bx, by, bz = 2 * rho_01.real, 2 * rho_01.imag, rho_00 - rho_11
        
        bloch_vector.set_data([0, bx], [0, by])
        bloch_vector.set_3d_properties([0, bz])

    # 4. Update timeline markers on Laser Diagnostics tracking sheet
    laser1_line.set_data(times[:frame+1], stirap_data[:frame+1, 0, 1]) # Monitoring Atom 0 lasers
    laser2_line.set_data(times[:frame+1], stirap_data[:frame+1, 0, 3])
    time_marker.set_xdata([times[frame], times[frame]])

# --- Runtime Export Exec ---
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