import os
import struct
import sys
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import numpy as np
import pandas as pd

# HEADLESS BACKEND (must be before pyplot import)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.image as mpimg
from matplotlib.colors import ListedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --- Simulation Defaults ---
Nx, Ny = 256, 256
L = 100.0
fps = 30
save_frames = True
kwf_file = None

# --- Process file from command line, if given. If not, provide a FileDialog ---
if len(sys.argv) > 1:
    kwf_file = sys.argv[1]
else:
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

output_dir = Path(kwf_file).stem
if save_frames and not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- Extended .kwf File Parser (Perfectly mirrors RealImag C++ stream) ---
def load_kwf(path):
    with open(path, "rb") as f:
        # Check file magic headers to verify version integrity
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

        # Read frames until EOF is reached
        while True:
            t_raw = f.read(8)
            if not t_raw:
                break
            t = struct.unpack("<d", t_raw)[0]

            # Parse full complex state vector coefficients
            sv_bytes = state_dim * 2 * 8
            sv_raw = f.read(sv_bytes)
            if len(sv_raw) < sv_bytes: break
            sv_floats = struct.unpack(f"<{state_dim * 2}d", sv_raw)
            state_vector = np.array(sv_floats[0::2]) + 1j * np.array(sv_floats[1::2])

            # Extract variable length string annotations for individual qubits
            frame_qubit_titles = []
            for _ in range(num_qubits):
                raw_len = f.read(2)
                if not raw_len: break
                cap_len = struct.unpack("<H", raw_len)[0]
                caption = f.read(cap_len).decode("utf-8")
                frame_qubit_titles.append(caption)

            # Unpack Bloch components
            b_bytes = num_qubits * 4 * 8
            b_raw = f.read(b_bytes)
            if len(b_raw) < b_bytes: break
            b_vals =  np.array(struct.unpack(f"<{num_qubits * 4}d", b_raw)).reshape(num_qubits, 4)

            # Parse trace distance / state purity data
            p_bytes = num_qubits * 8
            p_raw = f.read(p_bytes)
            if len(p_raw) < p_bytes: break
            qubit_purities = np.array(struct.unpack(f"<{num_qubits}d", p_raw))

            # Unpack STIRAP laser diagnostic telemetry
            s_bytes = num_qubits * 4 * 8
            s_raw = f.read(s_bytes)
            if len(s_raw) < s_bytes: break
            s_vals = np.array(struct.unpack(f"<{num_qubits * 4}d", s_raw)).reshape(num_qubits, 4)

            # Read spatial 2D grid wavefunction payload
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

"""
# --- CSV Export Routine ---
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
"""

# Reconstruct complex spatial grids from interwoven Real/Imag components
reshaped_raw = raw.reshape(n_timesteps, num_qubits, Ny, Nx, 2)
psi_all_atoms = (reshaped_raw[..., 0] + 1j * reshaped_raw[..., 1]).astype(complex)

x_coord = np.linspace(-L, L, Nx)
y_coord = np.linspace(-L, L, Ny)

# --- Phase to RGB Conversion Map ---
def phase_to_rgb(psi):
    """Maps complex phase angle to Hue and absolute amplitude to value saturation."""
    phase = (np.angle(psi) + np.pi) / (2 * np.pi)
    phase = phase - np.floor(phase)
    amplitude = np.abs(psi) 
    if amplitude.max() > 0:
        amplitude /= amplitude.max()
        amplitude = amplitude**0.2

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
fig, axs = plt.subplots(nrows, ncols, figsize=(6 * ncols, 7.5), squeeze=False)
fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.32, wspace=0.15)

spatial_images = []
caption_boxes = []

# --- Custom High-Fidelity Bloch Sphere Insets per Atom Box ---
ax_blochs = []
bloch_vectors = []
bloch_warnings = []

def setup_bloch(axis):
    """Draws wireframe geometry and baseline basis state annotations for 3D state visualizations."""
    axis.set_axis_off()

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    axis.plot_wireframe(x, y, z, color="gray", alpha=0.4, linewidth=0.5)

    # Reference grid lines
    axis.plot([-1.1, 1.1], [0, 0], [0, 0], color="white", alpha=0.05)
    axis.plot([0, 0], [-1.1, 1.1], [0, 0], color="white", alpha=0.05)
    axis.plot([0, 0], [0, 0], [-1.1, 1.1], color="white", alpha=0.1)

    # Unit coordinate vectors
    axis.plot([0, 1], [0, 0], [0, 0], color="red", alpha=0.6)
    axis.plot([0, 0], [0, 1], [0, 0], color="green", alpha=0.6)
    axis.plot([0, 0], [0, 0], [0, 1], color="blue", alpha=0.6)

    # Standard basis state text layout
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
    ax_q.set_title(f"Qubit {q} Subspace", fontsize=11, weight='bold')

    # HUD text overlay for runtime diagnostics within the atom grid
    box = ax_q.text(
        0.02, 0.98, "", transform=ax_q.transAxes, fontsize=8, color="white",
        verticalalignment="top", multialignment="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.6, edgecolor="none")
    )
    caption_boxes.append(box)

    # Place the matching Bloch sphere overlay in the lower right corner of the spatial axes
    pos = ax_q.get_position()
    ax_b = fig.add_axes([pos.x1 - 0.18, pos.y0 + 0.01, 0.18, 0.25], projection='3d')
    ax_b.set_facecolor((0, 0, 0, 0))
    
    setup_bloch(ax_b)
    
    b_vec, = ax_b.plot([0, 0], [0, 0], [0, 0], color="cyan", linewidth=2.5, marker='o', markersize=4, zorder=10)
    
    # Text alert indicator for mixed / entangled states where purity drops
    b_warn = ax_b.text2D(
        0.5, 0.5, "MIXED\nSTATE", color="coral", ha="center", va="center", fontsize=8, weight="bold",
        transform=ax_b.transAxes, bbox=dict(facecolor='black', alpha=0.7, edgecolor='none', pad=2)
    )
    b_warn.set_visible(False)
    
    ax_blochs.append(ax_b)
    bloch_vectors.append(b_vec)
    bloch_warnings.append(b_warn)

ref_pos = 0.02 
last_pos = axs[0, -1].get_position()
total_width = last_pos.x1 - ref_pos

# --- Re-arranged Bottom Row with Logo and Text ---

# 1. Project Logo Panel (Far Left)
ax_logo = fig.add_axes([ref_pos, 0.11, total_width * 0.14, 0.15])
ax_logo.set_axis_off()
logo_path = os.path.join(os.path.dirname(__file__) if '__file__' in locals() else '.', '../doc/logo.png')
if os.path.exists(logo_path):
    logo_img = mpimg.imread(logo_path)
    ax_logo.imshow(logo_img)
else:
    ax_logo.text(0.5, 0.5, "[ Logo Placeholder ]", color="gray", ha="center", va="center")

# Metadata and experiment branding tags
fig.text(ref_pos, 0.1, "KetCat AB INITIO NEUTRAL ATOM\nQUANTUM COMPUTER SIMULATOR v3.0", 
         color="black", fontsize=6, weight="bold", ha="left", va="top")

fig.text(ref_pos, 0.07, "3-Qubit Greenberger–Horne–Zeilinger state demo\n(maximally entangled state)", 
         color="black", fontsize=7, weight="bold", style="italic", ha="left", va="top")

# 2. Shared Laser Diagnostics Panel (Next to Logo)
ax_stirap = fig.add_axes([ref_pos + 0.23, 0.06, 0.3, 0.18])
ax_stirap.set_facecolor((0, 0, 0, 0.8))
ax_stirap.set_yscale('log')
ax_stirap.tick_params(colors='black', labelsize=8)
for spine in ax_stirap.spines.values():
    spine.set_color('black')
ax_stirap.set_xlabel("Laser control  - Global time (sec)", color="black", fontsize=8)
ax_stirap.set_ylabel("Intensity (W/cm²)", color="black", fontsize=8)

laser1_line, = ax_stirap.plot([], [], color="lime", linewidth=1.5, label="Pump")
laser2_line, = ax_stirap.plot([], [], color="deepskyblue", linewidth=1.5, label="Stokes")
time_marker = ax_stirap.axvline(x=times[0], color="white", linestyle="--", alpha=0.5)

ax_stirap.set_xlim(times[0], times[-1])

max_laser_val = max(np.max(stirap_data[:, :, [1, 3]]), 200.0)
ax_stirap.set_ylim(0.1, max_laser_val * 1.5) 
ax_stirap.legend(loc="upper right", framealpha=0.6, facecolor="black", edgecolor="none", fontsize=7, labelcolor="white")

# 3. Shared Global Phase Colorbar (Perfectly in the middle)
colorbar_rgb = phase_to_rgb(np.exp(1j * np.linspace(-np.pi, np.pi, 256)).reshape(1, -1))
custom_cmap = ListedColormap(colorbar_rgb[0])
cax = fig.add_axes([ref_pos + total_width * 0.61, 0.06, total_width * 0.02, 0.18])
norm = Normalize(vmin=-np.pi, vmax=np.pi)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=custom_cmap), cax=cax)
cb.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
cb.set_ticklabels([r"$-\pi$", r"$-\pi$/2", r"$0$", r"$\pi$/2", r"$\pi$"])
cb.ax.tick_params(labelsize=8)
cb.set_label("Phase (rad)", color="black", fontsize=8, labelpad=5)

# 4. Qiskit-style Probability Histogram Panel (Right side of the bottom row)
ax_qiskit = fig.add_axes([ref_pos + total_width * 0.72, 0.06, total_width * 0.32, 0.18])
ax_qiskit.set_facecolor((0, 0, 0, 0.4))
ax_qiskit.set_title("Basis Probabilities & Phase Disks", color="black", fontsize=8, pad=4)
ax_qiskit.set_ylim(0, 1.25) 
bitstrings = [f"|{bin(i)[2:].zfill(num_qubits)}⟩" for i in range(state_dim)]

# Fixed initialization passing float probability magnitudes to bypass Matplotlib ComplexWarnings
bars = ax_qiskit.bar(bitstrings, np.abs(state_vectors[0]) ** 2, color="#6f42c1", edgecolor="#563d7c")
ax_qiskit.tick_params(colors='black', labelsize=7, axis='x', rotation=90 if num_qubits > 2 else 45)
for spine in ax_qiskit.spines.values():
    spine.set_color('black')

# --- Phase Disk Initialization atop the Bar Charts ---
phase_disks = []
phase_hands = []

# Disks are tracked at a static Y line offset to prevent layout jittering during updates
disk_y_pos = 0.5
max_disk_radius = 0.3  # Bound radius parameters to bar widths

for i in range(state_dim):
    # Base configuration boundary circle
    disk_outline = plt.Circle((i, disk_y_pos), max_disk_radius, color='darkgray', fill=False, alpha=0.5, linewidth=0.8)
    # Scalable inner circle tracker tied to amplitude changes
    disk_fill = plt.Circle((i, disk_y_pos), 0.0, color='cyan', alpha=0.6, zorder=4)
    # Polar phase dial needle tracker
    hand, = ax_qiskit.plot([i, i], [disk_y_pos, disk_y_pos], color='white', linewidth=1.3, zorder=5)
    
    ax_qiskit.add_patch(disk_outline)
    ax_qiskit.add_patch(disk_fill)
    
    phase_disks.append(disk_fill)
    phase_hands.append(hand)

current_limits = [35] * num_qubits

# --- Frame Renderer Pipeline ---
def update(frame):
    global current_limits
    current_probs = state_vectors[frame]
    
    for q in range(num_qubits):
        psi = psi_all_atoms[frame, q]
        spatial_images[q].set_data(phase_to_rgb(psi))
        
        # Substitute layout delimiters for line wrapping inside HUD fields
        formatted_caption = captions[frame][q].replace("|", "\n")
        caption_boxes[q].set_text(formatted_caption)

        # Dynamic, adaptive auto-zooming calculation tracking localized wave packet density expansion
        amplitude_sq = np.abs(psi) ** 2
        max_amp = amplitude_sq.max()
        target_limit = L
        if max_amp > 1e-8:
            threshold = 0.0001 * max_amp
            inside_indices = np.argwhere(amplitude_sq >= threshold)
            if len(inside_indices) > 0:
                y_indices, x_indices = inside_indices[:, 0], inside_indices[:, 1]
                max_bound = max(abs(x_coord[x_indices.min()]), abs(x_coord[x_indices.max()]),
                                abs(y_coord[y_indices.min()]), abs(y_coord[y_indices.max()]))
                target_limit = np.clip(max_bound * 1.25, 8.0, L)
        
        # Smooth asymptotic interpolation for plot boundary changes
        current_limits[q] = 0.9 * current_limits[q] + 0.1 * target_limit
        axs[0, q].set_xlim(-current_limits[q], current_limits[q])
        axs[0, q].set_ylim(-current_limits[q], current_limits[q])

        # Evaluate state purity bounds to toggle Bloch vs mixed system warnings
        purity_q = purities[frame][q]
        if num_qubits > 1 and purity_q < 0.90:
            bloch_vectors[q].set_visible(False)
            bloch_warnings[q].set_visible(True)
            ax_blochs[q].set_title(f"Purity={purity_q:.2f}", color="coral", fontsize=7, y=0.98, weight="bold")
        else:
            bloch_vectors[q].set_visible(True)
            bloch_warnings[q].set_visible(False)
            ax_blochs[q].set_title(f"Purity={purity_q:.2f}", color="darkgray", fontsize=7, y=0.98)
            
            # Reconstruct exact coordinates from mixed density matrices mapped via raw stream input
            bloch = bloch_data[frame][q]
            a = complex(bloch[0], bloch[1])
            b = complex(bloch[2], bloch[3])
            bx = 2 * (a.real * b.real + a.imag * b.imag)
            by = 2 * (a.real * b.imag - a.imag * b.real)
            bz = (np.abs(a)**2) - (np.abs(b)**2)
            bloch_vectors[q].set_data([0, bx], [0, by])
            bloch_vectors[q].set_3d_properties([0, bz])

    # Refresh probability state bar heights
    for bar, c_amp in zip(bars, current_probs):
        bar.set_height(np.abs(c_amp) ** 2)

    # Process phase disk values
    for i in range(state_dim):
        c_amp = state_vectors[frame][i]  # Direct complex amplitude from global vector stream
        amplitude = np.abs(c_amp)        # Amplitude magnitude determining radius scaling
        phase = np.angle(c_amp)          # Extracted complex wave phase angle in radians

        # Normalize scaling limits
        r_scale = np.clip(amplitude, 0.0, 1.0) * max_disk_radius
        phase_disks[i].set_radius(r_scale)
        
        # Map angular polar coordinates to 2D Cartesian display vectors
        dx = r_scale * np.cos(phase)
        dy = r_scale * np.sin(phase)
        
        # Set phase arrow data
        phase_hands[i].set_data([i, i + dx], [disk_y_pos, disk_y_pos + dy])

    # Update global diagnostic STIRAP history lines
    y1_vals = np.clip(np.max(stirap_data[:frame+1, :, 1], axis=1), 0.1, None)
    y2_vals = np.clip(np.max(stirap_data[:frame+1, :, 3], axis=1), 0.1, None)
    
    laser1_line.set_data(times[:frame+1], y1_vals)
    laser2_line.set_data(times[:frame+1], y2_vals)
    time_marker.set_xdata([times[frame], times[frame]])

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