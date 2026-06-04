import os
import struct
import sys
import tkinter as tk
from tkinter import filedialog

# HEADLESS BACKEND (must be before pyplot import)
import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import pandas as pd

output_dir = "frames"

# Hilbert space properties
Nx, Ny = 256, 256
L = 100.0

fps = 30
interval_ms = 10
save_frames = True

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

def load_kwf(path):
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic != b"KWF\x01":
            raise ValueError(f"Not a valid .kwf file (bad magic: {magic!r})")

        mode = struct.unpack("<B", f.read(1))[0]
        size = struct.unpack("<Q", f.read(8))[0]

        global Nx, Ny, L
        N = struct.unpack("<Q", f.read(8))[0]
        Nx = Ny = N
        L = struct.unpack("<d", f.read(8))[0]

        floats_per_state = 2 if mode == 0 else 1
        frame_floats = floats_per_state * size
        frame_bytes = frame_floats * 8

        captions, times, frames, qubits, stirap = [], [], [], [], []

        while True:
            raw_len = f.read(2)
            if not raw_len:
                break

            cap_len = struct.unpack("<H", raw_len)[0]
            caption = f.read(cap_len).decode("utf-8")

            t_raw = f.read(8)
            if not t_raw:
                break

            t = struct.unpack("<d", t_raw)[0]

            q_raw = f.read(32)
            if len(q_raw) < 32:
                break

            q_vals = struct.unpack("<4d", q_raw)
            alpha = complex(q_vals[0], q_vals[1])
            beta = complex(q_vals[2], q_vals[3])

            s_raw = f.read(32)
            if len(s_raw) < 32:
                break

            s_vals = struct.unpack("<4d", s_raw)

            payload = f.read(frame_bytes)
            if len(payload) < frame_bytes:
                break

            arr = np.frombuffer(payload, dtype="<f8")

            captions.append(caption)
            times.append(t)
            frames.append(arr)
            qubits.append((alpha, beta))
            stirap.append(s_vals)

    times = np.array(times)
    frames = np.array(frames)
    stirap = np.array(stirap)

    if mode == 0:
        data = (frames[:, 0::2] + 1j * frames[:, 1::2])
    else:
        data = frames

    return captions, times, data, qubits, mode, stirap

captions, times, raw, qubits, mode, stirap_data = load_kwf(kwf_file)
n_timesteps = raw.shape[0]

# Generate CSV with alpha/beta data
csv_data = []

for i in range(n_timesteps):
    a, b = qubits[i]
    csv_data.append([times[i], a.real, a.imag, b.real, b.imag])

df_export = pd.DataFrame(
    csv_data,
    columns=['Time', 'Alpha_Real', 'Alpha_Imag', 'Beta_Real', 'Beta_Imag']
)

df_export.to_csv("qubit_data.csv", index=False)
print("CSV exported to qubit_data.csv")

if mode == 0:
    psi_all = raw.reshape(n_timesteps, Ny, Nx)
else:
    psi_all = np.sqrt(raw).reshape(n_timesteps, Ny, Nx).astype(complex)

def phase_to_rgb(psi):
    phase = (np.angle(psi) + np.pi) / (2 * np.pi)
    phase = phase - np.floor(phase)

    amplitude = np.abs(psi) ** 2

    if amplitude.max() > 0:
        amplitude /= amplitude.max()

    amplitude = np.power(amplitude, 1.0)

    r = np.zeros_like(phase)
    g = np.zeros_like(phase)
    b = np.zeros_like(phase)

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

    return np.clip(
        np.stack((r * amplitude, g * amplitude, b * amplitude), axis=-1),
        0,
        1
    )

# -------------------- Auto-scaling extent additions (minimal) --------------------
# Auto-scaling extent params
auto_extent = True                 # turn on/off auto-scaling of the main plot extent
extent_mode = "global"             # "global" (stable) or "frame" (dynamic per frame)
threshold_rel = 0.01               # threshold vs max |psi|^2 to define "content"
padding_frac = 0.05                # padding around detected content box

# Precompute amplitude and grid spacings
ampl = (np.abs(psi_all) ** 2)
dx = 2 * L / Nx
dy = 2 * L / Ny

def compute_extent_for_mask(mask_2d):
    # mask_2d: shape (Ny, Nx), True where content exists
    if not mask_2d.any():
        return [-L, L, -L, L]

    rows = mask_2d.any(axis=1)  # Ny
    cols = mask_2d.any(axis=0)  # Nx

    y0 = int(np.argmax(rows))
    y1 = int(len(rows) - 1 - np.argmax(rows[::-1]))
    x0 = int(np.argmax(cols))
    x1 = int(len(cols) - 1 - np.argmax(cols[::-1]))

    # padding in pixels
    px = int(round(padding_frac * (x1 - x0 + 1)))
    py = int(round(padding_frac * (y1 - y0 + 1)))

    x0 = max(0, x0 - px)
    x1 = min(Nx - 1, x1 + px)
    y0 = max(0, y0 - py)
    y1 = min(Ny - 1, y1 + py)

    # convert pixel indices to edge coordinates
    xmin = -L + x0 * dx
    xmax = -L + (x1 + 1) * dx
    ymin = -L + y0 * dy
    ymax = -L + (y1 + 1) * dy
    return [xmin, xmax, ymin, ymax]

# Decide initial auto extent
if auto_extent and ampl.size and ampl.max() > 0:
    if extent_mode == "global":
        # union of content across all frames
        mask_global = (ampl >= (threshold_rel * ampl.max())).any(axis=0)  # (Ny, Nx)
        auto_extent_vals = compute_extent_for_mask(mask_global)
    else:
        # initialize from frame 0; update() will adjust each frame
        m0 = ampl[0] >= (threshold_rel * max(1e-16, ampl[0].max()))
        auto_extent_vals = compute_extent_for_mask(m0)
else:
    auto_extent_vals = [-L, L, -L, L]
# -----------------------------------------------------------------------------

# Setup plot
fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)

im = ax.imshow(
    phase_to_rgb(psi_all[0]),
    origin="lower",
    extent=auto_extent_vals,  # changed from [-L, L, -L, L]
    interpolation="bilinear"
)

# Keep axes in sync with the imshow extent
ax.set_xlim(auto_extent_vals[0], auto_extent_vals[1])
ax.set_ylim(auto_extent_vals[2], auto_extent_vals[3])

# Bloch-sphere Inset
ax_bloch = fig.add_axes([0.56, 0.01, 0.35, 0.35], projection='3d')
ax_bloch.set_facecolor((0, 0, 0, 0))

# STIRAP Intensity Plot
ax_stirap = fig.add_axes([0.16, 0.11, 0.3, 0.2])
ax_stirap.set_facecolor((0, 0, 0, 0.3))
ax_stirap.tick_params(colors='white', labelsize=8)

for spine in ax_stirap.spines.values():
    spine.set_color('white')

ax_stirap.set_xlabel("Time (sec)", color="white", fontsize=8)
ax_stirap.set_ylabel("Intensity (W/cm²)", color="white", fontsize=8)

l1_val = stirap_data[0, 0]
l2_val = stirap_data[0, 2]

text_l1 = ax_stirap.text(
    0.01, 0.9,
    f"{l1_val:.1f} nm",
    color="lime",
    transform=ax_stirap.transAxes,
    fontsize=8
)

text_l2 = ax_stirap.text(
    0.01, 0.8,
    f"{l2_val:.1f} nm",
    color="deepskyblue",
    transform=ax_stirap.transAxes,
    fontsize=8
)

laser1_line, = ax_stirap.plot([], [], color="lime", linewidth=1.5)
laser2_line, = ax_stirap.plot([], [], color="deepskyblue", linewidth=1.5)

time_marker = ax_stirap.axvline(
    x=times[0],
    color="white",
    linestyle="--",
    alpha=0.5
)

ax_stirap.set_xlim(times[0], times[-1])
ax_stirap.set_ylim(0, np.max(stirap_data[:, [1, 3]]) * 1.1)

def setup_bloch(axis):
    axis.set_axis_off()

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)

    axis.plot_wireframe(
        x, y, z,
        color="gray",
        alpha=0.5,
        linewidth=0.5
    )

    axis.plot([-1.1, 1.1], [0, 0], [0, 0], color="white", alpha=0.05)
    axis.plot([0, 0], [-1.1, 1.1], [0, 0], color="white", alpha=0.05)
    axis.plot([0, 0], [0, 0], [-1.1, 1.1], color="white", alpha=0.1)

    axis.plot([0,1],[0,0],[0,0], color="red")
    axis.plot([0,0],[0,1],[0,0], color="green")
    axis.plot([0,0],[0,0],[0,1], color="blue")

    o = 1.2

    axis.text(0, 0, o, r"$|0\rangle$", color="white",
              ha="center", va="bottom", fontsize=9)

    axis.text(0, 0, -o, r"$|1\rangle$", color="white",
              ha="center", va="top", fontsize=9)

    axis.text(o, 0, 0, r"$|+\rangle$", color="white",
              ha="left", va="center", fontsize=8)

    axis.text(-o, 0, 0, r"$|-\rangle$", color="white",
              ha="right", va="center", fontsize=8)

    axis.text(0, o, 0, r"$|+i\rangle$", color="white",
              ha="left", va="center", fontsize=8)

    axis.text(0, -o, 0, r"$|-i\rangle$", color="white",
              ha="right", va="center", fontsize=8)

    axis.view_init(elev=20, azim=45)

bar_samples = 256
phase_range = np.linspace(-np.pi, np.pi, bar_samples)

colorbar_input = np.exp(1j * phase_range).reshape(1, -1)
colorbar_rgb = phase_to_rgb(colorbar_input)

custom_cmap = ListedColormap(colorbar_rgb[0])

setup_bloch(ax_bloch)

bloch_vector, = ax_bloch.plot(
    [0, 0],
    [0, 0],
    [0, 0],
    color="cyan",
    linewidth=2.5,
    marker='o',
    markersize=4
)

divider = make_axes_locatable(ax)

cax = divider.append_axes("right", size="3%", pad=0.1)

norm = Normalize(vmin=-np.pi, vmax=np.pi)

cb = fig.colorbar(
    cm.ScalarMappable(norm=norm, cmap=custom_cmap),
    cax=cax
)

cb.set_label("Phase (radians)")
cb.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])

cb.set_ticklabels([
    r"$-\pi$",
    r"$-\pi/2$",
    r"$0$",
    r"$\pi/2$",
    r"$\pi$"
])

ax.set_xlabel("x (a.u.)")
ax.set_ylabel("y (a.u.)")
ax.set_title("")

caption_text = ax.text(
    0.01,
    0.99,
    captions[0].replace("|", "\n"),
    transform=ax.transAxes,
    fontsize=12,
    color="white",
    verticalalignment="top",
    horizontalalignment="left",
    multialignment="left",
    bbox=dict(
        boxstyle="round,pad=0.3",
        facecolor="black",
        alpha=0.45,
        edgecolor="none"
    ),
)

def update(frame):
    im.set_data(phase_to_rgb(psi_all[frame]))

    # Dynamic per-frame auto-extent (active only if extent_mode == "frame")
    if auto_extent and extent_mode == "frame":
        frame_max = ampl[frame].max()
        if frame_max > 0:
            mf = ampl[frame] >= (threshold_rel * frame_max)
            new_ext = compute_extent_for_mask(mf)
            im.set_extent(new_ext)
            ax.set_xlim(new_ext[0], new_ext[1])
            ax.set_ylim(new_ext[2], new_ext[3])

    caption_text.set_text(
        captions[frame].replace("|", "\n")
    )

    a, b = qubits[frame]

    bx = 2 * (a.real * b.real + a.imag * b.imag)
    by = 2 * (a.real * b.imag - a.imag * b.real)
    bz = (np.abs(a)**2) - (np.abs(b)**2)

    bloch_vector.set_data_3d(
        [0, bx],
        [0, by],
        [0, bz]
    )

    start = 0

    laser1_line.set_data(
        times[start:frame+1],
        stirap_data[start:frame+1, 1]
    )

    laser2_line.set_data(
        times[start:frame+1],
        stirap_data[start:frame+1, 3]
    )

    time_marker.set_xdata([times[frame], times[frame]])

    text_l1.set_text(
        f"Pump: {stirap_data[frame, 0]:.1f} nm"
    )

    text_l2.set_text(
        f"Stokes: {stirap_data[frame, 2]:.1f} nm"
    )

# HEADLESS EXPORT LOOP
print(f"Exporting {n_timesteps} frames...")

for frame in range(n_timesteps):
    update(frame)

    if save_frames:
        filename = os.path.join(output_dir, f"{frame:04d}.png")

        fig.savefig(
            filename,
            dpi=150
        )

        if frame % 20 == 0:
            print(f"Exporting: {filename}")

print("Done.")
plt.close(fig)
