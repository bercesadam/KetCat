import os
import struct
import sys
import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

output_dir = "frames2"

# Hilbert space properties
# Default values, to be overwritten by the loader
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
        frame_bytes  = frame_floats * 8
        captions, times, frames = [], [], []
        while True:
            raw = f.read(2)
            if len(raw) < 2:
                break
            cap_len = struct.unpack("<H", raw)[0]
            caption = f.read(cap_len).decode("utf-8")
            t_raw = f.read(8)
            if len(t_raw) < 8:
                break
            t = struct.unpack("<d", t_raw)[0]
            payload = f.read(frame_bytes)
            if len(payload) < frame_bytes:
                break
            arr = np.frombuffer(payload, dtype="<f8")
            captions.append(caption)
            times.append(t)
            frames.append(arr)
    times  = np.array(times)
    frames = np.array(frames)
    if mode == 0:
        data = (frames[:, 0::2] + 1j * frames[:, 1::2])
    else:
        data = frames
    return captions, times, data, mode, int(size)

captions, times, raw, mode, size = load_kwf(kwf_file)
n_timesteps = raw.shape[0]

if mode == 0:
    psi_all = raw.reshape(n_timesteps, Ny, Nx)
else:
    psi_all = np.sqrt(raw).reshape(n_timesteps, Ny, Nx).astype(complex)

def phase_to_rgb(psi):
    phase = (np.angle(psi) + np.pi) / (2 * np.pi)
    phase = phase - np.floor(phase)
    amplitude = np.abs(psi)**2
    if amplitude.max() > 0:
        amplitude /= amplitude.max()
    amplitude = np.power(amplitude, 0.6)
    r, g, b = np.zeros_like(phase), np.zeros_like(phase), np.zeros_like(phase)
    m1 = (phase < 0.25); f1 = phase[m1]/0.25
    r[m1], g[m1], b[m1] = 1.0, 0.0, f1
    m2 = (phase >= 0.25) & (phase < 0.50); f2 = (phase[m2]-0.25)/0.25
    r[m2], g[m2], b[m2] = 1.0-f2, 0.0, 1.0
    m3 = (phase >= 0.50) & (phase < 0.75); f3 = (phase[m3]-0.50)/0.25
    r[m3], g[m3], b[m3] = 0.0, f3, 1.0
    m4 = (phase >= 0.75); f4 = (phase[m4]-0.75)/0.25
    r[m4], g[m4], b[m4] = f4, 1.0-f4, 1.0-f4
    return np.clip(np.stack((r*amplitude, g*amplitude, b*amplitude), axis=-1), 0, 1)

bar_samples = 256
phase_range = np.linspace(-np.pi, np.pi, bar_samples)
colorbar_input = np.exp(1j * phase_range).reshape(1, -1)
colorbar_rgb = phase_to_rgb(colorbar_input)
custom_cmap = ListedColormap(colorbar_rgb[0])

fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)

im = ax.imshow(
    phase_to_rgb(psi_all[0]),
    origin="lower",
    extent=[-L, L, -L, L],
    interpolation="bilinear"
)
divider = make_axes_locatable(ax) 
cax = divider.append_axes("right", size="3%", pad=0.1)
norm = Normalize(vmin=-np.pi, vmax=np.pi)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=custom_cmap), cax=cax)
cb.set_label("Phase (radians)")
cb.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
cb.set_ticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
ax.set_xlabel("x (a.u.)")
ax.set_ylabel("y (a.u.)")
ax.set_title("")

caption_text = ax.text(
    0.01, 0.99,
    captions[0].replace("|", "\n"),
    transform=ax.transAxes,
    fontsize=12,
    color="white",
    verticalalignment="top",
    horizontalalignment="left",
    multialignment="left",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.45, edgecolor="none"),
)

def update(frame):
    im.set_data(phase_to_rgb(psi_all[frame]))
    caption_text.set_text(captions[frame].replace("|", "\n"))
    if save_frames:
        filename = os.path.join(output_dir, f"{frame:04d}.png")
        plt.savefig(filename, dpi=150)
        if frame % 20 == 0:
            print(f"Exporting: {filename}")
    return [im, caption_text]

ani = FuncAnimation(
    fig,
    update,
    frames=n_timesteps,
    interval=interval_ms,
    repeat=False,
    blit=False
)

plt.show(block=True)
