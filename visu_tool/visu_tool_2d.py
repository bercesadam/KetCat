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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
        captions, times, frames, qubits = [], [], [], []
        while True:
            # Caption length + payload
            raw_len = f.read(2)
            if not raw_len: break
            cap_len = struct.unpack("<H", raw_len)[0]
            caption = f.read(cap_len).decode("utf-8")
            
            # Timestamp
            t_raw = f.read(8)
            if not t_raw: break
            t = struct.unpack("<d", t_raw)[0]
            
            # Bloch vectors
            q_raw = f.read(32) 
            if len(q_raw) < 32: break
            q_vals = struct.unpack("<4d", q_raw)
            alpha = complex(q_vals[0], q_vals[1])
            beta = complex(q_vals[2], q_vals[3])
            
            # State vector
            payload = f.read(frame_bytes)
            if len(payload) < frame_bytes: break
            arr = np.frombuffer(payload, dtype="<f8")

            captions.append(caption)
            times.append(t)
            frames.append(arr)
            qubits.append((alpha, beta))

    times  = np.array(times)
    frames = np.array(frames)
    if mode == 0:
        data = (frames[:, 0::2] + 1j * frames[:, 1::2])
    else:
        data = frames
    return captions, times, data, qubits, mode,

captions, times, raw, qubits, mode = load_kwf(kwf_file)
n_timesteps = raw.shape[0]

if mode == 0:
    psi_all = raw.reshape(n_timesteps, Ny, Nx)
else:
    psi_all = np.sqrt(raw).reshape(n_timesteps, Ny, Nx).astype(complex)

def phase_to_rgb(psi):
    phase = (np.angle(psi) + np.pi) / (2 * np.pi)
    phase = phase - np.floor(phase)
    amplitude = np.abs(psi)**2
    if amplitude.max() > 0: amplitude /= amplitude.max()
    amplitude = np.power(amplitude, 0.6)
    r, g, b = np.zeros_like(phase), np.zeros_like(phase), np.zeros_like(phase)
    # Egyszerűsített színezés a fázishoz
    m1 = (phase < 0.25); f1 = phase[m1]/0.25
    r[m1], g[m1], b[m1] = 1.0, 0.0, f1
    m2 = (phase >= 0.25) & (phase < 0.50); f2 = (phase[m2]-0.25)/0.25
    r[m2], g[m2], b[m2] = 1.0-f2, 0.0, 1.0
    m3 = (phase >= 0.50) & (phase < 0.75); f3 = (phase[m3]-0.50)/0.25
    r[m3], g[m3], b[m3] = 0.0, f3, 1.0
    m4 = (phase >= 0.75); f4 = (phase[m4]-0.75)/0.25
    r[m4], g[m4], b[m4] = f4, 1.0-f4, 1.0-f4
    return np.clip(np.stack((r*amplitude, g*amplitude, b*amplitude), axis=-1), 0, 1)

# Setup plot
fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
im = ax.imshow(phase_to_rgb(psi_all[0]), origin="lower", extent=[-L, L, -L, L], interpolation="bilinear")

# Bloch-sphere Inset
ax_bloch = fig.add_axes([0.5, 0.00, 0.4, 0.4], projection='3d')
ax_bloch.set_facecolor((0, 0, 0, 0))  # Teljesen áttetsző háttér

def setup_bloch(axis):
    axis.set_axis_off()
    
    # Wireframe sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    axis.plot_wireframe(x, y, z, color="gray", alpha=0.3, linewidth=0.5)
    
    # Axes
    axis.plot([-1.1, 1.1], [0, 0], [0, 0], color="white", alpha=0.05) # X
    axis.plot([0, 0], [-1.1, 1.1], [0, 0], color="white", alpha=0.05) # Y
    axis.plot([0, 0], [0, 0], [-1.1, 1.1], color="white", alpha=0.1)  # Z

    # Label offset
    o = 1.2 
    
    # Z axis labels
    axis.text(0, 0, o, r"$|0\rangle$", color="white", ha="center", va="bottom", fontsize=9, alpha=0.6)
    axis.text(0, 0, -o, r"$|1\rangle$", color="white", ha="center", va="top", fontsize=9, alpha=0.6)
    
    # X axis labels
    axis.text(o, 0, 0, r"$|+\rangle$", color="white", ha="left", va="center", fontsize=8, alpha=0.4)
    axis.text(-o, 0, 0, r"$|-\rangle$", color="white", ha="right", va="center", fontsize=8, alpha=0.4)
    
    # Y axis labels
    axis.text(0, o, 0, r"$|i\rangle$", color="white", ha="left", va="center", fontsize=8, alpha=0.4)
    axis.text(0, -o, 0, r"$-|i\rangle$", color="white", ha="right", va="center", fontsize=8, alpha=0.4)
    
    # Set up POV
    axis.view_init(elev=20, azim=45)

setup_bloch(ax_bloch)

# Set up vector (will be updated in each cycle)
bloch_vector, = ax_bloch.plot([0, 0], [0, 0], [0, 0], color="cyan", linewidth=2.5, marker='o', markersize=4)

# Colorbar and captions
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.1)
cb = fig.colorbar(cm.ScalarMappable(norm=Normalize(-np.pi, np.pi), cmap=ListedColormap(phase_to_rgb(np.exp(1j*np.linspace(-np.pi, np.pi, 256)).reshape(1,-1))[0])), cax=cax)
cb.set_ticks([-np.pi, 0, np.pi])

caption_text = ax.text(0.01, 0.99, captions[0].replace("|", "\n"), transform=ax.transAxes, fontsize=10, color="white", verticalalignment="top", bbox=dict(facecolor="black", alpha=0.4))

def update(frame):
    im.set_data(phase_to_rgb(psi_all[frame]))
    caption_text.set_text(captions[frame].replace("|", "\n"))
    
    # Bloch vector calculation
    a, b = qubits[frame]
    bx = 2 * (a.real * b.real + a.imag * b.imag)
    by = 2 * (a.real * b.imag - a.imag * b.real)
    bz = (np.abs(a)**2) - (np.abs(b)**2)
    
    bloch_vector.set_data_3d([0, bx], [0, by], [0, bz])
    
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
