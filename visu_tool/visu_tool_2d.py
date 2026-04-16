import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.cm as cm

# ==========================
# USER SETTINGS
# ==========================
csv_file = "../out/build/x64-Release/simulation.csv"
output_dir = "frames2"

Nx, Ny = 256, 256
# Centered coordinates: -100 to +100 = 200 total size
L = 100.0 

fps = 30
interval_ms = 50
save_frames = True

if save_frames and not os.path.exists(output_dir):
    os.makedirs(output_dir)

# ==========================
# LOAD DATA
# ==========================
df = pd.read_csv(csv_file, index_col=False)
captions = df.iloc[:, 0].astype(str).values
times    = df.iloc[:, 1].values
raw      = df.iloc[:, 2:].values

n_timesteps = raw.shape[0]
psi_all = (raw[:, 0::2] + 1j * raw[:, 1::2]).reshape(n_timesteps, Ny, Nx)

# ==========================
# COLORMAP LOGIC
# ==========================
def phase_to_rgb(psi):
    phase = (np.angle(psi) + np.pi) / (2 * np.pi)
    phase = phase - np.floor(phase)
    
    amplitude = np.abs(psi)**2
    if amplitude.max() > 0:
        amplitude /= amplitude.max()
    amplitude = np.power(amplitude, 0.5)

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

# ==========================
# COLORBAR PREPARATION
# ==========================
bar_samples = 256
phase_range = np.linspace(-np.pi, np.pi, bar_samples)
colorbar_input = np.exp(1j * phase_range).reshape(1, -1)
colorbar_rgb = phase_to_rgb(colorbar_input)
custom_cmap = ListedColormap(colorbar_rgb[0])

# ==========================
# SETUP PLOT
# ==========================
fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(right=0.82) 

# extent=[-100, 100, -100, 100] shifts (0,0) to the center
im = ax.imshow(
    phase_to_rgb(psi_all[0]),
    origin="lower",
    extent=[-L, L, -L, L],
    interpolation="bilinear"
)

cax = fig.add_axes([0.85, 0.15, 0.03, 0.7]) 
norm = Normalize(vmin=-np.pi, vmax=np.pi)
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=custom_cmap), cax=cax)
cb.set_label('Phase (radians)')
cb.set_ticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
cb.set_ticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])

ax.set_xlabel("x (a.u.)")
ax.set_ylabel("y (a.u.)")
title_text = ax.set_title(f"{captions[0]}")

# ==========================
# ANIMATION & SAVING
# ==========================
def update(frame):
    im.set_data(phase_to_rgb(psi_all[frame]))
    title_text.set_text(f"{captions[frame]}")
    
    if save_frames:
        filename = os.path.join(output_dir, f"{frame:04d}.png")
        plt.savefig(filename, dpi=150)
        if frame % 20 == 0:
            print(f"Exporting: {filename}")

    return [im, title_text]

ani = FuncAnimation(
    fig, 
    update, 
    frames=n_timesteps, 
    interval=interval_ms, 
    repeat=False,
    blit=False
)

plt.show()