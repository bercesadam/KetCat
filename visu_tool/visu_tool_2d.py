import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import hsv_to_rgb

# ==========================
# USER SETTINGS
# ==========================
csv_file = "C:/Users/User/source/repos/bercesadam/QuantumCircuitsinCompiler/out/build/x64-Release/simulation.csv"

Nx = 256   # grid size in x
Ny = 256   # grid size in y

fps = 30
interval_ms = 50

# ==========================
# LOAD CSV
# ==========================
data = np.loadtxt(csv_file, delimiter=",", skiprows=1)

times = data[:, 0]
raw = data[:, 1:]

n_timesteps = raw.shape[0]
expected_cols = Nx * Ny * 2

if raw.shape[1] != expected_cols:
    raise ValueError(
        f"CSV shape mismatch: expected {expected_cols} columns "
        f"(for {Nx}x{Ny} complex grid), got {raw.shape[1]}"
    )

# ==========================
# RECONSTRUCT COMPLEX FIELD
# shape = (T, Ny, Nx)
# ==========================
psi_all = np.empty((n_timesteps, Ny, Nx), dtype=np.complex128)

for t in range(n_timesteps):
    row = raw[t]

    real_part = row[0::2]
    imag_part = row[1::2]

    psi = real_part + 1j * imag_part
    psi_all[t] = psi.reshape((Ny, Nx))


# ==========================
# COLOR MAPPING:
# phase -> hue
# density -> brightness
# ==========================
def complex_to_rgb(psi):
    density = np.abs(psi) ** 2
    phase = np.angle(psi)

    # normalize density frame-wise for visibility
    max_density = density.max()
    if max_density > 0:
        density = density / max_density

    hue = (phase + np.pi) / (2 * np.pi)   # map [-pi, pi] -> [0,1]
    saturation = np.ones_like(hue)
    value = density

    hsv = np.stack((hue, saturation, value), axis=-1)
    rgb = hsv_to_rgb(hsv)

    return rgb


# ==========================
# SETUP PLOT
# ==========================
fig, ax = plt.subplots(figsize=(8, 8))

rgb0 = complex_to_rgb(psi_all[0])

im = ax.imshow(
    rgb0,
    origin="lower",
    interpolation="nearest"
)

title = ax.set_title(f"t = {times[0]:.3f}")

ax.set_xlabel("x")
ax.set_ylabel("y")

# phase colorbar reference
phase_gradient = np.linspace(-np.pi, np.pi, 256)
phase_hsv = np.zeros((20, 256, 3))
phase_hsv[..., 0] = (phase_gradient + np.pi) / (2 * np.pi)
phase_hsv[..., 1] = 1.0
phase_hsv[..., 2] = 1.0

fig2, ax2 = plt.subplots(figsize=(8, 1.5))
ax2.imshow(hsv_to_rgb(phase_hsv), aspect="auto", extent=[-np.pi, np.pi, 0, 1])
ax2.set_title("Phase color scale")
ax2.set_yticks([])
ax2.set_xlabel("Phase [rad]")


# ==========================
# ANIMATION UPDATE
# ==========================
def update(frame):
    rgb = complex_to_rgb(psi_all[frame])
    im.set_data(rgb)
    title.set_text(f"t = {times[frame]:.3f}")
    return [im, title]


ani = FuncAnimation(
    fig,
    update,
    frames=n_timesteps,
    interval=interval_ms,
    blit=False
)

plt.show()