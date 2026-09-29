"""Reproduce the plots in the rotordynamics theory manual."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import numpy as np


# Stationary uniform-shaft reference.
length = 1.0
position = np.linspace(0.0, length, 401)

# Extended flexible-bearing Jeffcott reference.
young_modulus = 2.0e11
shaft_diameter = 0.02
disk_mass = 5.0
disk_diametral_inertia = 0.025
disk_polar_inertia = 0.05
bearing_stiffness = 1.0e6

second_moment = np.pi * shaft_diameter**4 / 64.0
flexural_rigidity = young_modulus * second_moment
translation_stiffness = 1.0 / (
    length**3 / (48.0 * flexural_rigidity)
    + 1.0 / (2.0 * bearing_stiffness)
)
rotation_stiffness = 1.0 / (
    length / (12.0 * flexural_rigidity)
    + 2.0 / (bearing_stiffness * length**2)
)

speeds_rpm = np.linspace(0.0, 6_000.0, 401)
speeds = speeds_rpm * 2.0 * np.pi / 60.0
translation_frequency = np.sqrt(translation_stiffness / disk_mass)
gyroscopic_term = disk_polar_inertia * speeds
discriminant = np.sqrt(
    gyroscopic_term**2
    + 4.0 * disk_diametral_inertia * rotation_stiffness
)
backward_conical = (
    discriminant - gyroscopic_term
) / (2.0 * disk_diametral_inertia)
forward_conical = (
    discriminant + gyroscopic_term
) / (2.0 * disk_diametral_inertia)

figure, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))

for mode in range(1, 4):
    fixed_free_shape = np.sin((2 * mode - 1) * np.pi * position / (2 * length))
    axes[0].plot(
        position,
        fixed_free_shape,
        label=f"Axial/torsional {mode}",
    )
axes[0].plot(
    position,
    np.sin(np.pi * position / length),
    color="black",
    linestyle="--",
    linewidth=1.8,
    label="Flexural 1",
)
axes[0].axhline(0.0, color="0.75", linewidth=0.8)
axes[0].set(
    title="Normalized stationary-shaft mode shapes",
    xlabel="Shaft coordinate z (m)",
    ylabel="Arbitrary normalized amplitude",
)
axes[0].legend(fontsize="small")

translation_hz = translation_frequency / (2.0 * np.pi)
axes[1].plot(
    speeds_rpm,
    np.full_like(speeds_rpm, translation_hz),
    label="Cylindrical pair",
)
axes[1].plot(
    speeds_rpm,
    backward_conical / (2.0 * np.pi),
    label="Backward conical",
)
axes[1].plot(
    speeds_rpm,
    forward_conical / (2.0 * np.pi),
    label="Forward conical",
)
axes[1].plot(
    speeds_rpm,
    speeds_rpm / 60.0,
    color="black",
    linestyle="--",
    label="1x synchronous",
)

cylindrical_critical = translation_frequency
conical_critical = np.sqrt(
    rotation_stiffness / (disk_diametral_inertia + disk_polar_inertia)
)
for critical_speed in (cylindrical_critical, conical_critical):
    critical_rpm = critical_speed * 60.0 / (2.0 * np.pi)
    axes[1].plot(critical_rpm, critical_rpm / 60.0, "o", color="tab:red")

axes[1].set(
    title="Extended-Jeffcott Campbell diagram",
    xlabel="Rotor speed (rpm)",
    ylabel="Vibration frequency (Hz)",
    xlim=(0.0, speeds_rpm[-1]),
    ylim=(0.0, None),
)
axes[1].legend(fontsize="small")

for axis in axes:
    axis.grid(alpha=0.25)

figure.tight_layout()
output = Path(__file__).resolve().parent / "images" / "analytical_rotordynamics.png"
figure.savefig(output, dpi=180, bbox_inches="tight")
print(f"Wrote {output}")
