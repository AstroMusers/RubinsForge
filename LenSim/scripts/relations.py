import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c, G
from astropy.cosmology import Planck18 as cosmo

# Constants for lensing calculations
c_km_s = c.to('km/s').value
G_kpc = G.to('kpc3 / (Msun s2)').value

# Example lens and source redshift
z_lens = 0.5
z_source = 2.0

# Compute lensing critical surface density (in M_sun / kpc^2)
D_d = cosmo.angular_diameter_distance(z_lens).value  # in Mpc
D_s = cosmo.angular_diameter_distance(z_source).value  # in Mpc
D_ds = cosmo.angular_diameter_distance_z1z2(z_lens, z_source).value  # in Mpc

Sigma_crit = (c_km_s**2 / (4 * np.pi * G_kpc)) * (D_s / (D_d * D_ds)) * 1e6  # convert Mpc to kpc (1e6)

# Define ranges for r_s and r_c in kpc
r_s_values = np.linspace(5, 50, 5)  # scale radius from 5 to 50 kpc (galaxy to group scale)
r_c_fractions = [0, 0.05, 0.1, 0.2, 0.5]  # rc as fractions of r_s

# Approximate Einstein radius scaling
fig, ax = plt.subplots(figsize=(10, 6))

for rc_frac in r_c_fractions:
    theta_Es = []
    for r_s in r_s_values:
        r_c = rc_frac * r_s

        # Mass reduction approximation: as core grows, less mass is enclosed centrally
        mass_reduction_factor = (1 - rc_frac * 0.7)  # slightly softened approximation
        mass_enclosed = r_s**3 * mass_reduction_factor

        # Theta_E scaling: proportional to sqrt(M_enclosed) / D_d (arbitrary scaling for visualization)
        theta_E = 1.5 * (mass_enclosed**0.5) / D_d
        theta_Es.append(theta_E)

    ax.plot(r_s_values, theta_Es, marker='o', label=f'rc = {rc_frac} × rs')

ax.set_xlabel("Scale Radius r_s (kpc)")
ax.set_ylabel("Relative θ_E (arbitrary units)")
ax.set_title("Impact of Core Radius on Einstein Radius vs. Scale Radius")
ax.grid(True)
ax.legend(title="Core fraction")
plt.show()
