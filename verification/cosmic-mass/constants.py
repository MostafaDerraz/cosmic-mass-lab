"""Constants used by verify_cmt.py for 'Mass and Charge as One General Number' (CODATA 2022, SI units).

Mohr, Newell, Taylor and Tiesinga, Rev. Mod. Phys. 97, 025002 (2025).
"""
import math

c = 299792458.0                     # speed of light, m/s (exact)
e = 1.602176634e-19                 # elementary charge, C (exact)
hbar = 1.054571817e-34              # reduced Planck constant, J s (exact to the digits shown)
G = 6.67430e-11                     # Newtonian constant of gravitation, m^3 kg^-1 s^-2, relative uncertainty 2.2e-5
G_REL_UNC = 0.00015 / 6.67430
mu0 = 1.25663706127e-6              # vacuum magnetic permeability, N A^-2
eps0 = 8.8541878188e-12             # vacuum electric permittivity, F/m
K = 1.0 / (4.0 * math.pi * eps0)    # Coulomb constant, N m^2 C^-2
alpha_fs = 7.2973525643e-3          # fine-structure constant
a0 = 5.29177210544e-11              # Bohr radius, m
Ry_eV = 13.605693122990             # Rydberg energy hc R_inf, eV
u = 1.66053906892e-27               # atomic mass constant (dalton), kg
me = 9.1093837139e-31               # electron mass, kg
mp = 1.67262192595e-27              # proton mass, kg
mn = 1.67492750056e-27              # neutron mass, kg
me_u = 5.485799090441e-4            # electron mass, u
mp_u = 1.0072764665789              # proton mass, u
mn_u = 1.00866491606                # neutron mass, u

GU = math.sqrt(K / G)               # |G|, the size of the Composer Unit, kg/C
mS = e * GU                         # the Stoney mass, kg
mP = math.sqrt(hbar * c / G)        # the Planck mass, kg

# The Axiom |G| = |cos p + sin p|, with the length read in units of 1e10 kg/C (Section 4).
ELL = GU / 1e10
ROOT1 = math.asin(ELL / math.sqrt(2.0)) - math.pi / 4.0     # arcsine root, radians (10.1395 degrees)
ROOT2 = math.pi / 2.0 - ROOT1                               # mirror root, radians (79.8605 degrees)
# The theory holds the mirror root at 79.86 degrees, as in the worked examples of the General Numbers
# paper; the root moves by about 0.0009 degree with one standard error of G. Every number of the paper
# that depends on the angle is computed at this held angle.
PHI2 = math.radians(79.86)
PHI1 = math.pi / 2.0 - PHI2
SIN, COS = math.sin(PHI2), math.cos(PHI2)
DEG = 180.0 / math.pi
LEAN_PER_DALTON = DEG * SIN * u / mS                        # degree per dalton per elementary charge


def Be_eV(Z):
    """Total binding energy of the Z electrons of a neutral atom, eV (Lunney, Pearson and Thibault 2003)."""
    return 14.4381 * Z ** 2.39 + 1.55468e-6 * Z ** 5.35
