"""Convert water-SPC/E parameters from the original paper to epsilon and sigma."""

# U = -(A/r)^6+(B/r)^12
# that makes
# sigma = B^2/A and epsilon = (1/4)*(A/B)^12

from decimal import *

A = Decimal(0.37122)  # ((kJ/mol)^1/6)nm
B = Decimal(0.3428)  # (kJ/mol)^1/12nm
kB = Decimal(1.380649 * 1e-23)  # J/K
NA = Decimal(6.02214076 * 1e23)  # 1/mol

sigma = Decimal(10) * B**2 / A  # Ang
epsilon_kJmol = Decimal(0.25) * (A / B) ** 12  # kJ/mol

epsilon_kcalmol = epsilon_kJmol / Decimal(4.184)

epsilon_K = Decimal(1000) * epsilon_kJmol / (NA * kB)


print("sigma is {} Angstrom".format(sigma))
print(
    "epsilon is {} kJ/mol = {} kcal/mol = {} K".format(
        epsilon_kJmol, epsilon_kcalmol, epsilon_K
    )
)
