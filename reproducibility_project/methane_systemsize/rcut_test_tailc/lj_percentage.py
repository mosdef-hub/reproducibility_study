"""Script to calculate LJ potential/epsilon at various rcuts."""
import matplotlib.pyplot as plt
import numpy as np


def get_lj(r, sigma=3.73, epsilon=148):
    """Return LJ potential given r, sigma, and potential."""
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)


rs = np.linspace(4, 20, 1000)

es = get_lj(rs)

plt.plot(rs, es)
plt.savefig("lj.png", dpi=600)


r = 3.73 * (2 ** (1 / 6))

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))


r = 10

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))


r = 14

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))


r = 18

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))


r = 22

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))


r = 24

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))


r = 25

e = get_lj(r)
epsilon = 148
print("Ratio (e at r = {} and epsilon) is {}".format(r, abs(e / epsilon)))
