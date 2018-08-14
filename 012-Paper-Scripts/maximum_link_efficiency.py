""""
09/08/2018
Author: Rohan Ramasamy

Script to model the maximum link efficiency achievable.
"""

import numpy as np
from matplotlib import pyplot as plt


def get_link_efficiecny(epsilon, receiver_sizes):
    plt.figure()

    for receiver_size in receiver_sizes:
        z = np.linspace(800e3, 5e6, 1000)
        f = (epsilon * z + receiver_size) ** 2 / receiver_size ** 2
        plt.plot(z * 1e-3, 100.0 / f, label="{}m".format(receiver_size))

    plt.ylabel("Link efficiency [%]")
    plt.xlabel("Range [km]")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    epsilon = 1e-6
    radii = [1.0, 10.0, 100.0]

    get_link_efficiecny(epsilon, radii)

