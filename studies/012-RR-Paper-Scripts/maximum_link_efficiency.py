""""
09/08/2018
Author: Rohan Ramasamy

Script to model the maximum link efficiency achievable.
"""

import numpy as np
from matplotlib import pyplot as plt


def get_link_efficiecny(names, receiver_sizes):
    plt.figure()

    z = np.linspace(800e3, 5e6, 1000)
    for i, receiver_size in enumerate(receiver_sizes):
        for epsilon in [1e-6, 1e-7]:
            print(receiver_size, names[i])
            radii = np.sqrt(receiver_size / np.pi)
            f = (epsilon * z + radii) ** 2 / radii ** 2
            style = "-" if epsilon == 1e-6 else "--"
            c = "r" if names[i] is "Sorato" else "b"
            plt.plot(z * 1e-3, 100.0 / f,
                     label="{} {}$\ \mu Rad$".format(names[i], round(epsilon * 1e6, 2)),
                     linestyle=style, c=c)

    plt.ylabel("Link efficiency [%]")
    plt.xlabel("Range [km]")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    names = ["Sorato", "AMALIA"]
    areas = np.asarray([0.079, 0.366])

    get_link_efficiecny(names, areas)

