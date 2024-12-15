"""Helper script to visualize the surface_source.h5 files created with this test.
"""

import h5py
import matplotlib.pyplot as plt


if __name__ == "__main__":

    # Select an option
    # "show": 3D visualization using matplotlib
    # "savefig": 2D representation using matplotlib and storing the fig under plot_2d.png
    option = "show"
    # option = "savefig"

    # Select the case from its folder name
    folder = "case-20"

    # Reading the surface source file
    with h5py.File(f"{folder}/surface_source_true.h5", "r") as fp:
        source_bank = fp["source_bank"][()]
        r_xs = source_bank["r"]["x"]
        r_ys = source_bank["r"]["y"]
        r_zs = source_bank["r"]["z"]

    print("Size of the source bank: ", len(source_bank))

    # Select data range to visualize
    idx_1 = 0
    idx_2 = -1

    # Show 3D representation
    if option == "show":

        fig = plt.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(projection="3d", proj_type="ortho")
        ax1.scatter(r_xs[idx_1:idx_2], r_ys[idx_1:idx_2], r_zs[idx_1:idx_2], marker=".")
        ax1.view_init(0, 0)
        ax1.xaxis.set_ticklabels([])
        ax1.set_ylabel("y-axis [cm]")
        ax1.set_zlabel("z-axis [cm]")
        ax1.set_aspect("equal", "box")

        plt.show()

    # Save 2D representations
    elif option == "savefig":

        fig = plt.figure(figsize=(14, 5))
        ax1 = fig.add_subplot(121, projection="3d", proj_type="ortho")
        ax1.scatter(r_xs[idx_1:idx_2], r_ys[idx_1:idx_2], r_zs[idx_1:idx_2], marker=".")
        ax1.view_init(0, 0)
        ax1.xaxis.set_ticklabels([])
        ax1.set_ylabel("y-axis [cm]")
        ax1.set_zlabel("z-axis [cm]")
        ax1.set_aspect("equal", "box")

        ax2 = fig.add_subplot(122, projection="3d", proj_type="ortho")
        ax2.scatter(r_xs[idx_1:idx_2], r_ys[idx_1:idx_2], r_zs[idx_1:idx_2], marker=".")
        ax2.view_init(90, -90)
        ax2.zaxis.set_ticklabels([])
        ax2.set_xlabel("x-axis [cm]")
        ax2.set_ylabel("y-axis [cm]")
        ax2.set_aspect("equal", "box")

        plt.savefig("plot_2d.png")
