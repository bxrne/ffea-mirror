# Compute and write node to node distances between two rods in a N^2 loop
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

try:
    import ffeatools.modules.FFEA_rod as rod
except ModuleNotFoundError as e:
    print("ModuleNotFoundError:", e)
    if sys.version_info[0] != 2:
        print("Script must be run with a Python 2 executable.")
    sys.exit()


def node_node_distance(node_pos_a, node_pos_b):
    """Arguments are [N, 3] numpy arrays. Return N^2 node-node distance array."""

    if node_pos_a.shape[1] != 3:
        raise Exception("Unexpected node position shape:" + str(node_pos_a.shape))

    result = np.zeros((node_pos_a.shape[0], node_pos_b.shape[0]))
    for i, ra in enumerate(node_pos_a):
        for j, rb in enumerate(node_pos_b):
            result[i, j] = np.linalg.norm(rb - ra)

    return result


def main():

    parser = argparse.ArgumentParser(description="Paths to two .rodtraj files.")
    parser.add_argument("-t1", "--rod_1_traj", type=str, required=True)
    parser.add_argument("-t2", "--rod_2_traj", type=str, required=True)
    parser.add_argument("-r", "--radius", type=float, required=True)
    parser.add_argument("-l", "--length", type=float, required=True)
    parser.add_argument("-o", "--out_dir", type=str, required=True)
    args = parser.parse_args()

    rod1 = rod.FFEA_rod(filename=args.rod_1_traj)
    rod2 = rod.FFEA_rod(filename=args.rod_2_traj)

    dist = node_node_distance(rod1.current_r[-1], rod2.current_r[-1])

    # Discrete colorbar
    cmap = plt.cm.RdBu
    bounds = np.array([0, args.radius * 1e9, args.radius * 2 * 1e9, dist.max() * 1e9])
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    plt.imshow(dist * 1e9, cmap=cmap, norm=norm, interpolation="nearest")
    plt.xlabel("Node index (rod 1)")
    plt.ylabel("Node index (rod 2)")
    plt.title("Node-node distance (nm)")
    plt.colorbar()
    name = args.rod_1_traj.split("/")[-1].split(".")[0][:-2]
    plt.savefig(args.out_dir + "/" + name + "_heatmap_distance.png", dpi=300)

    np.savetxt(args.out_dir + "/" + name + "_distance.txt", dist)


if __name__ == "__main__":
    main()
