# Compute and write node to node distances between two rods in a N^2 loop
import sys
import argparse
import numpy as np

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

    norm = np.zeros((node_pos_a.shape[0], node_pos_b.shape[0]))
    for i, ra in enumerate(rod_a_r):
        for j, rb in enumerate(rod_b_r):
            result[i, j] = np.linalg.norm(rb - ra)

    return norm


def main():

    parser = argparse.ArgumentParser(description="Paths to two .rodtraj files.")
    parser.add_argument("-a", "--rod_a_traj", type=str, required=True)
    parser.add_argument("-b", "--rod_b_traj", type=str, required=True)
    args = parser.parse_args()

    a = rod.FFEA_rod(filename=args.rod_a_traj)
    b = rod.FFEA_rod(filename=args.rod_b_traj)

    norm = distance(a.current_r, b.current_r)


if __name__ == "__main__":
    main()
