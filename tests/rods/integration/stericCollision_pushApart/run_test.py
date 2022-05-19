import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import ffeatools.modules.FFEA_rod as FFEA_rod


def check_node_distances():
    """
    Return true if the average node-node distance is greater than the radius
    sum of both rods for most of the trajectory duration, but within a sensible
    upper limit.
    """

    # load both rod trajectories
    rod1 = FFEA_rod.FFEA_rod(filename="outputs/1.rodtraj")
    rod2 = FFEA_rod.FFEA_rod(filename="outputs/2.rodtraj")

    print(rod1.get_trajectory_length())
    print(rod2.get_trajectory_length())

    # measure the node-node distance for like-numbered nodes (0-0, 1-1, etc)

    # compute the mean distance of each node pair with time
    mean_node_distances = np.zeros((20, 10))
    # plot (num_nodes lines, distance y-axis, time x-axis)
    for node, dist in enumerate(mean_node_distances):
        time = np.arange(0, dist.size)
        plt.plot(time, dist, label=node)

    plt.legend()
    plt.xlabel("Simulation step")
    plt.ylabel("Mean node-node distance")
    plt.savefig("meanDistance_vs_step.png", dpi=300)

    return False


def main():

    shutil.rmtree("outputs", ignore_errors=True)
    os.makedirs("outputs")

    py_return_status = subprocess.call(["python", "create_straight_rod.py"])

    if py_return_status != 0:
        sys.exit("Error calling create_straight_rod.py")

    if len(glob.glob("*.ffeatest")) > 1:
        sys.exit("Error: more than one .ffeatest file - cannot determine test "
                 "name")
    else:
        test_name = glob.glob("*.ffeatest")[0].split(".")[0]

    # test function located in: src/ffea_test.cpp
    ffea_return_status = subprocess.call(["ffea", test_name + ".ffeatest"])

    print("FFEA returned with error code " + str(ffea_return_status))

    check_node_distances()

    return ffea_return_status


if __name__ == "__main__":
    err = main()
    sys.exit(err)
