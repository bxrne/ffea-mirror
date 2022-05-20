import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')  # enable plotting over SSH
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

    # measure the node-node distance for like-numbered nodes (0-0, 1-1, etc)
    r1 = rod1.current_r
    r2 = rod2.current_r
    dr = np.linalg.norm(r2 - r1, axis=2)
    # time averaged
    means = np.mean(dr, axis=0)

    radius1 = rod1.material_params[0][0][2]
    radius2 = rod2.material_params[0][0][2]
    radius_sum = radius1 + radius2

    plt.plot(np.arange(means.size), means * 1e9, 'bo', label="Data")
    plt.plot(np.arange(means.size), np.ones(means.size) * radius_sum * 1e9,
             'r--', label="Radius sum")
    plt.plot(np.arange(means.size), np.zeros(means.size), 'k--', lw=0.5)

    plt.legend()
    plt.title("Simulation steps: " + str(r1.shape[0]))
    plt.xlabel("Node index")
    plt.ylabel("Time-averaged node-node distance (nm)")
    plt.savefig("meanDistance_vs_step.png", dpi=300)

    if (means >= radius_sum).all():
        return True

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

    if check_node_distances() and ffea_return_status == 0:
        return 0

    print("Some node-node distances are less than the sum of the rod radii")

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
