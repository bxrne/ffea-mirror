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


def check_node_distances(rod1, rod2):
    """
    Return true if the average node-node distance is greater than the radius
    sum of both rods for most of the trajectory duration, but within a sensible
    upper limit.
    """

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


def plot_displacement(rod1, rod2):

    # Displacement from initial position
    disp1 = rod1.current_r - rod1.current_r[0]
    disp2 = rod2.current_r - rod2.current_r[0]

    print("Plotting displacement...")
    for frame in range(rod1.num_frames):
        for node in range(rod1.num_elements):
            displacement_rod1 = np.linalg.norm(disp1[frame, node])
            displacement_rod2 = np.linalg.norm(disp2[frame, node])
            plt.plot(frame, displacement_rod1*1e9, marker='o', color='red', ms=0.8)
            plt.plot(frame, displacement_rod2*1e9, marker='o', color='blue', ms=0.8)

    plt.xlabel("Step")
    plt.ylabel("Distance of node from initial position (nm)")

    plt.savefig("Distance_vs_Time.png", dpi=300)
    print("Saved figure 'Distance_vs_Time.png'")
    plt.close()
    print("Done")


def plot_force(rod1, rod2):

    force1 = rod1.steric_force
    force2 = rod2.steric_force

    print("Plotting force...")
    for frame in range(rod1.num_frames):
        for node in range(rod1.num_elements):
            force_rod1 = np.linalg.norm(force1[frame, node])
            force_rod2 = np.linalg.norm(force2[frame, node])
            plt.plot(frame, force_rod1*1e12, marker='o', color='red', ms=0.8)
            plt.plot(frame, force_rod2*1e12, marker='o', color='blue', ms=0.8)

    plt.xlabel("Step")
    plt.ylabel("Steric force on node (pN)")

    plt.savefig("StericForce_vs_Time.png", dpi=300)
    print("Saved figure 'StericForce_vs_Time.png'")
    plt.close()
    print("Done")


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
    print("Calling FFEA from Python script")
    ffea_return_status = subprocess.call(["ffea", test_name + ".ffeatest"])

    # load both rod trajectories
    rod1 = FFEA_rod.FFEA_rod(filename="outputs/1.rodtraj")
    rod2 = FFEA_rod.FFEA_rod(filename="outputs/2.rodtraj")

    plot_displacement(rod1, rod2)
    plot_force(rod1, rod2)

    # if check_node_distances(rod1, rod2) and ffea_return_status == 0:
    #     return 0

    print("THIS TEST IS INCOMPLETE")

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
