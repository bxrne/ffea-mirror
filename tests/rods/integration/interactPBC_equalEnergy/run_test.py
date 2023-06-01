import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import ffeatools.ffea_rod as ffea_rod

def main():

    shutil.rmtree("outputs", ignore_errors=True)
    os.makedirs("outputs")

    py_return_status = subprocess.call(["python", "create_straight_rod.py"])

    if py_return_status != 0:
        sys.exit("Error calling create_straight_rod.py")

    if len(glob.glob("*.ffea")) > 1:
        sys.exit("Error: more than one .ffea file - cannot determine test "
                 "name")
    else:
        ffea_script_name = glob.glob("*.ffea")[0].split(".")[0]

    # test function located in: src/ffea_test.cpp
    print("Calling FFEA from Python script")
    subprocess.call(["ffea", ffea_script_name + ".ffea"])

    # load both rod trajectories
    print("Analysing output of FFEA simulation" + str("\n"))
    rod1a = ffea_rod.ffea_rod(filename="1A.rodtraj")
    rod1b = ffea_rod.ffea_rod(filename="1B.rodtraj")
    rod2a = ffea_rod.ffea_rod(filename="2A.rodtraj")
    rod2b = ffea_rod.ffea_rod(filename="2B.rodtraj")

    r = {
        "1a" : rod1a.current_r[-1],
        "1b" : rod1b.current_r[-1],
        "2a" : rod2a.current_r[-1],
        "2b" : rod2b.current_r[-1]
    }

    ax = plt.figure().add_subplot(projection='3d')
    ax.legend()
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.set_zlabel('z (nm)')

    for rod in r.keys():
        ax.plot(r[rod][:, 0]*1e9, r[rod][:, 1]*1e9, r[rod][:, 2]*1e9, '-', label=rod)

    ax.view_init(elev=20., azim=-35, roll=0)
    plt.legend()
    plt.show()
    # plt.savefig("position_preview.png", dpi=400)

    print("THIS TEST IS INCOMPLETE")

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)