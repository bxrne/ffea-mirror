import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import ffeatools.ffea_rod as ffea_rod

import matplotlib as mpl
mpl.use('Agg')

def main():

    # subprocess.call(["bash", "cleanup.sh"])
    # py_return_status = subprocess.call(["python", "create_straight_rod.py"])

    # if py_return_status != 0:
    #     sys.exit("Error calling create_straight_rod.py")

    # if len(glob.glob("*.ffea")) > 1:
    #     sys.exit("Error: more than one .ffea file - cannot determine test "
    #              "name")
    # else:
    #     ffea_script_name = glob.glob("*.ffea")[0].split(".")[0]

    # # test function located in: src/ffea_test.cpp
    # print("Calling FFEA from Python script")
    # subprocess.call(["ffea", ffea_script_name + ".ffea"])

    # # load both rod trajectories
    # print("Analysing output of FFEA simulation" + str("\n"))
    # rod1 = ffea_rod.ffea_rod(filename="1A.rodtraj")
    # rod2 = ffea_rod.ffea_rod(filename="1B.rodtraj")

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
