import os
import shutil
import sys
import glob
import subprocess
import numpy as np

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
    rod0 = ffea_rod.ffea_rod(filename="outputs/1.rodtraj")
    rod1 = ffea_rod.ffea_rod(filename="outputs/2.rodtraj")

    force_dim = 2.4364705882352941e-11
    print("Dividing by FFEA force unit: " + str(force_dim) + str("\n"))

    force0 = rod0.steric_force[1, :, :] / force_dim
    force1 = rod1.steric_force[1, :, :] / force_dim

    # force on nodes of central element
    force_02 = force0[2, :]
    force_03 = force0[3, :]
    force_12 = force1[2, :]
    force_13 = force1[3, :]

    # force occurs halfway along element
    print("# === INTRA-ELEMENT FORCES EQUAL? === #")
    if np.linalg.norm(force_02 - force_03) < 1e-6:
        print("Rod 0, nodes 2 and 3: force magnitudes are equal (PASS)")
        intra_element = True
    else:
        print("Rod 0, nodes 2 and 3: force magnitudes are NOT equal (FAIL)")
        intra_element = False

    print("node 2: ", force_02)
    print("node 3: ", force_03)
    print("delta: ", np.linalg.norm(force_02 - force_03))
    print("")

    if np.linalg.norm(force_12 - force_13) < 1e-6:
        print("Rod 1, nodes 2 and 3: force magnitudes are equal (PASS)")
        intra_element = True
    else:
        print("Rod 1, nodes 2 and 3: force magnitudes are NOT equal (FAIL)")
        intra_element = False

    print("node 2: ", force_12)
    print("node 3: ", force_13)
    print("delta: ", np.linalg.norm(force_12 - force_13))
    print("")

    print("# === INTER-ELEMENT FORCES SYMMETRIC? === #")
    # force_02 + force_03 = -(force_12 + force_13)
    if np.linalg.norm(force_02 + force_03 + force_12 + force_13) < 1e-6:
        print("Inter-element forces are symmetric (PASS)")
        inter_element = True
    else:
        print("Inter-element forces are NOT symmetric (FAIL)")
        inter_element = False

    print("rod 0, node 2: ", force_02)
    print("rod 0, node 3: ", force_03)
    print("rod 1, node 2: ", force_12)
    print("rod 1, node 3: ", force_13)
    print("delta: ", np.linalg.norm(force_02 + force_03 + force_12 + force_13))
    print("")

    print("# === FORCE ONLY ON CENTRAL ELEMENTS? === #")
    force0_outer = np.delete(force0, [6, 7, 8, 9, 10, 11])
    force1_outer = np.delete(force1, [6, 7, 8, 9, 10, 11])
    if np.any(force0_outer != 0):
        print("Rod 0, steric force is on NON-central elements (FAIL)")
        print("offending nodes: ", list(np.where(force0_outer != 0)[0] // 3))
        central_element = False
    else:
        print("Rod 0, steric force is on central element only (PASS)")
        central_element = True

    print("Rod 0 forces:")
    print(force0)
    print("")

    if np.any(force1_outer != 0):
        print("Rod 1, steric force is on NON-central elements (FAIL)")
        print("offending nodes: ", list(np.where(force1_outer != 0)[0] // 3))
        central_element = False
    else:
        print("Rod 1, steric force is on central element only (PASS)")
        central_element = True

    print("Rod 1 forces:")
    print(force1)
    print("")

    print("Rod 0 num neighbours:")
    print(rod0.num_neighbours[1, :])
    print("")

    print("Rod 1 num neighbours:")
    print(rod1.num_neighbours[1, :])
    print("")

    # print("# === FORCES MATCH EXPECTED RESULT? === #")
    match_boe = False

    if intra_element and inter_element and central_element and match_boe:
        return 0
    else:
        print("intra-element forces equal: ", intra_element)
        print("inter-element forces symmetric: ", inter_element)
        print("forces only on central elements: ", central_element)
        print("forces match BoE calculation: ", match_boe)

    print("THIS TEST IS INCOMPLETE")

    return 1


if __name__ == "__main__":
    err = main()
    sys.exit(err)
