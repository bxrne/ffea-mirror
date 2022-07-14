import os
import shutil
import sys
import glob
import subprocess
import numpy as np
import ffeatools.modules.FFEA_rod as FFEA_rod


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
    print("Analysing output of FFEA simulation")
    rod1 = FFEA_rod.FFEA_rod(filename="outputs/1.rodtraj")
    rod2 = FFEA_rod.FFEA_rod(filename="outputs/2.rodtraj")

    force1 = rod1.steric_force[1, :]
    force2 = rod2.steric_force[1, :]

    # force on central (3rd) element
    force_12 = force1[6:9]
    force_13 = force1[9:12]
    force_22 = force2[6:9]
    force_23 = force2[9:12]

    # check intra-element forces are equal (force is applied halfway along element)
    if np.linalg.norm(force_12 - force_13) < 1e-6:
        print("Rod 1, nodes 2 and 3: force magnitudes are equal")
        intra_element = True
    else:
        print("FAIL: Rod 1, nodes 2 and 3: force magnitudes are NOT equal")
        print("node 2: ", force_12)
        print("node 3: ", force_13)
        print("delta: ", np.linalg.norm(force_12 - force_13))
        intra_element = False

    if np.linalg.norm(force_22 - force_23) < 1e-6:
        print("Rod 2, nodes 2 and 3: force magnitudes are equal")
        intra_element = True
    else:
        print("FAIL: Rod 2, nodes 2 and 3: force magnitudes are NOT equal")
        print("node 2: ", force_22)
        print("node 3: ", force_23)
        print("error: ", np.linalg.norm(force_22 - force_23))
        intra_element = False

    # check inter-element forces are symmetric
    # force_12 + force_13 = -(force_22 + force_23)
    if np.linalg.norm(force_12 + force_13 + force_22 + force_23) < 1e-6:
        print("Inter-element forces are symmetric")
        inter_element = True
    else:
        print("FAIL: Inter-element forces are NOT symmetrice")
        print("rod 1, node 2: ", force_12)
        print("rod 1, node 3: ", force_13)
        print("rod 2, node 2: ", force_22)
        print("rod 2, node 3: ", force_23)
        print("delta: ", np.linalg.norm(force_12 + force_13 + force_22 + force_23))
        inter_element = False

    # check force only occurs on central elements
    force1_outer = np.delete(force1, [6, 7, 8, 9, 10, 11])
    force2_outer = np.delete(force2, [6, 7, 8, 9, 10, 11])
    if np.any(force1_outer != 0):
        print("FAIL: Forces applied to non-central nodes (rod 1)")
        print("non-central, non-zero forces: ", force1_outer[force1_outer!=0])
        print("on nodes: ", list(np.where(force1_outer!=0)[0]))
        central_element = False
    else:
        print("Only central nodes have steric force (rod 1)")
        central_element = True

    if np.any(force2_outer != 0):
        print("FAIL: Forces applied to non-central nodes (rod 2)")
        print("non-central, non-zero forces: ", force2_outer[force2_outer!=0])
        print("on nodes: ", list(np.where(force2_outer!=0)[0]))
        central_element = False
    else:
        print("Only central nodes have steric force (rod 2)")
        central_element = True

    # TODO: check forces match BoE calculation
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
