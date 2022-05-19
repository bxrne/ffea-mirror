import os
import shutil
import sys
import glob
import subprocess
import ffeatools.modules.FFEA_rod


def main():

    shutil.rmtree("outputs", ignore_errors=True)
    os.makedirs("outputs")
    
    py_return_status = subprocess.call(["python", "create_straight_rod.py"])
    
    if py_return_status != 0:
        sys.exit("Error calling create_straight_rod.py")

    if len(glob.glob("*.ffeatest")) > 1:
        sys.exit("Error: more than one .ffeatest file - cannot determine test name")
    else:
        test_name = glob.glob("*.ffeatest")[0].split(".")[0]

    ffea_return_status = subprocess.call(["ffea", test_name + ".ffeatest"])

    # if FFEA simulation gets to the end, then pass.
    #   search stdout for a raised exception
    #   search trajectory for expected final timestep
    #   else, throw an exception (unable to XYZ)

    print("FFEA returned with error code " + str(ffea_return_status))

    return ffea_return_status

if __name__ == "__main__":
    err = main()
    sys.exit(err)