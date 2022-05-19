import os
import sys
import glob
import subprocess
import ffeatools.modules.FFEA_rod

def create_rod():
    pass

def main():
    
    if len(glob.glob("*.ffeatest")) > 1:
        print("Error: more than one .ffeatest file - cannot determine test name")
        sys.exit()
    else:
        test_name = glob.glob("*.ffeatest")[0].split(".")[0]
    
    subprocess.call(["ffea", test_name + ".ffeatest"])
    
    test_pass = False
    
    if test_pass:
        return 0
    else:
        return 1

if __name__ == "__main__":
    error_status = main()
    sys.exit(error_status)