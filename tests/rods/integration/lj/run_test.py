import sys
import shutil
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
import numpy as np

def main():

    params = OmegaConf.load("params.yml")

    subprocess.run([f"python", "create_straight_rod.py"])

    ffea_files = glob.glob(f"*.ffea")
    count = 0
    num_configs = len(ffea_files)

    # ! - Tests assume steric_force_factor = 20 (rod_structure.h)

    for i, path in enumerate(sorted(ffea_files), 1):

        name = path.split("/")[-1].split(".")[0]

        print(f"{name:s}\t({i:d}/{num_configs:d})")
        print("FFEA simulation...")
        ffea_result = subprocess.run(
            ["ffea", f"{name:s}.ffea"],
            capture_output=True,
            text=True,
        )

        print("Writing stdout and stderr...")
        with open(f"{name:s}.stdout", "w") as f:
            f.write(ffea_result.stdout)
        with open(f"{name:s}.stderr", "w") as f:
            f.write(ffea_result.stderr)

        if False:
            count += 1
            print("Passed\n")
        else:
            print("Failed\n")

    print(f"Total passed: {count:d}/{num_configs:d}")

    if Path("FFEA_meta.json").exists():
        Path("FFEA_meta.json").unlink()

    if count == num_configs:
        return 0
    else:
        return 1


if __name__ == "__main__":
    main()
    sys.exit()
