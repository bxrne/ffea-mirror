import sys
import shutil
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
import numpy as np

import ffeatools.ffea_rod as ffea_rod


def analysis(rod_1_fname: str, rod_2_fname: str):

    rod1 = ffea_rod.ffea_rod(filename=rod_1_fname)
    rod2 = ffea_rod.ffea_rod(filename=rod_2_fname)

    print("ROD 1:")
    print(rod1.vdw_site_type)
    print(rod1.vdw_site_pos)

    print("ROD 2:")
    print(rod2.vdw_site_type)
    print(rod2.vdw_site_pos)


def main():

    params = OmegaConf.load("params.yml")

    subprocess.run([f"python", "create_straight_rod.py"])

    ffea_files = glob.glob(f"*.ffea")
    count = 0
    num_configs = len(ffea_files)

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

        analysis()

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
