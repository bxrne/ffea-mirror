import sys
import shutil
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
import numpy as np


def main():

    params = OmegaConf.load("params.yml")

    if Path(params["in_dir"]).exists():
        shutil.rmtree(params["in_dir"])
    Path(params["in_dir"]).mkdir()

    if Path(params["out_dir"]).exists():
        shutil.rmtree(params["out_dir"])
    Path(params["out_dir"]).mkdir()
    Path(f"{params['out_dir']:s}/delete").mkdir()

    # NOTE: THE 'python2' EXEC IS A HACK THAT RUNS ON MY HOME PC.
    # PLEASE UPDATE FFEATOOLS TO PYTHON 3 TO FIX THIS
    # Set up simulations
    subprocess.run(["python2", "create_straight_rod.py"])
    subprocess.run(["python", "ffea_from_template.py"])

    ffea_files = sorted(glob.glob(f"{params['in_dir']:s}/*.ffea"))
    count = 0
    num_configs = len(ffea_files)
    failed_configs = []

    for i, path in enumerate(ffea_files, 1):

        name = path.split("/")[-1].split(".")[0]

        print(f"Configuration:\t{name:s}\t({i:d}/{num_configs:d})")

        print("FFEA simulation...")
        result = subprocess.run(
            ["ffea", f"{name:s}.ffea"],
            capture_output=True,
            text=True,
            cwd=params["in_dir"],
        )

        with open(f"{params['out_dir']:s}/{name:s}.stdout", "w") as f:
            f.write(result.stdout)
        with open(f"{params['out_dir']:s}/{name:s}.stderr", "w") as f:
            f.write(result.stderr)

        if "thruFail" in name and result.returncode != 0:
            count += 1
            print("Passed\n")
            continue

        print("ffeatools analysis...")
        subprocess.run(
            [
                "python2",
                "node_node_distance.py",
                "--rod_1_traj",
                f"{params['out_dir']:s}/{name:s}_1.rodtraj",
                "--rod_2_traj",
                f"{params['out_dir']:s}/{name:s}_2.rodtraj",
                "--radius",
                str(params["radius"]),
                "--length",
                str(params["length"]),
                "--out_dir",
                params["out_dir"],
            ]
        )

        dist = np.loadtxt(f"{params['out_dir']:s}/{name:s}_distanceHeatmap.txt")
        if dist[dist <= 2 * params["radius"]].size == 0:
            count += 1
            print("Passed\n")
        else:
            print("Failed\n")
            failed_configs.append(name)

    print(f"Passed: {count:d}/{num_configs:d}")

    with open(f"{params['out_dir']:s}/failed.txt", "w") as f:
        for item in failed_configs:
            f.write(f"{item:s}\n")

    if Path("FFEA_meta.json").exists():
        Path("FFEA_meta.json").unlink()

    if count == num_configs:
        return 0
    else:
        return 1


if __name__ == "__main__":
    main()
    sys.exit()
