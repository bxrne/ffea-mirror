import sys
import shutil
import glob
import subprocess
from pathlib import Path
from omegaconf import OmegaConf
import numpy as np


def copy_failed(new_dir: str, failed_path: str, input_dir: str):

    Path(new_dir).mkdir()

    with open(failed_path, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            lines[i] = line[:-1] + ".ffea"
            shutil.copyfile(f"{input_dir:s}/{lines[i]:s}", f"{new_dir:s}/{lines[i]:s}")

    check_copied = len(glob.glob(f"{new_dir:s}/*.ffea"))
    check_failed = len(lines)

    if check_copied != check_failed:
        raise Exception(
            f"Inconsistent copying of failed FFEA files: {check_failed:d} failed, {check_copied:d} copied to dir {new_dir:s}"
        )


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
    subprocess.run([f"{params['py2_exe_path']:s}", "create_straight_rod.py"])
    subprocess.run(["python", "ffea_from_template.py"])

    if params["run_failed_only"]:
        with open(params["fail_path"], "r") as f:
            ffea_files = f.readlines()

        for i, file in enumerate(ffea_files):
            ffea_files[i] = f"{params['in_dir']:s}/{file[:-1]}.ffea"

        print("Running failed configurations only\n")
    else:
        ffea_files = glob.glob(f"{params['in_dir']:s}/*.ffea")

    count = 0
    num_configs = len(ffea_files)
    failed_configs = []

    for i, path in enumerate(sorted(ffea_files), 1):

        name = path.split("/")[-1].split(".")[0]

        print(f"{name:s}\t({i:d}/{num_configs:d})")

        print("FFEA simulation...")
        ffea_result = subprocess.run(
            ["ffea", f"{name:s}.ffea"],
            capture_output=True,
            text=True,
            cwd=params["in_dir"],
        )

        with open(f"{params['out_dir']:s}/{name:s}.stdout", "w") as f:
            f.write(ffea_result.stdout)
        with open(f"{params['out_dir']:s}/{name:s}.stderr", "w") as f:
            f.write(ffea_result.stderr)

        if "thruFail" in name and ffea_result.returncode != 0:
            count += 1
            print("Passed\n")
            continue
        elif "thruFail" not in name and ffea_result.returncode != 0:
            print("Failed\n")
            failed_configs.append(name)
            continue

        print("ffeatools analysis...")
        subprocess.run(
            [
                f"{params['py2_exe_path']:s}",
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

    print(f"Total passed: {count:d}/{num_configs:d}")

    with open(f"{params['out_dir']:s}/failed.txt", "w") as f:
        for item in failed_configs:
            f.write(f"{item:s}\n")

    if Path("FFEA_meta.json").exists():
        Path("FFEA_meta.json").unlink()

    if len(failed_configs) > 0:
        fail_dir = f"{params['in_dir']:s}_failed"
        if Path(fail_dir).exists():
            shutil.rmtree(fail_dir)
        copy_failed(
            new_dir=fail_dir,
            failed_path=params["fail_path"],
            input_dir=params["in_dir"],
        )

    if count == num_configs:
        return 0
    else:
        return 1


if __name__ == "__main__":
    main()
    sys.exit()
