# Ryan Cocking 2023
from pathlib import Path
import shutil
import sys
from omegaconf import OmegaConf


def copy_template(filename: str):
    template = Path("template.ffea")
    if not template.exists():
        raise Exception(f"{template} not found")
    shutil.copy(template.name, filename)
    if not Path(filename).exists():
        raise Exception(f"Could not copy {template.name} to {filename}")
    return Path(filename)


def write_test_config_yaml(r: float, l: float, config_filename: str):
    config = {
        "parallel": {
            "+x": [r, 0, 0, 0, 0, 0],
            "-x": [-r, 0, 0, 0, 0, 0],
            "+y": [0, r, 0, 0, 0, 0],
            "-y": [0, -r, 0, 0, 0, 0],
            "+z": [0, 0, l + r, 0, 0, 0],
            "-z": [0, 0, -l - r, 0, 0, 0],
        },
        # "perp": {
        # },
        # "oblique":{
        # }
    }

    conf = OmegaConf.create(config)
    yaml_str = OmegaConf.to_yaml(conf)
    with open(config_filename, "w") as f:
        f.write(yaml_str)

    return config


def edit_ffea_field(file_line: str, template_val: str, new_val: str) -> str:
    """Edits a field string, e.g. "<foo = bar>", read from a .ffea file."""
    whitespace = file_line.split("<")[0]
    field = file_line.split("<")[-1].split(">")[0]  # assumes no comments in file
    field.replace(template_val, new_val)
    return f"{whitespace:s}<{field:s}>\n"


def edit_ffea_file(ffea_file_path: str, trans: tuple[float], rot: tuple[float]):
    """ "Edits rod fields in a .ffea file in-place. Assumes template has a certain format."""
    if len(trans) and len(rot) != 3:
        raise Exception("Translation and rotation vectors must be length 3")

    count = 1
    name = ffea_file_path.split("/")[-1].split(".")[0]

    with open(ffea_file_path, "r") as fr:
        lines = fr.readlines()

    for i, line in enumerate(lines):
        if "output" and "rodtraj" in line:
            lines[i] = edit_ffea_field(line, f"{count:d}.rodtraj", f"{name:s}.rodtraj")
            count += 1

        if "centroid_pos" in line:
            lines[i] = edit_ffea_field(
                line, "(0.0,0.0,0.0)", f"({trans[0]:.2f},{trans[1]:.2f},{trans[2]:.2f})"
            )

        if "rotation" in line:
            lines[i] = edit_ffea_field(
                line, "(0.0,0.0,0.0)", f"({rot[0]:.2f},{rot[1]:.2f},{rot[2]:.2f})"
            )

    with open(ffea_file_path, "w") as fw:
        fr.writelines(lines)

    return lines


def main():
    # define rod configurations by hand, check in PyMOL

    params = OmegaConf.load("params.yml")
    radius = params["radius"]
    length = params["length"]
    test_config = write_test_config_yaml(
        radius, length, "test_configdihaiadshadjada.yml"
    )

    for key1, val1 in test_config.items():
        for key2, val2 in val1.items():

            path = params["in_dir"] + f"/{key1}_{key2}.ffea"
            copy_template(path)

            # in-line replace new .ffea file according to the config
    pass


if __name__ == "__main__":
    main()
