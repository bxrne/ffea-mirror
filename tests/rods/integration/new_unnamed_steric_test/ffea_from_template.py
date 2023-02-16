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


def write_rod_test_config(r: float, l: float, config_filename: str):
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


def main():
    # define rod configurations by hand, check in PyMOL

    params = OmegaConf.load("params.yml")
    radius = params["radius"]
    length = params["length"]
    test_config = write_rod_test_config(radius, length, "test_config.yml")

    for key1, val1 in test_config.items():
        for key2, val2 in val1.items():

            path = params["in_dir"] + f"/{key1}_{key2}.ffea"
            copy_template(path)

            # in-line replace new .ffea file according to the config
    pass


if __name__ == "__main__":
    main()
