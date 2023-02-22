import sys

try:
    import ffeatools.modules.FFEA_rod as rod
    from ffeatools.modules.FFEA_rod import rod_creator as rc
except ModuleNotFoundError as e:
    print("ModuleNotFoundError:", e)
    if sys.version_info[0] != 2:
        print("Script must be run with a Python 2 executable.")
    sys.exit()


def yaml_to_dict(file_name):
    """Replicate reading a dictionary from a yaml file with OmegaConf"""
    with open(file_name, "r") as f:
        lines = f.readlines()

    params = {}
    for line in lines:
        key = line.split(":")[0]
        value = line.split(":")[-1].strip()

        try:
            if value.isdigit():
                params[key] = int(value)
            else:
                params[key] = float(value)
        except ValueError:
            params[key] = str(value).strip('"')

    return params


# Default rod parameters
params = yaml_to_dict("params.yml")
stretch = 3.5e-11  # N
twist = 5e-29  # N.m^2
bend = 3.5e-29  # m^4.Pa
diameter = 2 * params["radius"]

# Straight line in z-axis
def x_func(t):
    return 0


def y_func(t):
    return 0


def z_func(t):
    return t


def main():

    # Blank rod
    my_rod = rod.FFEA_rod(num_elements=params["num_nodes"])

    my_nodes = rc.create_rod_parametric(
        x_func, y_func, z_func, 0, params["length"], params["num_nodes"]
    )

    my_rod.current_r[0] = my_nodes
    my_rod.equil_r[0] = my_nodes

    my_rod.current_m[0], my_rod.equil_m[0] = rc.create_material_frame(my_rod)

    rc.set_params(
        rod=my_rod,
        stretch_constant=stretch,
        torsion_constant=twist,
        radius=params["radius"],
        bending_modulus=bend,
    )

    my_rod.write_rod(params["in_dir"] + "/z-axis.rod")


if __name__ == "__main__":
    main()
