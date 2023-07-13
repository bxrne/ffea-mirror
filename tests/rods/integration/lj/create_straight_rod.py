# Create a FFEA / KOBRA rod that is straight in the z-axis.
from omegaconf import OmegaConf
import ffeatools.ffea_rod as ffea_rod
from ffeatools.rod.rod_creator import rod_creator as rc
import ffeatools.ffea_lj as ffea_lj

# Default rod parameters
params = OmegaConf.load("params.yml")
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

    ffea_input_script = params["ffea_input_script"]

    # Blank rod
    my_rod = ffea_rod.ffea_rod(num_elements=params["num_nodes"])

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

    # Since no blobs have been parameterised, we must write the Lennard-Jones matrix
    # vdw_pairs : tList[tuple[int, int, float, float]]
    rc.write_lj_matrix(ffea_input_script, [(0,0,1e-20,1e-9)])

    my_rod.write_rod(params["in_dir"] + "/z-axis.rod")


if __name__ == "__main__":
    main()
