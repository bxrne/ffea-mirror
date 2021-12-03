# Create four rods for the rod_neighbour_list_construction unit test
# Rods are straight lines in x, y, z
# Rods 0, 1, 2 are all close to eachother, but rod 3 is isolated
#
# neighbours:
# rod 0: 1 and 2 (two neighbouring elements per element)
# rod 1: 0 (one neighbouring element per element)
# rod 2: 0 "         "          "            "
# rod 3: none

import ffeatools.modules.FFEA_rod as rod
import numpy as np

def r_func(t):
    return t

num_nodes = 4
rods = []

# Create rods
for i in range(4):
    my_rod = rod.FFEA_rod(num_elements=num_nodes)
    
    my_nodes = rod.rod_creator.create_rod_parametric(r_func, r_func, r_func, 0, 1e-8, num_nodes)
    
    my_rod.current_r[0] = my_nodes
    my_rod.equil_r[0] = my_nodes
    
    # 5.77 nm
    xx = my_rod.current_r[0][0][0] - my_rod.current_r[0][1][0]
    yy = my_rod.current_r[0][0][1] - my_rod.current_r[0][1][0]
    zz = my_rod.current_r[0][0][2] - my_rod.current_r[0][1][2]
    
    cutoff_radius = np.sqrt(xx**2 + yy**2 + zz**2)
    print("cutoff_radius = {0} nm".format(cutoff_radius/1e-9))
    
    my_rod.current_m[0], my_rod.equil_m[0] = rod.rod_creator.create_material_frame(my_rod)
    
    stretch_constant = 3.5e-11 # N
    twist_constant = 5e-29 # N.m^2
    radius = 1e-9 # m
    bend_constant = 3.5e-29 # m^4.Pa
    rod.rod_creator.set_params(my_rod, stretch_constant, twist_constant, radius, bending_modulus=bend_constant)
    
    rods.append(my_rod)
    
# 0, 4.3, -2.9, 58 nm
rod_shift = [0, cutoff_radius*0.75, -cutoff_radius*0.5, cutoff_radius*10]
    
# Shift rods and write to file
for i, my_rod in enumerate(rods):
    my_rod.current_r[0] += rod_shift[i]
    my_rod.write_rod("collider_{0}.rod".format(i))