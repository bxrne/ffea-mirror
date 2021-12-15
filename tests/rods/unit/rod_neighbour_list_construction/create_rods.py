# Create four rods for the rod_neighbour_list_construction unit test
# Rods are straight lines in x with some rotation applied.
# Rods 0, 1, 2 are all close to eachother, but rod 3 is isolated
#
# neighbours:
# rod 0: 1 and 2 (two neighbouring elements per element)
# rod 1: 0 (one neighbouring element per element)
# rod 2: 0 "         "          "            "
# rod 3: none

import ffeatools.modules.FFEA_rod as rod
import numpy as np
import matplotlib.pyplot as plt

def straight(t):
    return t

def zero(t):
    return 0

num_rods = 4
num_nodes = 4
skew_deg = [0, 10, -10 , 5]
rods = []

# Create rods
for i in range(num_rods):
    my_rod = rod.FFEA_rod(num_elements=num_nodes)
    
    my_nodes = rod.rod_creator.create_rod_parametric(straight, zero, zero, 0, 1e-8, num_nodes)
    
    my_rod.current_r[0] = my_nodes
    my_rod.equil_r[0] = my_nodes
    
    # Rotation
    for j in range(0, len(my_rod.current_r[0])):
        theta = skew_deg[i]*np.pi/180
        x = my_rod.current_r[0, j, 0]
        y = my_rod.current_r[0, j, 1]
        my_rod.current_r[0, j, 0] = x*np.cos(theta) - y*np.sin(theta)
        my_rod.current_r[0, j, 1] = x*np.sin(theta) + y*np.cos(theta)
    
    # first element
    x0 = my_rod.current_r[0][0][0]
    x1 = my_rod.current_r[0][1][0]
    y0 = my_rod.current_r[0][0][1]
    y1 = my_rod.current_r[0][1][1]
    z0 = my_rod.current_r[0][0][2]
    z1 = my_rod.current_r[0][1][2]
    p = np.array([x1 - x0, y1 - y0, z1 - z0])
    mag_p = abs(np.linalg.norm(p))
    
    print("Rod", i)
    print(skew_deg[i], "deg")
    print("r0", x0, y0, z0)
    print("r1", x1, y1, z1)
    print("p", p)
    print("|p| = {0:.3f} nm".format(mag_p/1e-9))
    print("r", my_rod.current_r[0])
    print("")
    
    my_rod.current_m[0], my_rod.equil_m[0] = rod.rod_creator.create_material_frame(my_rod)
    
    stretch_constant = 3.5e-11 # N
    twist_constant = 5e-29 # N.m^2
    radius = 1e-9 # m
    bend_constant = 3.5e-29 # m^4.Pa
    rod.rod_creator.set_params(my_rod, stretch_constant, twist_constant, radius, bending_modulus=bend_constant)
    
    rods.append(my_rod)
    

# 0, 4.3, -2.9, 58 nm
rod_shift = [0, mag_p*0.75, -mag_p*0.5, mag_p*10]
    
# Shift rods and write to file
for i, my_rod in enumerate(rods):
    my_rod.current_r[0, :, 1] += rod_shift[i]
    plt.plot(my_rod.current_r[0, :, 0], my_rod.current_r[0, :, 1], 'o', label=i)
    my_rod.write_rod("collider_{0}.rod".format(i))
    
plt.legend()
plt.savefig("Rods.png")