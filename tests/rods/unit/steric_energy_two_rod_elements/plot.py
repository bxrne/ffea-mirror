#!/usr/bin/env python2
# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
# 

import matplotlib
matplotlib.use('Agg')
import os
import sys
import matplotlib.pyplot as plt
import subprocess
import numpy as np

# comma-separated vectors, tab-separated columns
# r_a_0 [f3]    r_a_1 [f3]   r_b_0 [f3]   r_b_1 [f3]   |c_ab| [f]   U_pos_0  [f]  U_neg_0 [f]   U_pos_1 [f]   U_neg_1 [f]
def plot():
    with open("log.txt", "r") as f:
        data = np.loadtxt("log.txt")
        ra0 = data[:, 0:3]
        ra1 = data[:, 3:6]
        rb0 = data[:, 6:9]
        rb1 = data[:, 9:12]
        cab = data[:, 12]
        up0 = data[:, 13:16]
        un0 = data[:, 16:19]
        up1 = data[:, 19:22]
        un1 = data[:, 22:25]

        if(up1.shape != un1.shape):
            print("Error reading file - array shape mismatch")
            print(up1.shape)
            print(un1.shape)
            return 1
    
        m=1
        # node 0
        plt.plot(cab[::5]*1e9, up0[:, 0::3][::5], "ks", ms=m, label="$U_{0,+x}$")
        plt.plot(cab[::5]*1e9, up0[:, 1::3][::5], "ys", ms=m, label="$U_{0,+y}$")
        plt.plot(cab[::5]*1e9, up0[:, 2::3][::5], "gs", ms=m, label="$U_{0,+z}$")
        plt.plot(cab[::5]*1e9, un0[:, 0::3][::5], "bs", ms=m, label="$U_{0,-x}$")
        plt.plot(cab[::5]*1e9, un0[:, 1::3][::5], "cs", ms=m, label="$U_{0,-y}$")
        plt.plot(cab[::5]*1e9, un0[:, 2::3][::5], "rs", ms=m, label="$U_{0,-z}$")
        
        # node 1
        plt.plot(cab[::5]*1e9, up1[:, 0::3][::5], "kd", ms=m, label="$U_{1,+x}$")
        plt.plot(cab[::5]*1e9, up1[:, 1::3][::5], "yd", ms=m, label="$U_{1,+y}$")
        plt.plot(cab[::5]*1e9, up1[:, 2::3][::5], "gd", ms=m, label="$U_{1,+z}$")
        plt.plot(cab[::5]*1e9, un1[:, 0::3][::5], "bd", ms=m, label="$U_{1,-x}$")
        plt.plot(cab[::5]*1e9, un1[:, 1::3][::5], "cd", ms=m, label="$U_{1,-y}$")
        plt.plot(cab[::5]*1e9, un1[:, 2::3][::5], "rd", ms=m, label="$U_{1,-z}$")
        
        plt.plot(cab[::5]*1e9, 1.38e-23*300*np.ones(len(cab[::5])), 'k--', lw=1, label="kT")

        line = f.readlines()[0]
        radius_a = float(line.split(" ")[2])
        radius_b = float(line.split(" ")[5])
        plt.title("radius a: {0:.1f} nm, radius_b: {1:.1f} nm".format(radius_a, radius_b))	
       
        up0m = np.mean(up0, axis=1)
        un0m = np.mean(un0, axis=1)
        up1m = np.mean(up1, axis=1)
        un1m = np.mean(un1, axis=1)
        mean0 = np.mean(np.array([up0m, un0m]), axis=0)
        mean1 = np.mean(np.array([up1m, un1m]), axis=0)
        m0, c0 = np.polyfit(cab, mean0, 1)
        m1, c1 = np.polyfit(cab, mean1, 1)
        # TODO: match this gradient to a result from the code, to pass the test
        print("Mean gradient on node 0: {0:e} J / m".format(m0))
        print("Mean gradient on node 1: {0:e} J / m".format(m1))
 
        plt.xlabel("Rod-rod distance [nm]")
        plt.ylabel("Node energy [J]")
        plt.legend()
        plt.savefig("distance_vs_energy.png", dpi=300)


def auto():
    test_name = os.getcwd().split("/")[-1]
    # The FFEA test is called from within python, to allow plotting of results
    subprocess.call(["ffea", test_name + ".ffeatest"])
    plot()

    return 1
    
if __name__ == "__main__":
    error_status = auto()
    sys.exit(error_status)
