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

# This script runs the FFEA test, plots some stuff, and passes judgement (0|1)
# based     on analysis of the data

# Stop matplotlib complaining about $DISPLAY not being set
import matplotlib
matplotlib.use('Agg')

import os
import sys
import matplotlib.pyplot as plt
import subprocess
import numpy as np

kT = 1.38064852e-23 * 300

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
    
        # === distance vs. energy ===#
        m=1.5
        # node 0
        plt.plot(cab[::5]*1e9, up0[:, 0::3][::5], "ms", ms=m, label="$U_{0,+x}$")
        plt.plot(cab[::5]*1e9, up0[:, 1::3][::5], "ys", ms=m, label="$U_{0,+y}$")
        plt.plot(cab[::5]*1e9, up0[:, 2::3][::5], "gs", ms=m, label="$U_{0,+z}$")
        plt.plot(cab[::5]*1e9, un0[:, 0::3][::5], "bs", ms=m, label="$U_{0,-x}$")
        plt.plot(cab[::5]*1e9, un0[:, 1::3][::5], "cs", ms=m, label="$U_{0,-y}$")
        plt.plot(cab[::5]*1e9, un0[:, 2::3][::5], "rs", ms=m, label="$U_{0,-z}$")
        # node 1
        plt.plot(cab[::5]*1e9, up1[:, 0::3][::5], "md", ms=m, label="$U_{1,+x}$")
        plt.plot(cab[::5]*1e9, up1[:, 1::3][::5], "yd", ms=m, label="$U_{1,+y}$")
        plt.plot(cab[::5]*1e9, up1[:, 2::3][::5], "gd", ms=m, label="$U_{1,+z}$")
        plt.plot(cab[::5]*1e9, un1[:, 0::3][::5], "bd", ms=m, label="$U_{1,-x}$")
        plt.plot(cab[::5]*1e9, un1[:, 1::3][::5], "cd", ms=m, label="$U_{1,-y}$")
        plt.plot(cab[::5]*1e9, un1[:, 2::3][::5], "rd", ms=m, label="$U_{1,-z}$")
        # kT
        plt.plot(cab[::5]*1e9, kT*np.ones(len(cab[::5])), 'k--', lw=1, label="kT")
        plt.plot(cab[::5]*1e9, 5*kT*np.ones(len(cab[::5])), 'k:', lw=1, label="5kT")

        line = f.readlines()[0]
        radius_a = float(line.split(" ")[2])
        radius_b = float(line.split(" ")[5])
        radius_title = "$R_a$: {0:.1f}, $R_b$: {1:.1f} [nm]".format(radius_a, radius_b)
        plt.title(radius_title)
        plt.xlabel("Rod-rod distance [nm]")
        plt.ylabel("Node energy [J]")
        plt.legend()
        plt.savefig("distance_vs_energy.png", dpi=300)
        plt.yscale("log")
        plt.savefig("distance_vs_logEnergy.png", dpi=300)
        plt.close()

        # gradient
        up0m = np.mean(up0, axis=1)
        un0m = np.mean(un0, axis=1)
        up1m = np.mean(up1, axis=1)
        un1m = np.mean(un1, axis=1)
        mean0 = np.mean(np.array([up0m, un0m]), axis=0)
        mean1 = np.mean(np.array([up1m, un1m]), axis=0)
        m0, c0 = np.polyfit(cab, mean0, 1)
        m1, c1 = np.polyfit(cab, mean1, 1)
        
        # === step vs. energy === #
        steps = np.arange(0, len(up0[:, 0::3]))
        # node 0
        plt.plot(steps[::5], up0[:, 0::3][::5], "ms", ms=m, label="$U_{0,+x}$")
        plt.plot(steps[::5], up0[:, 1::3][::5], "ys", ms=m, label="$U_{0,+y}$")
        plt.plot(steps[::5], up0[:, 2::3][::5], "gs", ms=m, label="$U_{0,+z}$")
        plt.plot(steps[::5], un0[:, 0::3][::5], "bs", ms=m, label="$U_{0,-x}$")
        plt.plot(steps[::5], un0[:, 1::3][::5], "cs", ms=m, label="$U_{0,-y}$")
        plt.plot(steps[::5], un0[:, 2::3][::5], "rs", ms=m, label="$U_{0,-z}$")  
        # node 1
        plt.plot(steps[::5], up1[:, 0::3][::5], "md", ms=m, label="$U_{1,+x}$")
        plt.plot(steps[::5], up1[:, 1::3][::5], "yd", ms=m, label="$U_{1,+y}$")
        plt.plot(steps[::5], up1[:, 2::3][::5], "gd", ms=m, label="$U_{1,+z}$")
        plt.plot(steps[::5], un1[:, 0::3][::5], "bd", ms=m, label="$U_{1,-x}$")
        plt.plot(steps[::5], un1[:, 1::3][::5], "cd", ms=m, label="$U_{1,-y}$")
        plt.plot(steps[::5], un1[:, 2::3][::5], "rd", ms=m, label="$U_{1,-z}$")
        # kT
        plt.plot(steps[::5], kT*np.ones(len(cab[::5])), 'k--', lw=1, label="kT")
        plt.plot(steps[::5], 5*kT*np.ones(len(cab[::5])), 'k:', lw=1, label="5kT")
        plt.xlabel("Step")
        plt.ylabel("Node energy [J]")
        plt.legend()
        plt.savefig("step_vs_energy.png", dpi=300)
        plt.yscale("log")
        plt.savefig("step_vs_logEnergy.png", dpi=300)
        plt.close()

        # === step vs. distance === #
        plt.plot(steps[::5], cab[::5]*1e9, "b-", label="$|c_{ab}|$")
        plt.plot(steps[::5], (radius_a + radius_b)*np.ones(len(cab[::5])), 'k--', label="$R_a+R_b$")
        plt.title(radius_title)
        plt.xlabel("Step")
        plt.ylabel("Rod-rod distance [nm]")
        plt.legend()
        plt.savefig("step_vs_distance.png", dpi=300)

        return [ra0, ra1, rb0, rb1, cab, up0, un0, up1, un1], np.array([m0, m1])


def auto():
    test_name = os.getcwd().split("/")[-1]
    # The FFEA test is called from within python, to allow plotting of results
    subprocess.call(["ffea", test_name + ".ffeatest"])
    read_data, slopes = plot()
    u = np.array([read_data[5], read_data[6], read_data[7], read_data[8]])
    u_tolerance = 5*kT
    test_pass = True
    
    if (u > u_tolerance).any():
        print("FAIL - Steric repulsion energy should be smaller than {0:.1f} kT".format(u_tolerance/kT))
        test_pass = False
    elif (u < 0).any():
        print("FAIL - Steric repulsion energy should be positive")
        test_pass = False

    if (slopes > 0).any():
        print("FAIL - Gradient of rod-rod distance vs. energy plot should be negative")
        test_pass = False

    print("\tMaximum energy: {0:e} J ({1:.1f} kT)".format(u.max(), u.max()/kT))
    print("\tMean gradient on node 0: {0:e} J/m".format(slopes[0]))
    print("\tMean gradient on node 1: {0:e} J/m".format(slopes[1]))

    if test_pass:
        return 0
    else:
        print("\tView plots for further info")
        return 1

    
if __name__ == "__main__":
    error_status = auto()
    sys.exit(error_status)
