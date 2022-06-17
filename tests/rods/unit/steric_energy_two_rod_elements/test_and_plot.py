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
#from site import enablerlcompleter
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
        ep = data[:, 13:16]
        en = data[:, 16:19]

        if(ep.shape != en.shape):
            print("Error reading file - array shape mismatch")
            print(ep.shape)
            print(en.shape)
            return 1

        # === distance vs. element energy ===#
        m=1.5
        plt.plot(cab[::5]*1e9, ep[:, 0::3][::5], "ms", ms=m, label="$U_{0,+x}$")
        plt.plot(cab[::5]*1e9, ep[:, 1::3][::5], "ys", ms=m, label="$U_{0,+y}$")
        plt.plot(cab[::5]*1e9, ep[:, 2::3][::5], "gs", ms=m, label="$U_{0,+z}$")
        plt.plot(cab[::5]*1e9, en[:, 0::3][::5], "bs", ms=m, label="$U_{0,-x}$")
        plt.plot(cab[::5]*1e9, en[:, 1::3][::5], "cs", ms=m, label="$U_{0,-y}$")
        plt.plot(cab[::5]*1e9, en[:, 2::3][::5], "rs", ms=m, label="$U_{0,-z}$")
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
        epm = np.mean(ep, axis=1)
        enm = np.mean(en, axis=1)
        mean = np.mean(np.array([epm, enm]), axis=0)
        m0, _ = np.polyfit(cab, mean, 1)

        # === step vs. element energy === #
        steps = np.arange(0, len(ep[:, 0::3]))
        # node 0
        plt.plot(steps[::5], ep[:, 0::3][::5], "ms", ms=m, label="$U_{0,+x}$")
        plt.plot(steps[::5], ep[:, 1::3][::5], "ys", ms=m, label="$U_{0,+y}$")
        plt.plot(steps[::5], ep[:, 2::3][::5], "gs", ms=m, label="$U_{0,+z}$")
        plt.plot(steps[::5], en[:, 0::3][::5], "bs", ms=m, label="$U_{0,-x}$")
        plt.plot(steps[::5], en[:, 1::3][::5], "cs", ms=m, label="$U_{0,-y}$")
        plt.plot(steps[::5], en[:, 2::3][::5], "rs", ms=m, label="$U_{0,-z}$")
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

        return [ra0, ra1, rb0, rb1, cab, ep, en], m0


def auto():
    test_name = os.getcwd().split("/")[-1]
    # The FFEA test is called from within python, to allow plotting of results
    subprocess.call(["ffea", test_name + ".ffeatest"])
    read_data, slope = plot()
    u = np.array([read_data[5], read_data[6]])
    u_tolerance = 50*kT
    test_pass = True

    if (u > u_tolerance).any():
        print("FAIL - Steric repulsion energy should be smaller than {0:.1f} kT".format(u_tolerance/kT))
        test_pass = False
    elif (u < 0).any():
        print("FAIL - Steric repulsion energy should be positive")
        test_pass = False

    if (slope > 0):
        print("FAIL - Gradient of rod-rod distance vs. energy plot should be negative")
        test_pass = False

    print("\tMaximum energy: {0:e} J | {1:.2f} kT".format(u.max(), u.max()/kT))
    print("\tMean gradient on element: {0:e} J/m | {1:.2f} pN | {2:.2f} kT / nm".format(slope, slope*1e12, slope/(kT/1e-9)))

    if test_pass:
        return 0
    else:
        print("\tView plots for further info")
        return 1


if __name__ == "__main__":
    error_status = auto()
    sys.exit(error_status)
