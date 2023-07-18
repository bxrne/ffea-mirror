//
//  This file is part of the FFEA simulation package
//
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file.
//
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
//
//  To help us fund FFEA development, we humbly ask that you cite
//  the research papers on the package.
//

/*
 *      rod_vdw.cpp
 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#include "rod_vdw.h"

namespace rod
{

    float vdw_energy_6_12(float r_mag_inv, float eps, float sig)
    {
        return 4 * eps * (std::pow(sig * r_mag_inv, 12) - std::pow(sig * r_mag_inv, 6));
    }

    float vdw_force_6_12(float r_mag_inv, float eps, float sig)
    {
        return 24 * eps * r_mag_inv * (2 * std::pow(sig * r_mag_inv, 12) - std::pow(sig * r_mag_inv, 6));
    }

    // r_min_inv = 1 / (2^1/6 * sigma)
    float vdw_energy_interp(float r_mag, float eps, float r_min_inv)
    {
        return eps * (2 * std::pow(r_mag * r_min_inv, 3) - 3 * std::pow(r_mag * r_min_inv, 2));
    }

    // r_min_inv = 1 / (2^1/6 * sigma)
    float vdw_force_interp(float r_mag, float eps, float r_min_inv)
    {
        return 6 * eps * r_min_inv * (std::pow(r_mag * r_min_inv, 2) - r_mag * r_min_inv);
    }

    float force_between_sites(SSINT_matrix *ssint_matrix)
    {
        map<string, scalar> pmap;
        pmap = ssint_matrix->get_SSINT_params(vdw_type_i, vdw_type_j);

        float eps = pmap["Emin"];
        float sigma = pmap["Rmin"];  // mislabelled in LJ_matrix code
    }
}