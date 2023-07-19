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
 *      rod_vdw.h
 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#ifndef ROD_VDW
#define ROD_VDW

#include "rod_interactions.h"
#include "LJ_matrix.h"

namespace rod
{
    float vdw_energy_6_12(float r_mag_inv, float eps, float sig);
    float vdw_force_6_12(float r_mag_inv, float eps, float sig);
    float vdw_energy_interp(float r_mag, float eps, float r_min_inv);
    float vdw_force_interp(float r_mag, float eps, float r_min_inv);

    struct VDWSite
    {
        int rod_id;
        int elem_id;
        int site_id;
        int vdw_type;
        float L_rod;   // position of site along rod as a fraction of its length, 0 < L < 1
        float L_elem;

        VDWSite(int rodid, int elemid, int siteid, int vdwtype, float lrod, float lelem);

        void get_position(float r[3], float p[3], OUT float r_site[3]);
    };
}