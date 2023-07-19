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

    VDWSite::VDWSite(int rodid, int elemid, int siteid, int vdwtype, float lrod, float lelem);
    {
        int rod_id = rodid;
        int elem_id = elemid;
        int site_id = siteid;
        int vdw_type = vdwtype;
        float L_rod = lrod;
        float L_elem = lelem;
    };

    void VDWSite::get_position(float r[3], float p[3], OUT float r_site[3])
    {
        float p_hat[3] = {0};
        rod::normalize(p, p_hat);
        vec3d(n){r_site[n] = L_elem * p_hat;}
    };

    set_vdw_nbrs(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float p_a[3], float p_b[3], float r_a[3], float r_b[3],
    float radius_a, float radius_b, std::vector<InteractionData> &neighbours_a,
    std::vector<InteractionData> &neighbours_b, bool periodic, std::vector<float> box_dim)
{
    float mid_a[3] = {0};
    float mid_b[3] = {0};
    float mid_ab[3] = {0};
    float c_a[3] = {0};
    float c_b[3] = {0};
    float c_ab[3] = {0};
    std::vector<int> img = {0, 0, 0};
    float shift[3] = {0};
    float cutoff = 0;

    // PBC correction
    if (periodic)
    {
        img = nearest_periodic_image(c_a, c_b, box_dim);
        vec3d(n){shift[n] = box_dim[n] * img.at(n);}
        vec3d(n){r_b[n] -= shift[n];}
    }

    if (rod::dbg_print)
    {
        std::cout << "cull distant elements:\n";
        std::printf("  elem index a :  %d\n", elem_id_a);
        std::printf("  elem index b :  %d\n", elem_id_b);
        std::printf("  |p_a| :         %.3f\n", rod::absolute(p_a));
        std::printf("  |p_b| :         %.3f\n", rod::absolute(p_b));
        std::printf("  |midpoint ab| : %.3f\n", rod::absolute(mid_ab));
        std::printf("  cutoff     :    %.3f\n", vdw_cutoff);
    }

    if (rod::absolute(mid_ab) < std::max(rod::absolute(p_a), rod::absolute(p_b)) + radius_a + radius_b)
    {
        rod::element_minimum_displacement(p_a, p_b, r_a, r_b, c_a, c_b);
        vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }

        if (rod::dbg_print)
        {
            std::cout << "rod-rod distance:\n";
            printf("  |c_ab| : %.3e\n", rod::absolute(c_ab));
        }

        if (rod::absolute(c_ab) > 1e-5 and rod::absolute(c_ab) < (radius_a + radius_b))
        {
            if (rod::dbg_print)
            {
                std::cout << "  neighbour pair: " << rod_id_a << "|" << elem_id_a
                << " and " << rod_id_b << "|" << elem_id_b << "\n\n";
                std::cout << "  generating interaction structs\n";
            }

            // ! Save modified r_a, r_b, image ID [i,j,k], image origin coordinates [x,y,z]
            InteractionData stericDataA(
                rod_id_a,
                rod_id_b,
                elem_id_a,
                elem_id_b,
                radius_a,
                radius_b,
                c_a,
                c_b,
                shift,
                r_a,
                r_b);
            neighbours_a.push_back(stericDataA);

            InteractionData stericDataB(
                rod_id_b,
                rod_id_a,
                elem_id_b,
                elem_id_a,
                radius_b,
                radius_a,
                c_b,
                c_a,
                shift,
                r_b,
                r_a);
            neighbours_b.push_back(stericDataB);
        }
        else if (rod::absolute(c_ab) <= 1e-5)
        {
            std::cout << "ERROR!\n";
            throw std::runtime_error("Rod elements have fully overlapped.");
        }
    }
    else if (rod::dbg_print)
    {
        std::cout << "  culled\n";
    }
}

    // float force_between_sites(SSINT_matrix *ssint_matrix)
    // {
    //     map<string, scalar> pmap;
    //     pmap = ssint_matrix->get_SSINT_params(vdw_type_i, vdw_type_j);

    //     float eps = pmap["Emin"];
    //     float sigma = pmap["Rmin"];  // mislabelled in LJ_matrix code


    // }
}