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

    VDWSite::VDWSite(int rodid, int elemid, int siteid, int vdwtype, float lrod, float lelem)
    {
        int rod_id = rodid;
        int elem_id = elemid;
        int site_id = siteid;
        int vdw_type = vdwtype;
        float L_rod = lrod;
        float L_elem = lelem;
    };

    // Position of VDW site is dynamic, due to rod element stretching
    void VDWSite::get_position(float r[3], float p[3], OUT float r_site[3])
    {
        float p_hat[3] = {0};
        rod::normalize(p, p_hat);
        vec3d(n){r_site[n] = L_elem * p_hat[n];}
    }

    void VDWSite::print_info()
    {
        std::printf("rod_id: %d\n",rod_id);
        std::printf("elem_id: %d\n", elem_id);
        std::printf("site_id: %d\n", site_id);
        std::printf("vdw_type: %d\n", vdw_type);
        std::printf("L_rod: %.3f\n", L_rod);
        std::printf("L_elem: %.3f\n", L_elem);
    }

    void VDWSite::print_info(float r[3], float p[3])
    {
        std::printf("rod_id: %d\n", rod_id);
        std::printf("elem_id: %d\n", elem_id);
        std::printf("site_id: %d\n", site_id);
        std::printf("vdw_type: %d\n", vdw_type);
        std::printf("L_rod: %.3f\n", L_rod);
        std::printf("L_elem: %.3f\n", L_elem);
        float r_site[3] = { 0 };
        get_position(r, p, r_site);
        rod::print_array("r_site", r_site, 3);
    }

    void set_vdw_nbrs(VDWSite site_a, VDWSite site_b, float p_a[3], float p_b[3],
        float r_a[3], float r_b[3], float radius_a, float radius_b,
        std::vector<InteractionData>& nbr_a, std::vector<InteractionData>& nbr_b,
        bool periodic, std::vector<float> box_dim, float vdw_cutoff,
        float epsilon, float sigma)
    {
        float c_a[3] = { 0 };
        float c_b[3] = { 0 };
        float c_ab[3] = { 0 };
        std::vector<int> img = { 0, 0, 0 };
        float shift[3] = { 0 };
        float mag = 0;

        site_a.get_position(r_a, p_a, c_a);
        site_b.get_position(r_b, p_b, c_b);

        if (periodic)
        {
            img = nearest_periodic_image(c_a, c_b, box_dim);
            vec3d(n) { shift[n] = box_dim[n] * img.at(n); }
            vec3d(n) { c_b[n] -= shift[n]; }
        }

        vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
        mag = rod::absolute(c_ab);

        if (rod::dbg_print)
            printf("  |c_ab| : %.3e\n", mag);

        if (mag > radius_a + radius_b and mag < vdw_cutoff)
        {
            InteractionData vdwDataA(
                site_a.rod_id,
                site_b.rod_id,
                site_a.elem_id,
                site_b.elem_id,
                c_a,
                c_b,
                shift,
                epsilon,
                sigma);
            nbr_a.push_back(vdwDataA);

            InteractionData vdwDataB(
                site_b.rod_id,
                site_a.rod_id,
                site_b.elem_id,
                site_a.elem_id,
                c_b,
                c_a,
                shift,
                epsilon,
                sigma);
            nbr_b.push_back(vdwDataB);
        }
        else if (rod::dbg_print)
        {
            std::cout << "  culled\n";
        }
    }
}