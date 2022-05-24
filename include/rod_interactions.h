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
 *      rod_interactions.h
 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#ifndef ROD_INTERACTIONS
#define ROD_INTERACTIONS

#include "rod_math_v9.h"

namespace rod
{

void rod_distance_correction(
    float c_a[3],
    float c_b[3],
    float r_a[3],
    float r_b[3],
    float p_a[3],
    float p_b[3],
    OUT float c_a_out[3],
    float c_b_out[3]);

void get_shortest_distance_to_rod(
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    OUT float c_a[3],
    float c_b[3]);

void assign_neighbours_to_elements(
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    float radius_a,
    float radius_b,
    std::vector<float> &element_a_neighbours,
    std::vector<float> &element_b_neighbours);

void get_steric_perturbation_energy(
    float perturbation_amount,
    int perturbation_dimension,
    float force_constant,
    float r_a[3],
    float p_a[3],
    float c_a[3],
    float c_b[3],
    float radius_a,
    float radius_b,
    OUT float energies[2]);

float steric_energy_linear(
    float force_scaling_factor,
    float intersect_distance,
    float radius_sum);

float element_energy_from_perturbation(
    int perturb_dim,
    float perturb_delta,
    float force_scaling_factor,
    float contact_a[3],
    float contact_b[3],
    float radius_sum);

std::array<float, 6> node_steric_force_interpolation(
    float contact[3],
    float node_1[3],
    float element_length,
    float element_force[3]);
}
#endif
