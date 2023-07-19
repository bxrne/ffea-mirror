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

struct InteractionData
{
    int rod_id_self;
    int rod_id_nbr;
    int elem_id_self;
    int elem_id_nbr;
    float radius_self;
    float radius_nbr;
    float contact_self[3];
    float contact_nbr[3];
    float c_ab[3];
    float img_shift[3];
    float r_self[3];
    float r_nbr[3];
    float epsilon;
    float sigma;
    float r_min;
    float r_min_inv;

    InteractionData(
        int rod_id_a,
        int rod_id_b,
        int elem_id_a,
        int elem_id_b,
        float rad_a,
        float rad_b,
        float c_a[3],
        float c_b[3],
        float shift[3],
        float r_a[3],
        float r_b[3]);

    InteractionData(
        int rod_id_a,
        int rod_id_b,
        int elem_id_a,
        int elem_id_b,
        float c_a[3],
        float c_b[3],
        float shift[3],
        float eps,
        float sig);

    bool elements_intersect();
};


// previously: snap_to_nodes
void finite_length_correction(
    float (&c)[3],
    float r[3],
    float p[3]);

void nearest_node_correction(
    float (&c_a)[3],
    float c_b[3],
    float r_a[3],
    float r_b[3],
    float p_a[3],
    float p_b[3]);

void element_minimum_displacement(
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    float c_a[3],
    float c_b[3]);

std::vector<int> nearest_periodic_image(
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    std::vector<float> box_dim);

std::vector<int> nearest_periodic_image(
    float displacement[3],
    std::vector<float> box_dim);

std::vector<int> nearest_periodic_image(
    float a[3],
    float b[3],
    std::vector<float> box_dim);

void set_steric_nbrs(
    int rod_id_a,
    int rod_id_b,
    int elem_id_a,
    int elem_id_b,
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    float radius_a,
    float radius_b,
    std::vector<InteractionData> &neighbours_a,
    std::vector<InteractionData> &neighbours_b,
    bool periodic,
    std::vector<float> box_dim);

float steric_energy_linear(
    float force_scaling_factor,
    float intersect_distance);

float steric_energy_squared(
    float force_scaling_factor,
    float intersect_distance);

// Perturbation in +x/-x/+y/-y/+z/-z
float intersection_distance(
    int dim,
    float delta,
    float c_a[3],
    float c_b[3],
    float radius_sum);

// Perturbation in +c_ab/-c_ab
float intersection_distance(
    float delta,
    float c_a[3],
    float c_b[3],
    float radius_sum);

// No perturbation
float intersection_distance(
    float c_a[3],
    float c_b[3],
    float radius_sum);

std::vector<float> element_steric_energy(
    float delta,
    float force_strength,
    float radius_sum,
    float c_a[3],
    float c_b[3]);

std::vector<float> node_force_interpolation(
    float contact[3],
    float node_1[3],
    float element_length,
    const std::vector<float> &element_force);
}
#endif
