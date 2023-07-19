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

#include "rod_interactions.h"

namespace rod
{

// Steric interactions
InteractionData::InteractionData(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float radius_a, float radius_b, float c_a[3], float c_b[3],
    float shift[3], float r_a[3], float r_b[3])
{
    rod_id_self = rod_id_a;
    rod_id_nbr = rod_id_b;
    elem_id_self = elem_id_a;
    elem_id_nbr = elem_id_b;
    radius_self = radius_a;
    radius_nbr = radius_b;
    vec3d(n){contact_self[n] = c_a[n];}
    vec3d(n){contact_nbr[n] = c_b[n];}
    vec3d(n){c_ab[n] = c_b[n] - c_a[n];}
    vec3d(n){img_shift[n] = shift[n];}
    vec3d(n){r_self[n] = r_a[n];}
    vec3d(n){r_nbr[n] = r_b[n];}

    epsilon = 0;
    sigma = 0;
    r_min = 0;
    r_min_inv = 0;
};

// VDW interactions
InteractionData::InteractionData(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float c_a[3], float c_b[3], float shift[3], float eps,
    float sig)
{
    rod_id_self = rod_id_a;
    rod_id_nbr = rod_id_b;
    elem_id_self = elem_id_a;
    elem_id_nbr = elem_id_b;
    vec3d(n){contact_self[n] = c_a[n];}
    vec3d(n){contact_nbr[n] = c_b[n];}
    vec3d(n){c_ab[n] = c_b[n] - c_a[n];}
    vec3d(n){img_shift[n] = shift[n];}
    epsilon = eps;
    sigma = sig;
    r_min = std::pow(2, 1/6) * sigma;
    r_min_inv = 1 / r_min;

    radius_self = 0;
    radius_nbr = 0;
    vec3d(n){r_self[n] = 0;}
    vec3d(n){r_nbr[n] = 0;}
};

bool InteractionData::elements_intersect()
{
    float c_ab[3] = {
        this->contact_nbr[0] - this->contact_self[0],
        this->contact_nbr[1] - this->contact_self[1],
        this->contact_nbr[2] - this->contact_self[2]
    };

    if (rod::dbg_print)
    {
        std::cout << "|c_ab|: " << rod::absolute(c_ab) << "\n";
        std::cout << "radius sum: " << this->radius_self + this->radius_nbr << "\n";
    }

    if (rod::absolute(c_ab) < this->radius_self + this->radius_nbr)
        return true;
    else
        return false;
}


/**
 * If an interaction point (c) lies outside the finite length of its rod element,
 * replace its value with the nearest node. Required for certain situations, e.g.
 * nearly-parallel rods.
 */
void finite_length_correction(float c[3], float r[3], float p[3])
{
    float r_c[3] = { 0 };
    float dot = 0;

    vec3d(n) { r_c[n] = c[n] - r[n]; }
    dot = dot_product_3x1(p, r_c);

    if (dot <= 0)
    {
        vec3d(n) { c[n] = r[n]; }
    }
    else if (dot >= rod::absolute(p) * rod::absolute(p))
    {
        vec3d(n) { c[n] = r[n] + p[n]; }
    }

    if (rod::dbg_print)
    {
        std::cout << "finite_length_correction():\n";
        printf("  p.(c - r) : %.3e\n", dot);
    }
}


/** Compare the centreline displacement, c_ab, to the following
 * rod-rod displacements:
 *
 * d1 = c_b - r_a1
 * d2 = c_b - r_a2
 * d3 = c_a - r_b1
 * d4 = c_a - r_b2
 *
 * Find the smallest displacement and, if needed, adjust c_ab such that it
 * becomes the smallest displacement.
*/

void nearest_node_correction(float c_a[3], float c_b[3],
    float r_a[3], float r_b[3], float p_a[3], float p_b[3])
{
    float displ[3][3] = { 0 };
    float r_a2[3] = { 0 };
    float r_b2[3] = { 0 };
    std::vector<float> dist(3, 0);

    vec3d(n) { r_a2[n] = r_a[n] + p_a[n]; }
    vec3d(n) { r_b2[n] = r_b[n] + p_b[n]; }
    vec3d(n) { displ[0][n] = c_b[n] - c_a[n]; }
    vec3d(n) { displ[1][n] = c_b[n] - r_a[n]; }
    vec3d(n) { displ[2][n] = c_b[n] - r_a2[n]; }

    // index corresponding to the smallest distance
    vec3d(n) { dist[n] = rod::absolute(displ[n]); }
    auto iter = std::min_element(dist.begin(), dist.end());
    int index_min = std::distance(dist.begin(), iter);

    if (index_min == 1)
        vec3d(n) { c_a[n] = r_a[n]; }
    else if (index_min == 2)
        vec3d(n) { c_a[n] = r_a2[n]; }
    else if (index_min != 0)
        throw std::out_of_range("Invalid min index in rod::nearest_node_connection()");

    if (rod::dbg_print)
    {
        std::cout << "nearest_node_correction():\n";
        print_array("  r_a1", r_a, 3);
        print_array("  r_a2", r_a2, 3);
        printf("  |c_ab| :       %.3f\n", dist[0]);
        printf("  |c_b - r_a1| : %.3f\n", dist[1]);
        printf("  |c_b - r_a2| : %.3f\n", dist[2]);
        print_array("  c_a (corrected)", c_a, 3);
        std::cout << "\n";
    }
}

/** @brief Compute the two points, c_a and c_b, that form the centreline
 * displacement (interaction vector) joining two rod elements together, where
 * c_a sits on the element p_a.
 *
 * The displacement for oblique rod elements is calculated using an infinite
 * line assumption, so a correction is applied to account for the finite length
 * of elements.
*/
void element_minimum_displacement(float p_a[3], float p_b[3], float r_a[3],
    float r_b[3], float c_a[3], float c_b[3])
{
    float l_a[3] = {0}; // l_a = p_a / |p_a|
    float l_b[3] = {0};
    float l_a_cross_l_b[3] = {0};
    float n_a[3] = {0};
    float n_b[3] = {0};
    float r_ab[3] = {0};
    float r_ba[3] = {0};

    normalize(p_a, l_a);
    normalize(p_b, l_b);
    cross_product(l_a, l_b, l_a_cross_l_b);

    // oblique / perpendicular elements
    if (rod::absolute(l_a_cross_l_b) > 1e-5)
    {
        cross_product(l_a, l_a_cross_l_b, n_a);
        cross_product(l_b, l_a_cross_l_b, n_b);

        vec3d(n) { r_ab[n] = r_b[n] - r_a[n]; }
        vec3d(n) { r_ba[n] = -r_ab[n]; }

        vec3d(n)
        {
            c_a[n] = r_a[n] +
                dot_product_3x1(r_ab, n_b) / dot_product_3x1(l_a, n_b) * l_a[n];
        }
        vec3d(n)
        {
            c_b[n] = r_b[n] +
                dot_product_3x1(r_ba, n_a) / dot_product_3x1(l_b, n_a) * l_b[n];
        }
    }
    // parallel elements
    else
    {
        vec3d(n) { c_a[n] = r_a[n] + 0.5 * p_a[n]; }
        vec3d(n) { c_b[n] = r_b[n] + 0.5 * p_b[n]; }
    }

    if (rod::dbg_print)
    {
        std::cout << "minimum distance between rod elements:\n";
        print_array("  p_a", p_a, 3);
        print_array("  p_b", p_b, 3);
        print_array("  l_a", l_a, 3);
        print_array("  l_b", l_b, 3);
        print_array("  l_a x l_b", l_a_cross_l_b, 3);
        std::cout << "  |l_a x l_b|: " << rod::absolute(l_a_cross_l_b) << "\n";
        print_array("  n_b", n_a, 3);
        print_array("  n_b", n_b, 3);
        print_array("  r_a", r_a, 3);
        print_array("  r_b", r_b, 3);
        print_array("  r_ab", r_ab, 3);
        print_array("  r_ba", r_ba, 3);
        print_array("  c_a (initial)", c_a, 3);
        print_array("  c_b (initial)", c_b, 3);
        std::cout << "\n";
    }

    // enforce finite length of rod element
    finite_length_correction(c_a, r_a, p_a);
    finite_length_correction(c_b, r_b, p_b);

    // if a node is closer than c_a or c_b, use that instead
    nearest_node_correction(c_a, c_b, r_a, r_b, p_a, p_b);
    nearest_node_correction(c_b, c_a, r_b, r_a, p_b, p_a);
}

/*
Return the integer coordinates of the nearest image with which to correct the
displacement of an interacting element pair.
For displacment pointing from A to B, the returned image is that nearest to B,
so the corrected displacement would require subtracting a box length from B
*/
std::vector<int> nearest_periodic_image(float p_a[3], float p_b[3], float r_a[3],
    float r_b[3], std::vector<float> box_dim)
{
    float mid_a[3] = { 0 };
    float mid_b[3] = { 0 };
    float mid_ab[3] = { 0 };
    std::vector<int> img(3, 0);

    rod::get_element_midpoint(p_a, r_a, mid_a);
    rod::get_element_midpoint(p_b, r_b, mid_b);
    vec3d(n) { mid_ab[n] = mid_b[n] - mid_a[n]; }
    vec3d(n) { img.at(n) = std::floor((mid_ab[n] + 0.5 * box_dim[n]) / box_dim[n]); }

    return img;
}

std::vector<int> nearest_periodic_image(float displacement[3], std::vector<float> box_dim)
{
    std::vector<int> img(3, 0);
    vec3d(n) { img.at(n) = std::floor((displacement[n] + 0.5 * box_dim[n]) / box_dim[n]); }
    return img;
}

std::vector<int> nearest_periodic_image(float a[3], float b[3], std::vector<float> box_dim)
{
    std::vector<int> img(3, 0);
    vec3d(n) { img.at(n) = std::floor((b[n] - a[n] + 0.5 * box_dim[n]) / box_dim[n]); }
    return img;
}

/*
Check if two rod elements, a and b, interact by calculating the shortest
distance between them and comparing to the sum of their radii. If this
passes, the interaction information is added to both elements' neighbour
lists.
*/
void set_steric_nbrs(int rod_id_a, int rod_id_b, int elem_id_a,
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

    // PBC correction
    if (periodic)
    {
        img = nearest_periodic_image(p_a, p_b, r_a, r_b, box_dim);
        vec3d(n){shift[n] = box_dim[n] * img.at(n);}
        vec3d(n){r_b[n] -= shift[n];}
    }

    rod::get_element_midpoint(p_a, r_a, mid_a);
    rod::get_element_midpoint(p_b, r_b, mid_b);
    vec3d(n){mid_ab[n] = mid_b[n] - mid_a[n];}

    if (rod::dbg_print)
    {
        std::cout << "cull distant elements:\n";
        std::printf("  elem index a :  %d\n", elem_id_a);
        std::printf("  elem index b :  %d\n", elem_id_b);
        std::printf("  |p_a| :         %.3f\n", rod::absolute(p_a));
        std::printf("  |p_b| :         %.3f\n", rod::absolute(p_b));
        std::printf("  radius sum :    %.3f\n", radius_a + radius_b);
        std::printf("  |midpoint ab| : %.3f\n", rod::absolute(mid_ab));
        std::printf("  cutoff     :    %.3f\n", std::max(rod::absolute(p_a), rod::absolute(p_b)) + radius_a + radius_b);
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

float steric_energy_linear(float force_scaling_factor, float intersect_distance)
{
    return force_scaling_factor * intersect_distance;
}

float steric_energy_squared(float force_scaling_factor, float intersect_distance)
{
    return force_scaling_factor * intersect_distance * intersect_distance;
}

// Perturbation in +x/-x/+y/-y/+z/-z
float intersection_distance(int dim, float delta, float c_a[3], float c_b[3],
    float radius_sum)
{
    float c_ab[3] = { 0 };
    c_b[dim] += delta;
    vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }

    return rod::absolute(c_ab) - radius_sum;
}

// Perturbation in +c_ab/-c_ab, where c_ab points away from current element
float intersection_distance(float delta, float c_a[3], float c_b[3], float radius_sum)
{
    float c_ab[3] = { 0 };
    float c_ab_norm[3] = { 0 };
    float c_ab_ptb[3] = { 0 };

    vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
    rod::normalize(c_ab, c_ab_norm);
    vec3d(n) { c_ab_ptb[n] = c_ab[n] + delta * c_ab_norm[n]; }

    return rod::absolute(c_ab_ptb) - radius_sum;
}

// No perturbation
float intersection_distance(float c_a[3], float c_b[3], float radius_sum)
{
    float c_ab[3] = { 0 };
    vec3d(n) { c_ab[n] = c_a[n] - c_b[n]; }

    return rod::absolute(c_ab) - radius_sum;
}


// Return the steric collision energy of a rod element as a vector with the
// format [r+dr, r-dr, r]. Perturbation is applied along the intersection vector.
std::vector<float> element_steric_energy(float delta, float force_strength,
    float radius_sum, float c_a[3], float c_b[3])
{
    float dist[3] = {0};
    float energy[3] = {0};

    dist[0] = intersection_distance(0.5 * delta, c_a, c_b, radius_sum);
    dist[1] = intersection_distance(-0.5 * delta, c_a, c_b, radius_sum);
    dist[2] = intersection_distance(c_a, c_b, radius_sum);

    energy[0] = steric_energy_squared(force_strength, dist[0]);
    energy[1] = steric_energy_squared(force_strength, dist[1]);
    energy[2] = steric_energy_squared(force_strength, dist[2]);

    if (rod::dbg_print)
    {
        float c_ab[3] = { 0 };
        vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
        std::cout << "element steric energy:\n";
        std::cout << "  delta : " << delta << "\n";
        std::cout << "  force_strength : " << force_strength << "\n";
        std::cout << "  radius_sum : " << radius_sum << "\n";
        print_array("  c_a", c_a, 3);
        print_array("  c_b", c_b, 3);
        std::cout << "  |c_ab| : " << rod::absolute(c_ab) << "\n";
        printf("  perturbed distance +ve : %.3e\n", dist[0]);
        printf("  perturbed distance -ve : %.3e\n", dist[1]);
        printf("  distance : %.3e\n", dist[2]);
        printf("  perturbed energy +ve : %.3e\n", energy[0]);
        printf("  perturbed energy -ve : %.3e\n", energy[1]);
        printf("  energy : %.3e\n", energy[2]);
        std::cout << "\n";

        if (rod::absolute(c_ab) > radius_sum)
        {
            std::cout << "ERROR!\n";
            throw std::runtime_error("Centreline distance cannot be "
                "greater than radius sum when calculating steric energy.");
        }
    }

    return std::vector<float> {energy[0], energy[1], energy[2]};
}

// Interpolate the force applied to an element onto both nodes
std::vector<float> node_force_interpolation(float contact[3], float node_start[3],
    float element_length, const std::vector<float> &element_force)
{
    float displacement[3] = { 0 };
    float l1 = 0;
    float l2 = 0;
    std::vector<float> force(6, 0);  // x0, y0, z0, x1, y1, z1

    vec3d(n) { displacement[n] = contact[n] - node_start[n]; }

    l1 = rod::absolute(displacement);
    l2 = element_length - l1;

    if (l1 > element_length)
    {
        std::cout << "L1: " << l1 << "\n";
        std::cout << "element length: " << element_length << "\n";
        throw std::invalid_argument("Displacement along rod element, L1 = |c - r|, cannot be greater than the element length, |p|.");
    }

    vec3d(n) { force[n] = element_force[n] * (element_length - l1) / element_length; }
    vec3d(n) { force[n + 3] = element_force[n] * (element_length - l2) / element_length; }

    if (rod::dbg_print)
    {
        std::cout << "node force interpolation:\n";
        print_array("  contact", contact, 3);
        print_array("  node_start", node_start, 3);
        std::cout << "  element_length, L : " << element_length << "\n";
        print_vector("  element_force", element_force);
        print_array("  displacement", displacement, 3);
        std::cout << "  l1 : " << l1 << "\n";
        std::cout << "  l2 : " << l2 << "\n";
        std::cout << "  weight start node : " << (element_length - l1) / element_length << "\n";
        std::cout << "  weight end node : " << (element_length - l2) / element_length << "\n";
        print_vector("  force start node", std::vector<float>(force.begin(), std::next(force.begin(), 3)) );
        print_vector("  force end node", std::vector<float>(std::next(force.begin(), 3), force.end()) );
        std::cout << "\n";
    }

    return force;
}

//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND I WILL EAT THE
//   (_(_/-(_/   BUGS. WOOF WOOF!

} // namespace rod
