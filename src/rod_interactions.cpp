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

 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#include "rod_interactions.h"

namespace rod
{

InteractionData::InteractionData(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float radius_a, float radius_b, float c_a[3], float c_b[3])
{
    rod_id_self = rod_id_a;
    rod_id_nbr = rod_id_b;
    element_id_self = elem_id_a;
    element_id_nbr = elem_id_b;
    radius_self = radius_a;
    radius_nbr = radius_b;
    vec3d(n){contact_self[n] = c_a[n];}
    vec3d(n){contact_nbr[n] = c_b[n];}
};

/** 1) Check that the points of the rod interaction vector, c_a and c_b, lie
 * within their respective rod elements. Certain situations (e.g.
 * almost-parallel rods) will mean this correction can be poor (e.g. c being
 * corrected to completely the wrong end of the rod), so a secondary correction
 * is required.
 *
 * 2) Compare the rod interaction vector, c_ab, to distances measured from nodes
 * to c_a and c_b: d1 = c_b - r_a d2 = c_b - r_a2 d3 = c_a - r_b d4 = c_a - r_b2
 * Find the smallest vector from these and c_ab, and assign that to be the new
 * interaction vector.
*/
void rod_distance_correction(float c_a[3], float c_b[3], float r_a[3],
                                float r_b[3], float p_a[3], float p_b[3],
                                OUT float c_a_out[3], float c_b_out[3])
{
    float r_a2[3] = {0, 0, 0};
    float r_b2[3] = {0, 0, 0};
    float rc_a[3] = {0, 0, 0};
    float rc_b[3] = {0, 0, 0};
    float dot_a = 0;
    float dot_b = 0;
    float p_a_sq = 0;
    float p_b_sq = 0;

    float c_ab[3] = {0, 0, 0};
    float d1[3] = {0, 0, 0};
    float d2[3] = {0, 0, 0};
    float d3[3] = {0, 0, 0};
    float d4[3] = {0, 0, 0};
    float d1_mag = 0;
    float d2_mag = 0;
    float d3_mag = 0;
    float d4_mag = 0;

    vec3d(n) { r_a2[n] = r_a[n] + p_a[n]; }
    vec3d(n) { r_b2[n] = r_b[n] + p_b[n]; }

    // Ensure the points defining the vector lie on their respective finite rods.
    // This part can mis-correct if rods are almost parallel with some tiny angle
    // between them.
    vec3d(n) { rc_a[n] = c_a[n] - r_a[n]; }
    vec3d(n) { rc_b[n] = c_b[n] - r_b[n]; }

    dot_a = dot_product_3x1(p_a, rc_a);
    dot_b = dot_product_3x1(p_b, rc_b);

    p_a_sq = rod::absolute(p_a) * rod::absolute(p_a);
    p_b_sq = rod::absolute(p_b) * rod::absolute(p_b);

    if (dot_a <= 0)
    {
        vec3d(n) { c_a[n] = r_a[n]; }
    }
    else if (dot_a >= p_a_sq)
    {
        vec3d(n) { c_a[n] = r_a2[n]; }
    }

    if (dot_b <= 0)
    {
        vec3d(n) { c_b[n] = r_b[n]; }
    }
    else if (dot_a >= p_b_sq)
    {
        vec3d(n) { c_b[n] = r_b2[n]; }
    }

    // Compare c_ab to vectors pointing from the nodes on one rod to the
    // interaction point on the opposing rod.
    // This part accounts for the mis-correction of the previous section by
    // explicitly working out the shortest distance between the two rods.
    vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
    vec3d(n) { d1[n] = c_b[n] - r_a[n]; }
    vec3d(n) { d2[n] = c_b[n] - r_a2[n]; }
    vec3d(n) { d3[n] = c_a[n] - r_b[n]; }
    vec3d(n) { d4[n] = c_a[n] - r_b2[n]; }

    d1_mag = rod::absolute(d1);
    d2_mag = rod::absolute(d2);
    d3_mag = rod::absolute(d3);
    d4_mag = rod::absolute(d4);

    // TODO: use std::map trick here to bypass if statements (see Jacan)
    // Replace c_ab with the smallest vector. Do nothing if c_ab is already
    // the smallest.
    if (rod::absolute(c_ab) > 0.99 * std::min({d1_mag, d2_mag, d3_mag, d4_mag}))
    {
        if (d1_mag <= std::min({d2_mag, d3_mag, d4_mag}))
        {
            vec3d(n) { c_a[n] = r_a[n]; }
        }
        else if (d2_mag <= std::min(d3_mag, d4_mag))
        {
            vec3d(n) { c_a[n] = r_a2[n]; }
        }
        else if (d3_mag <= d4_mag)
        {
            vec3d(n) { c_b[n] = r_b[n]; }
        }
        else
        {
            vec3d(n) { c_b[n] = r_b2[n]; }
        }
    }

    vec3d(n) { c_a_out[n] = c_a[n]; }
    vec3d(n) { c_b_out[n] = c_b[n]; }

    if (rod::dbg_print)
    {
        std::cout << "correction to rod-rod distance" << std::endl;
        printf("\tp_a.(c_a - r_a) : %.3e\n", dot_a);
        printf("\tp_b.(c_b - r_b) : %.3e\n", dot_b);
        printf("\t|c_ab| : %.3e\n", rod::absolute(c_ab));
        printf("\t|d1| : %.3e\n", d1_mag);
        printf("\t|d2| : %.3e\n", d2_mag);
        printf("\t|d3| : %.3e\n", d3_mag);
        printf("\t|d4| : %.3e\n", d4_mag);
        print_array("\tc_a (corrected)", c_a_out, 3);
        print_array("\tc_b (corrected)", c_b_out, 3);
        std::cout << std::endl;
    }
}

/** Compute one of the two points, c_a and c_b, that form the interaction vector
joining two rod elements together, where c_a sits
on the element p_a. Element radii are not considered at this stage.

\f| \boldsymbol{c}_a = \boldsymbol{r}_a + \frac{(\boldsymbol{r}_b -
\boldsymbol{r}_a)\cdot\boldsymbol{n}_b^p}{\boldsymbol{l}_a\cdot\boldsymbol{n}_b^p}
\ \boldsymbol{l}_a \f|

To account for the infinite line assumption, a correction is applied.
*/
void get_shortest_distance_to_rod(float p_a[3], float p_b[3], float r_a[3],
                                    float r_b[3], OUT float c_a[3],
                                    float c_b[3])
{
    float l_a[3] = {0.0, 0.0, 0.0}; // l_a = p_a / |p_a|
    float l_b[3] = {0.0, 0.0, 0.0};
    float check[3] = {0, 0, 0};
    float l_a_cross_l_b[3] = {0.0, 0.0, 0.0};
    float n_a[3] = {0.0, 0.0, 0.0};
    float n_b[3] = {0.0, 0.0, 0.0};
    float r_ab[3] = {0.0, 0.0, 0.0};
    float r_ba[3] = {0.0, 0.0, 0.0};

    normalize(p_a, l_a);
    normalize(p_b, l_b);

    vec3d(n) { check[n] = l_a[n] - l_b[n]; }
    if (rod::absolute(check) < 1e-7)
    {
        print_array("l_a", l_a, 3);
        print_array("l_b", l_b, 3);
        throw std::runtime_error(
            "Parallel rods detected; distance function will return nan.");
    }

    cross_product(l_a, l_b, l_a_cross_l_b);
    cross_product(l_a, l_a_cross_l_b, n_a);
    cross_product(l_b, l_a_cross_l_b, n_b);

    vec3d(n) { r_ab[n] = r_b[n] - r_a[n]; }
    vec3d(n) { r_ba[n] = r_a[n] - r_b[n]; }

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

    if (rod::dbg_print)
    {
        std::cout << "shortest rod-rod distance" << std::endl;
        print_array("\tp_a", p_a, 3);
        print_array("\tp_b", p_b, 3);
        print_array("\tl_a", l_a, 3);
        print_array("\tl_b", l_b, 3);
        print_array("\tl_a x l_b", l_a_cross_l_b, 3);
        print_array("\tn_b", n_a, 3);
        print_array("\tn_b", n_b, 3);
        print_array("\tr_a", r_a, 3);
        print_array("\tr_b", r_b, 3);
        print_array("\tr_ab", r_ab, 3);
        print_array("\tr_ba", r_ba, 3);
        print_array("\tc_a (initial)", c_a, 3);
        print_array("\tc_b (initial)", c_b, 3);
        std::cout << std::endl;
    }

    // Apply corrections to c_a and/or c_b if appropriate
    rod_distance_correction(c_a, c_b, r_a, r_b, p_a, p_b, c_a, c_b);
}

/*
Check if two rod elements, a and b, interact by calculating the shortest
distance between them and comparing to the sum of their radii. If this
passes, the interaction information is added to both elements' neighbour
lists.
*/
void set_element_neighbours(int rod_id_a, int rod_id_b,
                            int elem_id_a, int elem_id_b,
                            float p_a[3], float p_b[3],
                            float r_a[3], float r_b[3],
                            float radius_a, float radius_b,
                            std::vector<InteractionData> &neighbours_a,
                            std::vector<InteractionData> &neighbours_b)
{

    float c_a[3] = {0};
    float c_b[3] = {0};
    float c_ab[3] = {0};

    rod::get_shortest_distance_to_rod(p_a, p_b, r_a, r_b, c_a, c_b);
    vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }

    if (rod::dbg_print)
    {
        std::cout << "|c_ab|: "
                    << rod::absolute(c_ab) * mesoDimensions::length * 1e9 << " nm"
                    << std::endl;
        std::cout << "radius_a + radius_b: "
                    << (radius_a + radius_b) * mesoDimensions::length * 1e9 << " nm"
                    << std::endl;
    }

    if (rod::absolute(c_ab) < (radius_a + radius_b))
    {
        if (rod::dbg_print)
        {
            std::cout << "assigning neighbours" << std::endl;
        }

        InteractionData stericDataA(
            rod_id_a,
            rod_id_b,
            elem_id_a,
            elem_id_b,
            radius_a,
            radius_b,
            c_a,
            c_b);
        neighbours_a.push_back(stericDataA);

        InteractionData stericDataB(
            rod_id_b,
            rod_id_a,
            elem_id_b,
            elem_id_a,
            radius_b,
            radius_a,
            c_b,
            c_a);
        neighbours_b.push_back(stericDataB);
    }
    else if (rod::dbg_print)
    {
        std::cout << "no interaction" << std::endl;
    }
}

/**
 * @brief Linear steric energy between two rod elements.
 *
 * Normalises the intersection distance such that the repulsive force is
 * equal to the force scaling factor when the elements fully intersect.
 */
float steric_energy_linear(float force_scaling_factor, float intersect_distance,
    float radius_sum)
{
    return force_scaling_factor * intersect_distance / radius_sum;
}

/**
 * @brief Apply a spatial perturbation to the centreline displacement between
 * two rod elements, and return the intersection distance.
*/
float perturbed_intersection_distance(int perturb_dim, float perturb_delta,
    float contact_self[3], float contact_neighb[3], float radius_sum)
{
    float displacement[3] = { 0 };

    contact_neighb[perturb_dim] += perturb_delta;

    vec3d(n) { displacement[n] = contact_neighb[n] - contact_self[n]; }

    return std::max(0.0f, radius_sum - rod::absolute(displacement));
}

/**
 * @brief Steric repulsive force on a rod element due to its collision
 * with a neighbouring element.
 */
std::vector<float> element_steric_force(float delta, float force_strength,
    float radius_sum, float contact_self[3], float contact_neighb[3])
{
    float distance[6] = { 0 };  // +x, -x, +y, -y, +z, -z
    float energy[6] = { 0 };
    std::vector<float> force(3, 0);   // x, y, z

    distance[0] = perturbed_intersection_distance(0,  delta, contact_self, contact_neighb, radius_sum);
    distance[1] = perturbed_intersection_distance(0, -delta, contact_self, contact_neighb, radius_sum);
    distance[2] = perturbed_intersection_distance(1,  delta, contact_self, contact_neighb, radius_sum);
    distance[3] = perturbed_intersection_distance(1, -delta, contact_self, contact_neighb, radius_sum);
    distance[4] = perturbed_intersection_distance(2,  delta, contact_self, contact_neighb, radius_sum);
    distance[5] = perturbed_intersection_distance(2, -delta, contact_self, contact_neighb, radius_sum);

    energy[0] = steric_energy_linear(force_strength, distance[0], radius_sum);
    energy[1] = steric_energy_linear(force_strength, distance[1], radius_sum);
    energy[2] = steric_energy_linear(force_strength, distance[2], radius_sum);
    energy[3] = steric_energy_linear(force_strength, distance[3], radius_sum);
    energy[4] = steric_energy_linear(force_strength, distance[4], radius_sum);
    energy[5] = steric_energy_linear(force_strength, distance[5], radius_sum);

    vec3d(n) { force[n] = -(energy[n] - energy[n + 1]) / delta; }

    if (rod::dbg_print)
    {
        std::cout << "element steric force:\n";
        std::cout << "  delta : " << delta << "\n";
        std::cout << "  force_strength : " << force_strength << "\n";
        std::cout << "  radius_sum : " << radius_sum << "\n";
        print_array("  contact_self", contact_self, 3);
        print_array("  contact_neighb", contact_neighb, 3);
        std::cout << "  -----\n";
        print_array("  perturbed distance", distance, 6);
        print_array("  element energy", energy, 6);
        std::cout << "  dimensions : [+x, -x, y, -y, z, -z]\n";
        print_vector("  element force", force);
        std::cout << "\n";
    }

    return force;
}

/**
 * @brief Interpolate the force applied to an element onto both nodes.
 */
std::vector<float> node_force_interpolation(float contact[3], float node_start[3],
    float element_length, const std::vector<float> &element_force)
{
    float displacement[3] = { 0 };
    float l1 = 0;
    float l2 = 0;
    std::vector<float> force(6, 0);  // x1, y1, z1, x2, y2, z2

    vec3d(n) { displacement[n] = contact[n] - node_start[n]; }

    l1 = rod::absolute(displacement);
    l2 = element_length - l1;

    vec3d(n) { force[n] = element_force[n] * (element_length - l1) / element_length; }
    vec3d(n) { force[n + 3] = element_force[n] * (element_length - l2) / element_length; }

    if (rod::dbg_print)
    {
        std::cout << "node force interpolation:\n";
        print_array("  contact", contact, 3);
        print_array("  node_start", node_start, 3);
        std::cout << "  element_length, L : " << element_length << "\n";
        print_vector("  element_force", element_force);
        std::cout << "  -----\n";
        print_array("  displacement", displacement, 3);
        std::cout << "  l1 : " << l1 << "\n";
        std::cout << "  l2 : " << l2 << "\n";
        std::cout << "  weight start node : " << (element_length - l1) / element_length << "\n";
        std::cout << "  weight end node : " << (element_length - l2) / element_length << "\n";
        print_vector("  node force", force);
        std::cout << "\n";
    }

    return force;
}

//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND I WILL EAT THE
//   (_(_/-(_/   BUGS. WOOF WOOF!

} // namespace rod
