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
     * on the element p_a. Element radii are not considered at this stage.
     \f| \boldsymbol{c}_a = \boldsymbol{r}_a + \frac{(\boldsymbol{r}_b -
     \boldsymbol{r}_a)\cdot\boldsymbol{n}_b^p}{\boldsymbol{l}_a\cdot\boldsymbol{n}_b^p}
     \ \boldsymbol{l}_a \f|
     * To account for some bad behaviour arising from the infinite line assumption
     of this equation, a correction function is also called.
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
            throw std::invalid_argument(
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
     Check if two rod elements interact by calculating the shortest distance between
     them and comparing to the sum of their radii. If this passes, the interaction
     vector and radius are appended to both elements' neighbour lists.

     Each interaction with a single neighbour takes up 7 places in
     std::vector<float> associated with a given element in a rod's neighbour list.
     Indices 0-2 and 3-5 are the points located on itself and the other rod,
     respectively, that describe the interaction vector, c_ab. Index 6 is the radius
     of the other element.
    */
    void assign_neighbours_to_elements(float p_a[3], float p_b[3], float r_a[3],
                                       float r_b[3], float radius_a, float radius_b,
                                       std::vector<float> &element_a_neighbours,
                                       std::vector<float> &element_b_neighbours)
    {

        float c_a[3] = {0, 0, 0};
        float c_b[3] = {0, 0, 0};
        float c_ab[3] = {0, 0, 0};

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

            // Increase vector capacity before assignment (optional, might help with
            // memory stuff)
            element_a_neighbours.reserve(element_a_neighbours.size() + 7);
            element_b_neighbours.reserve(element_b_neighbours.size() + 7);

            // Update both neighbour lists with the interaction vector and the radius of
            // the 'other' element
            vec3d(n) { element_a_neighbours.push_back(c_a[n]); }
            vec3d(n) { element_a_neighbours.push_back(c_b[n]); }
            element_a_neighbours.push_back(radius_b);

            vec3d(n) { element_b_neighbours.push_back(c_b[n]); }
            vec3d(n) { element_b_neighbours.push_back(c_a[n]); }
            element_b_neighbours.push_back(radius_a);
        }
        else if (rod::dbg_print)
        {
            std::cout << "no interaction" << std::endl;
        }
    }

    /* Perturb the separation between two sterically interacting rod elements in a
       specific degree of freedom to get the potential energy associated with rod a.

       \f| U_{int,ab} = \alpha[|\boldsymbol{c}_b - \boldsymbol{c}_a| - (R_a + R_b)]
       \f|

       Args:
       - perturbation_amount - the amount of perturbation to do in the numerical
       differentiation.
       - perturbation_dimension - which dimension to get dE/dr in (x,y,z)
       - force_constant - arbitrary coefficient used to scale the severity of the
       steric repulsion [force units]
       - node_a - the 'start' node of the current element on rod a
       - element_a -
       - c_a, c_b - the points forming a straight line between two rods, a and b
       Returns:
       - energies - a 2-element array for energies interpolated onto the start [0]
       and end [1] nodes of element_a
    */
    void get_steric_perturbation_energy(float perturbation_amount,
                                        int perturbation_dimension,
                                        float force_constant, float r_a[3],
                                        float p_a[3], float c_a[3], float c_b[3],
                                        float radius_a, float radius_b,
                                        OUT float energies[2])
    {

        float c_ab[3] = {0, 0, 0};
        float intersect_distance = 0;
        float energy = 0;
        float displacement[3] = {0, 0, 0};
        float weight_start_node = 0;
        float weight_end_node = 0;

        c_b[perturbation_dimension] += perturbation_amount;
        vec3d(n) { c_ab[n] = c_b[n] - c_a[n]; }
        intersect_distance =
            std::max(0.0f, radius_a + radius_b - rod::absolute(c_ab));
        energy = force_constant * intersect_distance;

        // Energy must be interpolated onto nodes of the element on rod a.
        // E.g. if c is located on the start node, r, then the weighting on
        // the start node is 1.
        vec3d(n) { displacement[n] = c_a[n] - r_a[n]; }
        weight_end_node = rod::absolute(displacement) / rod::absolute(p_a);
        weight_start_node = std::max(0.0f, 1 - weight_end_node);
        energies[0] = weight_start_node * energy;
        energies[1] = weight_end_node * energy;

        if (rod::dbg_print)
        {
            std::cout << "steric perturbation energy" << std::endl;
            std::cout << "\tamount: " << perturbation_amount << std::endl;
            std::cout << "\tdimension: " << perturbation_dimension << std::endl;
            std::cout << "\tforce_constant: " << force_constant << std::endl;
            rod::print_array("\tc_a", c_a, 3);
            rod::print_array("\tc_b perturbed", c_b, 3);
            std::cout << "\t|c_ab|: " << rod::absolute(c_ab) << std::endl;
            std::cout << "\tradius sum: " << radius_a + radius_b << std::endl;
            std::cout << "\tintersect_distance: " << intersect_distance << std::endl;
            std::cout << "\t|c_a - r_a|: " << rod::absolute(displacement) << std::endl;
            std::cout << "\t|p_a|: " << rod::absolute(p_a) << std::endl;
            std::cout << "\telement energy: " << energy << std::endl;
            std::cout << "\tnode weights: " << weight_start_node << ", "
                      << weight_end_node << std::endl;
            std::cout << "\tinterpolated node energies: " << energies[0] << ", "
                      << energies[1] << std::endl;
        }
    }

    // Sum steric energy from adjacent elements onto a given node
    void get_steric_energy_on_node(int node_index, int node_min, int node_max,
                                   float *energy_array_positive, // length 2L array
                                   float *energy_array_negative,
                                   OUT float energy_positive[3],
                                   float energy_negative[3])
    {

        int i = node_index * 6;

        if (node_index == node_min)
        {
            // End nodes: energy contributions from one element only
            vec3d(n) { energy_positive[n] = energy_array_positive[i + n]; }
            vec3d(n) { energy_negative[n] = energy_array_negative[i + n]; }
        }
        else if (node_index == node_max - 1)
        {
            vec3d(n) { energy_positive[n] = energy_array_positive[i - 3 + n]; }
            vec3d(n) { energy_negative[n] = energy_array_negative[i - 3 + n]; }
        }
        else
        {
            // Central nodes: energy contributions from two elements
            vec3d(n)
            {
                energy_positive[n] =
                    energy_array_positive[i - 3 + n] + energy_array_positive[i + n];
            }
            vec3d(n)
            {
                energy_negative[n] =
                    energy_array_negative[i - 3 + n] + energy_array_negative[i + n];
            }
        }
    }

    // Sum and normalise the two unit vectors from adjacent elements onto a shared
    // node
    void get_unit_vector_on_node(int node_index, int node_min, int node_max,
                                 float *unit_vector_array, // length L array
                                 OUT float unit_vector[3])
    {

        float sum[3] = {0, 0, 0};

        int i = node_index * 3;

        if (node_index == node_min)
        {
            vec3d(n) { unit_vector[n] = unit_vector_array[i + n]; }
        }
        else if (node_index == node_max - 1)
        {
            vec3d(n) { unit_vector[n] = unit_vector_array[i - 3 + n]; }
        }
        else
        {
            vec3d(n)
            {
                sum[n] = unit_vector_array[i - 3 + n] + unit_vector_array[i + n];
            }
            rod::normalize(sum, unit_vector);
        }
    }

    //    __      _
    //  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
    //   `_/      )  SOURCE CODE AND I WILL EAT THE
    //   (_(_/-(_/   BUGS. WOOF WOOF!

} // namespace rod
