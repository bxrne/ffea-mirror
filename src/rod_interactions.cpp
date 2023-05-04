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

#include "rod_interactions.h"

namespace rod
{

InteractionData::InteractionData(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float radius_a, float radius_b, float c_a[3], float c_b[3])
{
    rod_id_self = rod_id_a;
    rod_id_nbr = rod_id_b;
    elem_id_self = elem_id_a;
    elem_id_nbr = elem_id_b;
    radius_self = radius_a;
    radius_nbr = radius_b;
    vec3d(n){contact_self[n] = c_a[n];}
    vec3d(n){contact_nbr[n] = c_b[n];}
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
 * @brief If either interaction point, c_a or c_b, lies outside the finite
 * length of its respective rod element, replace its value with the nearest
 * node. This is necessary for certain situations, e.g. nearly-parallel rods.
 *
 * @return: 6-element std::vector<float> containing c_a and c_b.
 */
std::vector<float> snap_to_nodes(std::vector<float> c_ab, float r_a[3],
    float r_b[3], float p_a[3], float p_b[3])
{
    float rc_a[3] = {0};
    float rc_b[3] = {0};
    float dot_a = 0;
    float dot_b = 0;

    vec3d(n) { rc_a[n] = c_ab.at(n) - r_a[n]; }
    vec3d(n) { rc_b[n] = c_ab.at(n + 3) - r_b[n]; }
    dot_a = dot_product_3x1(p_a, rc_a);
    dot_b = dot_product_3x1(p_b, rc_b);

    if (dot_a <= 0)
    {
        vec3d(n) { c_ab.at(n) = r_a[n]; }
    }
    else if (dot_a >= rod::absolute(p_a) * rod::absolute(p_a))
    {
        vec3d(n) { c_ab.at(n) = r_a[n] + p_a[n]; }
    }

    if (dot_b <= 0)
    {
        vec3d(n) { c_ab.at(n + 3) = r_b[n]; }
    }
    else if (dot_b >= rod::absolute(p_b) * rod::absolute(p_b))
    {
        vec3d(n) { c_ab.at(n + 3) = r_b[n] + p_b[n]; }
    }

    if (rod::dbg_print)
    {
        std::cout << "snap to nodes:\n";
        printf("  p_a.(c_a - r_a) : %.3e\n", dot_a);
        printf("  p_b.(c_b - r_b) : %.3e\n", dot_b);
    }

    return c_ab;
}


/** @brief Compare the centreline displacement, c_ab, to the following
 * rod-rod displacements:
 *
 * d1 = c_b - r_a1
 * d2 = c_b - r_a2
 * d3 = c_a - r_b1
 * d4 = c_a - r_b2
 *
 * Find the smallest displacement and, if needed, adjust c_ab such that it is
 * the smallest displacement.
 *
 * @return: 6-element std::vector<float> containing c_a and c_b.
*/
std::vector<float> compare_node_distances(std::vector<float> c_ab, float r_a[3],
    float r_b[3], float p_a[3], float p_b[3])
{
    float displacements[5][3] = {0};
    float r_a2[3] = {0};
    float r_b2[3] = {0};
    std::vector<float> distances(5, 0);

    vec3d(n) { r_a2[n] = r_a[n] + p_a[n]; }
    vec3d(n) { r_b2[n] = r_b[n] + p_b[n]; }

    vec3d(n) { displacements[0][n] = c_ab.at(n + 3) - c_ab.at(n); }
    vec3d(n) { displacements[1][n] = c_ab.at(n + 3) - r_a[n]; }
    vec3d(n) { displacements[2][n] = c_ab.at(n + 3) - r_a2[n]; }
    vec3d(n) { displacements[3][n] = c_ab.at(n) - r_b[n]; }
    vec3d(n) { displacements[4][n] = c_ab.at(n) - r_b2[n]; }

    int i = 0;
    for (auto &d : distances)
        d = rod::absolute(displacements[i++]);

    // index corresponding to the smallest distance
    auto iter = std::min_element(distances.begin(), distances.end());
    int index_min = std::distance(distances.begin(), iter);

    switch (index_min)
    {
    case 0:
        break;
    case 1:
        vec3d(n) { c_ab.at(n) = r_a[n]; }
        break;
    case 2:
        vec3d(n) { c_ab.at(n) = r_a2[n]; }
    case 3:
        vec3d(n) { c_ab.at(n + 3) = r_b[n]; }
        break;
    case 4:
        vec3d(n) { c_ab.at(n + 3) = r_b2[n]; }
        break;
    default:
        std::cout << "ERROR!\n";
        throw std::out_of_range("Invalid index to vector of displacements in rod::compare_node_distances()");
    }

    if (rod::dbg_print)
    {
        std::cout << "node to interaction point distance comparison:\n";
        printf("  |c_ab| :       %.3f\n", distances[0]);
        printf("  |c_b - r_a1| : %.3f\n", distances[1]);
        printf("  |c_b - r_a2| : %.3f\n", distances[2]);
        printf("  |c_a - r_b1| : %.3f\n", distances[3]);
        printf("  |c_a - r_b2| : %.3f\n", distances[4]);
        print_vector("  c_a (corrected)", c_ab.begin(), std::next(c_ab.begin(), 3));
        print_vector("  c_b (corrected)", std::next(c_ab.begin(), 3), c_ab.end());
        std::cout << "\n";
    }

    return c_ab;
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
    bool is_parallel;

    normalize(p_a, l_a);
    normalize(p_b, l_b);
    cross_product(l_a, l_b, l_a_cross_l_b);

    // oblique / perpendicular elements
    if (rod::absolute(l_a_cross_l_b) > 1e-5)
    {
        is_parallel = false;
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
        is_parallel = true;
        vec3d(n) { c_a[n] = r_a[n] + 0.5 * p_a[n]; }
        vec3d(n) { c_b[n] = r_b[n] + 0.5 * p_b[n]; }
    }

    if (rod::dbg_print)
    {
        std::cout << "minimum distance between rod elements:\n";
        std::cout << "  rods parallel: " << is_parallel << "\n";
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

    // Apply corrections to infinite line assumption
    std::vector<float> c_ab{c_a[0], c_a[1], c_a[2], c_b[0], c_b[1], c_b[2]};
    c_ab = snap_to_nodes(c_ab, r_a, r_b, p_a, p_b);
    c_ab = compare_node_distances(c_ab, r_a, r_b, p_a, p_b);
    vec3d(n) { c_a[n] = c_ab.at(n); }
    vec3d(n) { c_b[n] = c_ab.at(n + 3); }
}

/*
Check if two rod elements, a and b, interact by calculating the shortest
distance between them and comparing to the sum of their radii. If this
passes, the interaction information is added to both elements' neighbour
lists.
*/
void set_element_neighbours(int rod_id_a, int rod_id_b, int elem_id_a,
    int elem_id_b, float p_a[3], float p_b[3], float r_a[3], float r_b[3],
    float radius_a, float radius_b, std::vector<InteractionData> &neighbours_a,
    std::vector<InteractionData> &neighbours_b)
{
    float mp_a[3] = {0};
    float mp_b[3] = {0};
    float mp_ab[3] = {0};
    float c_a[3] = {0};
    float c_b[3] = {0};
    float c_ab[3] = {0};

    rod::get_element_midpoint(p_a, r_a, mp_a);
    rod::get_element_midpoint(p_b, r_b, mp_b);
    vec3d(n){mp_ab[n] = mp_b[n] - mp_a[n];}

    if (rod::dbg_print)
    {
        std::cout << "cull distant elements:\n";
        std::printf("  |p_a| :         %.3f\n", rod::absolute(p_a));
        std::printf("  |p_b| :         %.3f\n", rod::absolute(p_b));
        std::printf("  radius sum :    %.3f\n", radius_a + radius_b);
        std::printf("  |midpoint ab| : %.3f\n", rod::absolute(mp_ab));
        std::printf("  cutoff     :    %.3f\n", std::max(rod::absolute(p_a), rod::absolute(p_b)) + radius_a + radius_b);
    }

    if (rod::absolute(mp_ab) < std::max(rod::absolute(p_a), rod::absolute(p_b)) + radius_a + radius_b)
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

// Distance between two rod elements, with a spatial perturbation applied to
// the centreline displacement
float intersection_distance(int perturb_dim, float perturb_delta,
    float contact_self[3], float contact_neighb[3], float radius_sum)
{
    float centreline_displacement[3] = { 0 };
    contact_neighb[perturb_dim] += perturb_delta;
    vec3d(n) { centreline_displacement[n] = contact_neighb[n] - contact_self[n]; }

    return rod::absolute(centreline_displacement) - radius_sum;
}

// No perturbation
float intersection_distance(float contact_self[3], float contact_neighb[3],
    float radius_sum)
{
    float centreline_displacement[3] = { 0 };
    vec3d(n) { centreline_displacement[n] = contact_neighb[n] - contact_self[n]; }

    return rod::absolute(centreline_displacement) - radius_sum;
}


/**
 * @brief Steric repulsive force on a rod element due to its collision
 * with a neighbouring element.
 */
std::vector<float> element_steric_force(float delta, float force_strength,
    float radius_sum, float contact_self[3], float contact_neighb[3])
{
    float distance[6] = { 0 };      // +x, -x, +y, -y, +z, -z
    float energy[6] = { 0 };        // +x, -x, +y, -y, +z, -z
    float force[3] = {0}; // x, y, z

    distance[0] = intersection_distance(0,  0.5 * delta, contact_self, contact_neighb, radius_sum);
    distance[1] = intersection_distance(0, -0.5 * delta, contact_self, contact_neighb, radius_sum);
    distance[2] = intersection_distance(1,  0.5 * delta, contact_self, contact_neighb, radius_sum);
    distance[3] = intersection_distance(1, -0.5 * delta, contact_self, contact_neighb, radius_sum);
    distance[4] = intersection_distance(2,  0.5 * delta, contact_self, contact_neighb, radius_sum);
    distance[5] = intersection_distance(2, -0.5 * delta, contact_self, contact_neighb, radius_sum);

    energy[0] = steric_energy_linear(force_strength, distance[0]);
    energy[1] = steric_energy_linear(force_strength, distance[1]);
    energy[2] = steric_energy_linear(force_strength, distance[2]);
    energy[3] = steric_energy_linear(force_strength, distance[3]);
    energy[4] = steric_energy_linear(force_strength, distance[4]);
    energy[5] = steric_energy_linear(force_strength, distance[5]);

    // ! gradient implicitly negative due to use of intersection distance, rather
    // ! than centreline distance. See notes 30/6/22.
    force[0] = (energy[0] - energy[1]) / delta;
    force[1] = (energy[2] - energy[3]) / delta;
    force[2] = (energy[4] - energy[5]) / delta;

    float c_ab[3] = {0};
    vec3d(n) { c_ab[n] = contact_neighb[n] - contact_self[n]; }

    if (rod::dbg_print)
    {
        std::cout << "element steric force:\n";
        std::cout << "  delta : " << delta << "\n";
        std::cout << "  force_strength : " << force_strength << "\n";
        std::cout << "  radius_sum : " << radius_sum << "\n";
        print_array("  contact_self", contact_self, 3);
        print_array("  contact_neighb", contact_neighb, 3);
        std::cout << "  centreline distance : " << rod::absolute(c_ab) << "\n";
        print_array("  perturbed distance x", distance, 0, 1);
        print_array("  perturbed distance y", distance, 2, 3);
        print_array("  perturbed distance z", distance, 4, 5);
        print_array("  element energy x", energy, 0, 1);
        print_array("  element energy y", energy, 2, 3);
        print_array("  element energy z", energy, 4, 5);
        print_array("  element force", force, 3);
        std::cout << "\n";
    }

    // sanity check #1
    if (rod::absolute(c_ab) > radius_sum)
    {
        std::cout << "ERROR!\n";
        throw std::runtime_error("Centreline distance cannot be "
            "greater than radius sum when calculating steric energy.");
    }

    // sanity check #2
    float force_norm[3] = {0};
    float c_ab_norm[3] = {0};

    normalize(force, force_norm);
    normalize(c_ab, c_ab_norm);
    float dot = rod::dot_product_3x1(force_norm, c_ab_norm);

    if (dot >= 0)
    {
        std::cout << "ERROR!\n";
        std::string msg = "Force on rod element A (self) should point away from rod element B (neighbour):\n"
            "  force :      (" + std::to_string(force[0]) + ", " + std::to_string(force[1]) + ", " + std::to_string(force[2]) + ")\n"
            "  c_ab :       (" + std::to_string(c_ab[0]) + ", " + std::to_string(c_ab[1]) + ", " + std::to_string(c_ab[2]) + ")\n"
            "  force_norm : (" + std::to_string(force_norm[0]) + ", " + std::to_string(force_norm[1]) + ", " + std::to_string(force_norm[2]) + ")\n"
            "  c_ab :       (" + std::to_string(c_ab_norm[0]) + ", " + std::to_string(c_ab_norm[1]) + ", " + std::to_string(c_ab_norm[2]) + ")\n"
            "  dot :         " + std::to_string(dot) + "\n";
        throw std::runtime_error(msg);
    }
    return std::vector<float> {force[0], force[1], force[2]};
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
    std::vector<float> force(6, 0);  // x0, y0, z0, x1, y1, z1

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
