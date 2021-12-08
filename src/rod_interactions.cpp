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

namespace rod {

/* See if two rod elements, p_i and p_j, are within some interacting distance, defined by a cutoff radius.
 * Measure this cutoff from the midpoint of p_i to both nodes of p_j. Return true if either node is within range.
 * 
 * With this method, we will get 'possible' interactions by applying a crude distance calculation + if statement to
 * every rod element in the system, to narrow down the number of elements that will be used in a more intensive force 
 * calculation later on.
*/
// TODO: Use this in World.cpp by looping over rod pairs, then the nodes within those rods
bool elements_within_cutoff(float r_1i[3], float p_i[3], float r_1j[3], float p_j[3], float cutoff){
    float p_i_mid[3] = {0, 0, 0};
    float d1[3] = {0, 0, 0};
    float d2[3] = {0, 0, 0};

    rod::get_p_midpoint(p_i, r_1i, p_i_mid);
    vec3d(n){d1[n] = p_i_mid[n] - r_1j[n];}
    vec3d(n){d2[n] = p_i_mid[n] - p_j[n] + r_1j[n];}

    if(rod::absolute(d1) < cutoff || rod::absolute(d2) < cutoff){
        return true;
    }
    else{
        return false;
    }
}


/* Construct steric interaction neighbour lists for two rods, a and b

 *   Arguments:
 *   - *rod_a, *rod_b : pointers to rod objects
 *   Affects:
 *   - rod_a,b->steric_interaction_coordinates[element_index]
 *
 *   Loops over every element of the two rods (O(N**2) loop) and checks if they might interact
 *   using a crude distance calculation. If this passes, then a more precise distance calculation
 *   is performed and the interaction coordinates are saved.
 * 
 *   Each element of the rods will be assigned a std::vector<float> containing the coordinate
 *   pairs describing steric interactions with other rod elements. These are stored in groups
 *   of 6, such that indices 0, 1, 2 are the coordinate on the current rod, and indices
 *   3, 4, 5 are the coordinate on the other rod, e.g. {ax, ay, az, bx, by, bz ...}. Hence for
 *   one element experiencing 5 interactions, this will create a vector containing 30 floats.
 */
void create_neighbour_list(rod::Rod *rod_a, rod::Rod *rod_b){
    float r_a[3] = {0, 0, 0};
    float r_b[3] = {0, 0, 0};
    float p_a[3] = {0, 0, 0};
    float p_b[3] = {0, 0, 0};
    bool within_cutoff = false;
    float l_a[3] = {0, 0, 0};
    float l_b[3] = {0, 0, 0};
    float l_a_cross_l_b[3] = {0, 0, 0};
    float c_a[3] = {0, 0, 0};
    float c_b[3] = {0, 0, 0};

    // Rather confusingly, num_elements actually refers to the number of NODES in the rod (8/12/21)
    for (int i=0; i<rod_a->num_elements-1; i++){
        for (int j=0; j<rod_b->num_elements-1; j++){

            // shift by i*3 and j*3 due to current_r being 1D
            vec3d(n){r_a[n] = rod_a->current_r[(i*3)+n];}
            vec3d(n){r_b[n] = rod_b->current_r[(j*3)+n];}
            rod_a->get_p(i, p_a, false);
            rod_b->get_p(j, p_b, false);

            within_cutoff = rod::elements_within_cutoff(r_a, p_a, r_b, p_b, rod::absolute(p_a));

            if(dbg_print){
                std::cout << "rod::create_neighbour_list()" << std::endl;
                std::cout << "  (a) rod " << rod_a->rod_no << ", element " << i << std::endl;
                std::cout << "  (b) rod " << rod_b->rod_no << ", element " << j << std::endl;
            }

            if(within_cutoff == true){
                if(dbg_print){std::cout << "  interaction!" << std::endl;}

                // Shortest distance between two rod elements
                rod::normalize(p_a, l_a);
                rod::normalize(p_b, l_b);
                rod::cross_product(l_a, l_b, l_a_cross_l_b);
                rod::get_point_on_connecting_line(p_a, p_b, l_a_cross_l_b, r_a, r_b, c_a);
                rod::get_point_on_connecting_line(p_b, p_a, l_a_cross_l_b, r_b, r_a, c_b);

// =============================================================================
                std::cout << "\nPOINTER TEST!! Assigning to rod " << rod_a->rod_no << std::endl;
                std::cout << "Num rod nodes: " << rod_a->num_elements << std::endl;
                std::cout << "Num rod elements: " << rod_a->num_elements-1 << std::endl;
                std::cout << "Num array elements in steric_interaction_coordinates: " << (rod_a->length/3)-1 << std::endl;

                for(int ele=0; ele<3; ele++){
                    vec3d(n){rod_a->steric_interaction_coordinates[ele].push_back(c_a[n]);}
                    vec3d(n){std::cout << "  element " << ele << ", c[" << n << "]: " << c_a[n] << std::endl;}
                }

                rod::print_vector("vector " + std::to_string(0), rod_a->steric_interaction_coordinates[0]);
                rod::print_vector("vector " + std::to_string(1), rod_a->steric_interaction_coordinates[1]);
                rod::print_vector("vector " + std::to_string(2), rod_a->steric_interaction_coordinates[2]);
                std::cout << "POINTER TEST OVER" << std::endl;
// =============================================================================                

                // TODO: seg fault occurs here
                // Update both rods with the interaction coordinate pair
                vec3d(n){rod_a->steric_interaction_coordinates[i].push_back(c_a[n]);}
                vec3d(n){rod_a->steric_interaction_coordinates[i].push_back(c_b[n]);}
                vec3d(n){rod_b->steric_interaction_coordinates[j].push_back(c_b[n]);}
                vec3d(n){rod_b->steric_interaction_coordinates[j].push_back(c_a[n]);}

                if(dbg_print){
                    rod::print_array("  c_a", c_a, 3);
                    rod::print_array("  c_b", c_b, 3);
                }
            }    
        }
    }
}



//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND IT WILL BECOME
//   (_(_/-(_/   A BUG-FREE ZONE. WOOF WOOF!

}//end namespace
