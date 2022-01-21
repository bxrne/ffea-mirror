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
 * 
 * NOTE: Not in use.
*/
bool elements_within_cutoff(float p_i[3], float p_j[3], float r_i[3], float r_j[3], float cutoff){
    float r_i_mid[3] = {0, 0, 0};
    float d1[3] = {0, 0, 0};
    float d2[3] = {0, 0, 0};

    rod::get_element_midpoint(p_i, r_i, r_i_mid);
    vec3d(n){d1[n] = r_i_mid[n] - r_j[n];}
    vec3d(n){d2[n] = r_i_mid[n] - p_j[n] + r_j[n];}

    if(rod::dbg_print){
        std::cout << "elements_within_cutoff()" << std::endl;
        std::cout << "  cutoff: " << cutoff << std::endl;
        std::cout << "  |d1|: " << rod::absolute(d1) << std::endl;
        std::cout << "  |d2|: " << rod::absolute(d2) << std::endl;
    }

    if(rod::absolute(d1) < cutoff || rod::absolute(d2) < cutoff){
        if(rod::dbg_print){std::cout << "  Within cutoff" << std::endl;}
        return true;
    }
    else{
        if(rod::dbg_print){std::cout << "  Outisde cutoff" << std::endl;}
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
 *   by calculating the shortest distance between two line elements, then comparing to the sum
 *   of the radii of both rods. If this passes, the interaction coordinates are appended to
 *   a vector.
 * 
 *   Each element of the rods will be assigned a std::vector<float> containing the coordinate
 *   pairs describing steric interactions with other rod elements. These are stored in groups
 *   of 6, such that indices 0, 1, 2 are the coordinate on the current rod, and indices
 *   3, 4, 5 are the coordinate on the other rod, with an additional element 6 to store the radius
 *   of the other rod, e.g. {ax, ay, az, bx, by, bz, R_b ...}. Hence for one element experiencing 5
 *   interactions, this will create a vector containing 35 floats.
 */
void create_neighbour_list(rod::Rod *rod_a, rod::Rod *rod_b){
    float r_a[3] = {0, 0, 0};
    float r_b[3] = {0, 0, 0};
    float p_a[3] = {0, 0, 0};
    float p_b[3] = {0, 0, 0};
    float c_a[3] = {0, 0, 0};
    float c_b[3] = {0, 0, 0};
    float c_ab[3] = {0, 0, 0};

    // Rather confusingly, num_elements actually refers to the number of NODES in the rod (8/12/21)
    for (int i=0; i<rod_a->num_elements-1; i++){
        for (int j=0; j<rod_b->num_elements-1; j++){

            // shift by i*3 and j*3 due to current_r being 1D
            vec3d(n){r_a[n] = rod_a->current_r[(i*3)+n];}
            vec3d(n){r_b[n] = rod_b->current_r[(j*3)+n];}
            rod_a->get_p(i, p_a, false);
            rod_b->get_p(j, p_b, false);

            // Shortest distance between two rod elements
            rod::get_interaction_vector(p_a, p_b, r_a, r_b, c_a, c_b);
            vec3d(n){c_ab[n] = c_b[n] - c_a[n];}

            if(rod::dbg_print){
                std::cout << "rod::create_neighbour_list()" << std::endl;
                std::cout << "  |c_ab|: " << rod::absolute(c_ab) << std::endl;
                std::cout << "  radii sum: " << rod_a->steric_radius + rod_b->steric_radius << std::endl;
            }

            if(rod::absolute(c_ab) < (rod_a->steric_radius + rod_b->steric_radius)){
                if(rod::dbg_print){std::cout << "  interaction - rod " << rod_a->rod_no << ", elem " << i << " | rod " << rod_b->rod_no << ", elem " << j << std::endl;}

                // Increase vector capacity by 6 (optional, might help with correct assignment to vector)
                rod_a->steric_interaction_coordinates.at(i).reserve(rod_a->steric_interaction_coordinates.at(i).size() + 7);
                rod_b->steric_interaction_coordinates.at(j).reserve(rod_b->steric_interaction_coordinates.at(j).size() + 7);

                // Update both rods with the interaction coordinate pair and the radius of the other rod
                vec3d(n){rod_a->steric_interaction_coordinates.at(i).push_back(c_a[n]);}
                vec3d(n){rod_a->steric_interaction_coordinates.at(i).push_back(c_b[n]);}
                rod_a->steric_interaction_coordinates.at(i).push_back(rod_b->steric_radius);

                vec3d(n){rod_b->steric_interaction_coordinates.at(j).push_back(c_b[n]);}
                vec3d(n){rod_b->steric_interaction_coordinates.at(j).push_back(c_a[n]);}
                rod_b->steric_interaction_coordinates.at(j).push_back(rod_a->steric_radius);
            }
            else{
                if(rod::dbg_print){std::cout << "  no interaction detected - rod " << rod_a->rod_no << ", elem " << i << " | rod " << rod_b->rod_no << ", elem " << j << std::endl;}
            }  
        }
    }
}



//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND IT WILL BECOME
//   (_(_/-(_/   A BUG-FREE ZONE. WOOF WOOF!

}//end namespace
