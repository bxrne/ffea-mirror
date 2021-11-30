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
 * With this method, we get 'possible' interactions by applying a crude distance calculation + if statement to
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
    vec3d(n){d2[n] = p_i_mid[n] - p_j[n] + r_1j;}

    if(rod::absolute(d1) < cutoff || rod::absolute(d2) < cutoff){
        return true;
    }
    return false;
}

/* Construct steric interaction neighbour list for a rod
 * 
 *   - i : id of the first rod
 *   - num_rods : the total number of rods in the simulation
 *   - rod_array : 1-D array containing pointers to all rod objects
*/
void create_neighbour_list(int i, int num_rods, rod::Rod **rod_array){
    float r_i[3] = {0, 0, 0};
    float r_j[3] = {0, 0, 0};
    float p_i[3] = {0, 0, 0};
    float p_j[3] = {0, 0, 0};
    bool within_cutoff = false;
    float l_i[3] = {0, 0, 0};
    float l_j[3] = {0, 0, 0};
    float l_i_cross_l_j[3] = {0, 0, 0};
    float c_i[3] = {0, 0, 0};
    float c_j[3] = {0, 0, 0};
    std::vector<float> rod_neighbour_list[rod_array[i]->num_elements];

    // other rod: j
    for (int j=i+1; j<num_rods; j++){
        // elements of rod i: mt
        for (int m=0; m<rod_array[i]->num_elements; m++){
            // neighbour list of unknown size for element m
            std::vector<float> element_neighbour_list;

            // elements of rod j: n
            for (int n=0; n<rod_array[j]->num_elements; n++){

                // current_r is 1D, so need to index with an integer shift
                vec3d(dim){r_i[dim] = rod_array[i]->current_r[(m*3)+dim];}
                vec3d(dim){r_j[dim] = rod_array[j]->current_r[(n*3)+dim];}

                rod_array[i]->get_p(m, p_i, false);
                rod_array[j]->get_p(n, p_j, false);

                within_cutoff = rod::elements_within_cutoff(r_i, p_i, r_j, p_j, rod::absolute(p_i));

                if(within_cutoff == true){
                    rod::normalize(p_i, l_i);
                    rod::normalize(p_j, l_j);
                    rod::cross_product(l_i, l_j, l_i_cross_l_j);

                    rod::get_point_on_connecting_line(p_i, p_j, l_i_cross_l_j, r_i, r_j, c_i);
                    rod::get_point_on_connecting_line(p_j, p_i, l_i_cross_l_j, r_j, r_i, c_j);

                    // Add c_i, c_j to the end of the rod neighbour list, adjacent
                    // to eachother {...c_ix, c_iy, c_iz, c_jx, c_jy, c_jz...}
                    //
                    // For rod element m, experiencing T collisions, there will 
                    // be 6*T floats in the neighbour list.
                    for(int dim=0; dim<3; dim++){
                        element_neighbour_list.push_back(c_i[dim]);
                    }
                    for(int dim=0; dim<3; dim++){
                        element_neighbour_list.push_back(c_j[dim]);
                    }
                 
                    rod::print_array("c_i", c_i, 3);
                    rod::print_array("c_j", c_j, 3);
                }
                
            }
            // consolidate dynamic neighbour list of element m into a static
            // 2D array that has rows equal to the number of elements, M, in
            // rod i. There will be 6*T*M floats in this array.
            rod_neighbour_list[m] = element_neighbour_list;
            delete element_neighbour_list;
        }
    }
}


//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND IT WILL BECOME
//   (_(_/-(_/   A BUG-FREE ZONE. WOOF WOOF!

}//end namespace
