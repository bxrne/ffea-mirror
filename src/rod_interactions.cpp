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

/* Construct steric interaction neighbour list for a rod
 * 
 *   - i : id of the first rod
 *   - num_rods : the total number of rods in the simulation
 *   - rod_array : 1-D array containing pointers to all rod objects
*/
void get_neighbour_list(int i, int num_rods, rod::Rod **rod_array){
    float r_i[3] = {0, 0, 0};
    float r_j[3] = {0, 0, 0};
    float p_i[3] = {0, 0, 0};
    float p_j[3] = {0, 0, 0};
    float l_i_cross_l_j[3] = {0, 0, 0};
    float c_i[3] = {0, 0, 0};
    float c_j[3] = {0, 0, 0};

    // neighbour_list = some kind of pointer?

    // other rod: j
    for (int j=i+1; j<num_rods; j++){
        // elements of rod i: m
        for (int m=0, m<rod_array[i]->num_elements; m++){
            // elements of rod j: n
            for (int n=0, n<rod_array[j]->num_elements; n++){
                rod::vec3d(ind){r_i[ind] = rod_array[i]->current_r[m][ind];}
                rod::vec3d(ind){r_j[ind] = rod_array[i]->current_r[m][ind];}

                rod::get_p_i(r_i, rod_array[i]->current_r[m+1], p_i);
                rod::get_p_i(r_j, rod_array[j]->current_r[n+1], p_j);

                rod::cross_product(p_i/rod::absolute(p_i), rod::absolute(p_j), l_i_cross_l_j);
                rod::get_point_on_connecting_line(p_i, p_j, l_i_cross_l_j, r_i, r_j, c_i);
                rod::get_point_on_connecting_line(p_j, p_i, l_i_cross_l_j, r_j, r_i, c_j);

                // add {c_i, c_j} to the rod neighbour list

                rod::print_array("c_i", c_i, 3);
                rod::print_array("c_j", c_j, 3);
            }
        } 
    }
}


//    __      _
//  o'')}____//  I AM DEBUG DOG. PUT ME IN YOUR
//   `_/      )  SOURCE CODE AND IT WILL BECOME
//   (_(_/-(_/   A BUG-FREE ZONE. WOOF WOOF!

}//end namespace