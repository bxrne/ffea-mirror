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

#ifndef ROD_INTERACTIONS
#define ROD_INTERACTIONS

#define _USE_MATH_DEFINES ///<  This has to come before including cmath

#include "rod_math_v9.h"
#include "rod_structure.h"

namespace rod {

void rod_distance_correction(float c_a[3], float c_b[3], float r_a[3], float r_b[3], float p_a[3], float p_b[3], OUT float c_a_out[3], float c_b_out[3]);
void get_shortest_distance_to_rod(float p_a[3], float p_b[3], float r_a[3], float r_b[3], OUT float c_a[3], float c_b[3]);
void create_neighbour_list(rod::Rod *rod_a, rod::Rod *rod_b);
// float get_spherical_volume_intersection(float separation, float radius_a, float radius_b);  // Not in use

}
#endif
