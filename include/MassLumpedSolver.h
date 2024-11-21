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

#ifndef MASSLUMPEDSOLVER_HPP_INCLUDED
#define MASSLUMPEDSOLVER_HPP_INCLUDED


#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseMatrixTypes.h"

class MassLumpedSolver : public Solver {
public:

    /** Builds the diagonal mass matrix and gets reciprocal of each value */
    void init(std::vector<mesh_node> &node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list) override;

    /** Applies inverse mass matrix (since diagonal: Mx = f => x = f_i/M_i */
    void solve(std::vector<arr3> &x) override;

    /** Applies the mass matrix to the given vector, 'in', putting the result in 'result'*/
    void apply_matrix(const std::vector<scalar> &in, std::vector<scalar> &result) override;

private:

    /** Number of rows in original matrix */
    int num_rows = 0;

    /** Mass matrix diagonal inversed */
    std::vector<scalar> inv_M;
};

#endif
