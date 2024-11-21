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

#ifndef SPARSESUBSTITUTIONSOLVER_H_INCLUDED
#define SPARSESUBSTITUTIONSOLVER_H_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseSubstitutionSolver.h"

#define INDEX(I, J) ((I) * num_rows + (J))

class SparseSubstitutionSolver : public Solver {
public:
    /** Builds the lower triangular Cholesky decomposed mass matrix */
    void init(std::vector<mesh_node >&node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list) override;

    /**
     * Solves the equation Ax = b for the unknown vector x, for
     * 3 right-hand-sides (b vectors) at once, using forward and backward substitution. This function
     * considers the case where A is sparse symmetric positive-definite and therefore the solver must
     * have constructed a variable band matrix containing the result of a Cholesky decomposition on
     * the original matrix A.
     *
     * This implementation is memory efficient as it modifies the vector x in place (no temporary vectors are needed).
     */
    void solve(std::vector<arr3> &x) override;

    void apply_matrix(const std::vector<scalar> &in, std::vector<scalar> &result) override;

private:

    /** Number of rows in this cholesky variable band matrix */
    int num_rows = 0;

    //@{
    /**
     * Array of size num_rows. Each element contains the index of the last non-zero
     * entry in the corresponding column. This is the end of the "band" of entries stored
     * in that column, that starts from the diagonal (travelling downwards or upwards for
     * the L and U matrices respectively).
     */
    std::vector<int> L_key, U_key;
    //@}

    /** Stores the inverse of the diagonal elements (need for LU solving) */
    std::vector<scalar> inverse_diag;

    //@{
    /**
     * Stores all the data from the bands of the band matrix (note that this may contain
     * some zeroes, but if it is a proper sparse band matrix there should be very few of these)
     * NOT INCLUDING THE DIAGONAL: only the inverse of the diagonal.
     *
     * The off-diagonal data for the lower and upper triangles are stored in L and U respectively.
     */
    std::vector<scalar> L, U;
    //@}
};

#endif
