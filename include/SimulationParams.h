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

#ifndef PARAMS_H_INCLUDED
#define PARAMS_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>
#include <map>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "FFEA_return_codes.h"
#include "FFEA_input_reader.h"
#include "FFEA_user_info.h"
#include "mat_vec_types.h"
#include "dimensions.h"

#define WALL_TYPE_PBC 0
#define WALL_TYPE_HARD 1
#define WALL_TYPE_STOP 2

#define UNSET 0
#define SET 1
#define DEFAULT 2

/**
 * Simulation parameters
 */

using namespace std;
namespace b_fs = boost::filesystem;

class SimulationParams {
public:
    scalar dt; ///< time step
    long long num_steps; ///< Number of time steps to run simulation for
    int check; ///< Every how many steps should the program 'check' the system i.e calculate energies, print snapshots etc.
    int num_blobs; ///< Number of blobs in the system
    int *num_conformations; ///< Number of conformations for each blob
    int *num_states; ///< Number of states for each blob
    int state_array_size;
    int conformation_array_size;
    int rng_seed; ///< Seed for random number generator

    scalar kT; ///< boltzmann's const times temperature

    int max_iterations_cg; ///< Max number of iterations when using conjugate gradient solver
    scalar epsilon2; ///< The tolerance threshold for CG solver (squared)

    int es_update; ///< Every how many steps should the electrostatic potential be recalculated
    int es_N_x; ///< X dimension of the 3D lookup grid (in number of cells)
    int es_N_y; ///< Y dimension of the 3D lookup grid (in number of cells)
    int es_N_z; ///< Z dimension of the 3D lookup grid (in number of cells)
    int move_into_box;  ///< If box is set, do we move the world to it's center at the start of a simulation? Default is yes!
    int restrict_motion[3]; ///< [x,y,z] array defining whether motion in the given direction should be nullified
    int num_dimensions;  /// Number of active dimensions after restricted_motion is applied
    scalar es_h; ///< Dimension of each cell in the lookup grid (in multiples of inverse kappa)

    scalar kappa; ///< Inverse Debye Screening length

    scalar epsilon_0; ///< Permittivity of free space
    scalar dielec_ext; ///< Exterior dielectric constant

    int restart; ///< Whether or not to restart the simulation from the last available time step

    int calc_vdw; ///< Whether or not to simulate van der waals interactions between surfaces
    int inc_self_vdw; ///< Whether or not to include van der Waals interactions derived from faces in the same blob.
    string vdw_type;  ///<Possible values: "lennard-jones" (default) or "steric".
    int calc_es; ///< Whether or not to simulate electrostatic interactions between proteins
    int calc_noise; ///< Whether or noise to simulate thermal noise for the system. Kind of the entire point of this simulation technique
    int calc_stokes; ///< Whether or not to include local action of the external fluid
    int calc_kinetics;  ///< Whether or not to calculate kinetic switching between different equilibrium states and binding sites  
    int calc_preComp; ///< Whether or not use preComputed potentials and forces   
    int calc_springs; ///< Whether or not to include the springs interactions defined in the springs block 
    int calc_ctforces; ///< Whether or not to include constant forces onto nodes defined in the ctforces block
    int force_pbc; ///< Whether or not to apply pbc to surface insteractions
    int kinetics_update; ///< How often to check for a state change. If rates are ~ >> dt then this can clearly be quite high   
    int wall_x_1;
    int wall_x_2;
    int wall_y_1;
    int wall_y_2;
    int wall_z_1;
    int wall_z_2;

    int sticky_wall_xz;

    scalar stokes_visc;

    scalar vdw_steric_factor; ///< Proportionality factor to the Steric repulsion.
    scalar vdw_cutoff; ///< Cutoff distance for the VdW interactions. 
    geoscalar vdw_steric_dr; ///< used to calculate the numerical derivative.

    string FFEA_script_filename;
    b_fs::path FFEA_script_path, FFEA_script_basename;
    string trajectory_out_fname;
    string kinetics_out_fname;
    string measurement_out_fname;
    string detailed_meas_out_fname;
    string vdw_in_fname;
    string bsite_in_fname;
    string icheckpoint_fname;  ///< Input Checkpoint file name
    string ocheckpoint_fname;  ///< Output Checkpoint file name
    string ctforces_fname; ///< Input file containing constant forces onto a list of nodes.
    string springs_fname; ///< Input file containing the springs details.
    string trajectory_beads_fname; ///< Output optional file.

    SimulationParams();

    ~SimulationParams();

    int validate(int sim_mode);

    /**
     * @brief DEPRECATED! Superseeded by: extract_params + assign.
     *
     * Expects a string of the form "lvalue = rvalue" where lvalue must be a recognised parameter name.
     */
    int parse_param_assignment(char *str);

    /**
     * Expects a string of the form "<system>\n<>\n<>.....\n<system>" where parameter data can be extracted from.
     */
    int extract_params(vector<string> script_vector);

    /** Expects a parameter label and value, which will be assigned if valid and rejected if not */
    int assign(string lvalue, string rvalue);

    /** Returns maximum number of states on any blob */
    int get_max_num_states();

    /** These set parameters are not private because the World needs them!! */
    int kinetics_out_fname_set;
    int trajbeads_fname_set; 

    /** Writes all params to params_out_fname for user's info */
    void write_to_file(FILE *fout);

private:
    int trajectory_out_fname_set;
    int measurement_out_fname_set;
    int icheckpoint_fname_set;
    int ocheckpoint_fname_set;
    int vdw_in_fname_set;
    int bsite_in_fname_set;

/** Check if the file oFile exists, and if so
  *     rename it to "__"+oFile+"__bckp.N",
  *     where N is an integer so that the resulting file is new.
  */
    int checkFileName(string oFile);

    string RemoveFileExtension(const string& FileName);
};
#endif
