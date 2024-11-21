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

#include "LJ_matrix.h"

void SSINT_matrix::init(const string &ssint_params_fname, const string &ssint_type, int calc_ssint, scalar ssint_cutoff) {
    // In that case we do not need an input parameter file,
    //      but just to initialise a number of values
    if (calc_ssint == 0) {
        init_steric();
    } else {
        init_ssint(ssint_params_fname, ssint_type, ssint_cutoff);
    }
}

void SSINT_matrix::init_steric() {

    num_ssint_face_types = 1;
    // Allocate the memory for that LJ pair
    try {
        params = std::vector<map<string, scalar>>(num_ssint_face_types * num_ssint_face_types);
    } catch (std::bad_alloc &) {
        throw FFEAException("Unable to allocate memory for LJ matrix.");
    }
    params[LJI(0, 0)]["Emin"] = 0;
    params[LJI(0, 0)]["Rmin"] = 0;
}


void SSINT_matrix::init_ssint(const string &ssint_params_fname, const string &ssint_type, scalar ssint_cutoff) {
    const int MAX_NUM_VARS = 3;
    string line;
    vector<string> head(MAX_NUM_VARS);
    vector<string> bufvec;
    vector<string> bufvec2;
    ifstream in (ssint_params_fname);
    bool is_rod = false;

    // What are the param possibilities (keys)?
    head[0] = "Emin";
    head[1] = "Rmin";
    head[2] = "k0";
    if (!in.is_open())
        throw FFEAFileException(ssint_params_fname);

    printf("\tReading in VDW forcefield parameters file: %s\n", ssint_params_fname.c_str());

    // first line should be the file type "ffea ssint forcefield params file" (.lj or .rodlj)
    getline(in, line);
    if (line != "ffea vdw forcefield params file" && line != "ffea ssint forcefield params file" && line != "ffea rod vdw forcefield params file")
    {
        in.close();
        throw FFEAException("This is not a 'ffea vdw forcefield params file' (read '%s').", line.c_str());
    }
    if (line.find("rod") != std::string::npos)
        is_rod = true;

    // read in the number of interaction face types in the file
    getline(in, line);
    boost::split(bufvec, line, boost::is_any_of(" "));
    cout << bufvec.at(0) << " " << atoi(bufvec.at(1).c_str()) << endl;
    num_ssint_face_types = atoi(bufvec.at(1).c_str());

    printf("\t\tNumber of VDW interaction types = %d\n", num_ssint_face_types);

    // Allocate the memory for all these LJ pairs
    try {
        params = std::vector<map<string, scalar>>(num_ssint_face_types * num_ssint_face_types);
    } catch (std::bad_alloc &) {
        in.close();
        throw FFEAException("Unable to allocate memory for LJ matrix.");
    }

    // Fill the matrix (with an as yet unknown number of parameters)
    for (int i = 0; i < num_ssint_face_types; i++)
    {
        // Use boost to get something we can deal with
	    getline(in, line);
        boost::split(bufvec, line, boost::is_any_of(")"));

        for (int j = 0; j < num_ssint_face_types; j++)
        {
            // Strip it to a comma delimited set
            string aset = bufvec.at(j);
            boost::trim(aset);
            aset = boost::erase_first_copy(aset, "(");
            boost::split(bufvec2, aset, boost::is_any_of(","));

            // First parameter is an energy one, second is a distance one, third is a curvature one. This may need further generalisation
            if (bufvec2.size() > 3)
                throw FFEAException("For SSINT_params set %d %d, currently unable to have more than 3 params. You specified %zu.", i, j, bufvec2.size());
            if (bufvec2.size() < 2)
                throw FFEAException("For SSINT_params set %d %d, currently must have at least 2 params (Emin, Rmin). You specified %zu.", i, j, bufvec2.size());

            for(int count = 0; count < bufvec2.size(); count++)
                params[LJI(i, j)][head[count]] = atof(bufvec2.at(count).c_str());

            if (params[LJI(i, j)]["Emin"] <= 0)
                throw FFEAException("Required: 'Emin' must be greater than 0 if you wish to use ssint (calc_ssint is 1).");
            if (params[LJI(i, j)]["Rmin"] <= 0)
                throw FFEAException("Required: 'Rmin' must be greater than 0 if you wish to use ssint (calc_ssint is 1).");

            // In gensoft, set the Rmin value to the cutoff distance for now, so avoid turning points
            if (ssint_type == "gensoft")
                params[LJI(i, j)]["Rmin"] = ssint_cutoff * mesoDimensions::length;

            // Defaulting the value if not found
            if(params[LJI(i, j)].count(head[2]) == 0)
                params[LJI(i, j)][head[2]] = 0.0;

            // ! - in the future, rods and blobs will use the same vdw file
            if (is_rod)
            {
                params[LJI(i, j)][head[0]] = params[LJI(i, j)][head[0]] / mesoDimensions::Energy ;
                params[LJI(i, j)][head[1]] = params[LJI(i, j)][head[1]] / mesoDimensions::length;
            }
            else
            {
                params[LJI(i, j)][head[0]] = params[LJI(i, j)][head[0]] * mesoDimensions::area * mesoDimensions::area / mesoDimensions::Energy ;
                params[LJI(i, j)][head[1]] = params[LJI(i, j)][head[1]] / mesoDimensions::length;
                params[LJI(i, j)][head[2]] = params[LJI(i, j)][head[2]] * mesoDimensions::area / mesoDimensions::Energy;
            }
	    }
    }

    in.close();

    printf("\t\tRead %d VDW forcefield parameter entries from %s\n", num_ssint_face_types * num_ssint_face_types, ssint_params_fname.c_str());
}

map<string, scalar> SSINT_matrix::get_SSINT_params(int type1, int type2) {
    if (type1 < 0 || type1 > num_ssint_face_types - 1) {
        printf("Frog1 %d %d\n", type1, num_ssint_face_types - 1);
    }
    if (type2 < 0 || type2 > num_ssint_face_types - 1) {
        printf("Frog2 %d %d\n", type2, num_ssint_face_types - 1);
    }

    return params[LJI(type1, type2)];
}

int SSINT_matrix::get_num_types() {
    return num_ssint_face_types;
}
