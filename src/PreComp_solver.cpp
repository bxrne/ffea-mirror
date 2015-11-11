/** 
 * \defgroup FMM Force Match
 * @{
 */

#include "PreComp_solver.h"

using namespace std;

PreComp_solver::PreComp_solver() { 
  nint = 0;
  msgc = 0;
  n_beads = 0;
} 

/** @brief destructor: deallocates pointers. */
PreComp_solver::~PreComp_solver() {
  delete U;
  delete F;
  delete b_types;
  delete b_elems; 
}


int PreComp_solver::msg(int whatever){
  cout << "--- PreComp " << msgc << ": " << whatever << endl;
  msgc += 1;
  return 0; 
} 

int PreComp_solver::msg(string whatever){
  cout << "--- PreComp " << msgc << ": " << whatever << endl;
  msgc += 1;
  return 0; 
} 


/** 
 * @brief Read input precomputed tables
 * @param[in] vector<string> types: types of beads present.
 * @param[in] int inputData: 1 means read .force and .pot files,
 *                 while 2 means read .pot and calculate the forces
 * @param[in] int approach: if .vdw -> use .vdw files and solve faces;
 *                          if .solid -> use .solid files and solve volumes.
 * @details read <type_i>-<type_j>.force and <type_i>-<type_j>.pot files
 *         for all the possible pairs, setting up the potential functions.
 *         All the .pot and .force files need to have the same x_range
 *         same Dx, and same number of points. It is allowed, however, that 
 *         they have a different number of lines at the beginning starting with "#". 
 */        
int PreComp_solver::init(PreComp_params *pc_params, SimulationParams *params, Blob **blob_array) {

   /* Firstly, we get the number of lines, x_range, 
    * and Dx from the first pair type potential file. 
    * We check that the x_range and Dx 
    * are the same for all the files. Secondly, 
    * we allocate the F and U arrays, and then
    * we store the y_values, re-reading the files.
    * The functions get_U and get_F will then be ready. 
    * Thirdly and fourthly, allocate and fill the 
    * elements and relative bead position arrays.
    * And finally, delete beads stuff from the blobs.
    */  
   
   if (pc_params->types.size() == 0) return 0;

   stringstream ssfile; 
   ifstream fin;
   string line; 
   vector<string> vec_line; 
   scalar x_0, x;
   scalar d2, d2_0;
   
   /*--------- FIRSTLY ----------*/ 
   /* get the number of interactions */
   for (int i=1; i<pc_params->types.size() + 1; i++){
     nint += i;
   } 
   msg(nint);
   // and the number of bead types: 
   ntypes = pc_params->types.size();

   /* read the first file, and get the number of lines
    * (n_values), Dx, and x_range */
   // the first file:
   ssfile << pc_params->folder << "/" << pc_params->types[0] << "-" << pc_params->types[0] << ".pot";
   fin.open(ssfile.str());
   // get the first line that does not start with "#"
   getline(fin, line);
   while (line.find("#", 0, 1) == 0) {
     getline(fin, line); 
   } 
   // get x_0, i. e., parse the line:
   boost::split( vec_line, line, boost::is_any_of(" \t"));
   x_0 = stod(vec_line[0]); 
   // and store it: 
   x_range[0] = x_0;
   getline(fin, line);
   boost::split( vec_line, line, boost::is_any_of(" \t"));
   x = stod(vec_line[0]); 
   Dx = x - x_0;
   n_values = 2;
   x_0 = x;
   while ( getline(fin, line) ) {
      n_values += 1;
      boost::split( vec_line, line, boost::is_any_of(" \t"));
      x = stod(vec_line[0]); 
      if (fabs(Dx - x + x_0) > 1e-6) { 
        FFEA_error_text();
        msg("Aborting: delta x was found to be non-uniform in:");
        msg(ssfile.str());
        return FFEA_ERROR;
      } 
      x_0 = x; 
   }  
   fin.close();

   /*--------- SECONDLY ----------*/ 
   // allocate:
   U = new scalar[n_values * nint];    
   F = new scalar[n_values * nint];    
      
   // and load potentials and forces:
   Dimensions dimens; 
   read_tabulated_values(*pc_params, "pot", U, pc_params->E_to_J / dimens.meso.Energy);
   if (pc_params->inputData == 1) {
     // scalar F_to_Jm = pc_params->E_to_J / pc_params->dist_to_m;
     scalar F_scale = ( pc_params->E_to_J / pc_params->dist_to_m ) / dimens.meso.force;
     read_tabulated_values(*pc_params, "force", F, F_scale);
   } else if (pc_params->inputData == 2) {
     calc_force_from_pot();
   } else {
     FFEA_error_text();
     msg("invalid value for precomp->inputData");
     return FFEA_ERROR;
   }
   Dx = Dx * pc_params->dist_to_m / dimens.meso.length ;
   x_range[0] /= dimens.meso.length; 
   x_range[1] = x_range[0] + Dx * n_values;



   /*------------ THIRDLY --------*/
   // allocate the list of elements that have a bead will use: 
   for (int i=0; i < params->num_blobs; i ++) {
     // we only consider one conformation per blob.
     if (params->num_conformations[i] > 1) {
       FFEA_error_text();
       msg("currently PreComp only deals with a single conformation per blob");
       return FFEA_ERROR;
     }
     n_beads += blob_array[i][0].get_num_beads();
   } 
   b_elems = new TELPtr[n_beads];
   // allocate the array that store the relative positions 
   //    of the beads to the elements where they belong to. 
   b_rel_pos = new scalar[n_beads*3];
   // allocate the array to compute the absolute positions of the beads:
   b_pos = new scalar[n_beads*3];
   // and allocate the bead types: 
   b_types = new int[n_beads]; 
   
   

   /*------------ FOURTHLY --------*/
   // get the elements of the list of "elements" that we will use: 
   vector3 u, v, w; 
   // vector3 s, e1, e2, e3;
   tetra_element_linear *e;
   matrix3 J, J_inv; // will hold the Jacobian for the current element. 
   scalar det;  // determinant for J.
   int m = 0;
   int n;
   // for each Blob: 
   for (int i=0; i < params->num_blobs; i ++) {
     // store the bead types: 
     n = blob_array[i][0].get_num_beads();
     memcpy(&b_types[m], blob_array[i][0].get_bead_type_ptr(), n*sizeof(int));

     // for each bead within this blob (remember that we only deal with conf 0):
     for (int j=0; j < n; j++) {
       v = blob_array[i][0].get_bead_position(j);
       d2_0 = 1e9;
       // get the closest node to this bead: 
       for (int k=0; k < blob_array[i][0].get_num_elements(); k++) { 
         e = blob_array[i][0].get_element(k);
         e->calc_centroid();
         d2 = (e->centroid.x - v.x)*(e->centroid.x - v.x) + 
              (e->centroid.y - v.y)*(e->centroid.y - v.y) + 
              (e->centroid.z - v.z)*(e->centroid.z - v.z);
       
         if (d2 < d2_0) {
           d2_0 = d2;
           b_elems[m+j] = e; 
         }
       } 
       // and get the relative coordinates within the element
       //   as a fraction of the basis vectors length.
       b_elems[m+j]->calculate_jacobian(J); 
       mat3_invert(J, J_inv, &det);
       vec3_vec3_subs(&v, &b_elems[j+m]->n[0]->pos, &w);
       vec3_mat3_mult(&w, J_inv, &u); 
       // now u has the relative coordinates, not under unit vectors
       //    but under full length vectors. And we store them:
       b_rel_pos[m+3*j] = u.x;
       b_rel_pos[m+3*j+1] = u.y;
       b_rel_pos[m+3*j+2] = u.z;
       
       /*
       //  prove it: v =? s
       vec3_vec3_subs(&b_elems[j]->n[1]->pos, &b_elems[j]->n[0]->pos, &e1);
       vec3_vec3_subs(&b_elems[j]->n[2]->pos, &b_elems[j]->n[0]->pos, &e2);
       vec3_vec3_subs(&b_elems[j]->n[3]->pos, &b_elems[j]->n[0]->pos, &e3);
       s.x = b_elems[j]->n[0]->pos.x + u.x*e1.x + u.y*e2.x + u.z*e3.x;
       s.y = b_elems[j]->n[0]->pos.y + u.x*e1.y + u.y*e2.y + u.z*e3.y;
       s.z = b_elems[j]->n[0]->pos.z + u.x*e1.z + u.y*e2.z + u.z*e3.z;
       print_vector3(&v);
       print_vector3(&s);
       s.x = b_elems[j]->n[0]->pos.x + u.x*J[0][0] + u.y*J[1][0] + u.z*J[2][0];
       s.y = b_elems[j]->n[0]->pos.y + u.x*J[0][1] + u.y*J[1][1] + u.z*J[2][1];
       s.z = b_elems[j]->n[0]->pos.z + u.x*J[0][2] + u.y*J[1][2] + u.z*J[2][2];
       print_vector3(&s);
       */

       /* the following is useless but useful while testing:
       b_pos[m+3*j] = v.x;
       b_pos[m+3*j+1] = v.y;
       b_pos[m+3*j+2] = v.z;
       */
       
       
     } 
     // and forget all about beads. 
     blob_array[i][0].forget_beads();
     m += n;
   } 
 
   cout << "done!" << endl;

   return FFEA_OK; 
}


int PreComp_solver::solve() {

    scalar d, f_ij, f_ijk_i, f_ijk_j; 
    vector3 dx, dtemp;
    int type_i; 
    scalar phi_i[4], phi_j[4];
    tetra_element_linear *e_i, *e_j;
    matrix3 Ji, Jj; // these are the jacobians for the elements i and j

    // 1 - Compute the position of the beads:
    compute_bead_positions();

    // 2 - Compute all the i-j forces:
    for (int i=0; i<n_beads; i++){ 
      type_i = b_types[i]; 
      phi_i[1] = b_rel_pos[3*i];
      phi_i[2] = b_rel_pos[3*i+1];
      phi_i[3] = b_rel_pos[3*i+2];
      phi_i[0] = 1 - phi_i[1] - phi_i[2] - phi_i[3];
      e_i = b_elems[i];  
      for (int j=i+1; j<n_beads; j++) {
        dx.x = (b_pos[3*j] - b_pos[3*i]);
        dx.y = (b_pos[3*j+1] - b_pos[3*i+1]);
        dx.z = (b_pos[3*j+2] - b_pos[3*i+2]);
        d = mag(&dx);
        dx.x = dx.x / d;
        dx.y = dx.y / d;
        dx.z = dx.z / d;
 
        f_ij = get_F(d, type_i, b_types[j]); 
        e_j = b_elems[j];
        vec3_scale(&dx, f_ij);
        dtemp = dx; 

        phi_j[1] = b_rel_pos[3*j];
        phi_j[2] = b_rel_pos[3*j+1];
        phi_j[3]= b_rel_pos[3*j+2];
        phi_j[0]= 1 - phi_j[1] - phi_j[2] - phi_j[3];
        // and apply the force to all the nodes in the elements i and j: 
        for (int k=0; k<4; k++) {
          // forces for e_i
          vec3_scale(&dx, -phi_i[k]);
          e_i->add_force_to_node(k, &dx);
          dx = dtemp; 
          // forces for e_j
          vec3_scale(&dx, phi_j[k]);
          e_j->add_force_to_node(k, &dx);
          dx = dtemp; 
        } 
      }
    }

  
    return FFEA_OK;
}


int PreComp_solver::compute_bead_positions() {

    matrix3 J; 
    for (int i=0; i<n_beads; i++){
       b_elems[i]->calculate_jacobian(J); 
       b_pos[3*i]   = b_elems[i]->n[0]->pos.x + b_rel_pos[3*i]*J[0][0] + b_rel_pos[3*i+1]*J[1][0] + b_rel_pos[3*i+2]*J[2][0];
       b_pos[3*i+1] = b_elems[i]->n[0]->pos.y + b_rel_pos[3*i]*J[0][1] + b_rel_pos[3*i+1]*J[1][1] + b_rel_pos[3*i+2]*J[2][1];
       b_pos[3*i+2] = b_elems[i]->n[0]->pos.z + b_rel_pos[3*i]*J[0][2] + b_rel_pos[3*i+1]*J[1][2] + b_rel_pos[3*i+2]*J[2][2];
    }
    return FFEA_OK;
}


/** Calculate F as the numerical derivative of U, instead of reading .force files 
  * @param[in] int total_int: total number of different interactions.
  */
int PreComp_solver::calc_force_from_pot() {

   scalar ym, y0, yM;
   scalar twoDx = 2*Dx;

   int index = -1;
   for (int i=0; i<nint; i++) {
     index += 1; 
     F[index] = (U[index] - U[index + 1]) / Dx;
     for (int j=1; j<n_values-1; j++) {
       index += 1;
       F[index] = (U[index -1] - U[index + 1]) / twoDx; 
     }
     index += 1;
     F[index] = 0; 
   } 
   return FFEA_OK;


}

/** Read either the .pot or .force files
  * and store its contents either in U or F */
int PreComp_solver::read_tabulated_values(PreComp_params &pc_params, string kind, scalar *Z, scalar scale_Z){

   stringstream ssfile; 
   ifstream fin;
   string line; 
   vector<string> vec_line; 
   scalar x_0, x;

   int m_values;
   int index = 0;
   // read the potentials/forces in while checking:
   for (int i=0; i<pc_params.types.size(); i++) {
     for (int j=i; j<pc_params.types.size(); j++) { 
       // open file i-j
       ssfile << pc_params.folder << "/" << pc_params.types[i] << "-" << pc_params.types[j] << "." << kind;
       fin.open(ssfile.str(), std::ifstream::in);
       // get the first line that does not start with "#"
       getline(fin, line);
       while (line.find("#", 0, 1) == 0) {
         getline(fin, line); 
       }  
       // get x_0, i. e., parse the line:
       boost::split( vec_line, line, boost::is_any_of(" \t"));
       x_0 = stod(vec_line[0]); 
       // check it: 
       if (x_range[0] != x_0){
          FFEA_error_text();
          msg("different starting point for file:");
          msg(ssfile.str());
          return FFEA_ERROR;
       }
       // and store the potential:
       Z[index] = stod(vec_line[1]) * scale_Z;
       index += 1;
       // get the next line, and check Dx:
       getline(fin, line);
       boost::split( vec_line, line, boost::is_any_of(" \t"));
       x = stod(vec_line[0]); 
       if (fabs(Dx - x + x_0) > 1e-6){
          FFEA_error_text();
          msg("different step size for x in file:");
          msg(ssfile.str());
          return FFEA_ERROR;
       }
       // and store the potential:
       Z[index] = stod(vec_line[1]) * scale_Z;
       index += 1;
       // Now, go for the rest of the file:
       m_values = 2;
       x_0 = x;
       while ( getline(fin, line) ) {
          m_values += 1;
          boost::split( vec_line, line, boost::is_any_of(" \t"));
          x = stod(vec_line[0]); 
          // check Dx at every line:
          if (fabs(Dx - x + x_0) > 1e-6) { 
            FFEA_error_text();
            msg("Aborting; delta x was found to be non-uniform for file:");
            msg(ssfile.str());
            return FFEA_ERROR;
          } 
          x_0 = x; 
          // and check that the file is not too long:
          if (m_values > n_values) { 
            FFEA_error_text();
            msg("Aborting; too many points for file:");
            msg(ssfile.str());
            return FFEA_ERROR;
          } 
          Z[index] = stod(vec_line[1]) * scale_Z;
          index += 1;
       }  
       // and check that it is not too short:
       if (m_values != n_values) { 
         FFEA_error_text();
         msg("Aborting; too many points for file:");
         msg(ssfile.str());
         return FFEA_ERROR;
       } 
       fin.close();
       ssfile.str("");
     }
   } 
     
   return 0;
}

/**
 * @brief Use finterpolate to return the value for U at x between types i and j.
 * 
 * @details Compute (through finterpolate) the value of the potential at x 
 * between bead types typei and typej. Kept just for clarity, may be removed
 * in the future.
 */
scalar PreComp_solver::get_U(scalar x, int typei, int typej) {
 
  return finterpolate(U, x, typei, typej); 

}


/**
 * @brief Use finterpolate to return the value for F at x between types i and j.
 *
 * @details Compute (through finterpolate) the value of the potential at x 
 * between bead types typei and typej. Kept just for clarity, may be removed
 * in the future.
 */
scalar PreComp_solver::get_F(scalar x, int typei, int typej) {
 
  return finterpolate(F, x, typei, typej); 

}

/** @brief Get the value Z value (either U or F) at x between types typei and typej.
  * @details The function is currently private and accessed through get_U and get_F.
  * For production purposes it may be reasonable to move it to public, and 
  * remove get_U and get_F
  */
scalar PreComp_solver::finterpolate(scalar *Z, scalar x, int typei, int typej){

   scalar y0, y1, x0, x1;
   int index = 0;
   int index_l = x/Dx;
   if (index_l < 0) 
     cout << "WTF?!" << endl; 

   // check that the index is not too high (all the tables are equally long): 
   if (index_l > n_values) {
      // cout << "returned zero for x: " << x << endl; 
      return 0.; 
   } 

   // sort so that typei <= typej
   if (typei > typej){
     int tmp = typei;
     typei = typej;
     typej = tmp;
   } 


   // get the index for the closest (bottom) value:
   for (int i=0; i<ntypes; i++){
     for (int j=i; j<ntypes; j++) {
       if ((typei == i) and (typej == j))
         goto end_loop;
       index += 1;
     }
   }
   end_loop:
   index = index*n_values;
   index += index_l;


   // interpolate:
   x0 = index_l * Dx; 
   x1 = x0 + Dx; 
   return Z[index] + (Z[index +1] - Z[index])*(x - x0)/Dx;

}  

/**@}*/
