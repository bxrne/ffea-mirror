#include "Steric_solver.h"


/** do_volumeExclusion calculates the force (and not the energy, yet) of two tetrahedra */
void Steric_solver::do_interaction(Face *f1, Face *f2){ 

    // First check two things (either of which results in not having to calculate anything):
    // Check that faces are facing each other, if not then they are not interacting
    if (dot(&f1->normal, &f2->normal) > ffea_const::zero) {
        return;
    }

    /* Robin suspects that this was leading to unstabilities... 
     *  but this steric solver is more stable than the LJ one. */
    // Check that faces are in front of each other
    vector3 sep = {f2->centroid.x - f1->centroid.x, f2->centroid.y - f1->centroid.y, f2->centroid.z - f1->centroid.z};
    if(dot(&sep, &f1->normal) < 0 && dot(&sep, &f2->normal) > 0) {
        return;
    }

    //  Firstly, check that no nodes are shared:
    if (f1->n[3] == f2->n[3]) {
          return;
    }
    for (int i=0; i<4; i++) {
        int in_i = f1->n[i]->index;
        for (int j=0; j<4; j++) {
           if (f2->n[j]->index == in_i){
               return;
           }
        }
    }

    //  Then, check whether the tetrahedra intersect.
    if (!f1->checkTetraIntersection(f2)) {
      return; 
    }

    //   and get the direction of the force for f1:
    /* TRIAL 2 */ 
    arr3 force1, force2, n1_b; 
    vec3Vec3SubsToArr3(f1->n[3]->pos, f2->n[3]->pos, force1);
    arr3Normalise<scalar,arr3>(force1); // that is the direction of the force for f1 (backwards). 

    /* TRIAL 1
    arr3 force1, force2, n1_b; 
    vec3ResizeToArr3(ffea_const::mOne, f1->normal, n1_b); 
    vec3Arr3AddToArr3(f2->normal, n1_b, force1); 
    arr3Normalise(force1); // that is the direction of the force for f1 (backwards). 
    */ 

    //////////////////////////////////////////////
    /// One more check ////// One more check /////
    arr3 inwards; 
    vec3Vec3SubsToArr3(f1->n[3]->pos, f1->centroid, inwards); 
    if (arr3arr3DotProduct<scalar,arr3>(inwards, force1) < ffea_const::zero) { 
      return; 
    } 
    // end checking //////////////////////////////
    
    // scalar vol_f = 1
    // Finally, get the intersection volume:
    scalar vol = f1->getTetraIntersectionVolume(f2); 
    /*if ((vol < 0) and (fabs(vol) > ffea_const::threeErr)) { 
      cout << "The IntersectoMetre had a precision problem: " << vol << endl;
      return; 
    }*/ 


    arr3Resize(vol, force1); 
    arr3Resize2(ffea_const::mOne, force1, force2);
    
    for (int j = 0; j < 3; j++) {
      f1->add_force_to_node(j, force1); 
      f1->add_bb_vdw_force_to_record(force1, f1->daddy_blob->blob_index);
      f2->add_force_to_node(j, force2); 
      f2->add_bb_vdw_force_to_record(force2, f1->daddy_blob->blob_index);
    } 


}
