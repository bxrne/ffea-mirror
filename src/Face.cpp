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

#include "Face.h"

Face::Face() {
    n[0] = NULL;
    n[1] = NULL;
    n[2] = NULL;
    n[3] = NULL;
    e = NULL;
    vdw_interaction_type = -1;
    area_0 = 0;
    zero_force();
    num_blobs = 0;
    vdw_xz_interaction_flag = false;
    vdw_bb_interaction_flag = NULL;
    kinetically_active = false;
    vdw_bb_force = NULL;
    vdw_bb_energy = NULL;
    vdw_xz_force = NULL;
    vdw_xz_energy = 0.0;
    daddy_blob = NULL;
}

Face::~Face() {
    n[0] = NULL;
    n[1] = NULL;
    n[2] = NULL;
    n[3] = NULL;
    e = NULL;
    vdw_interaction_type = -1;
    area_0 = 0;
    zero_force();
    num_blobs = 0;
    vdw_xz_interaction_flag = false;
    delete[] vdw_bb_interaction_flag;
    vdw_bb_interaction_flag = NULL;
    kinetically_active = false;
    delete[] vdw_bb_force;
    vdw_bb_force = NULL;
    delete[] vdw_bb_energy;
    vdw_bb_energy = NULL;
    delete[] vdw_xz_force;
    vdw_xz_force = NULL;
    vdw_xz_energy = 0.0;
    daddy_blob = NULL;
}

int Face::init(int index, tetra_element_linear *e, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *opposite, SecondOrderFunctions::stu centroid_stu, Blob *daddy_blob, SimulationParams *params) {

    this->index = index;
    this->e = e;
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
    n[3] = opposite;

    calc_area_normal_centroid();
    area_0 = area;

    this->centroid_stu.s = centroid_stu.s;
    this->centroid_stu.t = centroid_stu.t;
    this->centroid_stu.u = centroid_stu.u;

    this->num_blobs = params->num_blobs;
    vdw_bb_force = new vector3[num_blobs];
    vdw_bb_energy = new scalar[num_blobs];
    vdw_bb_interaction_flag = new bool[num_blobs];
    vdw_xz_force = new vector3;
    if (vdw_bb_force == NULL || vdw_bb_energy == NULL || vdw_bb_interaction_flag == NULL || vdw_xz_force == NULL) FFEA_ERROR_MESSG("Failed to store vectors in Face::init\n"); 

    for(int i = 0; i < num_blobs; ++i) {
	vector3_set_zero(&vdw_bb_force[i]);
        vdw_bb_energy[i] = 0.0;
	vdw_bb_interaction_flag[i] = false;
    }
    vector3_set_zero(vdw_xz_force);
    vdw_xz_energy = 0.0;


    this->daddy_blob = daddy_blob;

    return FFEA_OK;
}

int Face::init(int index, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *opposite, Blob *daddy_blob, SimulationParams *params) {

    this->index = index;
    this->e = NULL;
    n[0] = n0;
    n[1] = n1;
    n[2] = n2;
    n[3] = opposite;

    calc_area_normal_centroid();
    area_0 = area;

    this->centroid_stu.s = 0;
    this->centroid_stu.t = 0;
    this->centroid_stu.u = 0;

    this->num_blobs = params->num_blobs;
    vdw_bb_force = new vector3[num_blobs];
    vdw_bb_energy = new scalar[num_blobs];
    vdw_bb_interaction_flag = new bool[num_blobs];
    vdw_xz_force = new vector3;
    if (vdw_bb_force == NULL || vdw_bb_energy == NULL || vdw_bb_interaction_flag == NULL || vdw_xz_force == NULL) FFEA_ERROR_MESSG("Failed to store vectors in Face::init\n"); 

    for(int i = 0; i < num_blobs; ++i) {
        vector3_set_zero(&vdw_bb_force[i]);
        vdw_bb_energy[i] = 0.0;
        vdw_bb_interaction_flag[i] = false;
    }
    vector3_set_zero(vdw_xz_force);
    vdw_xz_energy = 0.0;

    zero_force();

    this->daddy_blob = daddy_blob;
}

void Face::set_vdw_interaction_type(int vdw_interaction_type) {
    this->vdw_interaction_type = vdw_interaction_type;
}

int Face::build_opposite_node() {

    // If opposite == NULL and VdW_type == steric, we can define a node specifically for this face, to make an 'element'
    // It will contain position data (maybe add some error checks in future)
    //   and same index as node 0, n[0]. It will be harmless, because this index
    //   is meant to add forces, but this is a static blob.

    if(n[3] == NULL) {
	n[3] = new mesh_node();
   if (n[3] == NULL) FFEA_ERROR_MESSG("Couldn't find memory for an opposite node\n"); 
	n[3]->index = n[0]->index;

	// Calculate where to stick it
	scalar length;
	vector3 in;

	// Get inward normal
	calc_area_normal_centroid();
	in.x = -1 * normal.x;
	in.y = -1 * normal.y;
	in.z = -1 * normal.z;

	// Now get a lengthscale
	length = sqrt(area);

	// Now place new node this far above the centroid
	n[3]->pos_0.x = centroid.x + length * in.x;
	n[3]->pos_0.y = centroid.y + length * in.y;
	n[3]->pos_0.z = centroid.z + length * in.z;
	n[3]->pos.x = n[3]->pos_0.x;
	n[3]->pos.y = n[3]->pos_0.y;
	n[3]->pos.z = n[3]->pos_0.z;
    }

   zero_force();
   
   return FFEA_OK;
}

void Face::set_kinetic_state(bool state) {
	this->kinetically_active = state;
}

void Face::calc_area_normal_centroid() {
    // (1/2) * |a x b|
    /* new_vector3:
    vector3 a, b;
    a.data = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z};
    b.data = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};
    */
    vector3 a = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z},
    b = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};
    normal.x = a.y * b.z - a.z * b.y;
    normal.y = a.z * b.x - a.x * b.z;
    normal.z = a.x * b.y - a.y * b.x;

    scalar normal_mag = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);

    area = .5 * normal_mag;

    // Normalise normal
    normal.x /= normal_mag;
    normal.y /= normal_mag;
    normal.z /= normal_mag;

    // Find centroid
    centroid.x = (1.0 / 3.0) * (n[0]->pos.x + n[1]->pos.x + n[2]->pos.x);
    centroid.y = (1.0 / 3.0) * (n[0]->pos.y + n[1]->pos.y + n[2]->pos.y);
    centroid.z = (1.0 / 3.0) * (n[0]->pos.z + n[1]->pos.z + n[2]->pos.z);
}

vector3 * Face::get_centroid() {
	centroid.x = (1.0 / 3.0) * (n[0]->pos.x + n[1]->pos.x + n[2]->pos.x);
	centroid.y = (1.0 / 3.0) * (n[0]->pos.y + n[1]->pos.y + n[2]->pos.y);
	centroid.z = (1.0 / 3.0) * (n[0]->pos.z + n[1]->pos.z + n[2]->pos.z);

	return &centroid;
}

void Face::print_centroid() {

	vector3 *c;
	c = get_centroid();
	fprintf(stderr, "Centroid: %f %f %f\n", c->x, c->y, c->z);
}

void Face::print_nodes() {

	for(int i = 0; i < 4; ++i) {
		fprintf(stderr, "Node %d: %f %f %f\n", i, n[i]->pos.x, n[i]->pos.y, n[i]->pos.z);
	}
	fprintf(stderr, "\n");
}

scalar Face::get_area() {

    vector3 temp;

    // (1/2) * |a x b|
    /* new_vector3:
    vector3 a, b;
    a.data = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z},
    b.data = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};
    */
    vector3 a = {n[1]->pos.x - n[0]->pos.x, n[1]->pos.y - n[0]->pos.y, n[1]->pos.z - n[0]->pos.z},
    b = {n[2]->pos.x - n[0]->pos.x, n[2]->pos.y - n[0]->pos.y, n[2]->pos.z - n[0]->pos.z};

    temp.x = a.y * b.z - a.z * b.y;
    temp.y = a.z * b.x - a.x * b.z;
    temp.z = a.x * b.y - a.y * b.x;

    scalar normal_mag = sqrt(temp.x * temp.x + temp.y * temp.y + temp.z * temp.z);

    area = .5 * normal_mag;
    return area;

}

/* Calculate the point p on this triangle given the barycentric coordinates b1, b2, b3 */
void Face::barycentric_calc_point(scalar b1, scalar b2, scalar b3, vector3 *p) {
    p->x = b1 * n[0]->pos.x + b2 * n[1]->pos.x + b3 * n[2]->pos.x;
    p->y = b1 * n[0]->pos.y + b2 * n[1]->pos.y + b3 * n[2]->pos.y;
    p->z = b1 * n[0]->pos.z + b2 * n[1]->pos.z + b3 * n[2]->pos.z;
}

void Face::barycentric_calc_point_f2(scalar b1, scalar b2, scalar b3, vector3 *p,scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index) {
    p->x = b1 * (n[0]->pos.x-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3])  + b2 * (n[1]->pos.x-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3]) + b3 * (n[2]->pos.x-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3]);
    p->y = b1 * (n[0]->pos.y-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3+1]) + b2 *(n[1]->pos.y-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3+1]) + b3 * (n[2]->pos.y-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3+1]);
    p->z = b1 * (n[0]->pos.z-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3+2]) + b2 * (n[1]->pos.z-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3+2]) + b3 * (n[2]->pos.z-blob_corr[f2_daddy_blob_index*(this->num_blobs)*3 + f1_daddy_blob_index*3+2]);
    //printf("blob_corr element for blob %d to blob %d is %f \n ",f1_daddy_blob_index,f2_daddy_blob_index,blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3]);

}

/* Returns the average electrostatic potential of this face */
scalar Face::average_phi() {
    return (1.0 / 3.0) * (n[0]->phi + n[1]->phi + n[2]->phi);
}

scalar Face::get_normal_flux() {
    vector3 dphi;
    e->get_grad_phi_at_stu(&dphi, centroid_stu.s, centroid_stu.t, centroid_stu.u);
    return dphi.x * normal.x + dphi.y * normal.y + dphi.z * normal.z;
}

/*
                scalar get_normal_flux()
                {
                        vector3 dphi =  {
                                        e->n[0]->phi * e->dpsi[0] + e->n[1]->phi * e->dpsi[1] + e->n[2]->phi * e->dpsi[2] + e->n[3]->phi * e->dpsi[3],
                                        e->n[0]->phi * e->dpsi[4] + e->n[1]->phi * e->dpsi[5] + e->n[2]->phi * e->dpsi[6] + e->n[3]->phi * e->dpsi[7],
                                        e->n[0]->phi * e->dpsi[8] + e->n[1]->phi * e->dpsi[9] + e->n[2]->phi * e->dpsi[10] + e->n[3]->phi * e->dpsi[11]
                                        };
                        return dphi.x * normal.x + dphi.y * normal.y + dphi.z * normal.z;
                }
 */

template <class brr3> void Face::add_force_to_node(int i, brr3 (&f)) {
    force[i].x += f[0];
    force[i].y += f[1];
    force[i].z += f[2];
}

/* void Face::add_force_to_node(int i, arr3 (&f)) {
    force[i].x += f[0];
    force[i].y += f[1];
    force[i].z += f[2];
}*/

void Face::add_force_to_node(int i, vector3 *f) {
    force[i].x += f->x;
    force[i].y += f->y;
    force[i].z += f->z;
}

void Face::add_force_to_node_atomic(int i, vector3 *f) {
    force[i].x += f->x;
    force[i].y += f->y;
    force[i].z += f->z;
}

void Face::add_bb_vdw_force_to_record(vector3 *f, int other_blob_index) {
    vdw_bb_force[other_blob_index].x += f->x;
    vdw_bb_force[other_blob_index].y += f->y;
    vdw_bb_force[other_blob_index].z += f->z;
}

template <class brr3> void Face::add_bb_vdw_force_to_record(brr3 &f, int other_blob_index) {
    vdw_bb_force[other_blob_index].x += f[0];
    vdw_bb_force[other_blob_index].y += f[1];
    vdw_bb_force[other_blob_index].z += f[2];
}

void Face::add_bb_vdw_energy_to_record(scalar energy, int other_blob_index) {
    vdw_bb_energy[other_blob_index] += energy;
}

void Face::add_xz_vdw_force_to_record(vector3 *f) {
    vdw_xz_force->x += f->x;
    vdw_xz_force->y += f->y;
    vdw_xz_force->z += f->z;
}

void Face::add_xz_vdw_energy_to_record(scalar energy) {
    vdw_xz_energy += energy;
}

void Face::zero_force() {
    for (int i = 0; i < 4; i++) {
        force[i].x = 0;
        force[i].y = 0;
        force[i].z = 0;
    }

}

void Face::zero_vdw_bb_measurement_data() {
    for (int i = 0; i < num_blobs; ++i) {
        vdw_bb_force[i].x = 0.0;
        vdw_bb_force[i].y = 0.0;
        vdw_bb_force[i].z = 0.0;
        vdw_bb_energy[i] = 0.0;
    }
}

void Face::zero_vdw_xz_measurement_data() {
    vdw_xz_force->x = 0.0;
    vdw_xz_force->y = 0.0;
    vdw_xz_force->z = 0.0;
    vdw_xz_energy = 0.0;
}

void Face::set_vdw_xz_interaction_flag(bool state) {
    vdw_xz_interaction_flag = state;
}

void Face::set_vdw_bb_interaction_flag(bool state, int other_blob_index) {
    vdw_bb_interaction_flag[other_blob_index] = state;
}

bool Face::is_vdw_active() {
    if (vdw_interaction_type == -1) {
        return false;
    } else {
        return true;
    }
}

bool Face::is_kinetic_active() {
	return kinetically_active;
}

bool Face::checkTetraIntersection(Face *f2) {

  // V1 = n
  // V2 = f2->n
  scalar tetA[4][3], tetB[4][3];
  for (int i=0; i<4; i++) {
     tetA[i][0] = this->n[i]->pos.x;
     tetA[i][1] = this->n[i]->pos.y;
     tetA[i][2] = this->n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  return (tet_a_tet(tetA, tetB));

}

scalar Face::getTetraIntersectionVolume(Face *f2){

  geoscalar tetA[4][3], tetB[4][3];

  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  return volumeIntersection<geoscalar,grr3>(tetA, tetB);


}

scalar Face::checkTetraIntersectionAndGetVolume(Face *f2){

  // V1 = n
  // V2 = f2->n
  scalar tetA[4][3], tetB[4][3];
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  if (!tet_a_tet(tetA, tetB)) return 0.0;
  return volumeIntersection<scalar,arr3>(tetA, tetB);



}

void Face::getTetraIntersectionVolumeAndArea(Face *f2, geoscalar &vol, geoscalar &area){

  geoscalar tetA[4][3], tetB[4][3];

  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  volumeAndAreaIntersection<geoscalar,grr3>(tetA, tetB, vol, area);

}

void Face::getTetraIntersectionVolumeAndGradient(Face *f2, grr3 &r, geoscalar &vol, geoscalar &dVdr){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3];
  geoscalar dr = 1e-3;
  geoscalar vol_m, vol_M;

  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x -dr*r[0];
     tetB[i][1] = f2->n[i]->pos.y -dr*r[1];
     tetB[i][2] = f2->n[i]->pos.z -dr*r[2];
     tetC[i][0] = f2->n[i]->pos.x + dr*r[0];
     tetC[i][1] = f2->n[i]->pos.y + dr*r[1];
     tetC[i][2] = f2->n[i]->pos.z + dr*r[2];
  }
  vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetB);
  vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);

  dVdr = (vol_M - vol_m)/(2*dr);
  vol  = (vol_m + vol_M)/2;

}

bool Face::getTetraIntersectionVolumeAndShapeFunctions(Face *f2, grr3 (&r), grr3 (&p), geoscalar &vol, grr4 (&phi1), grr4 (&phi2)){

  geoscalar tetA[4][3], tetB[4][3];
  grr3 ap1, ap2, p2, cm;

  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  vol = volumeIntersection<geoscalar,grr3>(tetA, tetB, cm);

  arr3arr3Add<geoscalar,grr3>(cm, r, p2);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetA[0], tetA[1], tetA[2], tetA[3], cm, phi1);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetB[0], tetB[1], tetB[2], tetB[3], cm, phi2);


  return true;

}


bool Face::getTetraIntersectionVolumeGradientAndShapeFunctions(Face *f2, grr3 (&r), geoscalar &vol, geoscalar &dVdr, grr4 (&phi1), grr4 (&phi2)){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3];
  grr3 ap1, ap2, cm1, cm2;
  geoscalar dr = 1e-3;
  geoscalar vol_m, vol_M;

  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x -dr*r[0];
     tetB[i][1] = f2->n[i]->pos.y -dr*r[1];
     tetB[i][2] = f2->n[i]->pos.z -dr*r[2];
     tetC[i][0] = f2->n[i]->pos.x + dr*r[0];
     tetC[i][1] = f2->n[i]->pos.y + dr*r[1];
     tetC[i][2] = f2->n[i]->pos.z + dr*r[2];
  }

  vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetB, cm1);
  vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC, cm2);

  dVdr = (vol_M - vol_m)/(2.*dr);
  vol = (vol_M + vol_m)/2.;

  arr3arr3Add<geoscalar,grr3>(cm1, cm2, cm1);
  arr3Resize<geoscalar,grr3>(0.5, cm1);

  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetA[0], tetA[1], tetA[2], tetA[3], cm1, phi1);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetB[0], tetB[1], tetB[2], tetB[3], cm1, phi2);

  return true;

}


bool Face::getTetraIntersectionVolumeGradientDirAndShapeFunctions(Face *f2, grr3 (&r), geoscalar &vol, geoscalar &dVdr, grr4 (&phi1), grr4 (&phi2)){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3], tetD[4][3];
  grr3 ap1, ap2, cm1;
  geoscalar dr = 1e-3;
  geoscalar vol_m, vol_M;

  // GET THE VOLUME AND THE DIRECTION OF THE GRADIENT:
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  // get the volume, the CM for the intersection, and the direction of the gradient:
  vol = volumeIntersection<geoscalar,grr3>(tetA, tetB, cm1, r);
  if (vol == 0) return false;

  // GET THE GRADIENT
  // 2nd Order:
  //   // a simple 1st order derivative was not working for the cubes & springs.
  // construct the two tetrahedra B, in the direction of the gradient
  for (int i=0; i<4; i++) {
     tetC[i][0] = f2->n[i]->pos.x + dr*r[0];
     tetC[i][1] = f2->n[i]->pos.y + dr*r[1];
     tetC[i][2] = f2->n[i]->pos.z + dr*r[2];
     tetD[i][0] = f2->n[i]->pos.x -dr*r[0];
     tetD[i][1] = f2->n[i]->pos.y -dr*r[1];
     tetD[i][2] = f2->n[i]->pos.z -dr*r[2];
  }

  vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetD);
  vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);


  dVdr = (vol_M - vol_m)/(2*dr);

  // 1st Order:
  //   // a simple 1st order derivative was not working for the cubes & springs.
  // construct the two tetrahedra B, in the direction of the gradient
  // for (int i=0; i<4; i++) {
  //    tetC[i][0] = f2->n[i]->pos.x + dr*r[0];
  //    tetC[i][1] = f2->n[i]->pos.y + dr*r[1];
  //    tetC[i][2] = f2->n[i]->pos.z + dr*r[2];
  // }

  // vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);
  // dVdr = (vol_M - vol)/dr;


  // GET THE LOCAL COORDINATES
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetA[0], tetA[1], tetA[2], tetA[3], cm1, phi1);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetB[0], tetB[1], tetB[2], tetB[3], cm1, phi2);

  return true;

}


bool Face::getTetraIntersectionVolumeGradientAndShapeFunctions(Face *f2, grr3 (&dVdr), geoscalar &vol, grr4 (&phi1), grr4 (&phi2)){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3], tetD[4][3];
  grr3 ap1, ap2, cm;
  geoscalar dr = 1e-2;
  geoscalar vol_m, vol_M;

  // GET THE VOLUME AND THE DIRECTION OF THE GRADIENT:
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x;
     tetB[i][1] = f2->n[i]->pos.y;
     tetB[i][2] = f2->n[i]->pos.z;
  }
  // get the volume, the CM for the intersection, and the direction of the gradient:
  vol = volumeIntersection<geoscalar,grr3>(tetA, tetB, cm); // , r);
  if (vol == 0) return false;

  // GET THE LOCAL COORDINATES where the force will be applied.
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetA[0], tetA[1], tetA[2], tetA[3], cm, phi1);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetB[0], tetB[1], tetB[2], tetB[3], cm, phi2);


  // GET THE GRADIENT
  // 1st Order:
  grr3 dx;
  for (int dir=0; dir<3; dir++) {
    arr3Initialise<grr3>(dx);
    dx[dir] = 1;
    for (int i=0; i<4; i++) {
       tetC[i][0] = f2->n[i]->pos.x + dr*dx[0];
       tetC[i][1] = f2->n[i]->pos.y + dr*dx[1];
       tetC[i][2] = f2->n[i]->pos.z + dr*dx[2];
    }
    vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);
    dVdr[dir] = (vol_M - vol)/dr;
  }

  /* // 2nd Order:
  grr3 dx;
  for (int dir=0; dir<3; dir++) {
    arr3Initialise<grr3>(dx);
    dx[dir] = 1;
    for (int i=0; i<4; i++) {
       tetC[i][0] = f2->n[i]->pos.x + dr*dx[0];
       tetC[i][1] = f2->n[i]->pos.y + dr*dx[1];
       tetC[i][2] = f2->n[i]->pos.z + dr*dx[2];
       tetD[i][0] = f2->n[i]->pos.x -dr*dx[0];
       tetD[i][1] = f2->n[i]->pos.y -dr*dx[1];
       tetD[i][2] = f2->n[i]->pos.z -dr*dx[2];
    }
    vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetD);
    vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);
    r[dir] = (vol_M - vol_m)/(2*dr);
  }
  dVdr = mag<geoscalar,grr3>(r);
  arr3Normalise<geoscalar,grr3>(r);  */


  return true;

}

bool Face::checkTetraIntersection(Face *f2,scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index) {

  // V1 = n
  // V2 = f2->n
  scalar tetA[4][3], tetB[4][3];
  double /*tempx=0, tempy=0,tempz=0,*/assx=0, assy=0, assz=0;
  for (int i=0; i<4; i++) {
     tetA[i][0] = this->n[i]->pos.x;
     tetA[i][1] = this->n[i]->pos.y;
     tetA[i][2] = this->n[i]->pos.z;

     tetB[i][0] = f2->n[i]->pos.x -blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     tetB[i][1] = f2->n[i]->pos.y -blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z -blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];
     //printf("blob_corr element for blob %d to blob %d is %f \n corrected location for f2.x is %f\n face %d x is at %f\n face %d x is at %f\n",f1_daddy_blob_index,f2_daddy_blob_index,blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3],tetB[i][0],this->index,this->n[i]->pos.x,f2->index,f2->n[i]->pos.x);
     //printf("blob_corr element for blob %d to blob %d is %F \n corrected location for f2.y is %F\n face %d y is at %F\n face %d y is at %F\n",f1_daddy_blob_index,f2_daddy_blob_index,blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1],tetB[i][1],this->index,this->n[i]->pos.y,f2->index,f2->n[i]->pos.y);
     //printf("blob_corr element for blob %d to blob %d is %F \n corrected location for f2.z is %F\n face %d z is at %F\n face %d z is at %F\n",f1_daddy_blob_index,f2_daddy_blob_index,blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2],tetB[i][2],this->index,this->n[i]->pos.z,f2->index,f2->n[i]->pos.z);
  }
  return (tet_a_tet(tetA, tetB));

}

scalar Face::getTetraIntersectionVolume(Face *f2, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){

  geoscalar tetA[4][3], tetB[4][3];
  double assx=0, assy=0, assz=0;
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;

     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];
  }
  return volumeIntersection<geoscalar,grr3>(tetA, tetB);


}

scalar Face::checkTetraIntersectionAndGetVolume(Face *f2, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){

  // V1 = n
  // V2 = f2->n
  scalar tetA[4][3], tetB[4][3];
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];
  }
  if (!tet_a_tet(tetA, tetB)) return 0.0;
  return volumeIntersection<scalar,arr3>(tetA, tetB);



}

void Face::getTetraIntersectionVolumeAndArea(Face *f2, geoscalar &vol, geoscalar &area,scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){

  geoscalar tetA[4][3], tetB[4][3];
  double assx=0, assy=0, assz=0;
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;

     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     //printf("Volume and intersection starts with %f\n f2.x is %f\nf1.x is %f\n blob_corr is %f\n",tetB[i][0],f2->n[i]->pos.x,this->n[i]->pos.x,blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3]);
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];
  }
  volumeAndAreaIntersection<geoscalar,grr3>(tetA, tetB, vol, area);
  //printf("Face %d and Face %d volume is %f area is %f\n",f2->index,this->index,vol,area);

}

template <class brr3> void Face::vec3Vec3SubsToArr3Mod(Face *f2, brr3 (&w),scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){
    //this->get_centroid();
     //f2->get_centroid();
    w[0] = f2->n[3]->pos.x-this->n[3]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
    w[1] = f2->n[3]->pos.y-this->n[3]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
    w[2] = f2->n[3]->pos.z-this->n[3]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];

    //printf("\n corr x is %f \n corr y is %f \ncorr z is %f \n f2.x  is %f\n f2.y  is %f\n f2.z  is %f\n f1.x  is %f\n f1.y  is %f\n f1.z  is %f\n",blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3],blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1],blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2],f2->n[3]->pos.x, f2->n[3]->pos.y, f2->n[3]->pos.z,this->n[3]->pos.x, this->n[3]->pos.y, this->n[3]->pos.z);
}


void Face::getTetraIntersectionVolumeAndGradient(Face *f2, grr3 &r, geoscalar &vol, geoscalar &dVdr, scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3];
  geoscalar dr = 1e-3;
  geoscalar vol_m, vol_M;

  double assx=0, assy=0, assz=0;
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;


     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];
     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3] -dr*r[0];
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1] -dr*r[1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2] -dr*r[2];
     tetC[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3] + dr*r[0];
     tetC[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1] + dr*r[1];
     tetC[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2] + dr*r[2];
  }
  vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetB);
  vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);

  dVdr = (vol_M - vol_m)/(2*dr);
  vol  = (vol_m + vol_M)/2;

}

bool Face::getTetraIntersectionVolumeGradientDirAndShapeFunctions(Face *f2, grr3 (&r), geoscalar &vol, geoscalar &dVdr, grr4 (&phi1), grr4 (&phi2), scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3], tetD[4][3];
  grr3 ap1, ap2, cm1;
  geoscalar dr = 1e-3;
  geoscalar vol_m, vol_M;

  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];
  }
  // get the volume, the CM for the intersection, and the direction of the gradient:
  vol = volumeIntersection<geoscalar,grr3>(tetA, tetB, cm1, r);
  if (vol == 0) return false;

  // construct the two tetrahedra B, in the direction of the gradient
  //   // a simple 1st order derivative was not working for the cubes & springs.
  for (int i=0; i<4; i++) {
     tetC[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3] + dr*r[0];
     tetC[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1] + dr*r[1];
     tetC[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2] + dr*r[2];
     tetD[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3] -dr*r[0];
     tetD[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1] -dr*r[1];
     tetD[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2] -dr*r[2];
  }

  vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetD);//cm1);
  vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);//cm2);


  dVdr = (vol_M - vol_m)/(2*dr);

  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetA[0], tetA[1], tetA[2], tetA[3], cm1, phi1);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetB[0], tetB[1], tetB[2], tetB[3], cm1, phi2);

  return true;

}

bool Face::getTetraIntersectionVolumeGradientAndShapeFunctions(Face *f2, grr3 (&dVdr), geoscalar &vol, grr4 (&phi1), grr4 (&phi2), scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index){

  geoscalar tetA[4][3], tetB[4][3], tetC[4][3], tetD[4][3];
  grr3 ap1, ap2, cm;
  geoscalar dr = 1e-2;
  geoscalar vol_m, vol_M;

  // GET THE VOLUME AND THE DIRECTION OF THE GRADIENT:
  for (int i=0; i<4; i++) {
     tetA[i][0] = n[i]->pos.x;
     tetA[i][1] = n[i]->pos.y;
     tetA[i][2] = n[i]->pos.z;
     tetB[i][0] = f2->n[i]->pos.x-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3];
     tetB[i][1] = f2->n[i]->pos.y-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1];
     tetB[i][2] = f2->n[i]->pos.z-blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2];

  }
  // get the volume, the CM for the intersection, and the direction of the gradient:
  vol = volumeIntersection<geoscalar,grr3>(tetA, tetB, cm); // , r);
  if (vol == 0) return false;

  // GET THE LOCAL COORDINATES where the force will be applied.
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetA[0], tetA[1], tetA[2], tetA[3], cm, phi1);
  getLocalCoordinatesForLinTet<geoscalar,grr3,grr4>(tetB[0], tetB[1], tetB[2], tetB[3], cm, phi2);


  // GET THE GRADIENT
  // 1st Order:
  grr3 dx;
  for (int dir=0; dir<3; dir++) {
    arr3Initialise<grr3>(dx);
    dx[dir] = 1;
    for (int i=0; i<4; i++) {
       tetC[i][0] = f2->n[i]->pos.x -blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3]+ dr*dx[0];
       tetC[i][1] = f2->n[i]->pos.y -blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+1]+ dr*dx[1];
       tetC[i][2] = f2->n[i]->pos.z -blob_corr[f1_daddy_blob_index*(this->num_blobs)*3 + f2_daddy_blob_index*3+2]+ dr*dx[2];
    }
    vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);
    dVdr[dir] = (vol_M - vol)/dr;
  }

  /* // 2nd Order:
  grr3 dx;
  for (int dir=0; dir<3; dir++) {
    arr3Initialise<grr3>(dx);
    dx[dir] = 1;
    for (int i=0; i<4; i++) {
       tetC[i][0] = f2->n[i]->pos.x + dr*dx[0];
       tetC[i][1] = f2->n[i]->pos.y + dr*dx[1];
       tetC[i][2] = f2->n[i]->pos.z + dr*dx[2];
       tetD[i][0] = f2->n[i]->pos.x -dr*dx[0];
       tetD[i][1] = f2->n[i]->pos.y -dr*dx[1];
       tetD[i][2] = f2->n[i]->pos.z -dr*dx[2];
    }
    vol_m = volumeIntersection<geoscalar,grr3>(tetA, tetD);
    vol_M = volumeIntersection<geoscalar,grr3>(tetA, tetC);
    r[dir] = (vol_M - vol_m)/(2*dr);
  }
  dVdr = mag<geoscalar,grr3>(r);
  arr3Normalise<geoscalar,grr3>(r);  */


  return true;

}
//////////////////////////////////////////////////
//// // // // Instantiate templates // // // // //
//////////////////////////////////////////////////
template void Face::add_force_to_node<arr3>(int i, arr3 (&f));
template void Face::add_bb_vdw_force_to_record<arr3>(arr3 &f, int other_blob_index);
template void Face::vec3Vec3SubsToArr3Mod<arr3>(Face *f2, arr3 (&w),scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);


#ifndef USE_DOUBLE
template void Face::add_force_to_node<grr3>(int i, grr3 (&f));
template void Face::add_bb_vdw_force_to_record<grr3>(grr3 &f, int other_blob_index);
template void Face::vec3Vec3SubsToArr3Mod<grr3>(Face *f2, grr3 (&w),scalar *blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);
#endif
