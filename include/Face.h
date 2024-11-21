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

#ifndef FACE_H_INCLUDED
#define FACE_H_INCLUDED

#include <cmath>

#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SecondOrderFunctions.h"
#include "SimulationParams.h"
#include "CheckTetrahedraOverlap.h"
#include "VolumeIntersection.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"

class Face {
public:
    Face();

    ~Face();

    /** The first 3 pointers are to those surface nodes that define this face,
      * while the last one belongs to the opposite node of the linear tetrahedron **/
    std::array<mesh_node*, 4> n;

    /** Index of this face */
    int index;

    /** Pointer to the element this is a face of */
    tetra_element_linear *e;

    /** Van der Waals interaction type **/
    int ssint_interaction_type;

    /** Initial, equilibrium area of this face **/
    scalar area_0;

    /** Last calculated area (by most recent call to calc_area_normal_centroid()) **/
    scalar area;

    /** Last calculated normal (by most recent call to calc_area_normal_centroid()) **/
    arr3 normal;

    /** Last calculated centroid (by most recent call to calc_area_normal_centroid()) **/
    arr3 centroid;

    /** Stores the current force applied to this face.
      * The first 3 "vectors" are to those surface nodes that define this face,
      * while the last one belongs to the opposite node of the linear tetrahedron **/
    std::array<arr3, 4> force;


    /** Stores the natural (shape function) coords of the centroid of this face in the parent element **/
    SecondOrderFunctions::stu centroid_stu;

    bool ssint_xz_interaction_flag;
    bool kinetically_active;

    //@{
    /** vdw measurement info **/
    int num_blobs;
    //@}

    /** Check whether the tetrahedron formed by this face an the opposite
      *   linear node does intersect with the corresponding tetrahedron in f2
      * Returns true if there is intersection.
      * Uses the "Fast Tetrahedron-Tetrahedron Overlap Algorithm".
      * It calls some private functions.
      **/
    bool checkTetraIntersection(const Face *f2);
    bool checkTetraIntersection(const Face *f2, const std::vector<scalar> &blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    /** Get the volume that the enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      * It calls volumeIntersection, at volumeIntersection.h
      **/
    scalar getTetraIntersectionVolume(const Face *f2);
    scalar getTetraIntersectionVolume(const Face *f2, const std::vector<scalar> &blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    /** Check whether the tetrahedron formed by this face an the opposite
      *   linear node does intersect with the corresponding tetrahedron in f2.
      * If so, return the overlapping volume, otherwise return 0.
      * Uses the "Fast Tetrahedron-Tetrahedron Overlap Algorithm" for checking if
      *  interaction occurs.
      **/
    scalar checkTetraIntersectionAndGetVolume(const Face *f2);
    scalar checkTetraIntersectionAndGetVolume(const Face *f2, const std::vector<scalar> &blob_corr, int f1_daddy_blob_index,int f2_daddy_blob_index);


    /** Get the volume and area that the enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      * It calls volumeIntersection, at volumeIntersection.h
      **/
    void getTetraIntersectionVolumeAndArea(const Face *f2, geoscalar &vol, geoscalar &area);
    void getTetraIntersectionVolumeAndArea(const Face *f2, geoscalar &vol, geoscalar &area, const std::vector<scalar> &blob_corr, int f1_daddy_blob_index,int f2_daddy_blob_index);


     /** Get the volume that enclose the intersection
      *   of the tetrahedron formed by this face an the opposite
      *   linear node with the corresponding tetrahedron in f2.
      *   In addition, return the gradient of this volume,
      *     calculated as dV/dx,dV/dy,dV/dz, and the internal coordinates
      *     of the point where the force is applied.
      * It calls 4 times volumeIntersection.
      * This version works well for the double loop i<j. 
      **/
    bool getTetraIntersectionVolumeTotalGradientAndShapeFunctions(const Face *f2, geoscalar dr, grr3 &dVdr, geoscalar &vol, grr4 &phi1, grr4 &phi2);
    bool getTetraIntersectionVolumeTotalGradientAndShapeFunctions(const Face *f2, geoscalar dr, grr3 &dVdr, geoscalar &vol, grr4 &phi1, grr4 &phi2, const std::vector<scalar> &blob_corr, int f1_daddy_blob_index, int f2_daddy_blob_index);


    Blob *daddy_blob;

    void init(int index, tetra_element_linear *e, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *oposite, SecondOrderFunctions::stu &_centroid_stu, Blob *daddy_blob, const SimulationParams &params);

    void init(int index, mesh_node *n0, mesh_node *n1, mesh_node *n2, mesh_node *opposite, Blob *daddy_blob, const SimulationParams &params);

    void set_ssint_interaction_type(int _ssint_interaction_type);

    void build_opposite_node();

    void calc_area_normal_centroid();
    void set_kinetic_state(bool state);

    arr3 &get_centroid();
    void print_centroid();
    void print_nodes();
    scalar get_area();

    /** Calculate the point p on this triangle given the barycentric coordinates b1, b2, b3. Altered to act periodically around box boundaries. **/
    void barycentric_calc_point(scalar b1, scalar b2, scalar b3, arr3 &p);
    void barycentric_calc_point_f2(scalar b1, scalar b2, scalar b3, arr3 &p, const std::vector<scalar> &blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    /** Returns the average electrostatic potential of this face **/
    scalar average_phi();

    scalar get_normal_flux();

    void add_force_to_node(int i, arr3 &f);
    template <class brr3> void add_force_to_node(int i, brr3 (&f));

    /**
     * @note Implementation matches add_force_to_node(int i, arr3 &f) ??
     */
    void add_force_to_node_atomic(int i, arr3 &f);

    void set_ssint_xz_interaction_flag(bool state);

    template <class brr3>
    void vec3Vec3SubsToArr3Mod(Face *f2, brr3 &w, const std::vector<scalar> &blob_corr,int f1_daddy_blob_index,int f2_daddy_blob_index);

    bool is_ssint_active();
    bool is_kinetic_active();

    scalar length_of_longest_edge();

    void zero_force() { initialise(force); }

private:
    int stuff;
    bool dealloc_n3;

};

#endif
