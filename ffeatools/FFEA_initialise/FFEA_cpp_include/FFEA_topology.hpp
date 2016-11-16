#ifndef FFEA_TOPOLOGY_INCLUDED
#define FFEA_TOPOLOGY_INCLUDED


#include <iostream>
using namespace std;

// Define an enum for element types
enum elementType {LINEAR, SECONDARY};

class FFEA_element {

	public:

		// Constructor / Destructor
		FFEA_element();
		~FFEA_element();
		
	protected:

		int *index;
		float **pos;
		elementType eltype;
};


// We don't need constructors! Only the assignment functions need overloading...
class FFEA_linear_element: public FFEA_element {

	public:

		// Assignment functions
		int set_indices(int *n);
		int set_indices(int n0, int n1, int n2, int n3);

		int set_pos(float **pos);
		//int set_pos(float *pos0, *pos1, *pos2, *pos3);
};


class FFEA_secondary_element: public FFEA_element {

	public:

		// Assignment functions
		int set_indices(int *n);
		int set_indices(int n0, int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8, int n9);

		int set_pos(float **pos);
		//int set_pos(float *pos0, *pos1, *pos2, *pos3, *pos, *pos5, *pos6, *pos7, *pos8, *pos9);	
};

class FFEA_topology {

	public:

		// Constructor / Destructor
		FFEA_topology();
		~FFEA_topology();

		// Loading Functions
		int load(string fname);
		int load_from_top(string fname);
		int load_from_vol(string fname);

		// Data access functions
		int get_num_elements();

		FFEA_element *element;

	private:

		int num_elements;
};

#endif
