#include "SparseMatrixUnknownPattern.h"

SparseMatrixUnknownPattern::SparseMatrixUnknownPattern() {
    num_rows = 0;
    row = NULL;
    diagonal = NULL;
}

SparseMatrixUnknownPattern::~SparseMatrixUnknownPattern() {
    num_rows = 0;
    delete[] row;
    delete[] diagonal;
}

int SparseMatrixUnknownPattern::init(int num_rows, int suggested_initial_size_for_row_vectors) {

    this->num_rows = num_rows;

    // Allocate the array of row vectors
    row = new vector<sparse_entry>[num_rows];

    if (row == NULL) {
        return FFEA_ERROR;
    }

    // Reserve a certain amount of space (for the inevitable large amount of growth that will come)
    for (int i = 0; i < num_rows; i++) {
        row[i].reserve(suggested_initial_size_for_row_vectors);
    }

    // Allocate the array for the diagonal elements (which are treated differently)
    diagonal = new scalar[num_rows];

    if (diagonal == NULL) {
        return FFEA_ERROR;
    }

    // Initialise diagonal to zero
    for (int i = 0; i < num_rows; i++) {
        diagonal[i] = 0;
    }

    return FFEA_OK;
}

void SparseMatrixUnknownPattern::add_off_diagonal_element(int row_index, int column_index, scalar val) {
    sparse_entry new_entry = {column_index, val};
    row[row_index].push_back(new_entry);
}

void SparseMatrixUnknownPattern::set_diagonal_element(int row_index, scalar val) {
    diagonal[row_index] = val;
}

void SparseMatrixUnknownPattern::calc_inverse_diagonal(scalar *inv_D) {
    for (int i = 0; i < num_rows; i++) {
        inv_D[i] = 1.0 / diagonal[i];
    }
}

void SparseMatrixUnknownPattern::zero() {
    for (int i = 0; i < num_rows; i++) {
        row[i].clear();
        diagonal[i] = 0;
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result' */
void SparseMatrixUnknownPattern::apply(scalar *in, scalar *result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] = diagonal[i] * in[i];
        int row_size = row[i].size();
        for (int j = 0; j < row_size; j++) {
            result[i] += row[i][j].val * in[row[i][j].column_index];
        }
    }
}

void SparseMatrixUnknownPattern::print() {
    int num_nonzero_elements = 0;
    for (int i = 0; i < num_rows; i++) {
        printf("Row %d:  ", i);
        printf("[%d %e] ", i, diagonal[i]);
        for (vector<sparse_entry>::iterator it = row[i].begin(); it != row[i].end(); it++) {
            printf("[%d %e] ", (*it).column_index, (*it).val);
        }
        printf("\n");
        num_nonzero_elements += row[i].size();
    }
    num_nonzero_elements += num_rows;

    printf("Num nonzero elements = %d, occupation = %f%%\n", num_nonzero_elements, (100.0 * num_nonzero_elements) / (num_rows * num_rows));
}
