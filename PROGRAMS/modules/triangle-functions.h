#include "Matrix.h"
#include <algorithm>
#include "helperFunctions.h"
#ifndef TRIANGLE_FUNCTIONS
#define TRIANGLE_FUNCTIONS

/**
 * @brief project a point c onto pq
 * 
 * @param c 
 * @param p 
 * @param q 
 * @return Matrix 
 */
Matrix project_point_on_segment(Matrix c, Matrix p, Matrix q) {
    double lambda = ((c + -1 * p).transpose() * (q + -1 * p)).matrixArray[0] /
                    ((q + -1 * p).transpose() * (q + -1 * p)).matrixArray[0];
    double lambda_seg = max(0, min(lambda, 1));
    return p + lambda_seg * (q + -1 * p);
}

/**
 * @brief Get the closest point on a triangle to a point 
 * 
 * @param point 
 * @param triangle_p 
 * @param triangle_q 
 * @param triangle_r 
 * @return Matrix 
 */
Matrix get_closest_point_on_triangle(Matrix point, Matrix triangle_p, Matrix triangle_q, Matrix triangle_r) {
    vector<Matrix> A_matrix_cols;
    A_matrix_cols.push_back(triangle_q + -1 * triangle_p);
    A_matrix_cols.push_back(triangle_r + -1 * triangle_p);
    Matrix A = Matrix(A_matrix_cols);
    Matrix result = (A.transpose()*A).inverse()*A.transpose()*(point + -1 * triangle_p);
    double lambda = result.matrixArray[0];
    double u = result.matrixArray[1];
    Matrix c = triangle_p+lambda*(triangle_q + -1 * triangle_p) + u * (triangle_r + -1 * triangle_p);
    if(lambda < 0) {
        c = project_point_on_segment(c, triangle_r, triangle_p); 
    }
    if(u < 0) {
        c = project_point_on_segment(c, triangle_p, triangle_q);
    }
    if(lambda + u > 1) {
        c = project_point_on_segment(c, triangle_q, triangle_r);
    }
    return c;
}

#endif