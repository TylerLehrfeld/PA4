

#include "Matrix.h"

#ifndef BOUNDING_SPHERE
#define BOUNDING_SPHERE

/**
 * @brief a class that creates a bounding sphere around a triangle
 * 
 */
class bounding_sphere {
private:
public:
    Matrix center;
    double radius;
    int* triangle;
    /**
     * @brief Construct a new bounding sphere object
     * 
     * @param triangle_vertex_indeces 
     * @param vertices 
     */
    bounding_sphere(int* triangle_vertex_indeces, vector<Matrix>& vertices) {
        triangle = triangle_vertex_indeces;
        //get longest edge
        double sizeAB = (vertices[triangle_vertex_indeces[0]] + -1 * vertices[triangle_vertex_indeces[1]]).magnitude();
        double sizeBC = (vertices[triangle_vertex_indeces[1]] + -1 * vertices[triangle_vertex_indeces[2]]).magnitude();
        double sizeAC = (vertices[triangle_vertex_indeces[0]] + -1 * vertices[triangle_vertex_indeces[2]]).magnitude();
        Matrix a, b, c;
        if(sizeAB > sizeBC && sizeAB > sizeBC) {
            a = vertices[triangle_vertex_indeces[0]];
            b = vertices[triangle_vertex_indeces[1]];
            c = vertices[triangle_vertex_indeces[2]];
        } else if(sizeBC > sizeAB && sizeBC > sizeAC) {
            a = vertices[triangle_vertex_indeces[1]];
            b = vertices[triangle_vertex_indeces[2]];
            c = vertices[triangle_vertex_indeces[0]];
        } else {
            a = vertices[triangle_vertex_indeces[0]];
            b = vertices[triangle_vertex_indeces[2]];
            c = vertices[triangle_vertex_indeces[1]];
        }
        Matrix f = (a + b)* .5;
        Matrix u = a + -1 *f;
        Matrix v = c + -1 * f;
        Matrix d = (u.cross(v)).cross(u);
        double y = ((v.transpose()*v).matrixArray[0] - (u.transpose()*u).matrixArray[0]) / ((2*d).transpose() * (v + -1*u)).matrixArray[0];
        y = y < 0 ? 0: y;
        center = f + y* d;
        radius = (a + -1 * center).magnitude();
    }
    ~bounding_sphere() {

    };
};

#endif