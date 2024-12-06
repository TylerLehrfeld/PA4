#include "bounding_sphere.h"
#include "triangle-functions.h"
#include <malloc.h>
#include <vector>
#include "Matrix.h"
#include <assert.h>
#include <iostream>
using std::malloc;
using std::vector;

#ifndef tri
#define tri

void test_get_closest_point() {
    Matrix point(3,1,{.5,.5,1});
    Matrix p(3,1,{0,0,0}); 
    Matrix q(3,1,{0,1,0});
    Matrix r(3,1,{1,0,0});
    assert(get_closest_point_on_triangle(point, p, q,r) == Matrix(3,1,{.5,.5,0}));
    Matrix point1(3,1,{1,1,1});
    Matrix p1(3,1,{0,0,0}); 
    Matrix q1(3,1,{0,1,0});
    Matrix r1(3,1,{1,0,0});
    assert(get_closest_point_on_triangle(point1, p1, q1,r1) == Matrix(3,1,{.5,.5,0}));
    Matrix point2(3,1,{-1,-1,1});
    Matrix p2(3,1,{0,0,0}); 
    Matrix q2(3,1,{0,1,0});
    Matrix r2(3,1,{1,0,0});
    assert(get_closest_point_on_triangle(point2, p2, q2,r2) == Matrix(3,1,{0,0,0}));
}

void test_bounding_sphere() {
    int* triangle = (int*)malloc(3*sizeof(int));
    triangle[0] = 0;
    triangle[1] = 1;
    triangle[2] = 2;
    vector<Matrix> vertices;
    vertices.push_back(Matrix(3,1,{-1,0,0}));
    vertices.push_back(Matrix(3,1,{1,0,0}));
    vertices.push_back(Matrix(3,1,{0,1,0}));
    bounding_sphere one(triangle, vertices);
    assert(one.center == Matrix(3,1,{0,0,0}));
    assert(one.radius == 1);
    vertices[2].matrixArray[1] = .5;
    bounding_sphere two(triangle, vertices);
    assert(two.center == Matrix(3,1,{0,0,0}));
    assert(two.radius == 1);
    vertices[2].matrixArray[0] = 1;
    vertices[2].matrixArray[1] = 1;
    bounding_sphere three(triangle, vertices);
    assert(three.center == Matrix(3,1,{0,.5,0}));
    assert(three.radius == sqrt(5)/2);
    
    free(triangle);
    std::cout <<"Bounding sphere success" << std::endl;
}



void testTriangleFunctions() {
    test_get_closest_point();
    test_bounding_sphere();
}

#endif