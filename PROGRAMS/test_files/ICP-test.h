#include <vector>
#include "../modules/Matrix.h"
#include "../modules/helperFunctions.h"
#include "../modules/oct_tree_node.h"
#include "../modules/main-helpers.h"
#include <stdlib.h>
#include <iostream>


/**
 * @brief generate a surface to test ICP on
 * 
 * @param points 
 * @param triangles 
 * @param meshcolumns 
 * @param meshrows 
 */
void generateSurface(vector<Matrix>& points, int** triangles, int meshcolumns, int meshrows) {
    int x = 0; int y = 0;
    for(int x = 0; x < meshcolumns; x++) {
        for(int y = 0; y < meshrows; y ++) {
            double dx = 2*(double)rand() / RAND_MAX - 1;
            double dy = 2*(double)rand() / RAND_MAX - 1;
            double dz = 2*(double)rand() / RAND_MAX - 1;
            //introduce randomness into mesh
            points.push_back(Matrix(3,1,{5*x+dx,5*y+dy,dz}));
        }
    }
    int c = 0;
    for(int i = 0; i < meshcolumns - 1; i++) {
        for(int j = 0; j < meshrows - 1; j++) {
            triangles[c][0] = meshcolumns*j + i;
            triangles[c][1] = meshcolumns*j + i + meshcolumns;
            triangles[c][2] = meshcolumns*j + i + meshcolumns + 1;
            triangles[c + 1][0] = meshcolumns*j + i;
            triangles[c + 1][1] = meshcolumns*j + i + 1;
            triangles[c + 1][2] = meshcolumns*j + i + meshcolumns + 1;
            c+=2;
        }
    }
}

/**
 * @brief test the ICP function implemented in get_F_reg
 * 
 */
void testICP() {
    
    int meshrows = 10;
    int meshcolumns = 5;
    int numTriangles = 2 * (meshcolumns - 1) * (meshrows - 1);
    vector<Matrix> points;
    int** triangles =(int**)malloc(numTriangles*sizeof(int*));
    for(int i = 0 ; i < numTriangles; i++) {
        triangles[i] = (int*)malloc(3*sizeof(int));
    }
    generateSurface(points, triangles, meshcolumns, meshrows);
    oct_tree_node* oct_tree = new oct_tree_node(points, triangles, numTriangles);
    
    double dx = 3*(double)rand() / RAND_MAX - 1;
    double dy = 3*(double)rand() / RAND_MAX - 1;
    double dz = 3*(double)rand() / RAND_MAX - 1;
            
    Transform randomSmallT(generateSmallRotation(), Matrix(3,1,{dx,dy,dz}));
    vector<Matrix> newPoints;
    for(int i = 0; i < points.size(); i++) {
        newPoints.push_back(randomSmallT*points[i]);
    }
    Transform Freg = get_F_reg(newPoints, points, triangles, oct_tree);
    randomSmallT = randomSmallT.inverse();
    assert((Freg.R_AB + -1 * randomSmallT.R_AB).overallsize() < 0.0001 && (Freg.p_AB + -1 * randomSmallT.p_AB).overallsize() < 0.005);
    std::cout << "successful testICP"  << std::endl;
    for(int i = 0 ; i < numTriangles; i++) {
        free(triangles[i]);
    }
    
    free(triangles);


}