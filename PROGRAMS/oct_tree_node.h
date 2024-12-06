#include <stdlib.h>

#include <vector>

#include "Matrix.h"
#include "bounding_sphere.h"
#include "helperFunctions.h"
#include "triangle-functions.h"
using std::malloc;

#ifndef OCTTREE
#define OCTTREE

const double min_diag = .1;
const double min_num_spheres = 20;

/**
 * @brief oct tree of bounding spheres to find closest point on mesh
 *
 */
class oct_tree_node {
private:
    /**
     * @brief finds the closest point on a mesh by recursively searching through
     * subtrees
     *
     * @param point
     * @param bound
     * @param closest
     * @param vertices
     */
    void find_closest(Matrix& point, double& bound, Matrix& closest,
                      vector<Matrix>& vertices) {
        double dist = bound + max_radius;
        if(point.matrixArray[0] > upper_bound_point.matrixArray[0] + dist) {
            return;
        }
        if(point.matrixArray[1] > upper_bound_point.matrixArray[1] + dist) {
            return;
        }
        if(point.matrixArray[2] > upper_bound_point.matrixArray[2] + dist) {
            return;
        }
        if(point.matrixArray[0] < lower_bound_point.matrixArray[0] - dist) {
            return;
        }
        if(point.matrixArray[1] < lower_bound_point.matrixArray[1] - dist) {
            return;
        }
        if(point.matrixArray[2] < lower_bound_point.matrixArray[2] - dist) {
            return;
        }
        if(has_subtrees) {
            for(int i = 0; i < num_subtrees; i++) {
                subtrees[i].find_closest(point, bound, closest, vertices);
            }
        } else {
            for(int i = 0; i < num_spheres; i++) {
                update_closest(spheres[i], point, bound, closest, vertices);
            }
        }
    }
    /**
     * @brief a helper method that updates the closest point found so far
     * 
     * @param bsphere_i 
     * @param point 
     * @param bound 
     * @param closest 
     * @param vertices 
     */
    void update_closest(bounding_sphere* bsphere_i, Matrix& point,
                        double& bound, Matrix& closest,
                        vector<Matrix>& vertices) {
        double dist = (point + -1 * bsphere_i->center).magnitude();
        if(dist - bsphere_i->radius > bound) {
            return;
        }
        Matrix cp = get_closest_point_on_triangle(
            point, vertices[bsphere_i->triangle[0]],
            vertices[bsphere_i->triangle[1]], vertices[bsphere_i->triangle[2]]);
        dist = (cp + -1 * point).magnitude();
        if(dist < bound) {
            bound = dist;
            closest = cp;
        }
    }
    /**
     * @brief Get the bounding info of the spheres. includes lower bound, upper
     * bound and max_radius
     *
     */
    void get_bounding_info() {
        Matrix centroid = origin;
        max_radius = spheres[0]->radius;
        upper_bound_point = spheres[0]->center;
        lower_bound_point = spheres[0]->center;
        for(int i = 0; i < num_spheres; i++) {
            if(spheres[i]->center.matrixArray[0] >
               upper_bound_point.matrixArray[0])
                upper_bound_point.matrixArray[0] =
                    spheres[i]->center.matrixArray[0];
            if(spheres[i]->center.matrixArray[1] >
               upper_bound_point.matrixArray[1])
                upper_bound_point.matrixArray[1] =
                    spheres[i]->center.matrixArray[1];
            if(spheres[i]->center.matrixArray[2] >
               upper_bound_point.matrixArray[2])
                upper_bound_point.matrixArray[2] =
                    spheres[i]->center.matrixArray[2];
            if(spheres[i]->center.matrixArray[0] <
               lower_bound_point.matrixArray[0])
                lower_bound_point.matrixArray[0] =
                    spheres[i]->center.matrixArray[0];
            if(spheres[i]->center.matrixArray[1] <
               lower_bound_point.matrixArray[1])
                lower_bound_point.matrixArray[1] =
                    spheres[i]->center.matrixArray[1];
            if(spheres[i]->center.matrixArray[2] <
               lower_bound_point.matrixArray[2])
                lower_bound_point.matrixArray[2] =
                    spheres[i]->center.matrixArray[2];
            if(spheres[i]->radius > max_radius)
                max_radius = spheres[i]->radius;
            centroid = centroid + spheres[i]->center;
        }
        center = (double)(1 / (double)num_spheres) * centroid;
    }

    /**
     * @brief construct the subtrees of the octotree
     *
     */
    void construct_subtrees() {
        if(num_spheres <= min_num_spheres ||
           (upper_bound_point + -1 * lower_bound_point).magnitude() <=
               min_diag) {
            has_subtrees = false;
            return;
        }
        // vector of length 8; num_spheres_in_each_box[0] = nnn and
        // num_spheres_in_each_box[7] = ppp where i can be represented by binary
        // with n = 0 and p = 1.
        vector<int> num_spheres_in_each_box;
        Splitsort(center, spheres, num_spheres_in_each_box, num_spheres);
        int num_viable_boxes = 0;
        for(int num : num_spheres_in_each_box) {
            if(num > 0) {
                num_viable_boxes++;
            }
        }
        num_subtrees = num_viable_boxes;
        int ind = 0;
        int subtree_ind = 0;
        for(int i = 0; i < 8; i++) {
            if(num_spheres_in_each_box[i] != 0) {
                subtrees.push_back(
                    oct_tree_node(spheres + ind, num_spheres_in_each_box[i]));
                subtree_ind++;
            }
            ind += num_spheres_in_each_box[i];
        }
        assert(subtree_ind == num_subtrees);
        assert(ind == num_spheres);
        has_subtrees = true;
    }
    /**
     * @brief Sorts the spheres by the origin into boxes
     *
     * @param center
     * @param spheres
     * @param num_spheres_in_each_box
     * @param num_spheres
     */
    void Splitsort(Matrix& center, bounding_sphere** spheres,
                   vector<int>& num_spheres_in_each_box, int num_spheres) {
        int ind = 0;
        int cum = 0;
        // for each nnn through ppp
        for(int j = 0; j < 8; j++) {
            for(int i = ind; i < num_spheres; i++) {
                Matrix dif = (spheres[i]->center + -1 * center);
                // if  j = xyz is and x = 1, we get pyz so x coordinate must be
                // positive, so we multiply it by negative one to get it less
                // than 0 and true.
                int x_is_opposite = (j & 4) > 0 ? -1 : 1;
                int y_is_opposite = (j & 2) > 0 ? -1 : 1;
                int z_is_opposite = (j & 1) > 0 ? -1 : 1;

                if(x_is_opposite * dif.matrixArray[0] < 0 &&
                   y_is_opposite * dif.matrixArray[1] < 0 &&
                   z_is_opposite * dif.matrixArray[2] < 0) {
                    if(dif.matrixArray[0] > 0) {
                        assert(j >= 4);
                    } else {
                        assert(j < 4);
                    }
                    if(dif.matrixArray[1] > 0) {
                        assert(j == 2 || j == 3 || j == 6 || j == 7);
                    }
                    if(dif.matrixArray[2] > 0) {
                        assert(j == 1 || j == 3 || j == 5 || j == 7);
                    }

                    swap(ind, i, spheres);
                    ind++;
                }
            }
            num_spheres_in_each_box.push_back(ind - cum);
            cum = ind;
        }
    }

    /**
     * @brief helper to swap in place bounding spheres
     * 
     * @param ind1 
     * @param ind2 
     * @param spheres 
     */
    void swap(int ind1, int ind2, bounding_sphere** spheres) {
        bounding_sphere* temp = spheres[ind1];
        spheres[ind1] = spheres[ind2];
        spheres[ind2] = temp;
    }

public:
    Matrix center;
    Matrix upper_bound_point;
    Matrix lower_bound_point;
    bounding_sphere** spheres;
    vector<oct_tree_node> subtrees;
    int num_spheres;
    bool has_subtrees;
    int num_subtrees;
    double max_radius;

    /**
     * @brief free allocated data
     * 
     */
    void free_oct_tree() {
        for(int i = 0; i < num_spheres; i++) {
            delete(spheres[i]);
        }
        free(spheres);
    }

    /**
     * @brief Construct a new oct tree node object
     * 
     * @param vertex_coordinates 
     * @param triangle_indeces 
     * @param num_triangles 
     */
    oct_tree_node(vector<Matrix>& vertex_coordinates, int** triangle_indeces,
                  int num_triangles) {
        spheres =
            (bounding_sphere**)malloc(sizeof(bounding_sphere*) * num_triangles);
        for(int i = 0; i < num_triangles; i++) {
            spheres[i] =
                new bounding_sphere(triangle_indeces[i], vertex_coordinates);
        }
        init_oct_tree_node(spheres, num_triangles);
    }

    /**
     * @brief Construct a new oct tree node object with bounding spheres
     * 
     * @param _spheres 
     * @param number_of_spheres 
     */
    oct_tree_node(bounding_sphere** _spheres, int number_of_spheres) {
        spheres = _spheres;
        num_spheres = number_of_spheres;
        get_bounding_info();
        construct_subtrees();
    }
    /**
     * @brief Construct a new oct tree node object with bounding spheres
     * 
     * @param _spheres 
     * @param number_of_spheres 
     */
    void init_oct_tree_node(bounding_sphere** _spheres, int number_of_spheres) {
        spheres = _spheres;
        num_spheres = number_of_spheres;
        get_bounding_info();
        construct_subtrees();
    }
    /**
     * @brief find the closest point on a mesh
     * 
     * @param point 
     * @param vertices 
     * @return Matrix 
     */
    Matrix find_closest(Matrix point, vector<Matrix>& vertices) {
        Matrix closest = get_closest_point_on_triangle(
            point, vertices[spheres[0]->triangle[0]],
            vertices[spheres[0]->triangle[1]],
            vertices[spheres[0]->triangle[2]]);
        double bound = (point + -1 * closest).magnitude();
        find_closest(point, bound, closest, vertices);
        return closest;
    }
};

#endif