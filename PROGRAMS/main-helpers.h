#include <vector>
#include <filesystem>
#include <fstream>
#include "Matrix.h"
#include "Transform.h"
#include "data_reader.h"
#include "helperFunctions.h"
#include "PointCloudTransform.h"
#include "triangle-functions.h"
#include "oct_tree_node.h"

#ifndef MAIN_HELPER
#define MAIN_HELPER

float D_MAX = .5;

/**
 * @brief Get the F_d by combining F_A and F_B
 * 
 * @param Body_A 
 * @param Body_B 
 * @param OPT_Body_A_point_cloud 
 * @param OPT_Body_B_point_cloud 
 * @return Transform 
 */
Transform get_F_d(vector<Matrix> Body_A, vector<Matrix> Body_B, vector<Matrix> OPT_Body_A_point_cloud, vector<Matrix> OPT_Body_B_point_cloud) {
    PointCloudTransform T;
    Transform F_OA = T.compute(Body_A, OPT_Body_A_point_cloud);
    Transform F_OB = T.compute(Body_B, OPT_Body_B_point_cloud);
    return F_OB.inverse() * F_OA;
};


/**
 * @brief use oct_tree to find closest point
 * 
 * @param s_k 
 * @param vertex_coordinates 
 * @param triangle_indeces 
 * @param num_triangles 
 * @param oct_tree 
 * @return Matrix 
 */
Matrix closest_point_on_mesh_fast(Matrix s_k, vector<Matrix>& vertex_coordinates, int** triangle_indeces, int num_triangles, oct_tree_node* oct_tree) {
    return oct_tree->find_closest(s_k, vertex_coordinates);
}


/**
 * @brief Get the F_reg transform
 * 
 * @param d_ks 
 * @param input 
 * @param vertex_coordinates_list 
 * @param triangle_indeces 
 * @param oct_tree 
 * @return Transform 
 */
Transform get_F_reg(vector<Matrix> d_ks, data_reader input, vector<Matrix>& vertex_coordinates_list, int** triangle_indeces, oct_tree_node* oct_tree) {
    //initial guess
    Transform T(I, origin);
    PointCloudTransform Transform_computer;
    for(int i = 0; i < 25; i++) {
        vector<Matrix> m_is;
        vector<Matrix> matched_d_ks;
        for(Matrix d_k : d_ks) {
            Matrix m_i = closest_point_on_mesh_fast(T*d_k, vertex_coordinates_list, triangle_indeces, input.get_num_triangles(), oct_tree);
            m_is.push_back(m_i);
            //matched_d_ks.push_back(d_k);
        }    
        //D_MAX *= .95;
        T = Transform_computer.compute(m_is, d_ks);
    }
    return T;
};

/**
 * @brief use linear search to find closest point
 * 
 * @param s_k 
 * @param vertex_coordinates 
 * @param triangle_indeces 
 * @param num_triangles 
 * @return Matrix 
 */
Matrix closest_point_on_mesh_slow(Matrix s_k, vector<Matrix> vertex_coordinates, int** triangle_indeces, int num_triangles) {
    Matrix closest_point = get_closest_point_on_triangle(s_k, vertex_coordinates[triangle_indeces[0][0]], vertex_coordinates[triangle_indeces[0][1]], vertex_coordinates[triangle_indeces[0][2]]);
    double smallest_dist_mag = (s_k + -1 * closest_point).magnitude();
    double smallest_dist_ind = 0;
    for(int i = 1; i < num_triangles; i++) {
        Matrix new_closest_point = get_closest_point_on_triangle(
            s_k, vertex_coordinates[triangle_indeces[i][0]],
            vertex_coordinates[triangle_indeces[i][1]],
            vertex_coordinates[triangle_indeces[i][2]]);
        double mag = (s_k + -1 *  new_closest_point).magnitude();
        if(mag < smallest_dist_mag) {
            smallest_dist_mag = mag;
            smallest_dist_ind = i;
            closest_point = new_closest_point;
        }
    }
    return closest_point;
}

/**
 * @brief Create an output file 
 * 
 * @param d_ks 
 * @param c_ks 
 * @param FILE_NAME 
 */
void create_output_file(vector<Matrix> d_ks, vector<Matrix> c_ks, string FILE_NAME) {
    FILE_NAME = FILE_NAME.substr(FILE_NAME.find("PA4"), FILE_NAME.find("Sample") - FILE_NAME.find("PA4"));
    FILE_NAME = "../OUTPUT/"+FILE_NAME+"Output.txt";
    std::ofstream outfile(FILE_NAME);
    if(!outfile.is_open()) {
        throw std::logic_error("The file wasn't opened");
    }
    int num_points = d_ks.size();

    outfile << num_points << " " << FILE_NAME << " " << 0 << std::endl;
    for(int i = 0; i < num_points; i++) {
        Matrix diff = c_ks[i] + (-1 * d_ks[i]);
      outfile << " " << d_ks[i].matrixArray[0] << " " << d_ks[i].matrixArray[1] << " "
              << d_ks[i].matrixArray[2] << "\t" << c_ks[i].matrixArray[0] << " "
              << c_ks[i].matrixArray[1] << " " << c_ks[i].matrixArray[2] << " "
              << diff.magnitude() << std::endl;
    }
};

/**
 * @brief compare output files with my results
 * 
 * @param input 
 * @param num_debug_tests 
 */
void compare_output_files(data_reader input, int num_debug_tests) {
    for(int i = 0; i < num_debug_tests; i++) {
        string filename = input.get_file_name(i);
        std::ifstream debug_file(filename);
        string new_filename = filename.substr(filename.find("PA4"), filename.find("Sample") - filename.find("PA4"));
        filename = "./DATA/"+new_filename+"Output.txt";
        new_filename = "../OUTPUT/"+new_filename+"Output.txt";
        std::ifstream outfile(new_filename);
        vector<Matrix> outputed_points = input.read_output_file(new_filename);
        vector<Matrix> debug_points = input.read_output_file(filename);
        double mean_squared_error_on_difference = 0;
        double msq_on_debug = 0;
        for(int i = 0; i < outputed_points.size(); i+=2) {
            mean_squared_error_on_difference+=pow((outputed_points[i] + -1*outputed_points[i+1]).magnitude(),2);
            msq_on_debug+=pow((debug_points[i] + -1*debug_points[i+1]).magnitude(),2);
        }
        std::cout << "mean squared error on my outputs: " << mean_squared_error_on_difference/(outputed_points.size()/2) << std::endl;
        std::cout << "mean squared error on given outputs: " << msq_on_debug/(debug_points.size()/2) << std::endl;
    }
}
        
#endif