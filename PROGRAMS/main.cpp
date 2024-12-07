/**
 * @file main.cpp
 * @author Tyler Lehrfeld
 * @brief This file will implement the main structure of the PA4 Program
 * @version 0.1
 * @date 2024-11-12
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "./modules/Matrix.h"
#include <Eigen/Dense>
#include <vector>
#include "./modules/Transform.h"
#include "./modules/data_reader.h"
#include "./modules/main-helpers.h"
#include <string>
using std::vector;
using std::string;
const string DATA_DIR = "./DATA";    
const int NUM_TESTS_WITH_DEBUG_DATA = 6;


/**
 * @brief The main function. It will do four things:
 *  input data
 *  create oct-tree
 *  calculate F_reg using ICP to find the closes point
 *  create output files
 * 
 * @return int 
 */
int main() {
    //input data
    data_reader input(DATA_DIR);

    //pointer point cloud
    vector<Matrix> Body_A = input.read_A_point_cloud_data();
    //rigid body point cloud
    vector<Matrix> Body_B = input.read_B_point_cloud_data();

    //mesh lists
    vector<Matrix> vertex_coordinates_list = input.read_mesh_vertex_coordinates();
    int** triangle_indeces = input.read_triangle_indeces();

    //create oct tree using mesh
    oct_tree_node* oct_tree = new oct_tree_node(vertex_coordinates_list, triangle_indeces, input.get_num_triangles());
    
    int num_debug_cases = input.get_num_cases();
    for(int debug_case_index = 0; debug_case_index < num_debug_cases; debug_case_index++) {
        
        vector<Matrix> d_ks;
        int num_sample_frames = input.get_num_sample_frames(debug_case_index);
        Matrix A_tip = input.get_A_tip();
        for(int sample_frame_index = 0; sample_frame_index < num_sample_frames; sample_frame_index++) {
            //read the point clouds on the pointer and rigid body in order to create a transform
            vector<Matrix> OPT_Body_A_point_cloud = input.read_Optical_A_point_cloud_body(debug_case_index, sample_frame_index);
            vector<Matrix> OPT_Body_B_point_cloud = input.read_Optical_B_point_cloud_body(debug_case_index, sample_frame_index);
            //F_d is the frame in which the d_i coordinates are measured
            Transform F_d = get_F_d(Body_A, Body_B, OPT_Body_A_point_cloud, OPT_Body_B_point_cloud);
            Matrix d_k = F_d * A_tip;
            d_ks.push_back(d_k);
        }

        //get F_reg using ICP
        Transform F_reg = get_F_reg(d_ks, vertex_coordinates_list, triangle_indeces, oct_tree);

        //the results of the transformed reading
        vector<Matrix> s_ks;

        //the closes points on the mesh post transformation
        vector<Matrix> c_ks;
        
        for(Matrix d_k : d_ks) {
            Matrix s_k = F_reg*d_k;
            Matrix c_k = closest_point_on_mesh_fast(s_k, vertex_coordinates_list, oct_tree);
            s_ks.push_back(s_k);
            c_ks.push_back(c_k);
        }

        //create the ouptut file for the debug case
        create_output_file(s_ks, c_ks, input.get_file_name(debug_case_index));
    }

    //outputs error for given debugs compared to my outputs I average the squared error in each file (s_k-c_k) to measure accuracy.
    compare_output_files(input, NUM_TESTS_WITH_DEBUG_DATA, num_debug_cases);
    input.free_data(triangle_indeces);
    oct_tree->free_oct_tree();
    free(oct_tree);
    return 0;
}