#include "./modules/Matrix.h"
#include <Eigen/Dense>
#include <vector>
#include "./modules/Transform.h"
#include "./modules/data_reader.h"
#include "./modules/main-helpers.h"
#include <string>
#include <chrono>
using std::vector;
using std::string;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

const string DATA_DIR = "./DATA";    
const int NUM_TESTS_WITH_DEBUG_DATA = 6;


/**
 * @brief The main function. It will do four things:
 *  input data
 *  calculate F_reg, 
 *  use ICP to find the closes point
 *  create output files
 * 
 * @return int 
 */
int main() {
    //input data
    data_reader input(DATA_DIR);

    vector<Matrix> Body_A = input.read_A_point_cloud_data();
    vector<Matrix> Body_B = input.read_B_point_cloud_data();
    vector<Matrix> vertex_coordinates_list = input.read_mesh_vertex_coordinates();
    int** triangle_indeces = input.read_triangle_indeces();
    oct_tree_node* oct_tree = new oct_tree_node(vertex_coordinates_list, triangle_indeces, input.get_num_triangles());
    
    int num_debug_cases = input.get_num_cases();

    // Timer accumulators for each find method
    auto total_slow_time = milliseconds(0);
    auto total_fast_time = milliseconds(0);

    for(int debug_case_index = 0; debug_case_index < num_debug_cases; debug_case_index++) {
        vector<Matrix> c_ks;
        vector<Matrix> d_ks;
        int num_sample_frames = input.get_num_sample_frames(debug_case_index);
        Matrix A_tip = input.get_A_tip();
        for(int sample_frame_index = 0; sample_frame_index < num_sample_frames; sample_frame_index++) {
            vector<Matrix> OPT_Body_A_point_cloud = input.read_Optical_A_point_cloud_body(debug_case_index, sample_frame_index);
            vector<Matrix> OPT_Body_B_point_cloud = input.read_Optical_B_point_cloud_body(debug_case_index, sample_frame_index);
            //F_d is the frame in which the d_i coordinates are measured
            Transform F_d = get_F_d(Body_A, Body_B, OPT_Body_A_point_cloud, OPT_Body_B_point_cloud);
            Matrix d_k = F_d * A_tip;
            d_ks.push_back(d_k);
            //I will be the initial guess of F_reg 
            Matrix s_k = Transform(I, origin) * d_k;
            
            auto slow_start = high_resolution_clock::now();
            closest_point_on_mesh_slow(s_k, vertex_coordinates_list, triangle_indeces, input.get_num_triangles());
            auto slow_end = high_resolution_clock::now();
            total_slow_time += duration_cast<milliseconds>(slow_end - slow_start);

            auto fast_start = high_resolution_clock::now();
            closest_point_on_mesh_fast(s_k, vertex_coordinates_list, oct_tree);
            auto fast_end = high_resolution_clock::now();
            total_fast_time += duration_cast<milliseconds>(fast_end - fast_start);
        }
    }
    std::cout << "total time using slow algorithm in ms: " << total_slow_time.count() << std::endl;
    std::cout << "total time using fast algorithm in ms: " << total_fast_time.count() << std::endl;
    input.free_data(triangle_indeces);
    return 0;
}