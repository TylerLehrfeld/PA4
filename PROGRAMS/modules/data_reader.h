#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include "Matrix.h"
#include "helperFunctions.h"
#include <iostream>
using std::string;
#ifndef DATA_READER
#define DATA_READER

/**
 * @brief Class to read and store data from the data directory
 *
 */
class data_reader {
public:
    /**
     * @brief Construct a new data reader object
     * 
     * @param data_dir 
     */
    data_reader(string data_dir) {  
        for(auto const& file : std::filesystem::directory_iterator(data_dir)) {
            string filename = file.path().generic_string();
            if(filename.find("SampleReadingsTest.txt") != -1)
                filenames.push_back(filename);
        }
        A_Body_IF_Stream = new std::ifstream("./DATA/Problem4-BodyA.txt");
        B_Body_IF_Stream = new std::ifstream("./DATA/Problem4-BodyB.txt");
        mesh_IF_Stream = new std::ifstream("./DATA/Problem4MeshFile.sur");
        A_tip = Matrix(3, 1, {0, 0, 0});
        for(string filename: filenames) {
            sample_readings_streams.push_back(new std::ifstream(filename)); 
        }
    }
    /**
     * @brief return the pcloud data for the A pointer
     * 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_A_point_cloud_data() {
        num_A_points = read_int(A_Body_IF_Stream);
        return read_pcloud_data(A_Body_IF_Stream, true, num_A_points);
    };

    /**
     * @brief return the pcloud data for the B rigid body
     * 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_B_point_cloud_data() {
        num_B_points = read_int(B_Body_IF_Stream);
        return read_pcloud_data(B_Body_IF_Stream, true, num_B_points);
    };

    /**
     * @brief returns a point cloud of vertices on the mesh
     * 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_mesh_vertex_coordinates() {
        int num_points = read_int(mesh_IF_Stream);
        return read_pcloud_data(mesh_IF_Stream, false, num_points);
    };

    /**
     * @brief returns the indeces of vertices for each triangle on the mesh
     * 
     * @return int** 
     */
    int** read_triangle_indeces() {
        num_triangles = read_int(mesh_IF_Stream);
        int** triangle_vertex_indeces = (int**)malloc(sizeof(int*) * num_triangles);
        for(int i = 0; i < num_triangles; i++) {
            triangle_vertex_indeces[i] = (int*)malloc(sizeof(int) * 3);
            triangle_vertex_indeces[i][0] = read_int(mesh_IF_Stream);
            triangle_vertex_indeces[i][1] = read_int(mesh_IF_Stream);
            triangle_vertex_indeces[i][2] = read_int(mesh_IF_Stream);
            //read -1s
            read_double(mesh_IF_Stream);
            read_double(mesh_IF_Stream);
            read_double(mesh_IF_Stream);
        }
        return triangle_vertex_indeces;
    };

    /**
     * @brief Get the number of triangles
     * 
     * @return int 
     */
    int get_num_triangles() {
        return num_triangles;
    }

    /**
     * @brief Get the number of test cases for the program to run
     * 
     * @return int 
     */
    int get_num_cases() {
        return filenames.size();
    }

    /**
     * @brief Get the number of sample frames for a debug case
     * 
     * @param debug_case_index 
     * @return int 
     */
    int get_num_sample_frames(int debug_case_index) {
        N_S_list.push_back(read_int(sample_readings_streams[debug_case_index]));
        int num_sample_frames = read_int(sample_readings_streams[debug_case_index]);
        string filename = read_string(sample_readings_streams[debug_case_index]);
        int dummynum = read_int(sample_readings_streams[debug_case_index]);
        return num_sample_frames;
    };

    /**
     * @brief Get A_tip from file
     * 
     * @return Matrix 
     */
    Matrix get_A_tip() {
        if(!(A_tip == origin)) {
            return A_tip;
        }
        double x = read_double(A_Body_IF_Stream);
        double y = read_double(A_Body_IF_Stream);
        double z = read_double(A_Body_IF_Stream);
        A_tip = Matrix(3, 1, {x, y, z});
        return A_tip;
    };
    /**
     * @brief read an optical frame of the pointer
     * 
     * @param debug_case_index 
     * @param sample_frame_index 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_Optical_A_point_cloud_body(int debug_case_index,
                                                   int sample_frame_index) {
        return read_pcloud_data(sample_readings_streams[debug_case_index], false, num_A_points);
    };
    /**
     * @brief read an optical frame of the reference body
     * 
     * @param debug_case_index 
     * @param sample_frame_index 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_Optical_B_point_cloud_body(int debug_case_index,
                                                   int sample_frame_index) {
        vector<Matrix> ret_cloud = read_pcloud_data(sample_readings_streams[debug_case_index], false, num_B_points);
        //read dummy points
        read_pcloud_data(sample_readings_streams[debug_case_index], false, N_S_list[debug_case_index] - num_A_points - num_B_points);
        return ret_cloud;
    };

    /**
     * @brief free allocated data
     * 
     * @param triangle_vertext_indeces 
     */
    void free_data(int** triangle_vertext_indeces) {
        A_Body_IF_Stream->close();
        delete(A_Body_IF_Stream);
        B_Body_IF_Stream->close();
        delete(B_Body_IF_Stream);
        mesh_IF_Stream->close();
        delete(mesh_IF_Stream);
        for(std::ifstream* stream: sample_readings_streams) {
            stream->close();
            delete(stream);
        }

        for(int i = 0; i < num_triangles; i++) {
            free(triangle_vertext_indeces[i]);
        }
        free(triangle_vertext_indeces);
    }

    /**
     * @brief Get the file name for a debug case
     * 
     * @param debug_case_index 
     * @return string 
     */
    string get_file_name(int debug_case_index) {
        return filenames[debug_case_index];
    }

    /**
     * @brief read an output file and return its points.
     * 
     * @param filename 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_output_file(string filename) {
        std::ifstream* output = new std::ifstream(filename);
        int num_points = read_int(output);
        filename = read_string(output);
        int dummy_int = read_int(output);
        vector<Matrix> pts;
        double x, y, z, err;
        for(int i = 0; i < num_points * 2; i++) {
            x = read_double(output);
            y = read_double(output);
            z = read_double(output);
            pts.push_back(Matrix(3, 1, {x, y, z}));
            if(i % 2 == 1)
                err = read_double(output);
        }
        output->close();
        delete(output);
        return pts;
    }
private:
    vector<string> filenames;
    Matrix A_tip;
    std::ifstream* A_Body_IF_Stream;
    std::ifstream* B_Body_IF_Stream;
    std::ifstream* mesh_IF_Stream;
    vector<std::ifstream*> sample_readings_streams;
    vector<int> N_S_list;
    int num_A_points;
    int num_B_points;
    int num_triangles;

    /**
     * @brief get the next string of a filestream
     *
     * @param fileStream
     * @return string
     */
    string read_string(std::ifstream* fileStream) {
        string nextStr = "";
        *fileStream >> nextStr;
        return nextStr;
    }

    /**
     * @brief read an int from a filestream
     *
     * @param fileStream
     * @return int
     */
    int read_int(std::ifstream* fileStream) {
        string nextStr = "";
        *fileStream >> nextStr;
        return std::stoi(nextStr);
    }

    /**
     * @brief read a double from a filestream
     *
     * @param fileStream
     * @return double
     */
    double read_double(std::ifstream* fileStream) {
        string nextStr = "";
        *fileStream >> nextStr;
        return std::stof(nextStr);
    }

    /**
     * @brief generic pcloud reading function
     * 
     * @param BodyFile 
     * @param has_filename 
     * @param num_points 
     * @return vector<Matrix> 
     */
    vector<Matrix> read_pcloud_data(std::ifstream* BodyFile, bool has_filename, int num_points) {
        vector<Matrix> pcloud;
        if(has_filename)
            string filename = read_string(BodyFile);
        for(int i = 0; i < num_points; i++) {
            double x = read_double(BodyFile);
            double y = read_double(BodyFile);
            double z = read_double(BodyFile);
            pcloud.push_back(Matrix(3,1,{x,y,z}));
        }
        return pcloud;
    };
    
};

#endif