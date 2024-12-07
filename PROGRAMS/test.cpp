#include <iostream>
#include <string>
#include <cmath>
#include "./test_files/Transform-test.h"
#include "./test_files/PointCloudTest.h"
#include "./test_files/Matrix-test.h"
#include "./test_files/triangle-tests.h"
#include "./test_files/ICP-test.h"

using std::cout;
using std::endl;
using std::string;

/**
 * @brief runs all of tests in which the appropriate args are given
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, char *argv[]) {
    for(int i = 0; i < argc; i++) {
        string argument = argv[i];
        if(argument == "matrix") {
            testMatrixClass();
        }
        if(argument == "transform") {
            testTransformClass();
        }
        if(argument == "point-cloud") {
            testPointCloudClasses();
        }
        if(argument == "triangles") {
            testTriangleFunctions();
        }
        if(argument == "ICP") {
            testICP();
        }
        
    }    
}