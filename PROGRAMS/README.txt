bounding_sphere.h: defines the class for bounding spheres
data_reader.h: defines the class for reading and writing data
helperFunctions.h: suite of miscelaneous helper functions
Makefile: the file that allows compilation with make and make test
main.cpp: the main file that creates the appropriate output files
Matrix.cpp: the file that implements the matrix class
Matrix.h: the file that defines the matrix class
oct_tree_node.h: the file that defines and implements an oct_tree for searching for a closest point
performance.cpp: the testing file that runs the brute force vs more efficient closest point search
PointCloudTransform.h: the file that defines the point cloud registration class
Transform.cpp: the file that implements the transform class
Transform.h: the file that defines the transform class
triangle-functions.h: this file houses funcitions that the closest point searcher uses related to triangles
main-helpers.h: the file that contains the helper files for the main functions
data-reader.h: the file that defines the data reader class
instructions.txt: instructions to make and run the executables

Testing files:
test.cpp: the testing file that runs all of the unit tests
Matrix-test.h: tests matrix functions
PointCloudTest.h: definines test for point cloud registration
Transform-test.h: defines tests for Frame transformations
triangle-testsh.h: 

Executalbles
main: executable
test: executable (with parameters, see instructions)
performance: executable that tests the two methods against each other

Eigen/Dense: in dependencies, there is a folder called eigen-3.4.0. It contains a library that I use for SVD. Here is the link: https://eigen.tuxfamily.org/dox/GettingStarted.html
