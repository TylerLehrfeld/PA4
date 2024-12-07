#include <iostream>
#include "../modules/PointCloudTransform.h"
#include "../modules/helperFunctions.h"
using std::cout;
using std::endl;

/**
 * @brief test the point cloud registration class using point cloud
 * generator data
 * 
 */
void testPointCloudClasses() {
    cout << "testing generator" << endl;
    PointCloudTransform T = PointCloudTransform();
    
    for(int i = 0; i < 30; i++) {
        vector<Matrix> PcloudA, PcloudB;
        Transform F_BA = generateRandomTransform();
        // create a list of points in 3D space (3x1 matrices);
        // then calculate its corresponding point in 3D space after transformation.
        int numPoints = (int)(randomdouble() * 20 + 10);
        for(int i = 0; i < numPoints; i++) {
            PcloudA.push_back(
                Matrix(3, 1,
                   {randomdouble() * 20 - 10, randomdouble() * 20 - 10,
                    randomdouble() * 20 - 10}));
            PcloudB.push_back(F_BA * PcloudA[i]);
        }
        Transform F_BA_Computed = T.compute(PcloudA, PcloudB);
        Transform IdentitiyTransform =
            F_BA * F_BA_Computed.inverse();
        assert(IdentitiyTransform.R_AB == I);
        assert(IdentitiyTransform.p_AB == origin);
        assert(F_BA * PcloudA[0] == F_BA_Computed * PcloudA[0]);
    }
    cout << "All point clouds are alligned: test point cloud succeeded" << endl;
}
