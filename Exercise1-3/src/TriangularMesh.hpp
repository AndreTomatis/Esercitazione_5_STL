#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace TriangularLibrary {

struct TriangularMesh
{
    unsigned int NumCell0Ds;
    unsigned int NumCell1Ds;
    unsigned int NumCell2Ds;

    vector<unsigned int> Cell0DsId;
    vector<unsigned int> Cell1DsId;
    vector<unsigned int> Cell2DsId;

    Eigen::MatrixXd Cell0DsCoordinates;
    Eigen::MatrixXi Cell1DsExtrema;
    vector<array<unsigned int, 3>> Cell2DsVertices;
    vector<array<unsigned int, 3>> Cell2DsEdges;
};

}
