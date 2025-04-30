#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
// ***************************************************************************

void importMarker(map<unsigned int, list<unsigned int>>& m, int marker, int id){
    if(marker != 0)
    {
        // Check if marker exists in the map
        if (m.find(marker) == m.end()) {
            // If not, insert an empty list for this marker
            m[marker] = list<unsigned int>();
        }

        // Append the id to the list corresponding to the marker
        m[marker].push_back(id);  
    }
}


bool ImportCell0Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell0Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        char delimiter;

        converter >>  id >> delimiter >> marker >> delimiter >> mesh.Cell0DsCoordinates(0, id) >> delimiter >> mesh.Cell0DsCoordinates(1, id);

        mesh.Cell0DsId.push_back(id);

        /// Store the markers
        map<unsigned int, list<unsigned int>>& m = mesh.Cell0DsMarkers;
        importMarker(m, marker, id);
        

    }

    return true;
}
// ***************************************************************************
bool ImportCell1Ds(PolygonalMesh& mesh)
{
    ifstream file("./Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        char delimiter;

        converter >> id >> delimiter >> marker >>  delimiter >> mesh.Cell1DsExtrema(0, id) >> delimiter >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);

        /// Store the markers
        map<unsigned int, list<unsigned int>>& m = mesh.Cell1DsMarkers;
        importMarker(m, marker, id);

        //Check if the edge has a non null length
        int& start = mesh.Cell1DsExtrema(0, id);
        int& end = mesh.Cell1DsExtrema(1, id);

        if(start == end)
        {
            cerr<<id<<" has length equal to zero" << endl;
        }
    }

    return true;
}
// ***************************************************************************
vector<string> SplitCSVLine(const string& line) {
    vector<string> tokens;
    string token;
    istringstream stream(line);
    while (getline(stream, token, ';')) {
        tokens.push_back(token);
    }
    return tokens;
}


bool ImportCell2Ds(PolygonalMesh& mesh)
{
    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);
    
    for (const string& line : listLines)
    {
        
        vector<string> tokens = SplitCSVLine(line);
        size_t idx = 0;

        if (tokens.size() < 6) // id, marker, nVertices, at least 1 vertex, nEdges, at least 1 edge
        {
            std::cerr << "Malformed line in Cell2Ds.csv: " << line << std::endl;
        }

        unsigned int id = stoi(tokens[idx++]);
        unsigned int marker = stoi(tokens[idx++]);

        unsigned int nVertices = stoi(tokens[idx++]);
        vector<unsigned int> vertices(nVertices);

        
        for (unsigned int i = 0; i < nVertices; ++i)
            vertices[i] = stoi(tokens[idx++]);

        unsigned int nEdges = stoi(tokens[idx++]);
        vector<unsigned int> edges(nEdges);
        for (unsigned int i = 0; i < nEdges; ++i)
            edges[i] = stoi(tokens[idx++]);
        
        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
    
        /// Store the markers
        map<unsigned int, list<unsigned int>>& m = mesh.Cell2DsMarkers;
        importMarker(m, marker, id);


        // Ensure non-degenerate polygons (non-zero area)
        vector<unsigned int>& polygon_vertices = mesh.Cell2DsVertices[id];
        const unsigned int vertex_count = polygon_vertices.size();

        double signed_area = 0.0;
        const MatrixXd& coords = mesh.Cell0DsCoordinates;

        for (size_t v = 0; v < vertex_count; ++v)
        {
            unsigned int idx_a = polygon_vertices[v];
            unsigned int idx_b = polygon_vertices[(v + 1) % vertex_count]; // wrap-around to first

            double xa = coords(0, idx_a);
            double ya = coords(1, idx_a);
            double xb = coords(0, idx_b);
            double yb = coords(1, idx_b);
            signed_area += xa * yb - xb * ya;
        }

        signed_area = fabs(0.5 * signed_area);

        if (signed_area < 1e-16)
        {
            std::cerr << "[ERROR] Degenerate polygon detected (ID = " << id << ")" << std::endl;
        }

    }

    return true;
}

}
