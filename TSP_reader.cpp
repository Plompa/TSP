//
// Created by lucmi on 21/02/2025.
//

#include "TSP_reader.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <array>
#include <cmath>

#include "Graph.h"

TSP_reader::TSP_reader(const std::string& filename) {
    std::string line;
    std::ifstream file(filename);

    if (not file) {
        throw std::runtime_error("Cannot open file / path not found");
    }

    //ich lese die datei ein indem ich immer auf ein keyword warte und dann den danebenstehenden parameter einlese

    await(file, line, "NAME");
    name = getArgument(line);

    await(file, line, "TYPE");
    sym = getArgument(line) == "TSP";

    await(file, line, "DIMENSION");
    numNodes = std::stoi(getArgument(line));
    distmatrix = std::vector<std::vector<int>>(numNodes, std::vector<int>(numNodes, 0));

    await(file, line, "EDGE_WEIGHT_TYPE");
    edgeWeightType = getArgument(line);

    if(edgeWeightType == "EXPLICIT") {
        await(file, line, "EDGE_WEIGHT_FORMAT");
        edgeWeightFormat = getArgument(line);

        await(file, line, "EDGE_WEIGHT_SECTION");

        std::vector<int> edgeWeights = getEdgeWeights(file, line); //liest erstmal alle kantengewichte ein

        //die verschiedenen interpretationen der kanten
        if(edgeWeightFormat == "FULL_MATRIX") {
            for (int x = 0; x < numNodes; x++) {
                for (int y = 0; y < numNodes; y++) {
                    distmatrix[x][y] = edgeWeights[x * numNodes + y];
                }
            }
        } else if (edgeWeightFormat == "LOWER_DIAG_ROW") {
            int counter = 0;
            for (int x = 0; x < numNodes; x++) {
                for (int y = 0; y <= x; y++) {
                    distmatrix[x][y] = edgeWeights[counter];
                    distmatrix[y][x] = edgeWeights[counter];
                    counter++;
                }
            }
        } else if (edgeWeightFormat == "UPPER_DIAG_ROW") {
            int counter = 0;
            for (int x = 0; x < numNodes; x++) {
                for (int y = x; y <= numNodes; y++) {
                    distmatrix[x][y] = edgeWeights[counter];
                    distmatrix[y][x] = edgeWeights[counter];
                    counter++;
                }
            }
        } else if (edgeWeightFormat == "UPPER_ROW") {
            int counter = 0;
            for (int x = 0; x < numNodes; x++) {
                for (int y = x; y < numNodes; y++) {
                    distmatrix[x][y] = edgeWeights[counter];
                    distmatrix[y][x] = edgeWeights[counter];
                    counter++;
                }
            }
        }
        node_ids = std::vector<int>(numNodes);
        for (int i = 0; i < numNodes; i++) {
            node_ids[i] = i;
        }
    } else {
        await(file, line, "NODE_COORD_SECTION");

        std::vector<std::vector<double>> coords_with_ids = getCoords_with_ids(file, line);
        std::vector<std::vector<double>> coords(numNodes);
        std::vector<int> node_ids(numNodes);
        for (int i = 0; i < numNodes; ++i) {
            coords[i] = std::vector<double>{coords_with_ids[i][1], coords_with_ids[i][2]};
            node_ids[i] = static_cast<int>(coords_with_ids[i][0]);
        }

        Graph graph(coords, edgeWeightType, node_ids); //für koordinaten besitzt der graph einen eigenen konstruktor
        distmatrix = graph.getDistmatrix();
        this -> node_ids = node_ids;
    }

    file.close();
}

std::vector<std::vector<double>> TSP_reader::getCoords_with_ids(std::ifstream &file, std::string &line) {
    std::vector<std::vector<double>> coords;
    while (getline(file, line)) {
        line = line.substr(line.find_first_not_of(' '), line.size()); //entfernt leerzeichen am anfang
        if(line == "EOF") {
            return coords;
        }
        std::vector<double> coord(3);
        coord[0] = std::stod(line.substr(0, line.find_first_of(' ')));

        line = line.substr(line.find_first_of(' '), line.size()); //entfernt das 1. argument
        line = line.substr(line.find_first_not_of(' '), line.size()); //entfernt leerzeichen
        coord[1] = std::stod(line.substr(0, line.find_first_of(' ')));

        line = line.substr(line.find_first_of(' '), line.size()); //entfernt das 2. argument
        line = line.substr(line.find_first_not_of(' '), line.size()); //entfernt leerzeichen
        coord[2] = std::stod(line.substr(0, line.find_first_of(' ')));

        coords.push_back(coord);
    }
    return coords;
}

std::vector<int> TSP_reader::getEdgeWeights(std::ifstream &file, std::string &line) {
    std::vector<int> edgeWeights;
    while (getline(file, line, ' ')) {
        if(line.empty()) {
            continue;
        }

        //line kann von der form int_a + \n + int_b sein
        std::string a = line;
        std::string b;
        if(line.find('\n') != std::string::npos) {
            a = line.substr(0, line.find_first_of('\n'));
            b = line.substr(line.find_first_of('\n') + 1, line.size());
        }

        if(!is_number(a)) {
            return edgeWeights;
        }
        edgeWeights.push_back(std::stoi(a));
        if(is_number(b)) {
            edgeWeights.push_back(std::stoi(b));
        } else if(!b.empty() && b != "\n") {
            return edgeWeights;
        }
    }
    return edgeWeights;
}

bool TSP_reader::is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

std::string TSP_reader::getArgument(std::string &str) {
    bool firstSpaceFound = false;
    int index = 0;
    for (auto it = str.begin(); it != str.end(); it++) {
        if (firstSpaceFound) {
            if (*it != ' ' and *it != ':') {
                break;
            }
        }
        if (*it == ' ') {
            firstSpaceFound = true;
        }
        index++;
    }
    std::string ret = str.substr(index, str.size());
    std::string::iterator end_pos = std::remove(ret.begin(), ret.end(), ' ');
    ret.erase(end_pos, ret.end());

    end_pos = std::remove(ret.begin(), ret.end(), '\n');
    ret.erase(end_pos, ret.end());

    ret.erase(remove_if(ret.begin(),ret.end(), [] (char c) -> bool {return !(c>=0 && c <128);}), ret.end());
    return ret;
}

void TSP_reader::await(std::ifstream &file, std::string &line, const std::string& cmp) {
    while (getline (file, line)) {
        std::string lineStart = line.substr(0, cmp.size());
        if(lineStart == cmp) {
            return;
        }
    }
    std::cout << "Couldnt find argument: " << cmp << std::endl;
}

std::vector<std::vector<int> > TSP_reader::getDistmatrix() {
    return distmatrix;
}

int TSP_reader::getNumNodes() const {
    return numNodes;
}

std::string TSP_reader::getName() const{
    return name;
}

std::vector<int> TSP_reader::getNode_ids() const {
    return node_ids;
}


