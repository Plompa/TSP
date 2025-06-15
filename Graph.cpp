//
// Created by lucmi on 24/02/2025.
//

#include "Graph.h"

#include <array>
#include <string>
#include <cmath>
#include <random>

Graph::Graph(int numNodes) : numNodes(numNodes) {
    distmatrix = std::vector< std::vector<int> >(numNodes, std::vector<int>(numNodes));
}

Graph::Graph(int numNodes, std::vector<std::vector<int>> distmatrix) : numNodes(numNodes), distmatrix(distmatrix) {
}

Graph::Graph(std::vector<std::array<double, 2>> coords, std::string norm) {
    numNodes = coords.size();
    distmatrix = std::vector<std::vector<int>>(numNodes, std::vector<int>(numNodes));
    for (auto & coord : coords) {
        int x = static_cast<int>(coord[0]);
        int y = static_cast<int>(coord[1]);
        int_coords.push_back(std::array<int, 2>{x, y});
    }

    calcDistmatrix(coords, norm);
}

void Graph::random2d(std::string norm, int max) {
    //std::random_device rd; // obtain a random number from hardware
    //std::mt19937 gen(rd()); // seed the generator
    std::mt19937 gen(1); // seed the generator
    std::uniform_int_distribution<> distr(0, max); // define the range

    std::vector<std::array<double, 2>> coords;

    for(int i = 0; i < Graph::numNodes; i++) {
        std::array<double, 2> coord{};
        coord[0] = distr(gen);
        coord[1] = distr(gen);
        coords.push_back(coord);
    }
    for (auto & coord : coords) {
        int x = static_cast<int>(coord[0]);
        int y = static_cast<int>(coord[1]);
        int_coords.push_back(std::array<int, 2>{x, y});
    }

    calcDistmatrix(coords, norm);
}

void Graph::calcDistmatrix(std::vector<std::array<double, 2>> coords, std::string norm) {
    //so wie ich das verstanden habe sollen wir immer ganzzahlen verwenden, daher runde ich einfach
    if(norm == "EUC_2D") {
        for (int x = 0; x < numNodes; x++) {
            for (int y = 0; y < numNodes; y++) {
                double dx = coords[x][0] - coords[y][0];
                double dy = coords[x][1] - coords[y][1];
                int dist = static_cast<int>(round(sqrt(dx * dx + dy * dy)));
                distmatrix[x][y] = dist;
            }
        }
    } else { //ceil
        for (int x = 0; x < numNodes; x++) {
            for (int y = 0; y < numNodes; y++) {
                double dx = coords[x][0] - coords[y][0];
                double dy = coords[x][1] - coords[y][1];
                int dist = static_cast<int>(ceil(sqrt(dx * dx + dy * dy)));
                distmatrix[x][y] = dist;
            }
        }
    }
}

std::vector<std::vector<int> > Graph::getDistmatrix() {
    return distmatrix;
}

int Graph::getNumNodes() const {
    return numNodes;
}

void Graph::setName(std::string name) {
    this -> name = name;
}

std::string Graph::getName() const {
    return name;
}

std::vector<std::array<int, 2>> Graph::get_int_coords() {
    return int_coords;
}
