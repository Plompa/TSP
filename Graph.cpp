//
// Created by lucmi on 24/02/2025.
//

#include "Graph.h"

#include <array>
#include <string>
#include <cmath>
#include <random>
#include <utility>

Graph::Graph(const int numNodes) : numNodes(numNodes) {
    distmatrix = std::vector< std::vector<int> >(numNodes, std::vector<int>(numNodes));
    node_ids = std::vector<int>(numNodes);
    for (int i = 0; i < numNodes; i++) {
        node_ids[i] = i;
    }
}

Graph::Graph(const int numNodes, std::vector<std::vector<int>> distmatrix, const std::vector<int> &node_ids) :
    numNodes(numNodes),
    distmatrix(std::move(distmatrix)) {
    this -> node_ids = node_ids;
}

Graph::Graph(std::vector<std::vector<double>> coords, const std::string& norm, const std::vector<int> &node_ids) {
    this -> node_ids = node_ids;
    numNodes = coords.size();
    distmatrix = std::vector<std::vector<int>>(numNodes, std::vector<int>(numNodes));
    for (auto & coord : coords) {
        int x = static_cast<int>(coord[0]);
        int y = static_cast<int>(coord[1]);
        int_coords.push_back(std::vector<int>{x, y});
    }

    calcDistmatrix(coords, norm);
}

void Graph::random2d(const std::string& norm, int max) {
    //std::random_device rd; // obtain a random number from hardware
    //std::mt19937 gen(rd()); // seed the generator
    std::mt19937 gen(1); // seed the generator //1
    std::uniform_int_distribution<> distr(0, max); // define the range

    std::vector<std::vector<double>> coords;

    for(int i = 0; i < Graph::numNodes; i++) {
        std::vector<double> coord(2);
        coord[0] = distr(gen);
        coord[1] = distr(gen);
        coords.push_back(coord);
    }
    for (auto & coord : coords) {
        int x = static_cast<int>(coord[0]);
        int y = static_cast<int>(coord[1]);
        int_coords.push_back(std::vector<int>{x, y});
    }

    calcDistmatrix(coords, norm);
}

void Graph::calcDistmatrix(const std::vector<std::vector<double>>& coords, const std::string& norm) {
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

std::vector<std::vector<int>> Graph::getDistmatrix() {
    return distmatrix;
}

int Graph::getNumNodes() const {
    return numNodes;
}

void Graph::setName(const std::string &name) {
    this -> name = name;
}

std::string Graph::getName() const {
    return name;
}

std::vector<std::vector<int>> Graph::get_int_coords() const {
    return int_coords;
}

std::vector<int> Graph::get_node_ids() const {
    return node_ids;
}

