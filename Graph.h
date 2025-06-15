//
// Created by lucmi on 24/02/2025.
//

#ifndef GRAPH_H
#define GRAPH_H
#include <string>
#include <vector>
#include <array>


class Graph {
public:
    Graph(int numNodes, std::vector<std::vector<int>> distmatrix);
    Graph(std::vector<std::vector<double>> coords, std::string norm);
    explicit Graph(int numNodes);

    void random2d(std::string norm, int max);

    std::vector<std::vector<int>> getDistmatrix();
    int getNumNodes() const;
    void setName(std::string name);
    std::string getName() const;
    std::vector<std::vector<int>> get_int_coords();

private:
    int numNodes;
    std::vector<std::vector<int>> distmatrix;
    std::string name = "unnamed";
    std::vector<std::vector<int>> int_coords;

    void calcDistmatrix(std::vector<std::vector<double>> coords, std::string norm);
};



#endif //GRAPH_H;
