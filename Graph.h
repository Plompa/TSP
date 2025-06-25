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
    Graph(int numNodes, std::vector<std::vector<int>> distmatrix, const std::vector<int>& node_ids);
    Graph(std::vector<std::vector<double>> coords, const std::string& norm, const std::vector<int>& node_ids);
    explicit Graph(int numNodes);

    void random2d(const std::string& norm, int max);

    std::vector<std::vector<int>> getDistmatrix();
    [[nodiscard]] int getNumNodes() const;
    [[nodiscard]] std::string getName() const;
    [[nodiscard]] std::vector<std::vector<int>> get_int_coords() const;
    [[nodiscard]] std::vector<int> get_node_ids() const;
    void setName(const std::string &name);


private:
    int numNodes;
    std::vector<std::vector<int>> distmatrix;
    std::string name = "unnamed";
    std::vector<std::vector<int>> int_coords;
    std::vector<int> node_ids;

    void calcDistmatrix(const std::vector<std::vector<double>>& coords, const std::string& norm);
};



#endif //GRAPH_H;
