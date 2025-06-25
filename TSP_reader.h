//
// Created by lucmi on 21/02/2025.
//

#ifndef TSP_READER_H
#define TSP_READER_H
#include <string>
#include <vector>

class TSP_reader {
public:
    explicit TSP_reader(const std::string& filename);
    std::vector<std::vector<int>> getDistmatrix();
    [[nodiscard]] int getNumNodes() const;
    std::string getName() const;
    std::vector<int> getNode_ids() const;

private:
    int numNodes;
    bool sym;
    std::string name;
    std::vector<std::vector<int>> distmatrix;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    std::vector<int> node_ids;

    void await(std::ifstream &file, std::string &line, const std::string& cmp);
    std::string getArgument(std::string &str);

    std::vector<int> getEdgeWeights(std::ifstream &file, std::string &line);

    bool is_number(const std::string &s);

    std::vector<std::vector<double>> getCoords_with_ids(std::ifstream &file, std::string &line);
};



#endif //TSP_READER_H
