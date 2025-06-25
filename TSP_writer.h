//
// Created by lucmi on 25/04/2025.
//

#ifndef TSP_WRITER_H
#define TSP_WRITER_H
#include <map>
#include <string>
#include <vector>


class TSP_writer {
public:
    TSP_writer();
    void save(const std::vector<int> &tour, int tour_len, const std::string &name, const std::vector<int> &node_ids);
    void save(const std::vector<int> &tour, int tour_len, const std::string &name, const std::vector<int> &node_ids, const std::map<
              std::string, std::vector<int>> &display_data);
};



#endif //TSP_WRITER_H
