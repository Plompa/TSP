//
// Created by lucmi on 22/02/2025.
//

#ifndef TSP_ENUMERATION_H
#define TSP_ENUMERATION_H
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "Graph.h"


class TSP_enumeration {
public:
    explicit TSP_enumeration(Graph graph);
private:
    const int MAXCOST = 100000000;
    std::vector<std::vector<int>> distmatrix;
    std::vector<int> tour, besttour;
    int bestlen = 1000000000;
    int ncount;
    std::unordered_map<std::vector<bool>, int> mst_map;
    std::unordered_map<std::vector<bool>, int> shortest_path_map;
    int map_uses = 0;
    int shortest_path_map_uses = 0;
    std::vector<int> delta;
    bool tourFound = false;
    std::vector<int> k_neighbour_max_dist;
    std::map<std::string, std::vector<int>> display_data;

    std::vector<int> optimal_tour = {33, 10, 16, 17, 11, 23, 34, 13, 32, 29, 20, 21, 25, 24, 1, 15, 9, 8, 0, 3, 18, 26, 19, 31, 7, 30, 14, 5, 12, 6, 4, 27, 2, 28, 22, 35};

    int dist(int i, int j) const;

    int unmodified_dist(int i, int j) const;

    int tour_length() const;

    int subtour_length(int s, int t) const;

    int best_tour_length() const;

    void tour_swap(int i, int j);

    void print_path_map() const;

    void permute(int k, int tourlen, bool has_zero);

    int mst(int count);

    void k_neigbour_init(int k);

    bool dynamic_shortest_path(int s, int t, int tourlen);

    int shortest_path(int k);

    int path_permute(int k, int t, int pathlen, int best_pathlen);

    int hk_bound();

    int one_tree_bound();

    std::vector<int> max_one_tree();

    void held_karp_tuning();

    void held_karp_tuning_2();

    int w_delta(const std::vector<int>& one_tree) const;

    int w_delta_2(const std::vector<int> &deg) const;

    double step(int m, int M, double step1);

    std::vector<int> held_karp_one_tree(int start);

    std::vector<int> held_karp_mst(int count) const;

    void twoOptInit();

    int two_opt(int cur_len);

    int nntour(int start);

    bool is_sym() const;
};



#endif //TSP_ENUMERATION_H
