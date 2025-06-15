//
// Created by lucmi on 25/04/2025.
//
#include <iostream>
#include <fstream>

#include "TSP_writer.h"

#include <bits/fs_fwd.h>
#include <bits/fs_path.h>

TSP_writer::TSP_writer() {
}

void TSP_writer::save(const std::vector<int> &tour, int tour_len, std::string name) {
    //hier passiert nichts spanndendes
    std::ofstream file(name + ".opt.tour");

    file << "NAME: " << name << std::endl;
    file << "COMMENT: Optimal Solution for " << name << "(" << tour_len << ")" << std::endl;
    file << "TYPE: TOUR" << std::endl;
    file << "DIMENSION: " << tour.size() << std::endl;
    file << "TOUR SECTION" << std::endl;

    for (int i = 0; i < tour.size(); i++) {
        file << tour[i] << " ";
    }
    file << std::endl;

    file << "EOF" << std::endl;
    file.close();
}

void TSP_writer::save(const std::vector<int>& tour, int tour_len, const std::map<std::string, std::vector<int>> &display_data, std::string name) {
    //hier passiert nichts spanndendes
    std::ofstream file(name + ".opt.tour");

    file << "NAME: " << name << std::endl;
    file << "COMMENT: Optimal Solution for " << name << "(" << tour_len << ")" << std::endl;
    file << "TYPE: TOUR" << std::endl;
    file << "DIMENSION: " << tour.size() << std::endl;
    file << "TOUR SECTION" << std::endl;

    for (int i = 0; i < tour.size(); i++) {
        file << tour[i] << " ";
    }
    file << std::endl;

    file << "DISPLAY_DATA SECTION" << std::endl;
    for (auto const&[descriptor, vector]: display_data) {
        file << descriptor << " ";
        for (auto val: vector) {
            file << val << " ";
        }
        file << std::endl;
    }

    file << "EOF" << std::endl;
    file.close();

    //std::cout << std::filesystem::current_path() << std::endl;
}
