#include <fstream>
#include <iostream>

#include "Einstiegsaufgabe.h"
#include "Graph.h"
#include "TSP_enumeration.h"
#include "TSP_reader.h"

int main(int argc, char * argv[]) {
    std::string nrw = "C:\\Users\\lucmi\\Downloads\\TSPLIB_ alle symmetrischen TSP Instanzen\\nrw1379.tsp\\nrw1379.tsp";
    std::string dantzig = "C:\\Users\\lucmi\\Downloads\\TSPLIB_ alle symmetrischen TSP Instanzen\\dantzig42.tsp\\dantzig42.tsp";
    std::string groetschel = "C:\\Users\\lucmi\\Downloads\\TSPLIB_ alle symmetrischen TSP Instanzen\\gr120.tsp\\gr120.tsp";
    std::string other = "C:\\Users\\lucmi\\Downloads\\TSPLIB_ alle symmetrischen TSP Instanzen\\eil76.tsp\\eil76.tsp";

    if(argc < 2) {
        std::cout << "please enter arguments" << std::endl;
        return 1;
    }

    Graph graph(0);

    std::string arg1(argv[1]);

    if(arg1 ==  "-path") {
        if(argc < 3) {
            std::cout << "please enter path" << std::endl;
            return 1;
        }
        std::string arg2(argv[2]);
        TSP_reader tsp_read(arg2);
        graph = Graph(tsp_read.getNumNodes(), tsp_read.getDistmatrix(), tsp_read.getNode_ids());
        graph.setName(tsp_read.getName());
    }
    else if(arg1 ==  "-rand") {
        //seed in graph

        int ncount = 0;

        if(argc >= 3) {
            ncount = atoi(argv[2]);
        } else {
            std::cout << "please enter city count" << std::endl;
            return 1;
        }

        graph =  Graph(ncount);
        graph.random2d("EUC_2D", 1000);
    } else {
        std::cout << "invalid arguments" << std::endl;
        return 1;
    }

    if(argc >= 4) {
        std::string arg3(argv[3]);
        if(arg3 == "-out") {
            if(argc >= 5) {
                graph.setName(argv[4]);
            } else {
                std::cout << "please enter output file name" << std::endl;
            }
        }
    }

    std::cout << "Number of nodes: " << graph.getNumNodes() << std::endl;
    std::cout << std::endl;

    std::cout << "Einstiegsaufgabe: " << std::endl;
    Einstiegsaufgabe aufg(graph);
    std::cout << std::endl;

    std::cout << "Enumeration algorithm: " << std::endl;
    TSP_enumeration tsp(graph);

    /*for (int i = 100; i < 1000; ++i) {
        std::cout << i << " cities" << std::endl;
        Graph graph2(i);
        graph2.random2d("EUC_2D", 1000);
        TSP_enumeration tsp2(graph2);
        std::cout << "-------------" << std::endl;
    }*/

    return 0;
}
