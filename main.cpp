#include <fstream>
#include <iostream>

#include "Einstiegsaufgabe.h"
#include "Graph.h"
#include "TSP_enumeration.h"
#include "TSP_reader.h"

int main(int argc, char * argv[]) {
    int out_offset = 0; //offsetet die -out parameter da <seed> optional ist

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
    } else if(arg1 ==  "-rand") {
        //seed in graph

        int ncount = 0;

        if(argc >= 3) {
            ncount = atoi(argv[2]); //anzahl knoten
        } else {
            std::cout << "please enter city count" << std::endl;
            return 1;
        }

        int seed = 1; //1 als seed klappt gut :)
        if(argc >= 4) {
            std::string arg3(argv[3]);
            if(arg3 != "-out") {
                seed = atoi(argv[3]); //seed
                out_offset = 1; //falls ein seed angegeben wurde werden die restlichen parameter um 1 verschoben
            }
        }

        graph =  Graph(ncount);
        graph.random2d("EUC_2D", 1000, seed);
    } else {
        std::cout << "invalid arguments" << std::endl;
        return 1;
    }

    if(argc >= out_offset + 4) {
        std::string arg3(argv[out_offset + 3]);
        if(arg3 == "-out") {
            if(argc >= out_offset + 5) {
                graph.setName(argv[out_offset + 4]);
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

    return 0;
}
