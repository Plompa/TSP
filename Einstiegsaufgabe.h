//
// Created by lucmi on 25/04/2025.
//

#ifndef EINSTIEGSAUFGABE_H
#define EINSTIEGSAUFGABE_H
#include "Graph.h"


class Einstiegsaufgabe {
public:
    Einstiegsaufgabe(Graph graph);
private:
    struct Edge {
        int head, tail, weight;
        Edge (int head, int tail, int weight) : head(head), tail(tail), weight(weight) {};
    };
    std::vector<Edge> edges;
};



#endif //EINSTIEGSAUFGABE_H
