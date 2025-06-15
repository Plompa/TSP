//
// Created by lucmi on 25/04/2025.
//

#include "Einstiegsaufgabe.h"

#include <algorithm>
#include <iostream>
#include <ostream>

#include "TSP_writer.h"

Einstiegsaufgabe::Einstiegsaufgabe(Graph graph) {
    for (int i = 0; i < graph.getNumNodes(); ++i) {
        for (int j = 0; j < graph.getNumNodes(); ++j) {
            edges.push_back(Edge(i, j, graph.getDistmatrix()[i][j]));
        } //der graph ist im wesentlichen eine adjazenzmatrix, um die kanten zu sortieren werden ersteinmal hilfsobjekte erstellt
    }
    std::sort(edges.begin(), edges.end(), //sortiert die kanten nach gewicht
        [](const Edge &e1, const Edge &e2) {return e1.weight < e2.weight;});

    int tour_len = 0;

    //in einer tour wird jeder knoten genau einmal als head und einmal als tail verwendet, mit diesen listen wird auf dieses limit geachtet
    std::vector<bool> used_as_head(graph.getNumNodes(), false);
    std::vector<bool> used_as_tail(graph.getNumNodes(), false);

    /*der kern des algos: während des einfügen der kanten bilden sich teiltouren - ketten
     * um sicherzugehen das keine kette vorzeitig zu einem kreis geschlossen wird, muss buch über die zusammenhangskomponenten geführt werden
     * used_as_head/tail stellen sicher, dass kein innerer knoten der kette betrachtet werden muss
     * also muss sich nur um die enden gekümmert werden; sie müssen einander kennen um in konstanter laufzeit einen potenziellen kreisschluss zu sehen
     * other_end bildet die endstücke aufeinander ab
     * insbesondere hat dieser einfügealgo konstante zeit, unionfind hätte log
     */
    std::vector<int> other_end(graph.getNumNodes(), -1);
    for (int i = 0; i < graph.getNumNodes(); ++i) other_end[i] = i; //anfangs gibt es keine kette, also ist other_end die identität

    //speichert die verwendeten kanten zur rekonstruktion der tour
    std::vector<Edge> tour_edges;

    for (auto cur : edges) {
        if(!used_as_head[cur.head] && !used_as_tail[cur.tail]) { //innere knoten von ketten sind nicht mehr zu betrachten
            if(tour_edges.size() < graph.getNumNodes() - 1) { //die letzte kante darf einen kreis schliessen
                if(other_end[cur.head] == cur.tail) { //falls das andere ende des heads der tail ist, wird ein kreis geschlossen -> nicht einfügen
                    continue;
                }
            }

            tour_edges.push_back(cur); //füge die kante ein
            tour_len += cur.weight; //aktualisiere tour länge

            int other_end_tail = other_end[cur.tail]; //aktualisiere other_end; die neuen enden der kette sind other_end von beiden enden der kante
            int other_end_head = other_end[cur.head]; //ich empfehle ein bild zu malen
            other_end[other_end_tail] = other_end_head; //nur die enden müssen aktualisiert werden
            other_end[other_end_head] = other_end_tail; //die köpfe der kette werden auf die neuen enden gesetzt

            used_as_head[cur.head] = true; //aktualisiere die used_as_head/tail
            used_as_tail[cur.tail] = true;
        }
    }

    std::cout << "NN-Tour length: " << tour_len << std::endl;
    std::cout << "Tour: " << std::endl;

    std::vector<int> tour;
    /* rekonstruiere tour:
     * starte bei 0, suche die kante (0, w) die in 0 anfängt, wiederhole mit w usw
     */
    int from = 0;
    for (int i = 0; i < graph.getNumNodes(); ++i) {
        std::cout << from << " ";
        tour.push_back(from);
        for (auto edge: tour_edges) {
            if(edge.tail == from) {
                from = edge.head;
                break;
            }
        }
    }
    std::cout << std::endl;

    TSP_writer writer;
    writer.save(tour, tour_len, graph.getName() + "_einst");
}

