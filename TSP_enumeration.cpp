//
// Created by lucmi on 22/02/2025.
//

#include "TSP_enumeration.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <ostream>
#include <queue>
#include <unordered_map>
#include <bits/ranges_algo.h>

#include "Graph.h"
#include "TSP_reader.h"
#include "TSP_writer.h"

TSP_enumeration::TSP_enumeration(Graph graph) {
    ncount = graph.getNumNodes();
    distmatrix = graph.getDistmatrix();
    tour = std::vector<int>(ncount);
    besttour = std::vector<int>(ncount);
    mst_map = std::unordered_map<std::vector<bool>, int>();
    shortest_path_map = std::unordered_map<std::vector<bool>, int>();
    delta = std::vector<int>(ncount);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //k_neigbour_init(10); //TODO vor hk tuning //7984

    if(!is_sym()) {
        std::cout << "nicht symmetrisch, two opt klappt nicht mehr(?)" << std::endl;
        return;
    }

    for (int i = 0; i < ncount; i++) tour[i] = i; //initialisieren für permutation über alle touren

    /*std::vector<int> test = held_karp_mst(ncount);
    std::cout << "w: " << w_delta(test) << " mst inzidenzen: ";
    for (int i = 0; i < ncount; i++) std::cout << test[i] << " ";
    std::cout << std::endl;*/

    held_karp_tuning();
    display_data["held_karp_delta"] = delta;

    /*test = held_karp_mst(ncount);
    std::cout << "w: " << w_delta(test) << " held-karp mst inzidenzen: ";
    for (int i = 0; i < ncount; i++) std::cout << test[i] << " ";
    std::cout << std::endl;
    std::cout << "HK-bound: " << one_tree_bound() << std::endl;*/


    twoOptInit();
    display_data["two_opt_tour"] = std::vector<int>(besttour);
    int two_opt_bound = bestlen;

    int one_tree_bound_val = one_tree_bound();
    std::cout << "1-tree bound: " << one_tree_bound_val << std::endl;

    int hk_bound_val = hk_bound();
    std::cout << "HK-bound: " << hk_bound_val << std::endl;

    bestlen = (one_tree_bound_val + two_opt_bound) / 2;
    //bestlen = one_tree_bound_val; //achtung tour found wird immer true gesetzt nach permute
    if (one_tree_bound_val > two_opt_bound) bestlen = two_opt_bound;
    //std::cout << "guess: " << bestlen << std::endl;
    //TODO dynamic shortest path darf nicht staged search speichern

    //ich berechne c_T falsch (muss mit originalen kosten)
    //bestlen = 2833; //20 knoten //findet nicht
    //bestlen = 2834; //20 knoten //findet

    while (!tourFound) {
        shortest_path_map = std::unordered_map<std::vector<bool>, int>(); //TODO das ist doof aber ohne klappts nicht

        //std::cout << "bound: " << bestlen << std::endl;

        for (int i = 0; i < ncount; i++) tour[i] = i;
        permute(ncount-1,0, tour[ncount - 1] == 0);
        //print_path_map();

        //tourFound = true;
        if(!tourFound) {
            bestlen = static_cast<int>(bestlen * 1.001) + 1;
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    //speichern von zusätzlichen daten in der opt.tour datei
    display_data["besttour"] = besttour;
    std::vector<int> coords_x{};
    std::vector<int> coords_y{};
    for (auto coord: graph.get_int_coords()) {
        coords_x.push_back(coord[0]);
        coords_y.push_back(coord[1]);
    }
    display_data["coords_x"] = coords_x;
    display_data["coords_y"] = coords_y;

    printf ("Modified Optimal Tour Length = %d\n", bestlen);
    printf ("Unmodified Optimal Tour Length = %d\n", best_tour_length());
    /*printf ("Optimal Tour: ");
    for (int i = 0; i < ncount; i++) printf ("%d ", besttour[i]);
    printf ("\n");*/

    std::cout << "mst map size: " << mst_map.size() << " uses: " << map_uses << std::endl;
    std::cout << "path map size: " << shortest_path_map.size() << " uses: " << shortest_path_map_uses << std::endl;

    std::cout << "calc time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms] = ";
    std::cout << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << "[s]" << std::endl;

    TSP_writer writer;
    writer.save(besttour, bestlen, display_data, graph.getName() + "_haupt");
}

void TSP_enumeration::permute(int k, int tourlen, bool has_zero) {

    if (tourlen+mst(k+1) >= bestlen) return; // >= two opt

    if(k <= ncount - 4) {
        if(!dynamic_shortest_path(ncount - 1, k, tourlen)) {
            return;
        }
    }

    if(!has_zero && tour[k] == 1) return;

    /*int cur_tourlen = tourlen;
    for (int s = ncount - 1; k <= s - 3; s--) {
        if(!dynamic_shortest_path(s, k, cur_tourlen)) {
            return;
        }
        cur_tourlen -= dist(tour[s-1], tour[s]);
    }*/

    int i;
    if (k == 1) {
        tourlen += (dist(tour[0],tour[1]) + dist(tour[ncount-1],tour[0]));
        if (tourlen < bestlen) { // < two opt
            tourFound = true;
            //std::cout << "found " << tourlen << std::endl;
            bestlen = tourlen;
            for (i = 0; i < ncount; i++) besttour[i] = tour[i];
        }
    } else {
        for (i = 0; i < k; i++) {
            tour_swap(i,k-1);
            permute(k-1, tourlen+dist(tour[k-1],tour[k]), has_zero || tour[k-1] == 0);
            tour_swap(i,k-1);
        }
    }
}

int TSP_enumeration::mst(int count) /* Adopted from Bentley, Unix Review 1996 */
//finds length of mst for first count cities
{
    std::vector<bool> key(ncount);
    for (int i = 0; i < count; i++) {
        key[tour[i]] = true;
    }
    auto findit = mst_map.find(key);
    if (findit != mst_map.end()) {
        map_uses++;
        return findit -> second;
    }

    int i, m, mini, newcity, mindist, thisdist, len = 0;
    int pcity[ncount], pdist[ncount];
    if (count <= 1) return 0;
    for (i = 0; i < count; i++) {
        pcity[i] = tour[i];
        pdist[i] = MAXCOST;
    }
    if (count != ncount) pcity[count++] = tour[ncount-1];
    newcity = pcity[count-1];
    for (m = count-1; m > 0; m--) {
        mindist = MAXCOST;
        for (i = 0; i < m; i++) {
            thisdist = dist(pcity[i],newcity);
            if (thisdist < pdist[i]) {
                pdist[i] = thisdist;
            }
            if (pdist[i] < mindist) {
                mindist = pdist[i];
                mini = i;
            }
        }
        newcity = pcity[mini];
        len += mindist;
        pcity[mini] = pcity[m-1];
        pdist[mini] = pdist[m-1];
    }

    mst_map[key] = len;
    return len;
}

/*eine methode die für jeden knoten die kanten der nächsten k städte behält, und für die restlichen
städte entfernt (das gewicht der kante wird auf maxcost gesetzt)
diese methode war ein test um für größere instanzen die kantenmenge zu reduzieren und somit eine
heuristik für die optimale tour*/
void TSP_enumeration::k_neigbour_init(int k) {
    if(k >= ncount) return;

    std::vector new_distmatrix(ncount, std::vector<int>(ncount, 100000));

    for (int i = 0; i < ncount; ++i) {
        std::vector<int> neighbour_dists;
        for (int j = 0; j < ncount; ++j) {
            if (i == j) {
                continue;
            }
            neighbour_dists.push_back(dist(i, j));
        }
        std::sort(neighbour_dists.begin(), neighbour_dists.end());
        int max_dist = neighbour_dists[k - 1];
        for (int j = 0; j < ncount; ++j) {
            if (dist(i, j) <= max_dist) {
                new_distmatrix[i][j] = distmatrix[i][j];
                new_distmatrix[j][i] = distmatrix[i][j];
            }
        }
    }
    distmatrix = new_distmatrix;
}

/* merkt sich für startknoten s, zielknoten t und knotenmenge S - die knoten auf dem s-t weg der tour -
 * die länge des weges. wurde bereits ein kürzerer s-t weg über S gefunden kann die tour nicht mehr
 * optimal sein, vergleiche aufgabe 8
 * es wird kein optimaler weg berechnet sondern nur der aktuell beste gemerkt!
 */
bool TSP_enumeration::dynamic_shortest_path(int s, int t, int tourlen) {
    std::vector<bool> key(3 * ncount);
    for (int i = t + 1; i < s; i++) {
        key[tour[i]] = true;
    }
    key[tour[s] + ncount] = true;
    key[tour[t] + 2 * ncount] = true;

    auto findit = shortest_path_map.find(key);
    if (findit != shortest_path_map.end()) {
        int best = findit -> second;
        if(s == ncount - 1) {
            if(tourlen < best) {
                shortest_path_map[key] = tourlen;
                return true;
            }
        } else {
            if(tourlen <= best) {
                shortest_path_map[key] = tourlen;
                return true;
            }
        }

        shortest_path_map_uses++;
        return false;
    } else {
        shortest_path_map[key] = tourlen;
        return true;
    }
}

/* berechnet den kürzesten weg vom startknoten zu k über die knotenmenge S auf dem weg
 * dynamisch mit hashtable um nicht denselben weg doppelt zu berechnen
 */
int TSP_enumeration::shortest_path(int k) {
    std::vector<bool> key(2 * ncount);
    for (int i = k; i < ncount; i++) {
        key[tour[i]] = true;
    }
    key[tour[k] + ncount] = true;

    auto findit = shortest_path_map.find(key);
    if (findit != shortest_path_map.end()) {
        shortest_path_map_uses++;
        return findit -> second;
    }

    int len = path_permute(ncount - 1, k, 0,MAXCOST);

    shortest_path_map[key] = len;
    return len;
}

/* berechnet den kürzesten weg von stadt k zu stadt t über städte auf dem aktuellen k-t weg
 */
int TSP_enumeration::path_permute(int k, int t, int pathlen, int best_pathlen) {
    //std::cout << std::string((ncount - 1 - k) * 3, ' ') << pathlen << std::endl;
    if (k == t + 1) {
        //std::cout << "path length = " << pathlen << std::endl;
        //std::cout << "+ " << k << "-" << k - 1 << std::endl;
        return std::min(pathlen + dist(tour[k-1],tour[k]), best_pathlen);
    }
    for (int i = t + 1; i < k; i++) {
        tour_swap(i,k-1);
        //std::cout << "+ " << k << "-" << k - 1 << std::endl;
        best_pathlen = path_permute(k - 1, t, pathlen + dist(tour[k-1],tour[k]), best_pathlen);
        tour_swap(i,k-1);
    }
    return best_pathlen;
}

int TSP_enumeration::hk_bound() {
    int best = MAXCOST;
    for (int start = 0; start < ncount; ++start) {
        std::vector<int> deg = held_karp_one_tree(start);
        int len = w_delta_2(deg);
        best = std::min(len, best);
    }
    return best;
}

int TSP_enumeration::one_tree_bound() {
    int best = MAXCOST;
    for (int start = 0; start < ncount; ++start) {
        std::vector<int> deg = held_karp_one_tree(start);
        int len = deg[ncount];
        best = std::min(len, best);
    }
    return best;
}

std::vector<int> TSP_enumeration::max_one_tree() {
    std::vector<int> ret;
    int best = 0;
    for (int start = 0; start < ncount; ++start) {
        std::vector<int> deg = held_karp_one_tree(start);
        int len = deg[ncount];
        if(len > best) {
            best = len;
            ret = deg;
        }
    }
    return ret;
}

void TSP_enumeration::held_karp_tuning() {
    std::vector<int> delta_sum(ncount);
    for (int start = 0; start < ncount; ++start) {
        for (int k = 0; k < ncount; ++k) {
            delta[k] = 0;
        }

        std::vector<int> deg = held_karp_one_tree(start);
        std::vector<int> deg_prev(ncount);
        for (int k = 0; k < ncount; ++k) deg_prev[k] = deg[k];
        std::vector<int> best_delta(ncount);

        int min_w = w_delta(deg);
        int len = deg[ncount];
        int step1 = len / (2 * ncount);
        int M = static_cast<int>(ncount * ncount / 50. + 0.5) + ncount + 15;

        for (int m = 0; m < M; ++m) {
            //std::cout << step(m, M, step1) << std::endl;
            int cur_w = w_delta(deg);
            if(cur_w < min_w) {
                min_w = cur_w;
                step1 = deg[ncount] / (2 * ncount);
                //for (int k = 0; k < ncount; ++k) best_delta[k] = delta[k];
            }
            for (int k = 0; k < ncount; ++k) {
                delta[k] += static_cast<int>(step(m, M, step1) * (deg[k] - 2));
                //Volgenant and Jonkler
                /*if(deg[k] == 2) continue;
                delta[k] += (0.6 * step(m, M, step1) * (deg[k] - 2)
                + 0.4 * step(m, M, step1) * (deg_prev[k] - 2));*/
            }
            for (int k = 0; k < ncount; ++k) deg_prev[k] = deg[k]; //update prev
            deg = held_karp_one_tree(start);
        }

        for (int k = 0; k < ncount; ++k) {
            delta_sum[k] += delta[k];
            //delta_sum[k] += best_delta[k];
        }
    }
    for (int k = 0; k < ncount; ++k) {
        delta[k] = delta_sum[k] / ncount;
    }
}

//TODO immer den kleinsten (oder größten?) one_tree benutzen und nur eine schleife statt average
void TSP_enumeration::held_karp_tuning_2() {
    for (int k = 0; k < ncount; ++k) {
        delta[k] = 0;
    }

    std::vector<int> deg = max_one_tree();

    std::vector<int> deg_prev(ncount);
    for (int k = 0; k < ncount; ++k) deg_prev[k] = deg[k];
    std::vector<int> best_delta(ncount);

    int min_w = w_delta_2(deg);
    int len = deg[ncount];
    int step1 = len / (2 * ncount);
    int M = static_cast<int>(ncount * ncount / 50. + 0.5) + ncount + 15;

    for (int m = 0; m < M; ++m) {
        //std::cout << step(m, M, step1) << std::endl;
        int cur_w = w_delta_2(deg);
        if(cur_w < min_w) {
            min_w = cur_w;
            step1 = deg[ncount] / (2 * ncount);
            //for (int k = 0; k < ncount; ++k) best_delta[k] = delta[k];
        }
        for (int k = 0; k < ncount; ++k) {
            delta[k] += static_cast<int>(step(m, M, step1) * (deg[k] - 2));
            //Volgenant and Jonkler
            /*if(deg[k] == 2) continue;
            delta[k] += (0.6 * step(m, M, step1) * (deg[k] - 2)
            + 0.4 * step(m, M, step1) * (deg_prev[k] - 2));*/
        }
        for (int k = 0; k < ncount; ++k) deg_prev[k] = deg[k]; //update prev

        deg = max_one_tree();
    }
}

int TSP_enumeration::w_delta(const std::vector<int>& deg) const {
    int sum = deg[ncount];
    for (int k = 0; k < ncount; ++k) {
        sum += delta[k] * (deg[k] - 2);
    }
    return sum;
}

int TSP_enumeration::w_delta_2(const std::vector<int>& deg) const {
    int sum = deg[ncount + 1];
    for (int k = 0; k < ncount; ++k) {
        sum += delta[k] * (deg[k] - 2);
    }
    return sum;
}

double TSP_enumeration::step(const int m, const int M, const double step1) {
    double step = (1.0 * (m - 1) * (2 * M - 5) / (2 * (M - 1))) * step1
    - (m - 2) * step1
    + (0.5 * (m - 1) * (m - 2) / ((M - 1) * (M - 2))) * step1;
    return step;
}

std::vector<int> TSP_enumeration::held_karp_one_tree(const int start) {
    tour_swap(ncount - 1, start);

    std::vector<int> incidentCityCount = held_karp_mst(ncount - 1);
    int len = incidentCityCount[incidentCityCount.size() - 1];
    incidentCityCount.pop_back();

    int mindist_1 = MAXCOST;
    int mindist_2 = MAXCOST;
    int mini_1 = 0;
    int mini_2 = 0;
    for (int i = 0; i < ncount - 1; ++i) {
        int curdist = dist(tour[i],start);

        if(curdist < mindist_1) {
            mindist_1 = curdist;
            mini_1 = i;
        } else if (curdist < mindist_2) {
            mindist_2 = curdist;
            mini_2 = i;
        }
    }
    incidentCityCount[tour[mini_1]]++;
    incidentCityCount[tour[mini_2]]++;
    incidentCityCount[start] = 2;

    tour_swap(ncount - 1, start);

    incidentCityCount.push_back(len);
    return incidentCityCount;
}

std::vector<int> TSP_enumeration::held_karp_mst(int count) const
/* Adopted from Bentley, Unix Review 1996 */
//finds length of mst for first count cities
{
    int i, m, mini, newcity, mindist, thisdist, len = 0, orig_len = 0;
    int pcity[ncount], pdist[ncount];
    int nearestCity[ncount];
    std::vector<int> incidentCityCount(ncount);
    if (count <= 1) return incidentCityCount;
    for (i = 0; i < count; i++) {
        pcity[i] = tour[i];
        pdist[i] = MAXCOST;
    }
    if (count != ncount) pcity[count++] = tour[ncount-1];
    newcity = pcity[count-1];
    for (m = count-1; m > 0; m--) {
        mindist = MAXCOST;
        for (i = 0; i < m; i++) {
            thisdist = dist(pcity[i],newcity);
            if (thisdist < pdist[i]) {
                pdist[i] = thisdist;
                nearestCity[i] = newcity;
            }
            if (pdist[i] < mindist) {
                mindist = pdist[i];
                mini = i;
            }
        }
        incidentCityCount[pcity[mini]]++;
        incidentCityCount[nearestCity[mini]]++;

        /*std::cout << "--------------------" << std::endl;
        int dist_a = dist(pcity[mini], nearestCity[mini]);
        int dist_b = dist(nearestCity[mini], pcity[mini]);
        if(dist_a != mindist) {
            std::cout << "mist " << dist_a << " " << dist_b << std::endl;
        }
        std::cout << "dist hat: " << dist(pcity[mini], nearestCity[mini]) << " dist soll: " << mindist << std::endl;
        std::cout << "kante: " << pcity[mini] << " " << nearestCity[mini] << std::endl;*/

        newcity = pcity[mini];
        len += mindist;
        orig_len += unmodified_dist(pcity[mini], nearestCity[mini]);
        pcity[mini] = pcity[m-1];
        pdist[mini] = pdist[m-1];
        nearestCity[mini] = nearestCity[m-1];
    }

    incidentCityCount.push_back(len);
    incidentCityCount.push_back(orig_len);
    return incidentCityCount;
}

void TSP_enumeration::twoOptInit() {
    int cur_len = 0;
    for (int i = 0; i < ncount; i++) {
        cur_len = nntour(i);
        //std::cout << "nn len: " << cur_len << std::endl;
        cur_len = two_opt(cur_len);
        //std::cout << "2opt len: " << cur_len << std::endl;
        if(cur_len < bestlen) {
            bestlen = cur_len;
            for (int j = 0; j < ncount; ++j) {
                besttour[j] = tour[j];
            }
        }
    }
    std::cout << "2opt len: " << bestlen << std::endl;
}

int TSP_enumeration::two_opt(int cur_len) {
    int len = cur_len;
    bool improved = true;
    while (improved) {
        improved = false;
        //for (int i = 0; i < ncount; i++) printf ("%d ", tour[i]);
        //printf ("\n");

        for (int i = 0; i < ncount; i++) {
            for (int j = i + 2; j < ncount; ++j) {
                if(i == 0 && j == ncount - 1) {
                    continue;
                }

                int dist_a = dist(tour[i],tour[(i + 1) % ncount]);
                int dist_b = dist(tour[j],tour[(j + 1) % ncount]);
                int dist_c = dist(tour[i],tour[j]);
                int dist_d = dist(tour[(i + 1) % ncount],tour[(j + 1) % ncount]);

                int curDist = dist_a + dist_b;
                int newDist = dist_c + dist_d;

                if (newDist < curDist) {
                    improved = true;
                    len = len - curDist + newDist;
                    //tour_swap(i, j);
                    //std::cout << tour[i] << " " << tour[j] << std::endl;

                    for (int k = 0; k < (j - i) / 2; ++k) {
                        tour_swap(i + k + 1, j - k);
                    }
                    int actual_len = tour_length();
                    if(len != actual_len) {
                        std::cout << "help" << std::endl;
                    }

                    i = ncount;
                    break;
                }
            }
        }
    }
    return len;
}

int TSP_enumeration::nntour(int start)
{
    int i, j, best, bestj, len = 0;
    for (i = 0; i < ncount; i++) tour[i] = (i + start) % ncount;
    for (i = 1; i < ncount; i++) {
        best = MAXCOST;
        for (j = i; j < ncount; j++) {
            if (dist(tour[i-1],tour[j]) < best) {
                best = dist(tour[i-1],tour[j]); bestj = j;
            }
        }
        len += best; tour_swap(i, bestj);
    }
    return len+dist(tour[ncount-1],tour[0]);
}

bool TSP_enumeration::is_sym() const {
    for (int i = 0; i < ncount; ++i) {
        for (int j = i + 1; j < ncount; ++j) {
            if(dist(i, j) != dist(j, i)) {
                return false;
            }
        }
    }
    return true;
}

int TSP_enumeration::dist(int i, int j) const {
    return distmatrix[i][j] + delta[i] + delta[j];
}

int TSP_enumeration::unmodified_dist(int i, int j) const {
    return distmatrix[i][j];
}

int TSP_enumeration::tour_length() const {
    int len = 0;
    for (int i = 1; i < ncount; i++) {
        len += dist(tour[i-1],tour[i]);
    }
    return len + dist(tour[ncount-1],tour[0]);
}

int TSP_enumeration::subtour_length(int s, int t) const {
    int len = 0;
    for (int i = s; i > t; i--) {
        len += dist(tour[i-1], tour[i]);
    }
    return len;
}

int TSP_enumeration::best_tour_length() const {
    int len = 0;
    for (int i = 1; i < ncount; i++) {
        len += unmodified_dist(besttour[i-1],besttour[i]);
    }
    return len + unmodified_dist(besttour[ncount-1],besttour[0]);
}

void TSP_enumeration::tour_swap(const int i, const int j)
{
    const int temp = tour[i];
    tour[i] = tour[j];
    tour[j] = temp;
}

void TSP_enumeration::print_path_map() const {
    for (const auto& i: shortest_path_map) {
        for (int j = ncount; j < 2 * ncount; ++j) {
            if (i.first[j]) {
                std::cout << j - ncount << " ";
            }
        }

        for (int j = 0; j < ncount; ++j) {
            if (i.first[j]) {
                std::cout << j << " ";
            }
        }

        for (int j = 2 * ncount; j < 3 * ncount; ++j) {
            if (i.first[j]) {
                std::cout << j - 2 * ncount << " ";
            }
        }
        std::cout << ": " << i.second << std::endl;
    }
}
