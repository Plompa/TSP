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

#include "Graph.h"
#include "TSP_writer.h"

TSP_enumeration::TSP_enumeration(Graph graph) {
    ncount = graph.getNumNodes();
    distmatrix = graph.getDistmatrix();
    tour = std::vector<int>(ncount);
    besttour = std::vector<int>(ncount);
    mst_map = std::unordered_map<std::vector<bool>, int>();
    mst_map_2 = std::unordered_map<std::bitset<300>, int>();
    shortest_path_map = std::unordered_map<std::vector<bool>, int>();
    pi = std::vector<int>(ncount);

    std::chrono::steady_clock::time_point total_begin = std::chrono::steady_clock::now();

    //k_neigbour_init(10);

    if(!is_sym()) {
        throw std::runtime_error("nicht symmetrisch, der Algo setzt symmetrie vorraus");
    }

    for (int i = 0; i < ncount; i++) tour[i] = i; //initialisieren für permutation über alle touren

    twoOptInit();
    display_data["two_opt_tour"] = std::vector<int>(besttour);

    held_karp_tuning();
    display_data["held_karp_pi"] = pi;

    bestlen = one_tree_bound();

    std::chrono::steady_clock::time_point staged_search_begin = std::chrono::steady_clock::now();
    //staged search, wenn der guess zu niedrig war um eine tour zu finden erhöhe solange bis eine gefunden wird
    while (!tourFound) {
        for (int i = 0; i < ncount; i++) tour[i] = i;
        permute(ncount-1,0, tour[ncount - 1] == 0);

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

    printf ("Optimal Tour Length = %d\n", best_tour_length());
    printf ("Optimal Tour: ");
    for (int i = 0; i < ncount; i++) printf ("%d ", graph.get_node_ids()[besttour[i]]);
    printf ("\n");

    std::cout << "mst map size: " << mst_map.size() << " uses: " << map_uses << std::endl;
    std::cout << "path map size: " << shortest_path_map.size() << " uses: " << shortest_path_map_uses << std::endl;

    std::cout << "total calc time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - total_begin).count() << "[ms] = ";
    std::cout << std::chrono::duration_cast<std::chrono::seconds> (end - total_begin).count() << "[s]" << std::endl;

    std::cout << "staged search time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - staged_search_begin).count() << "[ms] = ";
    std::cout << std::chrono::duration_cast<std::chrono::seconds> (end - staged_search_begin).count() << "[s]" << std::endl;

    TSP_writer writer;
    writer.save(besttour, best_tour_length(), graph.getName() + "_haupt", graph.get_node_ids(), display_data);
}

void TSP_enumeration::permute(int k, int tourlen, bool has_zero) {
    //die abbruchbedingungen sind so geordnet wie sie am effektivsten sind (habe alle permutationen getestet)

    //haupt abbruchbedingung; die länge der aktuellen tour plus die länge des mst der übrigen städte ist schon größer als die beste bekannte tour
    if(tourlen+mst(k+1) >= bestlen) return; // >= two opt

    //überprüfe ob es einen kürzeren bekannten weg gibt
    if(k <= ncount - 4) { //falls k > ncount - 4 ist der weg schon eindeutig
        if(!dynamic_shortest_path(ncount - 1, k, tourlen)) {
            return;
        }
    }

    //nur touren bei denen 0 vor 1 kommt werden betrachtet, da sonst beide richtungen einer tour betrachtet werden würden
    if(!has_zero && tour[k] == 1) return;

    int i;
    if (k == 1) {
        tourlen += (dist(tour[0],tour[1]) + dist(tour[ncount-1],tour[0]));
        if (tourlen < bestlen) {
            tourFound = true;
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
//gibt die länge des msts der ersten count nodes zurück, übernommen aus dem kapitel
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
    if (count != ncount) pcity[count++] = tour[ncount-1]; //das count++ ist suspekt. ohne ist der algo jedoch sehr langsam
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
    key[tour[t] + 2 * ncount] = true; //berechnet schlüssel der hashfunktion

    auto findit = shortest_path_map.find(key);
    if (findit != shortest_path_map.end()) {
        int best = findit -> second;

        if(tourlen <= best) { /* <= da bei < durch andere optimierungen eventuell keine küruzeste tour mehr
                * gefunden wird, zb staged search würde nur in einer stage den optimalen pfad entlanggehen aber
                * bei zu kleinem upperbound abbrechen und danach nie wieder den den pfad benutzen
                */
            shortest_path_map[key] = tourlen;
            return true;
        }

        shortest_path_map_uses++;
        return false;
    } else {
        shortest_path_map[key] = tourlen;
        return true;
    }
}

/* berechnet die hk bound nach valenzuela und jones
 * analog zu one_tree_bound mit unmodifizierten gewichten
 * gibt auch die knotengrade des maximalen 1 trees zurück
 * tatsächlich gibt die funktion nur die hk bound für ein festes pi und nicht die "echte"
 * bound welche nur approximiert werden kann
 */
std::vector<int> TSP_enumeration::hk_bound(){
    std::vector<int> best_deg {};
    int best = -MAXCOST;
    for (int start = 0; start < ncount; ++start) {
        std::vector<int> deg = one_tree(start); //berechnet die länge jedes 1 trees
        int len = w_pi_unmodified(deg);
        if(best < len) {
            best = len;
            best_deg = deg;
        }
    }
    best_deg.pop_back();
    best_deg.pop_back();
    best_deg.push_back(best);
    return best_deg;
}

//berechnet die länge des größten one trees, also eine lower bound. man beachte das die kosten modifiziert sind
int TSP_enumeration::one_tree_bound() {
    int best = -MAXCOST;
    for (int start = 0; start < ncount; ++start) {
        std::vector<int> deg = one_tree(start);
        int len = deg[ncount]; //die methode held_karp_one_tree gibt die länge des 1 trees im letzen eintrag des vektors
        //da jeder 1 tree eine untere Schranke für die Länge der optimalen Tour ist,
        //wird der maximale one tree gesucht, um möglichst nah an den optimalen wert zu kommen
        best = std::max(len, best);
    }
    return best;
}

/* optimierung der kantenlängen, so dass msts touren ähneln -> mst bound bricht besser ab, die wohl wichtigste optimierung
 * im buch wird dieses konzept nur sehr wage beschrieben, ich habe nach etwas googeln das paper von valenzuela und jones gefunden
 * und mich an ihrer ausführung orientiert
 * statt nur von einem knoten aus die 1 trees zu betrachten, betrachte ich die von jedem knoten einzelnd und mittele am ende über
 * die modifizierten gewichte, wie im kapitel beschrieben
 */
void TSP_enumeration::held_karp_tuning() {
    int best_hk_bound = -MAXCOST;

    std::vector<int> pi_sum(ncount);
    for (int start = 0; start < ncount; ++start) {
        for (int k = 0; k < ncount; ++k) {
            pi[k] = 0;
        }

        std::vector<int> deg = one_tree(start); //erhalte knotengrade bei aktuellem pi = 0
        std::vector<int> deg_prev(ncount);
        for (int k = 0; k < ncount; ++k) deg_prev[k] = deg[k];
        std::vector<int> best_pi(ncount);

        int min_w = w_pi(deg);
        const int len = deg[ncount];
        int step1 = len / (2 * ncount);
        int M = static_cast<int>(ncount * ncount / 50. + 0.5) + ncount + 15;

        for (int m = 0; m < M; ++m) {
            //std::cout << step(m, M, step1) << std::endl;
            int cur_w = w_pi(deg);
            if(cur_w < min_w) {
                min_w = cur_w;
                step1 = deg[ncount] / (2 * ncount);
            }
            for (int k = 0; k < ncount; ++k) {
                //passt pi an: ist der knotengrad höher als 2, so ist der knoten zu "beliebt" und die distanz wird erhöht
                //bei knotengrad kleiner als 2 wird die distanz verringert
                //bei knotengrad 2 passiert nichts
                //so "konvergieren" die knotengrade gegen 2
                pi[k] += static_cast<int>(step(m, M, step1) * (deg[k] - 2));
                //Volgenant and Jonkler bieten eine andere Funktion, war bei mir jedoch schlchter
                /*if(deg[k] == 2) continue;
                pi[k] += (0.6 * step(m, M, step1) * (deg[k] - 2)
                + 0.4 * step(m, M, step1) * (deg_prev[k] - 2));*/
            }
            for (int k = 0; k < ncount; ++k) deg_prev[k] = deg[k]; //update prev
            deg = one_tree(start); //aktualisiert die knotengrade nach verändertem pi
        }

        for (int k = 0; k < ncount; ++k) {
            pi_sum[k] += pi[k];
        }

        //die hk bound maximiert w(pi) über alle alle pi, da hier viele pis berechnet werden
        //verwende ich diese um die hk bound zu approximieren
        std::vector<int> hk_one_tree = hk_bound();
        if(hk_one_tree[ncount] > best_hk_bound) {
            best_hk_bound = hk_one_tree[ncount];
        }

    }
    for (int k = 0; k < ncount; ++k) {
        pi[k] = pi_sum[k] / ncount; //mittelt über alle pis
    }
    std::cout << "HK bound (lower bound): " << std::max(best_hk_bound, hk_bound()[ncount]) << std::endl;
}

//w(pi) im buch, pi ist die anzahl inzidenter kanten im mst
int TSP_enumeration::w_pi(const std::vector<int>& deg) const {
    int sum = deg[ncount];
    for (int k = 0; k < ncount; ++k) {
        sum += pi[k] * (deg[k] - 2);
    }
    return sum;
}

//w(pi) allerdings auf dem originales graph
int TSP_enumeration::w_pi_unmodified(const std::vector<int>& deg) const {
    int sum = deg[ncount + 1];
    for (int k = 0; k < ncount; ++k) {
        sum += pi[k] * (deg[k] - 2);
    }
    return sum;
}

//step function nach valenzuela und jones
double TSP_enumeration::step(const int m, const int M, const double step1) {
    const double step = (1.0 * (m - 1) * (2 * M - 5) / (2 * (M - 1))) * step1
                        - (m - 2) * step1
                        + (0.5 * (m - 1) * (m - 2) / ((M - 1) * (M - 2))) * step1;
    return step;
}

/* berechnet die länge des 1 trees, der startknoten ist derjenige welcher nicht im mst vorkommt
 * gibt außerdem die anzahl der inzidenten kanten im 1 tree für jeden knoten
 */
std::vector<int> TSP_enumeration::one_tree(const int start) {
    const int start_city = tour[start];

    tour_swap(ncount - 1, start);
    std::vector<int> incidentCityCount = held_karp_mst(ncount - 1);
    tour_swap(ncount - 1, start);

    int mindist_1 = MAXCOST; //berechnet die beiden kleinsten inzidenten kanten zu start
    int mindist_2 = MAXCOST;
    int unmodified_mindist_1 = MAXCOST;
    int unmodified_mindist_2 = MAXCOST;
    int mini_1 = 0;
    int mini_2 = 0;
    for (int i = 0; i < ncount; ++i) {
        if(tour[i] == start_city) {
            continue;
        }
        const int curdist = dist(tour[i],start_city);

        if(curdist < mindist_1) {
            mindist_2 = mindist_1;
            unmodified_mindist_2 = unmodified_mindist_1;
            mini_2 = mini_1;

            mindist_1 = curdist;
            unmodified_mindist_1 = unmodified_dist(tour[i],start_city);
            mini_1 = i;
        } else if (curdist < mindist_2) {
            mindist_2 = curdist;
            unmodified_mindist_2 = unmodified_dist(tour[i],start_city);
            mini_2 = i;
        }
    }
    incidentCityCount[tour[mini_1]]++; //passt die knotengrade an da von dem startknoten noch die 2 kleinsten kanten eingefügt werden
    incidentCityCount[tour[mini_2]]++;
    incidentCityCount[start_city] = 2;

    incidentCityCount[ncount] += mindist_1 + mindist_2; //aktualisiert die länge da die beiden kanten hinzugefügt werden
    incidentCityCount[ncount + 1] += unmodified_mindist_1 + unmodified_mindist_2;

    return incidentCityCount;
}

/* gibt die länge des mst der ersten count knoten zurück
 * gibt außerdem die anzahl der inzidenten kanten im mst für jeden knoten zurück, für die berechnung von w_pi
 */
std::vector<int> TSP_enumeration::held_karp_mst(int count) const
/* Adopted from Bentley, Unix Review 1996 */
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
    if (count != ncount) pcity[count + 1] = tour[ncount-1];
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
        incidentCityCount[pcity[mini]]++; //aktualisiert pi der beiden knoten
        incidentCityCount[nearestCity[mini]]++;

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

/* berechnet ausgehend von jedem knoten eine nearest-neighbour tour und davon ausgehend two opt schritte
 * bis keine verbesserung mehr gemacht wird
 */
void TSP_enumeration::twoOptInit() {
    int cur_len = 0;
    for (int i = 0; i < ncount; i++) {
        cur_len = nntour(i); //nearest neighbour tour
        cur_len = two_opt(cur_len); //anwenden der two opt schritte
        if(cur_len < bestlen) { //merken der minimalen tour
            bestlen = cur_len;
            for (int j = 0; j < ncount; ++j) {
                besttour[j] = tour[j];
            }
        }
    }
    std::cout << "2opt len (upper bound): " << bestlen << std::endl;
}

int TSP_enumeration::two_opt(const int cur_len) {
    int len = cur_len;
    bool improved = true;
    while (improved) {
        improved = false;

        //sucht nach knoten, deren kanten zu tauschen die länge verkürzt
        for (int i = 0; i < ncount; i++) {
            for (int j = i + 2; j < ncount; ++j) {
                if(i == 0 && j == ncount - 1) {
                    continue;
                }

                //unter dieser graph und tour struktur stellt sich der two opt swap schritt als gar nicht mal so einfach heraus

                int dist_a = dist(tour[i],tour[(i + 1) % ncount]);
                int dist_b = dist(tour[j],tour[(j + 1) % ncount]);
                int dist_c = dist(tour[i],tour[j]);
                int dist_d = dist(tour[(i + 1) % ncount],tour[(j + 1) % ncount]);

                int curDist = dist_a + dist_b;
                int newDist = dist_c + dist_d;

                if (newDist < curDist) {
                    improved = true;
                    len = len - curDist + newDist;

                    for (int k = 0; k < (j - i) / 2; ++k) {
                        tour_swap(i + k + 1, j - k); //dreht die reihenfolge der knoten bei swap der kanten, ein bild malen hilft
                    }
                    i = ncount; //um beide for schleifen zu beenden, hässlich ik
                    break;
                }
            }
        }
    }
    return len;
}

//gibt die länge der nearest neighbour tour ausgehend von start zurück
//übernommen aus dem kapitel
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

//überprüft ob der graph symmetrisch ist
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
    return distmatrix[i][j] + pi[i] + pi[j];
}

int TSP_enumeration::unmodified_dist(int i, int j) const {
    return distmatrix[i][j];
}

//übernommen aus dem kapitel
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

//übernommen aus dem kapitel
void TSP_enumeration::tour_swap(const int i, const int j)
{
    const int temp = tour[i];
    tour[i] = tour[j];
    tour[j] = temp;
}

//verworfene methoden

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

//debug methode
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