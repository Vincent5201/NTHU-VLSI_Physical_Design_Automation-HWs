#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <queue>
#include <cmath>

using namespace std;

struct Net {
    string name;
    int idx;
    int x0, y0, x1, y1;
    double maxcost;
    vector<pair<int, int>> paths;
    vector<pair<int, int>> ppaths;  // previous paths
    vector<pair<int, int>> bpaths;  // best paths
};

// used for A*-search routing
struct Node {
    int x, y;
    int px, py;                 // parent x/y
    double costo, costd;        // cost of overflow/distance
    
    Node(int xi, int yi, int pxi, int pyi, double costoi, double costdi): 
            x(xi), y(yi), px(pxi), py(pyi), costo(costoi), costd(costdi) {}
    bool operator>(const Node &other) const {
        return costo + costd > other.costo + other.costd;
    }
};


int gridx, gridy, capx, capy, numNets;
int bov, bWL, ov, WL;
vector<Net> nets;
vector<vector<pair<int, int>>> edgesUse;            // count of nets on edges
vector<vector<pair<double, double>>> edgesCost;     // cost od esges
vector<vector<double>> prefCol;                     // column prefix sum of edgesCost
vector<vector<double>> prefRow;                     // row prefix sum of edgesCost
vector<vector<Node*>> visited;                      // visited node in A*-search routing

int h = 50;
double k = 0.7;
int mybetta = 3;
double mygamma = 0.1;
double mylambda = 1;


void read_file(string &inputfile) {
    ifstream lef(inputfile);
    string line;
    int count = 0;
    
    while (getline(lef, line)) {
        istringstream iss(line);
        string token, tmp;
        iss >> token;

        if (token == "Net") {
            count++;
            Net net;
            net.idx = count;
            iss >> net.name;
            getline(lef, line);
            iss.clear();
            iss.str(line);
            iss >> tmp >> tmp;
            iss >> net.x0 >> net.y0;
            getline(lef, line);
            iss.clear();
            iss.str(line);
            iss >> tmp >> tmp;
            iss >> net.x1 >> net.y1;
            net.paths.reserve(4);
            net.ppaths.reserve(4);
            net.bpaths.reserve(4);
            nets[count] = move(net);

        } else if (token == "Grid") {
            iss >> gridx >> gridy;
        } else if (token == "Capacity") {
            iss >> capx >> capy;
        } else if (token == "NumNets") {
            iss >> numNets;
            nets = vector<Net>(numNets + 1);
        }
    }
    lef.close();    

    edgesUse = vector<vector<pair<int, int>>>(gridx, 
                vector<pair<int, int>>(gridy, {0, 0}));
    edgesCost = vector<vector<pair<double, double>>>(gridx, 
                vector<pair<double, double>>(gridy, {0, 0}));
    prefCol = vector<vector<double>>(gridx + 1, vector<double>(gridy + 1, 0));
    prefRow = vector<vector<double>>(gridx + 1, vector<double>(gridy + 1, 0));
}

// cost of use d, capacity s
inline double icost(double d, int s) {
    return 1 + h / (1 + exp(-k * (d - s)));
}

// sum of cost between (x1, y) (x2, y)
inline double sumx(int x1, int x2, int y) {
    if (x1 > x2)
        return prefRow[x1][y] - prefRow[x2][y];
    return prefRow[x2][y] - prefRow[x1][y];
}

// sum of cost between (x, y2) (x, y2)
inline double sumy(int x, int y1, int y2) {
    if (y1 > y2)
        return prefCol[x][y1] - prefCol[x][y2];
    return prefCol[x][y2] - prefCol[x][y1];
}

// for x in [xa, xb) do func(x)
template<typename Func>
inline void forXRange(int xa, int xb, Func func) {
    if (xa == xb)
        return;
    int x0 = min(xa, xb);
    int x1 = max(xa, xb);
    for (int x = x0; x < x1; x++)
        func(x);
}

// for y in [ya, yb) do func(y)
template<typename Func>
inline void forYRange(int ya, int yb, Func func) {
    if (ya == yb)
        return;
    int y0 = min(ya, yb);
    int y1 = max(ya, yb);
    for (int y = y0; y < y1; y++)
        func(y);
}

// update prefCol / prefRow
void update_pref(int x0, int y0, int x1, int y1) {
    if (x0 != x1 && y0 != y1) {
        cout << "error line2\n";
        return;
    }
    if (x0 == x1) {
        forYRange(y0, y1, [&](int y) { 
            edgesCost[x0][y].second = icost(edgesUse[x0][y].second, capy); });
        
        for (int j = 0; j < gridy; j++) {
            prefCol[x0][j+1] = prefCol[x0][j] + edgesCost[x0][j].second;
            if (prefCol[x0][j+1] < 0)
                cout << "error prefcol1\n";
        }
    } else {
        forXRange(x0, x1, [&](int x) { 
            edgesCost[x][y0].first = icost(edgesUse[x][y0].first, capx); });
        
        for (int i = 0; i < gridx; i++) {
            prefRow[i+1][y0] = prefRow[i][y0] + edgesCost[i][y0].first;
            if (prefRow[i+1][y0] < 0)
                cout << "error prefrow1\n";
        }
    }
}

// put a net on (x0, y0) to (x1, y1)
void put_line(int x0, int y0, int x1, int y1) {
    if (x0 != x1 && y0 != y1) {
        cout << "error line1\n";
        return;
    }
    if (x0 == x1) {
        forYRange(y0, y1, [&](int y) { edgesUse[x0][y].second++; });
    } else {
        forXRange(x0, x1, [&](int x) { edgesUse[x][y0].first++; });
    }
}

// calculate total overflow
int get_overflow(bool reset) {
    if (reset) {
        for (int x = 0; x < gridx; ++x)
            for (int y = 0; y < gridy; ++y)
                edgesUse[x][y] = {0, 0};
            
        for (int i = 0; i <= numNets; i++) {
            Net &net = nets[i];
            if (net.paths.size() == 0) {
                put_line(net.x0, net.y0, net.x1, net.y1);
            } else {
                vector<pair<int, int>> &paths = net.paths;
                put_line(net.x0, net.y0, paths[0].first, paths[0].second);
                put_line(net.x1, net.y1, 
                        paths[paths.size() - 1].first, paths[paths.size() - 1].second);
                for (int j = 1; j < (int)paths.size(); j++) {
                    put_line(paths[j - 1].first, paths[j - 1].second, 
                            paths[j].first, paths[j].second);
                }
            }
        }
    }

    int ans = 0;
    for (int i = 0; i < gridx; i++) {
        for (int j = 0; j < gridy; j++) {
            ans += max(0, edgesUse[i][j].first - capx);
            ans += max(0, edgesUse[i][j].second - capy);
        }
    }
    return ans;
}

// get WL on a straight net
int get_line_WL(int x0, int y0, int x1, int y1) {
    int ans = 0;
    if (x0 == x1) {
        ans += abs(y0 - y1);
    } else if (y0 == y1) {
        ans += abs(x0 - x1);
    } else {
        cout << "error line\n";
        cout << x0 << " " << x1 << " " << y0 << " " << y1 << endl;
    }
    return ans;
}

// get total WL
int get_WL() {
    int ans = 0;
    for (int i = 0; i <= numNets; i++) {
        Net &net = nets[i];
        if (net.paths.size() == 0) {
            ans += get_line_WL(net.x0, net.y0, net.x1, net.y1);
        } else {
            vector<pair<int, int>> &paths = net.paths;
            ans += get_line_WL(net.x0, net.y0, paths[0].first, paths[0].second);
            ans += get_line_WL(net.x1, net.y1, 
                    paths[paths.size() - 1].first, paths[paths.size() - 1].second);
            for (int j = 1; j < (int)paths.size(); j++) {
                ans += get_line_WL(paths[j - 1].first, paths[j - 1].second, 
                                    paths[j].first, paths[j].second);
            }
        }
    }
    return ans;
}

// check sum of nets is corrected
bool WL_checker() {
    int sum = 0;
    for (int i = 0; i < gridx; i++) {
        for (int j = 0; j < gridy; j++) {
            sum += edgesUse[i][j].first;
            sum += edgesUse[i][j].second;
        }
    }
    return sum == get_WL();
}


double dmax(double a, double b) {
    if (a >= b)
        return a;
    return b;
}

// get the egde has the max cost on a net
double get_maxcost(Net &net) {
    double ans = 0;
    int x0 = net.x0, x1 = net.x1, y0 = net.y0, y1 = net.y1;
    if (net.paths.size() == 0) {
        if (x0 == x1) {
            forYRange(y0, y1, [&](int y) { ans = dmax(ans, edgesCost[x0][y].second); });
        } else {
            forXRange(x0, x1, [&](int x) { ans = dmax(ans, edgesCost[x][y0].first); });
        }
    } else {
        vector<pair<int, int>> &paths = net.paths;
        int x2 = paths[0].first, y2 = paths[0].second;
        if (x0 == x2) {
            forYRange(y0, y2, [&](int y) { ans = dmax(ans, edgesCost[x2][y].second); });
        } else {
            forXRange(x0, x2, [&](int x) { ans = dmax(ans, edgesCost[x][y2].first); });
        }
        x2 = paths[paths.size() - 1].first;
        y2 = paths[paths.size() - 1].second;
        if (x1 == x2) {
            forYRange(y1, y2, [&](int y) { ans = dmax(ans, edgesCost[x2][y].second); });
        } else {
            forXRange(x1, x2, [&](int x) { ans = dmax(ans, edgesCost[x][y2].first); });
        } 
        for (int j = 1; j < (int)paths.size(); j++) {
            if (paths[j - 1].first == paths[j].first) {
                forYRange(paths[j - 1].second, paths[j].second, [&](int y) { 
                    ans = dmax(ans, edgesCost[paths[j].first][y].second); });
            } else {
                forXRange(paths[j - 1].first, paths[j].first, [&](int x) {
                    ans = dmax(ans, edgesCost[x][paths[j].second].first); });
            } 
        }
    }
    return ans;
}

// rid up a straight line
void rid_up_line(int x0, int y0, int x1, int y1) {
    if (x0 != x1 && y0 != y1) {
        cout << "error line1\n";
        return;
    }
    if (x0 == x1) {
        forYRange(y0, y1, [&](int y) { edgesUse[x0][y].second--; });
    } else {
        forXRange(x0, x1, [&](int x) { edgesUse[x][y0].first--; });
    }
}

// reroute a L shape net
bool reroute_L(int x0, int y0, int x1, int y1, int x2, int y2) {
    if (x0 == x2 || y0 == y2) {
        cout << "error L shape1\n";
        return 0;
    }
    if ((x0 != x1 && y0 != y1) || (x2 != x1 && y2 != y1)) {
        cout << "error L shape2\n";
        return 0;
    }

    rid_up_line(x0, y0, x1, y1);
    update_pref(x0, y0, x1, y1);
    rid_up_line(x2, y2, x1, y1);
    update_pref(x2, y2, x1, y1);

    double cost1 = 0, cost3 = 0;
    int x3, y3;
    if (x0 == x1) {
        x3 = x2; y3 = y0;
        cost1 += sumy(x0, y0, y1) + sumx(x1, x2, y2);
        cost3 += sumy(x2, y0, y1) + sumx(x0, x2, y0);
    } else {
        x3 = x0; y3 = y2;
        cost1 += sumy(x1, y0, y1) + sumx(x0, x1, y0);
        cost3 += sumy(x0, y0, y1) + sumx(x0, x2, y2);
    }

    if (cost3 < cost1) {
        put_line(x0, y0, x3, y3);
        update_pref(x0, y0, x3, y3);
        put_line(x2, y2, x3, y3);
        update_pref(x2, y2, x3, y3);
        return true;
    } else {
        put_line(x0, y0, x1, y1);
        update_pref(x0, y0, x1, y1);
        put_line(x2, y2, x1, y1);
        update_pref(x2, y2, x1, y1);
        return false;
    }
    return false;
}


void initial_stage() {

    edgesUse.assign(gridx, vector<pair<int,int>>(gridy, {0, 0}));
    edgesCost.assign(gridx, vector<pair<double,double>>(gridy, {0.0, 0.0}));
    prefCol.assign(gridx + 1, vector<double>(gridy + 1, 0.0));
    prefRow.assign(gridx + 1, vector<double>(gridy + 1, 0.0));

    // add 0.5 prob in two L shape for all nets
    for (int i = 1; i <= numNets; i++) {
        int x0 = nets[i].x0, x1 = nets[i].x1, y0 = nets[i].y0, y1 = nets[i].y1;
        forXRange(x0, x1, [&](int x) { edgesUse[x][y0].first++; });
        forXRange(x0, x1, [&](int x) { edgesUse[x][y1].first++; });
        forYRange(y0, y1, [&](int y) { edgesUse[x0][y].second++; });
        forYRange(y0, y1, [&](int y) { edgesUse[x1][y].second++; });
    }
    
    // calculate cost
    for (int i = 0; i < gridx; i++) {
        for (int j = 0; j < gridy; j++) {
            edgesCost[i][j].first = icost(edgesUse[i][j].first * 0.5, capx);
            edgesCost[i][j].second = icost(edgesUse[i][j].second * 0.5, capy);
        }
    }

    // decide all nets' two L shape (select one)
    for (int i = 1; i <= numNets; i++) {
        int x0 = nets[i].x0, x1 = nets[i].x1, y0 = nets[i].y0, y1 = nets[i].y1;
        if (x0 == x1 || y0 == y1)
            continue;
        double cost1 = 0, cost2 = 0;
        cost1 += sumx(x0, x1, y0) + sumy(x1, y0, y1);
        cost2 += sumx(x0, x1, y1) + sumy(x0, y0, y1);
        if (cost1 < cost2) {
            nets[i].paths.push_back({x1, y0});
        } else {
            nets[i].paths.push_back({x0, y1});
        }
    }

    // reset use
    edgesUse.assign(gridx, vector<pair<int,int>>(gridy, {0,0}));
    for (int i = 1; i <= numNets; i++) {
        int x0 = nets[i].x0, x1 = nets[i].x1, y0 = nets[i].y0, y1 = nets[i].y1;
        if (x0 == x1) {
            forYRange(y0, y1, [&](int y) { edgesUse[x0][y].second++; });
        } else if (y0 == y1) {
            forXRange(x0, x1, [&](int x) { edgesUse[x][y0].first++; });
        } else {
            int x2 = nets[i].paths[0].first, y2 = nets[i].paths[0].second;
            forXRange(x0, x2, [&](int x) { edgesUse[x][y0].first++; });
            forXRange(x1, x2, [&](int x) { edgesUse[x][y1].first++; });
            forYRange(y1, y2, [&](int y) { edgesUse[x1][y].second++; });
            forYRange(y0, y2, [&](int y) { edgesUse[x0][y].second++; });
        }
    }

    // calculate cost
    for (int i = 0; i < gridx; i++) {
        for (int j = 0; j < gridy; j++) {
            edgesCost[i][j].first = icost(edgesUse[i][j].first, capx);
            edgesCost[i][j].second = icost(edgesUse[i][j].second, capy);
        }
    }
    
    // set cost prefix sum
    for (int i = 0; i < gridx; i++) {
        for (int j = 0; j < gridy; j++) {
            prefCol[i][j+1] = prefCol[i][j] + edgesCost[i][j].second;
            if (prefCol[i][j+1] < 0)
                cout << "error prefcol\n";
        }
    }
    for (int j = 0; j < gridy; j++) {
        for (int i = 0; i < gridx; i++) {
            prefRow[i+1][j] = prefRow[i][j] + edgesCost[i][j].first;
            if (prefRow[i+1][j] < 0)
                cout << "error prefrow\n";
        }
    }

    int over = get_overflow(false);
    //cout << "Origin overflow: " << over << endl;
    //cout << "Origin WL: " << get_WL() << endl;
    //cout << "check WL: " << WL_checker() << endl;
    
    // try to flip each L shape net
    int over2 = over;
    bool stop = false;
    
    for (int t = 0; t < 10; t++) {
        // the net has the highest cost first
        for (int i = 1; i <= numNets; i++)
            nets[i].maxcost = get_maxcost(nets[i]);
        sort(nets.begin() + 1, nets.end(), [](Net &a, Net &b) {
            return a.maxcost > b.maxcost;
        });

        for (int i = 1; i <= numNets; i++) {
            if (nets[i].maxcost < 1 + h/2)
                break;

            int x0 = nets[i].x0, x1 = nets[i].x1, y0 = nets[i].y0, y1 = nets[i].y1;
            if (x0 == x1 || y0 == y1)
                continue;
            if (nets[i].paths.size() != 1)
                cout << "init L error\n";

            int x2 = nets[i].paths[0].first, y2 = nets[i].paths[0].second;
            if (reroute_L(x0, y0, x2, y2, x1, y1)) {
                int x3, y3;
                if (x2 == x0) {
                    x3 = x1; y3 = y0;
                } else {
                    x3 = x0; y3 = y1;
                }
                nets[i].paths[0].first = x3;
                nets[i].paths[0].second = y3;
            }
        }
        
        over2 = get_overflow(true);
        //cout << "overflow: " << over2 << endl;
        //cout << "check WL: " << WL_checker() << endl;
        
        // stop condition
        if (over2 < over) {
            over = over2;
            if (stop)
                break;
        } else {
            over = over2;
            stop = true;
        }
    }
}

// update cost between a straight line (x0, y0) (x1, y1)
void update_cost(int x0, int y0, int x1, int y1) {
    if (x0 == x1) {
        forYRange(y0, y1, [&](int y) { 
            edgesCost[x0][y].second = icost(edgesUse[x0][y].second, capy); });
    } else {
        forXRange(x0, x1, [&](int x) { 
            edgesCost[x][y0].first = icost(edgesUse[x][y0].first, capx); });
    }
}

// rid up a net
void rid_up_net(Net &net) {
    if (net.paths.size() == 0) {
        rid_up_line(net.x0, net.y0, net.x1, net.y1);
        update_cost(net.x0, net.y0, net.x1, net.y1);
    } else {
        vector<pair<int, int>> &paths = net.paths;
        rid_up_line(net.x0, net.y0, paths[0].first, paths[0].second);
        update_cost(net.x0, net.y0, paths[0].first, paths[0].second);
        rid_up_line(net.x1, net.y1, 
                paths[paths.size() - 1].first, paths[paths.size() - 1].second);
        update_cost(net.x1, net.y1, 
                paths[paths.size() - 1].first, paths[paths.size() - 1].second);
        for (int j = 1; j < (int)paths.size(); j++) {
            rid_up_line(paths[j - 1].first, paths[j - 1].second, 
                        paths[j].first, paths[j].second);
            update_cost(paths[j - 1].first, paths[j - 1].second, 
                        paths[j].first, paths[j].second);
        }
    }
}

// put a net
void put_net(Net &net) {
    if (net.paths.size() == 0) {
        put_line(net.x0, net.y0, net.x1, net.y1);
        update_cost(net.x0, net.y0, net.x1, net.y1);
    } else {
        vector<pair<int, int>> &paths = net.paths;
        put_line(net.x0, net.y0, paths[0].first, paths[0].second);
        update_cost(net.x0, net.y0, paths[0].first, paths[0].second);
        put_line(net.x1, net.y1, 
                paths[paths.size() - 1].first, paths[paths.size() - 1].second);
        update_cost(net.x1, net.y1, 
                paths[paths.size() - 1].first, paths[paths.size() - 1].second);
        for (int j = 1; j < (int)paths.size(); j++) {
            put_line(paths[j - 1].first, paths[j - 1].second, 
                        paths[j].first, paths[j].second);
            update_cost(paths[j - 1].first, paths[j - 1].second, 
                        paths[j].first, paths[j].second);
        }
    }
}


bool inline in_grid(int x, int y) {
    return x >= 0 && y >= 0 && x < gridx && y < gridy;
}


int inline dist(int x0, int y0, int x1, int y1) {
    return max(abs(x0 - x1), abs(y0 - y1));
}

// trace back A*-search routing to set paths
void trace(Net &net, Node *node) {
    if (node->px != net.x0 || node->py != net.y0)
        trace(net, visited[node->px][node->py]);
    
    net.paths.push_back({node->x, node->y});
    
}

// reduce redundant node on a line 
void reduce_node(Net &net) {
    if (net.paths.empty())
        return;

    if (net.paths.size() == 1) {
        if (net.x0 == net.x1 || net.y0 == net.y1)
            net.paths.clear();
        return;
    }

    vector<pair<int,int>> result;
    result.reserve(net.paths.size());
    int x0 = net.x0;
    int y0 = net.y0;
    int x1 = net.paths[0].first;
    int y1 = net.paths[0].second;
    int x2 = net.paths[1].first;
    int y2 = net.paths[1].second;
    bool sameDir = (x0 == x1 && x1 == x2) || (y0 == y1 && y1 == y2);
    if (!sameDir) {
        result.push_back(net.paths[0]);
    }

    for (int i = 1; i < (int)net.paths.size() - 1; i++) {
        x0 = net.paths[i - 1].first;
        y0 = net.paths[i - 1].second;
        x1 = net.paths[i].first;
        y1 = net.paths[i].second;
        x2 = net.paths[i + 1].first;
        y2 = net.paths[i + 1].second;
        sameDir = (x0 == x1 && x1 == x2) || (y0 == y1 && y1 == y2);
        if (!sameDir) {
            result.push_back(net.paths[i]);
        }
    }
    x0 = net.paths[net.paths.size() - 2].first;
    y0 = net.paths[net.paths.size() - 2].second;
    x1 = net.paths[net.paths.size() - 1].first;
    y1 = net.paths[net.paths.size() - 1].second;
    x2 = net.x1;
    y2 = net.y1;
    sameDir = (x0 == x1 && x1 == x2) || (y0 == y1 && y1 == y2);
    if (!sameDir) {
        result.push_back(net.paths[net.paths.size() - 1]);
    }
    net.paths = result;
}

// A*-search routing to set a net's paths
void A_search_routing(Net &net) {

    priority_queue<Node, vector<Node>, greater<Node>> pq;
    for (int i = 0; i < gridx; i++) {
        for (int j = 0; j < gridy; j++) {
            if (visited[i][j] == NULL)
                continue;
            delete visited[i][j];
            visited[i][j] = NULL;
        }
    }

    net.paths.clear();
    visited[net.x0][net.y0] = new Node(net.x0, net.y0, -1, -1, 0, 0);

    if (in_grid(net.x0 + 1, net.y0) && visited[net.x0 + 1][net.y0] == NULL)
        pq.push(Node(net.x0 + 1, net.y0, net.x0, net.y0,
                edgesCost[net.x0][net.y0].first, mylambda * dist(net.x0 + 1, net.y0, net.x1, net.y1)));
    if (in_grid(net.x0, net.y0 + 1) && visited[net.x0][net.y0 + 1] == NULL)
        pq.push(Node(net.x0, net.y0 + 1, net.x0, net.y0,
                edgesCost[net.x0][net.y0].second, mylambda * dist(net.x0, net.y0 + 1, net.x1, net.y1)));
    if (in_grid(net.x0 - 1, net.y0) && visited[net.x0 - 1][net.y0] == NULL)
        pq.push(Node(net.x0 - 1, net.y0, net.x0, net.y0,
                edgesCost[net.x0 - 1][net.y0].first, mylambda * dist(net.x0 - 1, net.y0, net.x1, net.y1)));
    if (in_grid(net.x0, net.y0 - 1) && visited[net.x0][net.y0 - 1] == NULL)
        pq.push(Node(net.x0, net.y0 - 1, net.x0, net.y0,
                edgesCost[net.x0][net.y0 - 1].second, mylambda * dist(net.x0, net.y0 - 1, net.x1, net.y1)));
    
    while (!pq.empty()) {
        Node node = pq.top();
        pq.pop();
        if (visited[node.x][node.y] != NULL)
            continue;
    
        visited[node.x][node.y] = new Node(node.x, node.y, node.px, node.py, node.costo, node.costd);

        if (node.x == net.x1 && node.y == net.y1) {
            if (node.px != net.x0 || node.py != net.y0)
                trace(net, visited[node.px][node.py]);
            break;
        }
        int x0 = node.x, y0 = node.y;

        if (in_grid(x0 + 1, y0) && visited[x0 + 1][y0] == NULL)
            pq.push(Node(x0 + 1, y0, x0, y0,
                    edgesCost[x0][y0].first, mylambda * dist(x0 + 1, y0, net.x1, net.y1)));
        if (in_grid(x0, y0 + 1) && visited[x0][y0 + 1] == NULL)
            pq.push(Node(x0, y0 + 1, x0, y0,
                    edgesCost[x0][y0].second, mylambda * dist(x0, y0 + 1, net.x1, net.y1)));
        if (in_grid(x0 - 1, y0) && visited[x0 - 1][y0] == NULL)
            pq.push(Node(x0 - 1, y0, x0, y0,
                    edgesCost[x0 - 1][y0].first, mylambda * dist(x0 - 1, y0, net.x1, net.y1)));
        if (in_grid(x0, y0 - 1) && visited[x0][y0 - 1] == NULL)
            pq.push(Node(x0, y0 - 1, x0, y0,
                    edgesCost[x0][y0 - 1].second, mylambda * dist(x0, y0 - 1, net.x1, net.y1)));
    }
}


void main_stage() {
    
    visited = vector<vector<Node*>>(gridx, vector<Node*>(gridy, NULL));
    int over = 2147483647;
    int over2;
    for (int t = 0; t < 10; t++) {
        for (int i = 0; i < gridx; i++) {
            for (int j = 0; j < gridy; j++) {
                edgesCost[i][j].first = icost(edgesUse[i][j].first, capx);
                edgesCost[i][j].second = icost(edgesUse[i][j].second, capy);
            }
        }
        // rid up and reroute each net from the net with the highest cost on one edge
        for (int i = 1; i <= numNets; i++)
            nets[i].maxcost = get_maxcost(nets[i]);
        sort(nets.begin() + 1, nets.end(), [](Net &a, Net &b) {
            return a.maxcost > b.maxcost;
        });
        //cout << nets[1].maxcost << endl;
        
        double bound = (nets[1].maxcost - 1) * 0.6 + 1;
        int count = 0;
        for (int i = 1; i <= numNets; i++) {
            if (nets[i].maxcost < bound)
                break;
            
            count++;
            rid_up_net(nets[i]);
            A_search_routing(nets[i]);
            put_net(nets[i]);
        }
        over2 = get_overflow(true);

        // stop condition
        if (over2 == 0) {
            over = over2;
            break;
        } else if (over2 < over) {
            over = over2;
            for (int i = 1; i <= numNets; i++) {
                nets[i].ppaths.clear();
                nets[i].ppaths = nets[i].paths;
            }
        } else {
            for (int i = 1; i <= numNets; i++) {
                nets[i].paths = nets[i].ppaths;
            }
            over = get_overflow(true);
            break;
        }

        //cout << "reroute " << count << endl;
        //cout << "overflow: " << over << endl;
        //cout << "WL: " << get_WL() << endl;
    }
}

void write_file(string &outputfile) {
    ofstream out(outputfile);
    sort(nets.begin() + 1, nets.end(), [](Net &a, Net &b) {
        return a.idx < b.idx;
    });
    
    int WL = get_WL();
    out << "Wirelength " << WL << endl;
    for (int i = 1; i <= numNets; i++) {
        reduce_node(nets[i]);
        out << "Net " << nets[i].name << endl;
        if (nets[i].paths.empty()) {
            out << "Segment " << nets[i].x0 << " " << nets[i].y0 << " ";
            out << nets[i].x1 << " " << nets[i].y1 << endl;
        } else {
            out << "Segment " << nets[i].x0 << " " << nets[i].y0 << " ";
            out << nets[i].paths[0].first << " " << nets[i].paths[0].second << endl;
            for (int j = 1; j < (int)nets[i].paths.size(); j++) {
                out << "Segment " << nets[i].paths[j - 1].first << " " << nets[i].paths[j - 1].second << " ";
                out << nets[i].paths[j].first << " " << nets[i].paths[j].second << endl;
            }
            out << "Segment " << nets[i].paths[nets[i].paths.size() - 1].first << " ";
            out << nets[i].paths[nets[i].paths.size() - 1].second << " ";
            out << nets[i].x1 << " " << nets[i].y1 << endl;
        }
    }
}


int main(int argc, char* argv[]) {

    string inputfile = argv[1];
    string outputfile = argv[2];
    bWL = 2147483647;
    bov = 2147483647;

    read_file(inputfile);

    for (h = 30; h <= 70; h += 10) {
        for (k = 0.4; k <= 0.8; k += 0.1) {
            int bound = 4;
            if (numNets >= 200000)
                bound--;
            if (numNets >= 100000)
                bound--;
            for (mylambda = 1; mylambda <= bound; mylambda++) {
                initial_stage();
                cout << "initial stage end" << endl;

                main_stage();
                cout << "main stage end" << endl;

                ov = get_overflow(true);
                WL = get_WL();

                cout << "overflow: " << ov << endl;
                cout << "WL: " << WL << endl;

                if (ov < bov) {
                    bov = ov;
                    bWL = WL;
                    for (int i = 1; i <= numNets; i++) {
                        nets[i].bpaths.clear();
                        nets[i].bpaths = nets[i].paths;
                    }
                } else if (ov == bov) {
                    if (WL < bWL) {
                        bWL = WL;
                        for (int i = 1; i <= numNets; i++) {
                            nets[i].bpaths.clear();
                            nets[i].bpaths = nets[i].paths;
                        }
                    }
                }
                for (int i = 1; i <= numNets; i++) {
                    nets[i].paths.clear();
                    nets[i].ppaths.clear();
                }
            }
        }
    }

    for (int i = 1; i <= numNets; i++) {
        nets[i].paths = nets[i].bpaths;
    }

    ov = get_overflow(true);
    WL = get_WL();
    cout << "boverflow: " << ov << endl;
    cout << "bWL: " << WL << endl;

    write_file(outputfile);
    cout << "write end" << endl;

    return 0;
}