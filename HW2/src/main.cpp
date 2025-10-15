#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include <unordered_map>

using namespace std;


struct Cell {
    int id;
    int size;
    int gain = 0;
    int part = 0;                   // cell's group
    int parent1 = -1;               // for coarsening
    int parent2 = -1;               // for coarsening
    bool locked = false;
    unordered_set<int> pins;        // nets connected to
    Cell(int i, int s) : id(i), size(s) {}
};

struct BucketH {
    BucketH* down = NULL;
    BucketH* up = NULL;
    unordered_set<int> s;
    int size;
    BucketH(int size) : size(size) {}
};

int way, tgtway, cellNum, netNum, maxPin, totalSize, epochs, level;
vector<int> waySize;                            // size of each way
vector<vector<Cell*>> cells;                    // cells' infomation ([level][cellid])
vector<vector<BucketH*>> buckets;               // bucket list ([group][gain])
vector<BucketH*> maxCell;                       // max gain cell of each group
vector<vector<unordered_set<int>>> net_way;     // cells on the net ([netId][group])
vector<vector<bool>> netLock;                   // any locked cell on net
vector<int> netNums;


void read_input(string& inputFile) {

    ifstream infile(inputFile);
    if (!infile) {
        std::cerr << "Cannot open input file " << inputFile << "\n";
        return;
    }

    int state = 0, count = 0;
    string line;

    while (getline(infile, line)) {
        size_t pos = line.find(' ');
        string s = line.substr(0, pos);
        if (s == "NumCells") {
            state = 0;
            count = 0;
            s = line.substr(pos + 1);
            cellNum = stoi(s);
            cells = vector<vector<Cell*>>(10, vector<Cell*>(cellNum + 1, NULL));

        } else if (s == "Cell") {
            if (state == 0) {
                count++;
                s = line.substr(pos + 1);
                pos = s.find(' ');
                s = s.substr(pos + 1);
                int size = stoi(s);
                cells[0][count] = new Cell(count, size);

            } else if (state == 1) {
                s = line.substr(pos + 2);
                int cellId = stoi(s);
                cells[0][cellId]->pins.insert(count);
            }
        } else if (s == "NumNets") {
            state = 1;
            count = 0;
            s = line.substr(pos + 1);
            netNum = stoi(s);
            netLock = vector<vector<bool>>(netNum + 1);
            netNums = vector<int>(netNum + 1);
            maxPin = netNum;
        } else if (s == "Net") {
            count++;
            s = line.substr(pos + 1);
            pos = s.find(' ');
            s = s.substr(pos + 1);
            netNums[count] = stoi(s);
        }
    }
    infile.close();

    net_way = vector<vector<unordered_set<int>>>(netNum + 1, vector<unordered_set<int>>(way));
    buckets = vector<vector<BucketH*>>(way, vector<BucketH*>(2 * maxPin + 1, NULL));

    return;
}


void reset_net_way() {
    for (auto &vec : net_way)
        for (auto &uset : vec)
            uset.clear();
    
    for (int i = 1; i <= cellNum; i++) {
        if (!cells[level][i])
            break;
        for (int netId : cells[level][i]->pins)
            net_way[netId][cells[level][i]->part].insert(i);
    }
}


int calculate_cutsize(int p1, int p2, bool reset) {
    // 2-way version

    if (reset)
        reset_net_way();

    int ret = 0;
    for (int i = 1; i <= netNum; i++) {
        if (net_way[i][p1].size() > 0 && net_way[i][p2].size() > 0)
            ret++;
    }
    return ret;
}


int calculate_cutsize4(bool reset) {
    // 4-way version

    if (reset)
        reset_net_way();

    int ret = 0;
    for (int i = 1; i <= netNum; i++) {
        int count = 0;
        for (int j = 0; j < 4; j++) {
            if ((int)net_way[i][j].size() > 0)
                count++;
        }
        if (count > 1)
            ret++;
    }
    return ret;
}


void write_outupt(string& outFile) {

    int finalCut = 0;
    if (way == 2) {
        finalCut = calculate_cutsize(0, 1, true);
    } else  if (way == 4) {
        finalCut = calculate_cutsize4(true);
    }
    
    vector<set<int>> result(way);
    for (int i = 1; i <= cellNum; i++) {
        Cell *c = cells[0][i];
        result[c->part].insert(c->id);
    }

    ofstream file(outFile);
    file << "CutSize " << finalCut << endl;

    file << "GroupA " << result[0].size() << endl;
    for (int c : result[0])
        file << "C" << c << endl;
    
    file << "GroupB " << result[1].size() << endl;
    for (int c : result[1])
        file << "C" << c << endl;

    if (way == 4) {
        file << "GroupC " << result[2].size() << endl;
        for (int c : result[2])
            file << "C" << c << endl;
        
        file << "GroupD " << result[3].size() << endl;
        for (int c : result[3])
            file << "C" << c << endl;
    }
    file.close();
    return;
}


void initial_partition(int srcPart, int dstPart, bool reset) {
    // 2-way version

    for (auto &vec : net_way)
        for (auto &uset : vec)
            uset.clear();
    waySize = vector<int>(way, 0);
    totalSize = 0;

    // partition and setting
    int count = 1;
    for (int i = 1; i <= cellNum; i++) {
        Cell *c = cells[level][i];
        if (!c)
            break;

        if (reset) {
            if (c->part == srcPart) {
                c->gain = 0;
                count++;
                if (count % 2)
                    c->part = dstPart;

                totalSize += c->size;
                waySize[c->part] += c->size;
                for (int netId : c->pins)
                    net_way[netId][c->part].insert(i);
            }
        } else {
            c->gain = 0;
            totalSize += c->size;
            waySize[c->part] += c->size;
            for (int netId : c->pins)
                net_way[netId][c->part].insert(i);
        }
    }
    return;
}


void coarsening() {
    // edge coarsening

    // calculate times of coarsening
    int times = 0;
    int tmpCellNum = cellNum;
    while (tmpCellNum > 1000) {
        tmpCellNum /= 2;
        times++;
    }
    times = min(times, 9);

    while (times > 0) {
        times--;
        int nlevel = level + 1;
        int count = 1;
        vector<bool> seleted(cellNum + 1, false);               // whether the cell is merged
        vector<vector<int>> net_cells(netNum + 1);              // cells on the net
        unordered_map<long long, double> edgeWeight;            // edge weight between two cells
        vector<vector<pair<int, double>>> adj(cellNum + 1);     // adjacency list of cells

        // build net_cells
        for (int i = 1; i <= cellNum; i++) {
            if (!cells[level][i])
                break;
            for (int netId : cells[level][i]->pins)
                net_cells[netId].push_back(i);
        }

        // calculate edge weight net by net and ignore large nets
        for (int i = 1; i <= netNum; i++) {
            int size = net_cells[i].size();
            if (size > 100)
                continue;
            for (int a = 0; a < size; a++) {
                for (int b = a + 1; b < size; b++) {
                    long long key = ((long long)net_cells[i][a] << 32) | net_cells[i][b];
                    edgeWeight[key] += (double)1 / (size - 1);
                }
            }
        }

        // build adjacency list
        for (auto &u : edgeWeight) {
            int a = u.first >> 32;
            int b = u.first & 0xffffffff;
            adj[a].push_back({b, u.second});
            adj[b].push_back({a, u.second});
        }

        // merge cells
        for (int i = 1; i <= cellNum; i++) {
            if (!cells[level][i])
                break;
            if (seleted[i])
                continue;
            
            int tgt = -1;
            double maxw = -1;
            for (pair<int, double> &p : adj[i]) {
                if (seleted[p.first])
                    continue;

                if (p.second > maxw) {
                    maxw = p.second;
                    tgt = p.first;
                }
            }

            seleted[i] = true;
            if (tgt != -1) {
                // merge i and tgt
                seleted[tgt] = true;
                cells[nlevel][count] = new Cell(count, cells[level][i]->size + cells[level][tgt]->size);
                cells[nlevel][count]->parent1 = i;
                cells[nlevel][count]->parent2 = tgt;
                cells[nlevel][count]->pins.insert(cells[level][i]->pins.begin(),
                                                  cells[level][i]->pins.end());
                cells[nlevel][count]->pins.insert(cells[level][tgt]->pins.begin(),
                                                    cells[level][tgt]->pins.end());
            } else {
                // send i to next level directly
                cells[nlevel][count] = new Cell(count, cells[level][i]->size);
                cells[nlevel][count]->parent1 = i;
                cells[nlevel][count]->pins.insert(cells[level][i]->pins.begin(),
                                                  cells[level][i]->pins.end());
            }
            count++;
        }
        level++;
        cout << "level" << level << " left " << count - 1 << endl;
    }
    return;
}


void initial_buckets(int p1, int p2) {
    // 2-way version
    for (int p = 0; p < way; p++)
        for (int i = 0; i < 2 * maxPin + 1; i++)
            buckets[p][i] = new BucketH(0);
    
    maxCell = vector<BucketH*>(way, NULL);
    
    // calculate gain of all cells
    for (int i = 1; i <= netNum; i++) {
        if (net_way[i][p1].size() == 1) {
            int cellId = *net_way[i][p1].begin();
            cells[level][cellId]->gain++;
        } else if (net_way[i][p1].size() == 0) {
            for (int cellId : net_way[i][p2])
                cells[level][cellId]->gain--;
        }
        if (net_way[i][p2].size() == 1) {
            int cellId = *net_way[i][p2].begin();
            cells[level][cellId]->gain++;
        } else if (net_way[i][p2].size() == 0) {
            for (int cellId : net_way[i][p1])
                cells[level][cellId]->gain--;
        }
    }
    
    // put cells into buckets
    for (int i = 1; i <= cellNum; i++) {
        Cell* c = cells[level][i];
        if (!c)
            break;
        if (c->part == p1 || c->part == p2) {
            buckets[c->part][c->gain + maxPin]->s.insert(c->id);
            buckets[c->part][c->gain + maxPin]->size++;
        }
        
    }
    
    // build bucket list
    BucketH* tmp = NULL;
    for (int i = 0; i < (int)buckets[p1].size(); i++) {
        BucketH* b = buckets[p1][i];
        if (b->size) {
            b->down = tmp;
            if (tmp)
                tmp->up = b;
            tmp = b;
        }
    }
    maxCell[p1] = tmp;
    
    tmp = NULL;
    for (int i = 0; i < (int)buckets[p2].size(); i++) {
        BucketH* b = buckets[p2][i];
        if (b->size) {
            b->down = tmp;
            if (tmp)
                tmp->up = b;
            tmp = b;
        }
    }
    maxCell[p2] = tmp;

    return;
}


bool valid_size(int part) {
    // check whether the part size is valid

    double ratio = (double)waySize[part] / totalSize;
    if (tgtway == 2 && (ratio < 0.45 || ratio > 0.55))
        return false;
    if (tgtway == 4 && (ratio < 0.225 || ratio > 0.275))
        return false;

    return true;
}


pair<BucketH*, BucketH*> unlink_cell(int cellId, int direc) {
    // remove a cell from the bucket list and return its up and down bucket

    Cell *c = cells[level][cellId];
    BucketH *b = buckets[c->part][c->gain + maxPin];
    BucketH *bd = b->down, *bu = b->up;
    b->size--;

    b->s.erase(cellId);
    if (b->size > 0) {
        if (direc == 1) {
            bd = b;
        } else {
            bu = b;
        }
    } else {
        if (bu) {
            bu->down = bd;
        } else {
            maxCell[c->part] = bd;
        }
        if (bd)
            bd->up = bu;
    }

    if (direc == 0)
        return {NULL, NULL};

    return {bu, bd};
}


void insert_cell(int cellId, pair<BucketH*, BucketH*>& bb) {
    // insert a cell into bucket list

    Cell *c = cells[level][cellId];
    BucketH *b = buckets[c->part][c->gain + maxPin];
    (b->s).insert(cellId);
    b->size++;

    // the cell is the first one in the bucket, need to link the bucket
    if (b->size == 1) {
        BucketH *bu = bb.first;
        BucketH *bd = bb.second;
        if (bu == NULL && bd == NULL) {
            // find up and down bucket
            int idx = c->gain + maxPin + 1;
            while (idx < 2 * maxPin + 1 && buckets[c->part][idx]->size == 0)
                idx++;
            if (idx < 2 * maxPin + 1)
                bu = buckets[c->part][idx];

            idx = c->gain + maxPin - 1;
            while (idx >= 0 && buckets[c->part][idx]->size == 0)
                idx--;
            if (idx >= 0)
                bd = buckets[c->part][idx];
        }
        
        b->up = bu;
        b->down = bd;
        if (bu) {
            bu->down = b;
        } else {
            maxCell[c->part] = b;
        }
        if (bd)
            bd->up = b;
    }
    
    return;
}


void inc_gain(int cellId) {
    if (cells[level][cellId]->locked)
        return;
    pair<BucketH*, BucketH*> bb = unlink_cell(cellId, 1);
    cells[level][cellId]->gain++;
    insert_cell(cellId, bb);
    return;
}


void dec_gain(int cellId) {
    if (cells[level][cellId]->locked)
        return;
    pair<BucketH*, BucketH*> bb = unlink_cell(cellId, -1);
    cells[level][cellId]->gain--;
    insert_cell(cellId, bb);
    return;
}


int choose_cell(int src) {
    BucketH *b = maxCell[src];
    if (!b)
        return -1;
    int cellId = *(b->s).begin();
    return cellId;
}


void move_cell(int cellId, int src, int dst, pair<BucketH*, BucketH*>& bb) {
    // 2-way version
    // move a cell from src to dst and update gain of related cells
    
    Cell *c = cells[level][cellId];
    waySize[src] -= c->size;
    waySize[dst] += c->size;
    c->part = dst;
    c->gain = 0;
    
    for (int netId : c->pins) {
        net_way[netId][src].erase(c->id);
        netLock[netId][dst] = true;
        
        if (!netLock[netId][src] || !netLock[netId][dst]) {
            // update gain
            if (net_way[netId][src].size() == 1) {
                inc_gain(*net_way[netId][src].begin());
            } else if (net_way[netId][src].size() == 0) {
                for (int cellId : net_way[netId][dst])
                    dec_gain(cellId);
            }
            if (net_way[netId][dst].size() == 1) {
                dec_gain(*net_way[netId][dst].begin());
            } else if (net_way[netId][dst].size() == 0) {
                for (int cellId : net_way[netId][src])
                    inc_gain(cellId);
            }
        }
        net_way[netId][dst].insert(c->id);
    }
    return;
}


void reconstruct(vector<int> steps, int count, int p1, int p2) {
    // 2-way version
    // steps[count] is the best step, so reverse steps after count

    for (int i = count + 1; i < cellNum; i++) {
        if (steps[i] == -1)
            break;
        
        if (cells[level][steps[i]]->part == p1) {
            cells[level][steps[i]]->part = p2;
        } else if (cells[level][steps[i]]->part == p2){
            cells[level][steps[i]]->part = p1;
        }
    }

    // reconstruct buckets
    initial_partition(p1, p2, false);
    initial_buckets(p1, p2);
    return;
}


void balance_buckets(int p1, int p2, int max_iter) {
    // 2-way FM between group p1 and p2

    int cut = calculate_cutsize(p1, p2, true);
    int lastCut = -1, bestCut = -1, bestStep = -1;
    int iter = 0;
    vector<int> steps(cellNum);
    
    while (max_iter > iter) {
        int cellCount = 0;
        bestStep = -1;
        fill(netLock.begin(), netLock.end(), vector<bool>(way, false));
        fill(steps.begin(), steps.end(), -1);
        
        while (cellCount < cellNum - 1) {
            int src = p1, dst = p2;
            
            if (!maxCell[p1] && !maxCell[p2])
                break;

            if (!maxCell[p1]) {
                src = p2;
                dst = p1;
            } else if (!valid_size(p1) || !valid_size(p2)) {
                // move from large group to small group
                if (waySize[p1] < waySize[p2]) {
                    src = p2;
                    dst = p1;
                }
            } else if (maxCell[p2] && 
                cells[level][choose_cell(p1)]->gain < cells[level][choose_cell(p2)]->gain) {
                // choose the larger gain
                src = p2;
                dst = p1;
            }

            int cellId = choose_cell(src);
            if (cellId == -1)
                break;
            pair<BucketH*, BucketH*> bb = unlink_cell(cellId, 0);
            cut -= cells[level][cellId]->gain;
            cells[level][cellId]->locked = true;
            move_cell(cellId, src, dst, bb);
            steps[cellCount] = cellId;

            //record best step
            if (valid_size(p1) && valid_size(p2)) {
                if (bestCut == -1 || cut < bestCut) {
                    bestCut = cut;
                    bestStep = cellCount;
                }
            }
            cellCount++;
        }

        for (int i = 1; i <= cellNum; i++) {
            if (!cells[level][i])
                break;
            cells[level][i]->locked = false;
        }

        reconstruct(steps, bestStep, p1, p2);
        
        if (bestCut < lastCut) {
            cut = bestCut;
            lastCut = bestCut;
            bestStep = 0;
        } else {
            break;
        }
        iter++;
    }
    cout << p1 << " " << p2 << " best cut " << bestCut << endl;
    return;
}


void FM_2way(int p1, int p2, bool reset) {
    initial_partition(p1, p2, reset);
    initial_buckets(p1, p2);
    balance_buckets(p1, p2, epochs);
}


void uncoarsening() {
    
    // first FM
    tgtway = 2;
    FM_2way(0, 1, true);
    if (way == 4) {
        FM_2way(0, 2, true);
        FM_2way(1, 3, true);
        FM_2way(0, 3, false);
        FM_2way(1, 2, false);
        FM_2way(2, 3, false);
    }

    // uncoarsening and FM
    while (level > 0) {
        int plevel = level - 1;
        for (int i = 1; i <= cellNum; i++) {
            Cell *c = cells[level][i];
            if (!c)
                break;
            cells[plevel][c->parent1]->part = c->part;
            if (c->parent2 != -1)
                cells[plevel][c->parent2]->part = c->part;
        }
        level--;
        tgtway = way;

        if (way == 2) {
            FM_2way(0, 1, false);
        } else if (way == 4) {
            FM_2way(0, 2, false);
            FM_2way(1, 3, false);
            FM_2way(0, 3, false);
            FM_2way(1, 2, false);
            FM_2way(0, 1, false);
            FM_2way(2, 3, false);
        }
    }
    return;
}

int main(int argc, char* argv[]) {

    // read args and setting
    if (argc != 4) {
        std::cerr << "argc != 4";
        return 1;
    }
    string inputFile = argv[1];
    string outputFile = argv[2];
    way = stoi(argv[3]);
    tgtway = 2;
    maxPin = 0;
    totalSize = 0;
    epochs = 100;
    level = 0;

    read_input(inputFile);

    coarsening();

    uncoarsening();
    
    // check size
    tgtway = way;
    totalSize = 0;
    waySize = vector<int>(way, 0);
    for (int i = 1; i <= cellNum; i++) {
        totalSize += cells[0][i]->size;
        waySize[cells[0][i]->part] += cells[0][i]->size;
    }
    for (int i = 0; i < way; i++) {
        cout << "way" << i << " size " << waySize[i] << endl;
        if (!valid_size(i))
            cout << "error size " << waySize[i] << endl;
    }

    // check cell
    int tmpsum = 0;
    vector<int> cellCount(way, 0);
    for (int i = 1; i <= cellNum; i++)
        cellCount[cells[0][i]->part]++;
    
    for (int i = 0; i < way; i++)
        tmpsum += cellCount[i];
    
    if (tmpsum != cellNum)
        cout << "error cell"<< endl;

    // performance
    int finalCut = 0;
    if (way == 2) {
        finalCut = calculate_cutsize(0, 1, true);
    } else  if (way == 4) {
        finalCut = calculate_cutsize4(true);
    }
    cout << "final cutsize: " << finalCut << endl;

    
    write_outupt(outputFile);
    
    return 0;
}