#include "deal_datas.h"
#include "utils.h"

#include <algorithm>
#include <iostream>
#include <math.h>


void sortNet(Net &net) {
    sort(net.compsX.begin(), net.compsX.end(), [](Component *a, Component *b){
        if (a->x1 == b->x1)
            return a->y1 < b->y1;
        return a->x1 < b->x1;
    });

    sort(net.pinsX.begin(), net.pinsX.end(), [](Pin *a, Pin *b){
        if (a->x1 == b->x1)
            return a->y1 < b->y1;
        return a->x1 < b->x1;
    });

    sort(net.compsY.begin(), net.compsY.end(), [](Component *a, Component *b){
        if (a->y1 == b->y1)
            return a->x1 < b->x1;
        return a->y1 < b->y1;
    });

    sort(net.pinsY.begin(), net.pinsY.end(), [](Pin *a, Pin *b){
        if (a->y1 == b->y1)
            return a->x1 < b->x1;
        return a->y1 < b->y1;
    });
}

// find c or another component's x1 > c (comps is sorted)
vector<Component*>::iterator findc(vector<Component*>& comps, int tgt) {
    auto it = lower_bound(comps.begin(), comps.end(), tgt,
                            []( Component* c, int value) {
                                return  c->x1 < value;
                            });
    return it;
}

// insert c to comps (comps is sorted)
void insertc(vector<Component*>& comps, Component* c) {
    auto it = lower_bound(comps.begin(), comps.end(), c,
        [](Component* a, Component* b) {
            return a->x1 < b->x1;
        });
    comps.insert(it, c);
}


void removec(vector<Component*>& comps, Component* tgt) {
    auto it = find(comps.begin(), comps.end(), tgt);
    if (it != comps.end()) {
        comps.erase(it);
    } 
    /*else {
        cout << "no remove" << tgt->name << ytoRow[tgt->y1]->name << endl;
    }*/
}

// calculate WL of nets connected with comp
long long calculate_comp_WL(Component *comp) {
    
    long long ans = 0;
    for (Net *net : comp->nets) {
        sort(net->compsX.begin(), net->compsX.end(), [](Component *a, Component *b){
            if (a->x1 == b->x1)
                return a->y1 < b->y1;
            return a->x1 < b->x1;
        });
        sort(net->compsY.begin(), net->compsY.end(), [](Component *a, Component *b){
            if (a->y1 == b->y1)
                return a->x1 < b->x1;
            return a->y1 < b->y1;
        });
        pair<int, int> small, large;
        
        small = {INT32_MAX, INT32_MAX};
        large = {0, 0};
        
        if (!net->compsX.empty()) {
            small.first = min(small.first, net->compsX.front()->x1);
            small.second = min(small.second, net->compsY.front()->y1);
            large.first = max(large.first, net->compsX.back()->x1);
            large.second = max(large.second, net->compsY.back()->y1);
        }

        if (!net->pinsX.empty()) {
            small.first = min(small.first, net->pinsX.front()->x1);
            small.second = min(small.second, net->pinsY.front()->y1);
            large.first = max(large.first, net->pinsX.back()->x1);
            large.second = max(large.second, net->pinsY.back()->y1);
        }
        
        ans += large.first - small.first + large.second - small.second;
    }

    ans *= lefMicrons;
    ans /= defMicrons;
    return ans;
}

// calculate total WL
long long calculate_WL() {
    long long ans = 0;
    for (auto &net : nets) {
        sortNet(net.second);
        pair<int, int> small, large;
        
        small = {INT32_MAX, INT32_MAX};
        large = {0, 0};
        
        if (!net.second.compsX.empty()) {
            small.first = min(small.first, net.second.compsX.front()->x1);
            small.second = min(small.second, net.second.compsY.front()->y1);
            large.first = max(large.first, net.second.compsX.back()->x1);
            large.second = max(large.second, net.second.compsY.back()->y1);
        }

        if (!net.second.pinsX.empty()) {
            small.first = min(small.first, net.second.pinsX.front()->x1);
            small.second = min(small.second, net.second.pinsY.front()->y1);
            large.first = max(large.first, net.second.pinsX.back()->x1);
            large.second = max(large.second, net.second.pinsY.back()->y1);
        }
        
        ans += large.first - small.first + large.second - small.second;
    }

    ans *= lefMicrons;
    ans /= defMicrons;

    return ans;
}


Area find_opt(Component *c) {
    Area opt;
    vector<int> allx, ally;
    allx.reserve(c->nets.size() * 2);
    ally.reserve(c->nets.size() * 2);
    
    for (Net *net : c->nets) {
        int sx, sy, lx, ly;
        if (net->pinsX.empty()) {
            sx = net->compsX.front()->x1;
            sy = net->compsY.front()->y1;
            lx = net->compsX.back()->x1;
            ly = net->compsY.back()->y1;
        } else if (net->compsX.empty()) {
            sx = net->pinsX.front()->x1;
            sy = net->pinsY.front()->y1;
            lx = net->pinsX.back()->x1;
            ly = net->pinsY.back()->y1;
        } else {
            sx = min(net->pinsX.front()->x1, net->compsX.front()->x1);
            sy = min(net->pinsY.front()->y1, net->compsY.front()->y1);
            lx = max(net->pinsX.back()->x1, net->compsX.back()->x1);
            ly = max(net->pinsY.back()->y1, net->compsY.back()->y1);
        }
        allx.push_back(sx);
        allx.push_back(lx);
        ally.push_back(sy);
        ally.push_back(ly);
    }

    int mid = allx.size() / 2;
    nth_element(allx.begin(), allx.begin() + mid, allx.end());
    opt.x2 = allx[mid];

    nth_element(allx.begin(), allx.begin() + mid - 1, allx.begin() + mid);
    opt.x1 = allx[mid - 1];

    nth_element(ally.begin(), ally.begin() + mid, ally.end());
    opt.y2 = ally[mid];

    nth_element(ally.begin(), ally.begin() + mid - 1, ally.begin() + mid);
    opt.y1 = ally[mid - 1];

    return opt;
}

// check empty space on left/right of component (*it)
int find_space(Row *row, int direc, vector<Component*>::iterator it) {
    if (it == row->comps.end() || (*it)->status == "FIXED")
        return 0;

    vector<Component*>::iterator nit = it;
    if (direc == 1) {
        nit++;
        if (nit == row->comps.end())
            return row->numX - (*it)->sx2;
        return (*nit)->sx1 - (*it)->sx2;
    } else {
        if (nit == row->comps.begin())
            return (*nit)->sx1;
        nit--;
        return (*it)->sx1 - (*nit)->sx2;
    }
    return 0;
}

// put component on tgtrow's tx1
void put_comp(Component *c, Row *tgtrow, int tx1, bool insert) {
    c->y1 = tgtrow->y;
    c->y2 = tgtrow->y2;
    c->x1 = tx1 * tgtrow->stepX + tgtrow->x;
    c->x2 = c->x1 + macros[c->type].sizeX;
    c->sx1 = tx1;
    c->sx2 = (c->x2 - tgtrow->x) / tgtrow->stepX;
    c->ori = tgtrow->ori;
    if (insert)
        insertc(tgtrow->comps, c);
}

// put c to tgtrow's tx1 and check WL
long long check_WL_on_Row(Component* c, int tx1, Row *tgtrow, long long oWL) {
    int x1 = tx1 * tgtrow->stepX + tgtrow->x;
    int oy = ytoRow[c->y1]->y;
    swap(c->x1, x1);
    c->y1 = tgtrow->y;
    long long dWL = calculate_comp_WL(c) - oWL;
    swap(c->x1, x1);
    c->y1 = oy;
    return dWL;
}


