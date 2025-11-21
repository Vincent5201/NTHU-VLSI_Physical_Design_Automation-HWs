#include "deal_datas.h"
#include "utils.h"
#include "algos.h"

#include <iostream>
#include <math.h>
#include <algorithm>


void check_leg() {
    for (auto &it : rows) {
        Row *row = &it.second;
        for (int i = 1; i < (int)row->comps.size(); i++) {
            if (row->comps[i - 1]->x2 > row->comps[i]->x1) {
                cout << row->comps[i - 1]->name << "overlap" << row->comps[i]->name << endl;
                cout << row->comps[i - 1]->sx2 << "overlap" << row->comps[i]->sx1 << endl;
                cout << row->comps[i - 1]->x2 << "overlap" << row->comps[i]->x1 << endl;
            }
        }
    }
}

// find an empty space for c on tgtrow
pair<int, long long> find_empty_on_Row(Component* c, Row *tgtrow, long long oWL) {
    
    // orientation setting
    if (c->ori != tgtrow->ori) {
        if (macros[c->type].syms.first == false)
            return {-1, 1};
        if (c->ori != "N" && c->ori != "FS")
            return {-1, 1};
        if (tgtrow->ori != "N" && tgtrow->ori != "FS")
            return {-1, 1};
    }
    
    // set left bound and right bound
    int t1 = c->opt.x1, t2 = c->opt.x2;
    t2 = min(t2, tgtrow->x + tgtrow->numX * tgtrow->stepX);
    t1 = min(t1, t2);
    vector<Component*>::iterator rit = findc(tgtrow->comps, t2);
    vector<Component*>::iterator lit = findc(tgtrow->comps, t1);
    if (lit != tgtrow->comps.begin())
        lit--;

    int tgtsize = c->sx2 - c->sx1; 
    int tgtsx1 = -1;            // selected site
    long long bWL = 0;          // best WL
    Component *lc;
    while (rit != lit && lit != tgtrow->comps.end()) {
        lc = *lit;
        lit++;
        
        int esize;      // empty size on tgtrow
        if (lit == tgtrow->comps.end()) {
            esize = tgtrow->numX - lc->sx2;
        } else {
            esize = (*lit)->sx1 - lc->sx2;
        }
        
        // if space is large enough, check the WL benefit
        if (tgtsize - esize <= find_space(tgtrow, 1, lit)) {
            int tx1 = lc->sx2;
            long long dWL = check_WL_on_Row(c, tx1, tgtrow, oWL);
            dWL += max((tgtsize - esize), 0) * tgtrow->stepX;
            if (dWL < bWL) {
                bWL = dWL;
                tgtsx1 = tx1;
                return {tgtsx1, bWL};
            }
        }
    }

    return {tgtsx1, bWL};
}

// find a component to swap with c on tgtrow
pair<pair<Component*, long long>, pair<int, int>> find_swap_on_Row(Component* c, Row *tgtrow, long long oWL) {
    
    if (c->ori != tgtrow->ori) {
        if (macros[c->type].syms.first == false)
            return {{nullptr, 1}, {-1, -1}};
        if (c->ori != "N" && c->ori != "FS")
            return {{nullptr, 1}, {-1, -1}};
        if (tgtrow->ori != "N" && tgtrow->ori != "FS")
            return {{nullptr, 1}, {-1, -1}};
    }

    // find left/right bound on tgtrow
    int t1 = c->opt.x1, t2 = c->opt.x2;
    t2 = min(t2, tgtrow->x + tgtrow->numX * tgtrow->stepX);
    t1 = min(t1, t2);
    vector<Component*>::iterator trit = findc(tgtrow->comps, t2);
    vector<Component*>::iterator tlit = findc(tgtrow->comps, t1);
    if (tlit != tgtrow->comps.begin())
        tlit--;

    // find left/right component of c
    Row *orow = ytoRow[c->y1];
    vector<Component*>::iterator oit = findc(orow->comps, c->x1);
    vector<Component*>::iterator olit = oit;
    if (oit != orow->comps.begin())
        olit--;
    vector<Component*>::iterator orit = oit;
    if (oit != orow->comps.end()) {
        orit++;
    }
    /*else {
        cout<< "no find" << c->name << "on" << orow->y << "with" << c->y1 << endl;
    }*/
    

    // setting site and space size on c
    int osize, ox1;
    if (orow->comps.size() == 1) {
        osize = orow->numX;
        ox1 = (*oit)->sx1;
    } else if (oit == orow->comps.begin()) {
        osize = (*orit)->sx1;
        ox1 = 0;
    } else if (orit == orow->comps.end()) {
        osize = orow->numX - (*olit)->sx2;
        ox1 = (*olit)->sx2;
    } else {
        osize = (*orit)->sx1 - (*olit)->sx2;
        ox1 = (*olit)->sx2;
    }

    Component *bc = nullptr;                // component to swap
    long long tWL, otWL, toWL, bWL = 0;
    int finalox1 = -1, finaltx1 = -1;       // bc put on finalox1, c put on finaltx1
    while (tlit != trit && tlit != tgtrow->comps.end()) {
        if ((*tlit)->status == "FIXED") {
            tlit++;
            continue;
        }
        
        Component *tgtcomp = *tlit;
        tlit++;

        // tgtcomp's size > space on orow(c's row)
        if (tgtcomp->sx2 - tgtcomp->sx1 > osize || tgtcomp == c)
            continue;
        
        // setting site and space size on tgtrow
        int tx1, tsize;
        vector<Component*>::iterator tpit = tlit;
        if (tgtrow->comps.size() == 1) {
            tx1 = 0;
            tsize = tgtrow->numX;
        } else if (tlit == tgtrow->comps.end()) {
            tpit--;
            tpit--;
            tx1 = (*tpit)->sx2;
            tsize = tgtrow->numX - tx1;
        } else {
            tpit--;
            if (tpit == tgtrow->comps.begin()) {
                tx1 = 0;
                tsize = (*tlit)->sx1;
            } else {
                tpit--;
                tx1 = (*tpit)->sx2;
                tsize = (*tlit)->sx1 - tx1;
            }
        }

        // if space is large enough, check the WL benefit
        int rspace = find_space(tgtrow, 1, tlit);
        if ((c->sx2 - c->sx1) - tsize <= rspace) {

            tWL = calculate_comp_WL(tgtcomp);
            otWL = check_WL_on_Row(c, tx1, tgtrow, oWL);
            toWL = check_WL_on_Row(tgtcomp, ox1, orow, tWL);
            long long dWL = otWL + toWL;
            dWL += max(0, (c->sx2 - c->sx1) - tsize) * tgtrow->stepX;
            if (dWL < bWL) {
                bWL = dWL;
                bc = tgtcomp;
                finalox1 = ox1;
                finaltx1 = tx1;
                return {{bc, bWL}, {finalox1, finaltx1}};
            }
        }
    }

    return {{bc, bWL}, {finalox1, finaltx1}};
}

// for each component, find an empty space in optomal region to put it first
// if no proper empty space, find another component to swap with it
void global_swap() {

    for (auto &comp : components) {    
        if (comp.second.status == "FIXED")
            continue;

        // initial setting
        Component* c = &comp.second;
        Area opt = c->opt;
        if (c->y1 < opt.y2 && c->y1 >= opt.y1)
            continue;

        Row *orow = ytoRow[c->y1];
        long long oWL = calculate_comp_WL(c);
        long long bWL = 0;

        // find empty space
        auto it1 = ytoRow.lower_bound(opt.y1);
        auto it2 = ytoRow.lower_bound(opt.y2);
        int bsx1 = -1;
        Row *brow = nullptr;
        while (it1 != it2) {
            Row *tgtrow = it1->second;
            it1++;
            if (orow->name == tgtrow->name)
                continue;
            if (tgtrow->site != orow->site)
                continue;
            if (opt.x1 >= tgtrow->x + tgtrow->numX * tgtrow->stepX)
                continue;

            pair<int, long long> tmp = find_empty_on_Row(c, tgtrow, oWL);
            if (tmp.first == -1)
                continue;
            
            if (tmp.second < bWL) {
                bWL = tmp.second;
                bsx1 = tmp.first;
                brow = tgtrow;
                break;
            }
        }

        // put to empty space
        if (bsx1 != -1 && bsx1 + c->sx2 - c->sx1 < brow->numX) {
            bool check = true;
            vector<Component*>::iterator nc = findc(brow->comps, bsx1 * brow->stepX + brow->x);
            if (nc != brow->comps.end() && bsx1 + c->sx2 - c->sx1 > (*nc)->sx1) {
                int tmp = bsx1 + c->sx2 - c->sx1 - (*nc)->sx1;
                if (tmp + (*nc)->sx2 >= brow->numX)
                    check = false;
            }
            if (check) {
                removec(orow->comps, c);
                put_comp(c, brow, bsx1, true);
                vector<Component*>::iterator oc = findc(brow->comps, c->x1);
                oc++;
                if (oc != brow->comps.end() && (*oc)->sx1 < c->sx2)
                    put_comp(*oc, brow, c->sx2, false);
                continue;
            }
        }
        
        // find another component
        it1 = ytoRow.lower_bound(opt.y1);
        it2 = ytoRow.lower_bound(opt.y2);
        int finalox1, finaltx1;
        Component *finalcomp = nullptr;
        while (it1 != it2) {
            Row *tgtrow = it1->second;
            it1++;
            
            if (orow->name == tgtrow->name)
                continue;
            if (tgtrow->site != orow->site)
                continue;
            if (opt.x1 >= tgtrow->x + tgtrow->numX * tgtrow->stepX)
                continue;

            auto tmp = find_swap_on_Row(c, tgtrow, oWL);
            if (tmp.first.first == nullptr)
                continue;

            if (tmp.first.second < bWL) {
                finalcomp = tmp.first.first;
                bWL = tmp.first.second;            
                finalox1 = tmp.second.first;
                finaltx1 = tmp.second.second;
                break;
            }
        }

        // swap c and finalcomp
        if (finalcomp != nullptr) {
            Row *tgtrow = ytoRow[finalcomp->y1];
            if (finaltx1 + c->sx2 - c->sx1 >= tgtrow->numX)
                continue;
            if (finalox1 + finalcomp->sx2 - finalcomp->sx1 >= orow->numX)
                continue;

            bool check = true;
            vector<Component*>::iterator nc = findc(tgtrow->comps, finaltx1 * tgtrow->stepX + tgtrow->x);
            if (nc != tgtrow->comps.end() && finaltx1 + c->sx2 - c->sx1 > (*nc)->sx1) {
                int tmp = finaltx1 + c->sx2 - c->sx1 - (*nc)->sx1;
                if (tmp + (*nc)->sx2 >= tgtrow->numX)
                    check = false;
            }
            if (check) {
                removec(tgtrow->comps, finalcomp);
                removec(orow->comps, c);
                
                put_comp(c, tgtrow, finaltx1, true);
                put_comp(finalcomp, orow, finalox1, true);
                
                vector<Component*>::iterator oc = findc(tgtrow->comps, c->x1);
                oc++;
                if (oc != tgtrow->comps.end() && (*oc)->sx1 < c->sx2)
                    put_comp(*oc, tgtrow, c->sx2, false);
            }
        }
    }
}

// for each component, find an empty space on its upper/lower two row
// if no proper empty space, find another component to swap with it
void vertical_swap() {
    for (auto &comp : components) {
        if (comp.second.status == "FIXED")
            continue;

        // initial setting
        Component* c = &comp.second;
        Area opt = c->opt;
        if (c->y1 < opt.y2 && c->y1 >= opt.y1)
            continue;

        Row *orow = ytoRow[c->y1];
        long long oWL = calculate_comp_WL(c);
        long long bWL = 0;

        // find empty space
        auto it1 = ytoRow.lower_bound(opt.y2);
        auto it2 = ytoRow.lower_bound(opt.y1);
        if (c->y2 > opt.y2) {
            it2 = ytoRow.lower_bound(c->y2);
        } else if (c->y1 < opt.y1) {
            it1 = it2;
            it2 = ytoRow.lower_bound(c->y1);
        } else {
            continue;
        }
        int bsx1 = -1, count = 2;
        Row *brow = nullptr;
        while (it1 != it2 && count > 0) {
            count--;
            Row *tgtrow = it1->second;
            if (c->y2 > opt.y2) {
                it1++;
            } else if (c->y1 < opt.y1) {
                it1--;
            }

            if (orow->name == tgtrow->name)
                continue;
            if (tgtrow->site != orow->site)
                continue;
            if (opt.x1 >= tgtrow->x + tgtrow->numX * tgtrow->stepX) {
                if (c->x1 >= tgtrow->x + tgtrow->numX * tgtrow->stepX) {
                    continue;
                } else {
                    opt.x1 = (c->x1 + tgtrow->x + tgtrow->numX * tgtrow->stepX) / 2;
                }
            }
            pair<int, long long> tmp = find_empty_on_Row(c, tgtrow, oWL);
            if (tmp.first == -1)
                continue;
            
            if (tmp.second < bWL) {
                bWL = tmp.second;
                bsx1 = tmp.first;
                brow = tgtrow;
                break;
            }
        }

        // put to empty space
        if (bsx1 != -1 && bsx1 + c->sx2 - c->sx1 < brow->numX) {
            bool check = true;
             vector<Component*>::iterator nc = findc(brow->comps, bsx1 * brow->stepX + brow->x);
            if (nc != brow->comps.end() && bsx1 + c->sx2 - c->sx1 > (*nc)->sx1) {
                int tmp = bsx1 + c->sx2 - c->sx1 - (*nc)->sx1;
                if (tmp + (*nc)->sx2 >= brow->numX)
                    check = false;
            }
            if (check) {
                removec(orow->comps, c);
                put_comp(c, brow, bsx1, true);
                vector<Component*>::iterator oc = findc(brow->comps, c->x1);
                oc++;
                if (oc != brow->comps.end() && (*oc)->sx1 < c->sx2)
                    put_comp(*oc, brow, c->sx2, false);
                continue;
            }
        }

        // find another component
        it1 = ytoRow.lower_bound(opt.y2);
        it2 = ytoRow.lower_bound(opt.y1);
        if (c->y2 > opt.y2) {
            it2 = ytoRow.upper_bound(c->y2);
        } else if (c->y1 < opt.y1) {
            it1 = it2;
            it2 = ytoRow.lower_bound(c->y1);
        }
        count = 2;
        int finalox1, finaltx1;
        Component *finalcomp = nullptr;
        while (it1 != it2 && count > 0) {
            count--;
            Row *tgtrow = it1->second;
            if (c->y2 > opt.y2) {
                it1++;
            } else if (c->y1 < opt.y1) {
                it1--;
            }

            if (orow->name == tgtrow->name)
                continue;
            if (tgtrow->site != orow->site)
                continue;
            if (opt.x1 >= tgtrow->x + tgtrow->numX * tgtrow->stepX)
                continue;
            auto tmp = find_swap_on_Row(c, tgtrow, oWL);
            if (tmp.first.first == nullptr)
                continue;
            if (tmp.first.second < bWL) {
                finalcomp = tmp.first.first;
                bWL = tmp.first.second;            
                finalox1 = tmp.second.first;
                finaltx1 = tmp.second.second;
                break;
            }
            
        }

        // swap c and finalcomp
        if (finalcomp != nullptr) {
            Row *tgtrow = ytoRow[finalcomp->y1];

            if (finaltx1 + c->sx2 - c->sx1 >= tgtrow->numX)
                continue;
            if (finalox1 + finalcomp->sx2 - finalcomp->sx1 >= orow->numX)
                continue;

            bool check = true;
            vector<Component*>::iterator nc = findc(tgtrow->comps, finaltx1 * tgtrow->stepX + tgtrow->x);
            if (nc != tgtrow->comps.end() && finaltx1 + c->sx2 - c->sx1 > (*nc)->sx1) {
                int tmp = finaltx1 + c->sx2 - c->sx1 - (*nc)->sx1;
                if (tmp + (*nc)->sx2 >= tgtrow->numX)
                    check = false;
            }
            if (check) {
                removec(tgtrow->comps, finalcomp);
                removec(orow->comps, c);
                
                put_comp(c, tgtrow, finaltx1, true);
                put_comp(finalcomp, orow, finalox1, true);

                vector<Component*>::iterator oc = findc(tgtrow->comps, c->x1);
                oc++;
                if (oc != tgtrow->comps.end() && (*oc)->sx1 < c->sx2)
                    put_comp(*oc, tgtrow, c->sx2, false);
            }
        }
    }
}

// reordering with window size = 3
void local_reodering() {
    
    for (auto &net : nets)
        sortNet(net.second);

    for (auto &it : rows) {
        
        Row *row = &it.second;
        Component *c1 = nullptr, *c2 = nullptr, *c3 = nullptr;

        for (int i = 0; i < (int)row->comps.size(); i++) {
            if (row->comps[i]->status == "FIXED") {
                c1 = c2 = c3 = nullptr;
                continue;
            }
            
            if (c1 == nullptr) {
                c1 = row->comps[i];
            } else if (c2 == nullptr) {
                c2 = row->comps[i];
            } else {
                c3 = row->comps[i];
                
                int bx1 = c1->x1, bx2 = c2->x1, bx3 = c3->x1;
                long long bWL = 0, WL = 0;
                int rbound = c3->x2, lbound = c1->x1;
                int lc1 = c1->x2 - c1->x1, lc2 = c2->x2 - c2->x1, lc3 = c3->x2 - c3->x1;
                
                // 1 2 3
                bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                bWL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                
                // 1 2 3
                c2->x1 = c1->x1 + lc1;
                WL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                if (WL < bWL) {
                    bWL = WL;
                    bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                }

                // 1 3 2
                c1->x1 = lbound;
                c3->x1 = c1->x1 + lc1;
                c2->x1 = rbound - lc2;
                WL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                if (WL < bWL) {
                    bWL = WL;
                    bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                }
                
                // 2 3 1
                c2->x1 = lbound;
                c3->x1 = c2->x1 + lc2;
                c1->x1 = rbound - lc1;
                WL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                if (WL < bWL) {
                    bWL = WL;
                    bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                }
                
                // 2 1 3
                c2->x1 = lbound;
                c1->x1 = c2->x1 + lc2;
                c3->x1 = rbound - lc3;
                WL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                if (WL < bWL) {
                    bWL = WL;
                    bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                }
                
                // 3 1 2
                c3->x1 = lbound;
                c1->x1 = c3->x1 + lc3;
                c2->x1 = rbound - lc2;
                WL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                if (WL < bWL) {
                    bWL = WL;
                    bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                }
                // 3 2 1
                c3->x1 = lbound;
                c2->x1 = c3->x1 + lc3;
                c1->x1 = rbound - lc1;
                WL = calculate_comp_WL(c1) + calculate_comp_WL(c2) + calculate_comp_WL(c3);
                if (WL < bWL) {
                    bWL = WL;
                    bx1 = c1->x1; bx2 = c2->x1; bx3 = c3->x1;
                }
                
                c1->x1 = bx1; c2->x1 = bx2; c3->x1 = bx3;
                c1->x2 = c1->x1 + lc1; c2->x2 = c2->x1 + lc2; c3->x2 = c3->x1 + lc3;
                
                if (c1 ->x1 < c2->x1 && c1->x1 < c3->x1) {
                    if (c2->x1 < c3->x1) {
                        c1 = c2;
                        c2 = c3;
                    } else {
                        c1 = c3;
                        //c2 = c2;
                    }
                } else if (c2 ->x1 < c3->x1 && c2->x1 < c1->x1) {
                    if (c1->x1 < c3->x1) {
                        //c1 = c1;
                        c2 = c3;
                    } else {
                        c2 = c1;
                        c1 = c3;
                    }
                } else {
                    if (c1->x1 < c2->x1) {
                        //c1 = c1;
                        //c2 = c2;
                    } else {
                        swap(c1, c2);
                    }
                }
                c3 = nullptr;
            }
        }
        sort(row->comps.begin(), row->comps.end(), [](Component *a, Component *b){
            return a->x1 < b->x1;
        });
    }
    
    for (auto &it : components) {
        if (it.second.status == "FIXED")
            continue;
        Component *c = &it.second;
        Row *row = ytoRow[c->y1];
        c->sx1 = (c->x1 - row->x) / row->stepX;
        c->sx2 = (c->x2 - row->x) / row->stepX;
    }
}

// for yi = f(xi), i=[1,3], solve an quadratic equation
int solveEq(int x1, long long y1, int x2, long long y2, int x3, long long y3) {
    double tmp = (x1 - x2) * (x1 - x3) * (x2 - x3);
    double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / tmp;
    double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / tmp;
    double c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / tmp;
    
    double vertexX = -b / (2 * a);
    int bestX = round(vertexX);
    int bestY = round(a * bestX * bestX + b * bestX + c);
    int tmpx1 = bestX - 1;
    int tmpx2 = bestX + 1;
    
    // find the bestX between x1 and x3
    if (tmpx1 >= x1 && round(a * tmpx1 * tmpx1 + b * tmpx1 + c) < bestY)
        bestX = tmpx1;

    if (tmpx2 <= x3 && round(a * tmpx2 * tmpx2 + b * tmpx2 + c) < bestY)
        bestX = tmpx2;

    return bestX;
}

// a simplified version
// for all components, find their best position individualy
void single_segment() {
    for (auto &net : nets)
        sortNet(net.second);

    for (auto &it : rows) {
        Row *row = &it.second;
        Component *c1 = nullptr, *c2 = nullptr;
        int lbound = 0;
        for (int i = 0; i < (int)row->comps.size(); i++) {
            if (row->comps[i]->status == "FIXED") {
                lbound = row->comps[i]->sx2;
                c1 = c2 = nullptr;
                continue;
            }
            
            if (c1 == nullptr) {
                c1 = row->comps[i];
            } else {
                c2 = row->comps[i];
                int bc1 = c1->x1;
                long long bWL = calculate_comp_WL(c1);

                int rbound = c2->sx1 - (c1->sx2 - c1->sx1);
                int c1size = c1->x2 - c1->x1;
                int x1, x2, x3;
                long long y1, y2, y3;
                x1 = lbound;
                c1->x1 = x1 * row->stepX + row->x;
                y1 = calculate_comp_WL(c1);
                if (y1 < bWL)
                    bc1 = c1->x1;

                x3 = rbound;
                c1->x1 = x3 * row->stepX + row->x;
                y3 = calculate_comp_WL(c1);
                if (y3 < bWL)
                    bc1 = c1->x1;

                x2 = (x1 + x3) / 2;
                if (x2 != x1) {
                    c1->x1 = x2 * row->stepX + row->x;
                    y2 = calculate_comp_WL(c1);
                    int bestX = solveEq(x1, y1, x2, y2, x3, y3);
                    if (bestX > x1 && bestX < x3) {
                        c1->x1 = bestX * row->stepX + row->x;
                        if (calculate_comp_WL(c1) < bWL)
                            bc1 = c1->x1;
                    }
                }
                c1->x1 = bc1;
                c1->x2 = c1->x1 + c1size;
                lbound = (c1->x2 - row->x) / row->stepX;
                c1 = c2;
                c2 = nullptr;
            }
        }
    }

    for (auto &it : components) {
        if (it.second.status == "FIXED")
            continue;
        Component *c = &it.second;
        Row *row = ytoRow[c->y1];
        c->sx1 = (c->x1 - row->x) / row->stepX;
        c->sx2 = (c->x2 - row->x) / row->stepX;
    }
}


void FastDP(int iters) {

    long long oWL = calculate_WL();
    long long tWL = oWL;
    cout << "lWL: " << tWL << endl;

    do {
        for (auto &net : nets)
            sortNet(net.second);
        for (auto &comp : components) {
            if (comp.second.status != "FIXED")
                comp.second.opt = find_opt(&components[comp.first]);
        }
        oWL = tWL;
        global_swap();
        //check_leg();
        cout << "end global" << endl;
        
        for (auto &net : nets)
            sortNet(net.second);
        for (auto &comp : components) {
            if (comp.second.status != "FIXED")
                comp.second.opt = find_opt(&components[comp.first]);
        }
        
        vertical_swap();
        //check_leg();
        cout << "end vertical" << endl;

        local_reodering();
        //check_leg();
        tWL = calculate_WL();
        cout << "lWL: " << tWL << endl;
        iters--;
        
    } while (iters > 0);

    
    single_segment();
    //check_leg();
    tWL = calculate_WL();
    cout << "sWL: " << tWL << endl;
}

void reordering(int iters) {
    long long oWL = calculate_WL();
    long long tWL = oWL;

    do {
        oWL = tWL;
        local_reodering();
        tWL = calculate_WL();
        cout << "lWL: " << tWL << endl;
        iters--;
    } while (iters > 0);

    single_segment();
    tWL = calculate_WL();
    cout << "sWL: " << tWL << endl;
}