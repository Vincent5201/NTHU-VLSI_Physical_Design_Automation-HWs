#ifndef UTILS_H
#define UTILS_H

long long calculate_WL();

void sortNet(Net &net);

long long calculate_comp_WL(Component *comp);

long long calculate_WL();

Area find_opt(Component *c);

int find_space(Row *row, int direc, vector<Component*>::iterator it);

void put_comp(Component *c, Row *tgtrow, int tx1, bool insert);

long long check_WL_on_Row(Component* c, int tx1, Row *tgtrow, long long oWL);

vector<Component*>::iterator findc(vector<Component*>& comps, int tgt);

void insertc(vector<Component*>& comps, Component* c);

void removec(vector<Component*>& comps, Component* tgt);


#endif
