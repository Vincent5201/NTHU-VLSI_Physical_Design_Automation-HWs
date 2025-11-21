#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "deal_datas.h"

int lefMicrons;
unordered_map<string, Site> sites;
unordered_map<string, Macro> macros;

int defMicrons, compCount, pinCount, netCount;
Area dieArea;
unordered_map<string, Row> rows;
unordered_map<string, Component> components;
unordered_map<string, Pin> pins;
unordered_map<string, Net> nets;
map<int, Row*> ytoRow;


void read_LEF(string &inputLEF) {
    
    ifstream lef(inputLEF);
    string line;
    
    while (getline(lef, line)) {
        istringstream iss(line);
        string token, tmp;
        iss >> token;

        if (token == "DATABASE") {
            iss >> tmp >> lefMicrons;
        } else if (token == "SITE") {
            Site site;
            iss >> site.name;
            
            getline(lef, line);
            iss.clear();
            iss.str(line);
            iss >> tmp >> site.sClass;

            getline(lef, line);
            iss.clear();
            iss.str(line);
            iss >> tmp >> site.dsizeX >> tmp >> site.dsizeY;

            getline(lef, line);
            sites[site.name] = move(site);

        } else if (token == "MACRO") {
            Macro macro;
            iss >> macro.name;
            string tmp;

            getline(lef, line);
            iss.clear();
            iss.str(line);
            iss >> tmp >> macro.mClass;

            if (macro.mClass == "BLOCK") {
                getline(lef, line);
                getline(lef, line);

                getline(lef, line);
                iss.clear();
                iss.str(line);
                iss >> tmp >> macro.dsizeX >> tmp >> macro.dsizeY;

                getline(lef, line);
                iss.clear();
                iss.str(line);
                iss >> tmp;
                while (iss >> tmp) {
                    if (tmp == ";")
                        break;
                    if (tmp == "X") {
                        macro.syms.first = true;
                    } else if (tmp == "Y") {
                        macro.syms.second = true;
                    }
                }
            } else {
                getline(lef, line);
                getline(lef, line);

                getline(lef, line);
                iss.clear();
                iss.str(line);
                iss >> tmp >> macro.dsizeX >> tmp >> macro.dsizeY;

                getline(lef, line);
                iss.clear();
                iss.str(line);
                iss >> tmp;
                while (iss >> tmp) {
                    if (tmp == ";")
                        break;
                    if (tmp == "X") {
                        macro.syms.first = true;
                    } else if (tmp == "Y") {
                        macro.syms.second = true;
                    }
                }

                getline(lef, line);
                iss.clear();
                iss.str(line);
                iss >> tmp >> macro.siteName;
            }
            while (getline(lef, line)) {
                if (line.find(macro.name) != string::npos)
                    break;
            }
            macros[macro.name] = move(macro);
        }
    }
    lef.close();
    cout << "read LEF" << endl;
}


void read_DEF(string &inputDEF) {
    ifstream def(inputDEF);
    string line;

    while (getline(def, line)) {
        istringstream iss(line);
        string token, tmp;
        iss >> token;

        if (token == "UNITS") {
            iss >> tmp >> tmp >> defMicrons;
        } else if (token == "DIEAREA") {
            iss >> tmp >> dieArea.x1 >> dieArea.y1 >> tmp >> tmp >> dieArea.x2 >> dieArea.y2;

        } else if (token == "ROW") {
            Row row;
            iss >> row.name >> row.site >> row.x >> row.y >> row.ori;
            iss >> tmp >> row.numX >> tmp >> row.numY >> tmp >> row.stepX >> row.stepY;
            rows[row.name] = move(row);

        } else if (token == "COMPONENTS") {
            iss >> compCount;
            while (getline(def, line)) {
                if (line.find("END COMPONENTS") != string::npos)
                    break;

                if (line.find("-") != string::npos) {
                    iss.clear();
                    iss.str(line);
                    Component comp;
                    iss >> tmp >> comp.name >> comp.type;

                    if (line.find("+") == string::npos) {
                        getline(def, line);
                        iss.clear();
                        iss.str(line);
                    }
                    iss >> tmp >> comp.status >> tmp >> comp.x1 >> comp.y1>> tmp >> comp.ori;
                    components[comp.name] = move(comp);
                }
            }
        } else if (token == "PINS") {
            iss >> pinCount;
            while (getline(def, line)) {
                if (line.find("END PINS") != string::npos)
                    break;

                if (line.find("-") != string::npos) {
                    Pin pin;
                    iss.clear();
                    iss.str(line);
                    iss >> tmp >> pin.name >> tmp >> tmp >> pin.net;

                    getline(def, line);
                    iss.clear();
                    iss.str(line);
                    iss >> tmp >> tmp >> tmp >> tmp >> pin.x1 >> pin.y1 >> tmp;
                    iss >> tmp >> pin.x2 >> pin.y2;

                    getline(def, line);
                    iss.clear();
                    iss.str(line);
                    iss >> tmp >> tmp >> tmp >> pin.x >> pin.y >> tmp >> pin.ori;

                    pins[pin.name] = move(pin);
                }
            }
        } else if (token == "NETS") {
            iss >> netCount;
            while (getline(def, line)) {
                if (line.find("END NETS") != string::npos)
                    break;

                if (line.find("-") != string::npos) {
                    Net net;
                    iss.clear();
                    iss.str(line);
                    iss >> tmp >> net.name;

                    while (getline(def, line)) {
                        if (line.find(";") != string::npos)
                            break;

                        iss.clear();
                        iss.str(line);
                        while (iss >> tmp) {
                            string tmp1, tmp2;
                            iss >> tmp1 >> tmp2 >> tmp;
                            if (tmp1 == "PIN") {
                                net.pinsX.push_back(&pins[tmp2]);
                            } else {
                                net.compsX.push_back(&components[tmp1]);
                            }
                        }
                    }
                    nets[net.name] = move(net);
                }
            }
        }
    }
    def.close();
    cout << "read DEF" << endl;
}


void transfer_micron() {
    for (auto &it : sites) {
        it.second.sizeX = round(it.second.dsizeX * defMicrons);
        it.second.sizeY = round(it.second.dsizeY * defMicrons);
    }

    for (auto &it : macros) {
        it.second.sizeX = round(it.second.dsizeX * defMicrons);
        it.second.sizeY = round(it.second.dsizeY * defMicrons);
    }
}

// rotate pin and put them to correct position
void transfer_pins() {
    for (auto &pin : pins) {
        if (pin.second.ori == "S") {
            int tx1 = -pin.second.x2;
            int ty1 = -pin.second.y2;
            int tx2 = -pin.second.x1;
            int ty2 = -pin.second.y1;
            pin.second.x1 = tx1;
            pin.second.x2 = tx2;
            pin.second.y1 = ty1;
            pin.second.y2 = ty2;

        } else if (pin.second.ori == "W") {
            int tx1 = -pin.second.y2;
            int ty1 = -pin.second.x1;
            int tx2 = -pin.second.y1;
            int ty2 = pin.second.x2;
            pin.second.x1 = tx1;
            pin.second.x2 = tx2;
            pin.second.y1 = ty1;
            pin.second.y2 = ty2;

        } else if (pin.second.ori == "E") {
            int tx1 = pin.second.y1;
            int ty1 = -pin.second.x2;
            int tx2 = pin.second.y2;
            int ty2 = -pin.second.x1;
            pin.second.x1 = tx1;
            pin.second.x2 = tx2;
            pin.second.y1 = ty1;
            pin.second.y2 = ty2;
        }
        
        pin.second.x1 += pin.second.x;
        pin.second.x2 += pin.second.x;
        pin.second.y1 += pin.second.y;
        pin.second.y2 += pin.second.y;      
    }
}


void set_comps() {
    
    for (auto &comp : components) {
        comp.second.x2 = comp.second.x1 + macros[comp.second.type].sizeX;
        comp.second.y2 = comp.second.y1 + macros[comp.second.type].sizeY;
    }

    for (auto &net : nets) {
        for (Component *comp : net.second.compsX)
            comp->nets.push_back(&nets[net.first]);
        net.second.compsY = net.second.compsX;
        net.second.pinsY = net.second.pinsX;
    }
}

void set_rows() {
    
    for (auto &row : rows) {
        row.second.y2 = row.second.y + sites[row.second.site].sizeY;
        ytoRow[row.second.y] = &rows[row.first];
    }

    for (auto &it : components) {
        Component *c = &it.second;
        if (c->status == "FIXED") {
            for (auto &it : rows) {
                Row *row = &it.second;
                c->sx1 = (c->x1 - row->x) / row->stepX;
                c->sx2 = (c->x2 - row->x) / row->stepX;
                if ((c->x2 - row->x) % row->stepX != 0)
                    c->sx2++;
                if (row->y2 > c->y1 && row->y < c->y2)
                    row->comps.push_back(c);
            }
        } else {
            Row *row = ytoRow[c->y1];
            row->comps.push_back(c);
            c->sx1 = (c->x1 - row->x) / row->stepX;
            c->sx2 = (c->x2 - row->x) / row->stepX;
        }
    }

    for (auto &it : rows) {
        Row *row = &it.second;
        sort(row->comps.begin(), row->comps.end(), [](Component *a, Component *b){
            return a->x1 < b->x1;
        });
    }
}

void write_DEF(string &inputDEF, string &outputDEF) {
    ifstream def(inputDEF);
    ofstream out(outputDEF);
    string line;

    while (getline(def, line)) {
        istringstream iss(line);
        string token, tmp;
        iss >> token;

        if (token == "COMPONENTS") {
            out << "COMPONENTS " << compCount << " ;" << endl;
            while (getline(def, line)) {
                if (line.find("END COMPONENTS") != string::npos) {
                    out << line << endl;
                    break;
                }

                if (line.find("-") != string::npos) {
                    istringstream iss2(line);
                    string dash, name, type;
                    iss2 >> dash >> name >> type;

                    if (components.find(name) != components.end()) {
                        Component &c = components[name];
                        out << "- " << c.name << " " << c.type
                            << " + " << c.status << " ( " << c.x1 << " " << c.y1
                            << " ) " << c.ori << endl << " ;" << endl;
                    } else {
                        out << line << endl;
                    }
                }
            }

        } else {
            out << line << endl;
        }
    }
    def.close();
    out.close();
    cout << "write DEF" << endl;
}
