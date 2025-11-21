#ifndef READ_FILES_H
#define READ_FILES_H

#include <string>
#include <vector>
#include <unordered_map>
#include <map>

using namespace std;

struct Site {
    string name;
    string sClass;
    double dsizeX, dsizeY;      // size in micron
    int sizeX, sizeY;           // transfer micron to unit
};

struct Macro {
    string name;
    string mClass;
    double dsizeX, dsizeY;                      // size in micron
    pair<bool, bool> syms = {false, false};     // record "X" & "Y" exist or not
    string siteName; 
    int sizeX, sizeY;                           // transfer micron to unit
};

struct Area {
    int x1, y1, x2, y2;
};

struct Net;

struct Component {
    string name;
    string type;
    string status;
    int x1, y1, x2, y2;     // coordinates of four corners
    string ori;
    vector<Net*> nets;      // nets containing this component
    int sx1, sx2;           // transfer coordinate to align site
    Area opt;               // optimal region
};

struct Row {
    string name;
    string site;
    int x, y, y2;
    string ori;
    int numX, numY;
    int stepX, stepY;
    vector<Component*> comps;   // components on this row
};

struct Pin {
    string name;
    string net;
    int x1, y1, x2, y2, x, y;
    string ori;
};

struct Net {
    string name;
    vector<Component*> compsX;  // components on this net, sorted in coordinate x1
    vector<Component*> compsY;  // components on this net, sorted in coordinate y1
    vector<Pin*> pinsX;
    vector<Pin*> pinsY;
};

extern int lefMicrons;
extern unordered_map<string, Site> sites;
extern unordered_map<string, Macro> macros;

extern int defMicrons, compCount, pinCount, netCount;
extern Area dieArea;
extern unordered_map<string, Row> rows;
extern map<int, Row*> ytoRow;
extern unordered_map<string, Component> components;
extern unordered_map<string, Pin> pins;
extern unordered_map<string, Net> nets;

void read_LEF(string &inputLEF);

void read_DEF(string &inputDEF);

void transfer_micron();

void transfer_pins();

void set_rows();

void set_comps();

void write_DEF(string &inputDEF, string &outputDEF);


#endif
