#include <iostream>

#include "deal_datas.h"
#include "utils.h"
#include "algos.h"

int main(int argc, char* argv[]) {

    string inputLEF = argv[1];
    string inputDEF = argv[2];
    string outputDEF = argv[3];


    read_LEF(inputLEF);
    read_DEF(inputDEF);

    transfer_micron();
    transfer_pins();
    set_comps();
    set_rows();
    cout << "setting end" << endl;

    long long WL = calculate_WL();
    cout << "origin WL:" << WL << endl;

    if (components.size() < 10000) {
        FastDP(8);
        reordering(20);
    } else if (components.size() < 100000) {
        FastDP(10);
        reordering(10);
    } else if (components.size() < 600000) {
        FastDP(3);
        reordering(2);
    } else {
        FastDP(2);
        reordering(1);
    }
    
    write_DEF(inputDEF, outputDEF);

    return 0;
}