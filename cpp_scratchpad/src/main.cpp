#include <iostream>
#include <vector>
#include "utils.h"
#include "geometry.h"
#include "shadow_calculator.h"

using namespace std;

int main(){
    string data_path = "../data/MRO_highfidelity_";

    Part MRO(data_path);
    MRO.print_info();
    MRO.print_suface_data(300);

    std::vector<double> v({1.0, 1.0, 0.0});

    std::vector<double> f = compute_shadow(MRO, v);

    

    return 0;

}