#ifndef SHADOW_CALCULATOR_H
#define SHADOW_CALCULATOR_H

#include "utils.h"
#include "geometry.h"
#include <stdexcept>
#include <numeric>

class Part {
    public:

    std::vector<surface> surfaces;

    int surfaces_number;

    Part(const std::string& path);

    std::vector<surface> _load_surfaces(const std::string& path);
    std::vector<properties> _load_properties(const std::string& path);

    void print_info();
    void print_suface_data(int id_number);
   

};

std::vector<double> compute_shadow(const Part& part, const std::vector<double>& v);

bool first_discrimination(const surface& panel, const std::vector<double>& v);

bool second_discrimination(const surface& panel_shadowed, const surface& panel_shadowing, const std::vector<double>& v);


#endif