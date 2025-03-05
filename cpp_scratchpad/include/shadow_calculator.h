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
    // default constructor
    Part(const std::string& path);
    // copy rotation constructor
    Part(const Part& original, const Eigen::Matrix3d& R);

    std::vector<surface> _load_surfaces(const std::string& path);
    std::vector<properties> _load_properties(const std::string& path);

    void print_info();
    void print_suface_data(int id_number);
   

};

std::vector<double> compute_shadow(const Part& part, const Eigen::Matrix3d& R, const Eigen::Vector3d& v, const int& q_max, bool multi_thread = false);

bool first_discrimination(const surface& panel, const Eigen::Vector3d& v);

std::vector<std::vector<double>> second_discrimination(const surface& panel_shadowed, const surface& panel_shadowing, const Eigen::Vector3d& v);

bool first_test(const surface& panel_shadowed, const surface& panel_shadowing, const std::vector<double>& v);

bool second_test(const surface& panel_shadowed, const surface& panel_shadowing, const std::vector<double>& v);

std::vector<double> pixelation(const std::vector<std::vector<int>>& sigma, const Part& part, const int& q_max, 
    const std::vector<std::vector<projection>>& projections, const std::vector<int>& indexes);
std::vector<double> rp_coefficients(const Part& part, const Eigen::Matrix3d& R, const std::vector<double>& fraction, const Eigen::Vector3d& v);

#endif