#ifndef SHADOW_CALCULATOR_H
#define SHADOW_CALCULATOR_H

#include "utils.h"
#include "geometry.h"
#include <stdexcept>
#include <numeric>
#include <eigen/Eigen/Dense>

class Part {
    public:

    std::vector<surface> surfaces;
    int surfaces_number;
    Eigen::Matrix3d& frame;
    Eigen::Vector3d& origin;

    // empty constructor
    Part();
    // default constructor
    Part(std::string path, std::vector<std::vector<double>> prop,
         Eigen::Vector3d origin_ = Eigen::Vector3d(0.0, 0.0, 0.0),
         Eigen::Matrix3d frame_ = Eigen::Matrix3d::Identity());

};

class Assembly {
    public:

    Part bus;
    std::vector<Part> parts;
    int parts_number;

    // default constructor
    template<typename... Parts>
    Assembly(const Part& bus_, const Parts&... additionalParts)
        : bus(bus_), parts{additionalParts...}, parts_number(1 + sizeof...(additionalParts)) {};

    template<typename... R>
    Part combine_parts(const Eigen::Matrix3d& R_bus, const R&... rotation_matrixes) const {

        if (parts.size() != sizeof...(rotation_matrixes)) {
            throw std::runtime_error("Mismatch between Parts and rotation matrices.");
        }

        Part assembly = bus;
        std::vector<Eigen::Matrix3d> rotations = {rotation_matrixes...}; 
        // rotate bus
        for (auto& surf : assembly.surfaces) {
            surf.normal = (R_bus * surf.normal).normalized();
            surf.triangle.a = R_bus * surf.triangle.a;
            surf.triangle.b = R_bus * surf.triangle.b;
            surf.triangle.c = R_bus * surf.triangle.c;
            surf.l = (R_bus * surf.l).normalized();
            surf.m = (R_bus * surf.m).normalized();
        }
        for (int i=0; i<parts_number-1; i++) {
            Part rotated_part = parts[i];
            for (auto& surf : rotated_part.surfaces) {
                surf.normal = (rotations[i] * surf.normal).normalized();
                surf.triangle.a = rotations[i] * surf.triangle.a;
                surf.triangle.b = rotations[i] * surf.triangle.b;
                surf.triangle.c = rotations[i] * surf.triangle.c;
                surf.l = (rotations[i] * surf.l).normalized();
                surf.m = (rotations[i] * surf.m).normalized();
            }
            assembly.surfaces.insert(assembly.surfaces.begin(), rotated_part.surfaces.begin(), rotated_part.surfaces.end());
        }
        return assembly;
    }

};

class Part_v1 {
    public:

    std::vector<surface> surfaces;

    int surfaces_number;
    // default constructor
    Part_v1(const std::string& path);
    // copy rotation constructor
    Part_v1(const Part_v1& original, const Eigen::Matrix3d& R);

    std::vector<surface> _load_surfaces(const std::string& path);
    std::vector<std::vector<double>> _load_properties(const std::string& path);

    void print_info();
    void print_suface_data(int id_number);
   

};

std::vector<double> compute_shadow(const Part_v1& part, const Eigen::Matrix3d& R, const Eigen::Vector3d& v, const int& q_max, bool multi_thread = false);

bool first_discrimination(const surface& panel, const Eigen::Vector3d& v);

std::vector<std::vector<double>> second_discrimination(const surface& panel_shadowed, const surface& panel_shadowing, const Eigen::Vector3d& v);

bool first_test(const surface& panel_shadowed, const surface& panel_shadowing, const std::vector<double>& v);

bool second_test(const surface& panel_shadowed, const surface& panel_shadowing, const std::vector<double>& v);

std::vector<double> pixelation(const std::vector<std::vector<int>>& sigma, const Part_v1& Part_v1, const int& q_max, 
    const std::vector<std::vector<projection>>& projections, const std::vector<int>& indexes);
std::vector<double> rp_coefficients(const Part_v1& Part_v1, const Eigen::Matrix3d& R, const std::vector<double>& fraction, const Eigen::Vector3d& v);

std::vector<double> rp_coefficients_surface(const surface& surf, const double& fraction, const Eigen::Vector3d& v);

#endif