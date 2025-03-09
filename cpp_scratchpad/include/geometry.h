#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <algorithm>
#include <eigen/Eigen/Dense>

struct point2D {
    double x, y;

    point2D(double x_, double y_) : x(x_), y(y_) {};
};

struct point3D {
    double x, y, z;

    point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {};
};

struct triangle2D {
    point2D a, b, c;

    triangle2D(point2D a_, point2D b_, point2D c_) : a(a_), b(b_), c(c_) {};
};

struct triangle3D {
    Eigen::Vector3d a, b, c;

    triangle3D(Eigen::Vector3d a_, Eigen::Vector3d b_, Eigen::Vector3d c_) : a(a_), b(b_), c(c_) {};
};

struct properties {
    int id;
    std::vector<double> optical, infrared;

    properties(int id_,
        double optical_absorption_,
        double optical_diffuse_reflection_,
        double optical_specular_reflection_,
        double infrared_absorption_, 
        double infrared_diffuse_reflection_, 
        double infrared_specular_reflection_) : 

        id(id_),
        optical( {optical_absorption_, optical_diffuse_reflection_, optical_specular_reflection_} ),
        infrared( {infrared_absorption_, infrared_diffuse_reflection_, infrared_specular_reflection_} ) {};

};

struct surface {
    // surface geometric definition
    triangle3D triangle;
    Eigen::Vector3d normal;
    // surface properties
    double ca, cd, cs;
    double area;
    std::string material;
    // plane definition
    Eigen::Vector3d l, m;
    std::vector<point2D> triangle_2d;
    double l_min, l_max, m_min, m_max;

    surface(triangle3D triangle_, Eigen::Vector3d normal_, double ca_, double cd_, double cs_, std::string material_);

    // surface(const surface& s_): triangle(s_.triangle), normal(s_.normal), surface_properties(s_.surface_properties),
    // area(s_.area), l(s_.l), m(s_.m), triangle_2d(s_.triangle_2d), l_min(s_.l_min), l_max(s_.l_max), m_min(s_.m_min), m_max(s_.m_max) {};
};

struct projection {
    std::vector<point2D> triangle;
    double l_min, l_max, m_min, m_max;
    // default contrusctor
    projection() : l_min(0), l_max(0), m_min(0), m_max(0) {}
    // constructor
    projection(std::vector<double> l_coords, std::vector<double> m_coords) {

        point2D A(l_coords[0], m_coords[0]), B(l_coords[1], m_coords[1]), C(l_coords[2], m_coords[2]);

        triangle = {A, B, C};

        l_min = *std::min_element(l_coords.begin(), l_coords.end());
        l_max = *std::max_element(l_coords.begin(), l_coords.end());
        m_min = *std::min_element(m_coords.begin(), m_coords.end());
        m_max = *std::max_element(m_coords.begin(), m_coords.end());
    }
};

// PIP algorithms
std::vector<bool> PIP_raycast(const std::vector<point2D>& polygon, const std::vector<point2D>& points);
std::vector<bool> PIP_raycast_triangle(const std::vector<point2D>& triangle, const std::vector<point2D>& points);
std::vector<std::vector<int>> PIP_raycast_triangle_acc(const std::vector<point2D>& triangle,
    const std::vector<double>& x_coord, const std::vector<double>& y_coord);
void PIP_raycast_triangle_acc_shadowing(std::vector<std::vector<int>>& W, const std::vector<point2D>& triangle,
        const std::vector<double>& x_coord, const std::vector<double>& y_coord,
        const int& i_min, const int& i_max, const int& j_min, const int& j_max);
std::vector<bool> PIP_winding_triangle(const std::vector<point2D> triangle, const std::vector<point2D> points);

// basic geometric functions
double area_triangle3D(const triangle3D& triangle);

// projection functions 
std::vector<double> parallel_projection(const surface& surface, const Eigen::Vector3d& point, const Eigen::Vector3d& v);

Eigen::Matrix3d R_body_wind(const double& alpha, const double& beta);

#endif