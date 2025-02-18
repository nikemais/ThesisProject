#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

struct point2D {
    double x, y;

    point2D(double x_, double y_) : x(x_), y(y_) {};
};

struct point3D {
    double x, y, z;

    point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {};
};

struct triangle2D {
    point2D a, b;

    triangle2D(point2D a_, point2D b_) : a(a_), b(b_) {};
};

struct triangle3D {
    point3D a, b, c;

    triangle3D(point3D a_, point3D b_, point3D c_) : a(a_), b(b_), c(c_) {};
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
    std::vector<double> normal;
    // surface properties
    properties surface_properties;
    double area;
    // plane definition
    point3D center;
    std::vector<double> l, m;

    surface(triangle3D triangle_, std::vector<double> normal_, properties prop_);
};


bool PIP_raycast(const std::vector<point2D>& polygon,
    const point2D& point);

// basic geometric functions
double dist(const point3D& p1, const point3D& p2);
double area_triangle3D(const triangle3D& triangle);
point3D center_triangle3D(const triangle3D& triangle);
std::vector<double> cross_product(const std::vector<double>& a, const std::vector<double>& b);
std::vector<std::vector<double>> merge_vectors(const std::vector<std::vector<double>>& vectors);

// projection functions 
std::vector<double> parallel_projection(const surface& surface, const point3D& point, const std::vector<double>& v);


#endif