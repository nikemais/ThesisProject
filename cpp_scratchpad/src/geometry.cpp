#include <vector>
#include <iostream>
#include <math.h>
#include <numeric>
#include "geometry.h"
#include "utils.h"
#include <eigen/Eigen/Dense>

bool PIP_raycast(const std::vector<point2D>& polygon,
                      const point2D& point)
{
    int n = polygon.size();
    int count = 0;

    for (int i = 0; i < n; i++) {
        point2D p1 = polygon[i];
        point2D p2 = polygon[(i + 1) % n];

        if ((point.y > std::min(p1.y, p2.y))
            && (point.y <= std::max(p1.y, p2.y))
            && (point.x <= std::max(p1.x, p2.x))) {
    
            double xIntersect = (point.y - p1.y)
                                    * (p2.x - p1.x)
                                    / (p2.y - p1.y)
                                + p1.x;
            
            if (p1.x == p2.x || point.x <= xIntersect) {
                count++;
            }
        }
    }
    return count % 2 == 1;
}

surface::surface(triangle3D triangle_, std::vector<double> normal_, properties prop_) : 
    triangle(triangle_), normal(normal_), surface_properties(prop_),
     area(area_triangle3D(triangle_)), center(center_triangle3D(triangle_)) {
        // create the plane versors: l is one of the edges while m is the cross product between n and l
        double l_norm = dist(triangle.b, triangle.a);
        l = {(triangle.b.x - triangle.a.x)/l_norm,
             (triangle.b.y - triangle.a.y)/l_norm,
             (triangle.b.z - triangle.a.z)/l_norm};
        m = cross_product(normal, l);
        
     };


std::vector<double> parallel_projection(const surface& surface, const point3D& point, const std::vector<double>& v) {
    /**
     * Function used to compute the projection of a point onto a plane, it will yield the planar coordinates 
     * and the lambda coefficient. 
     * It solves the linear system [l m v] [{c1 c2 lambda}]' = P - C or A x = b
     */
    // second approach which does not use eigen/dense
    // std::vector<std::vector<double>> A = merge_vectors({surface.l, surface.m, v});
    // std::vector<double> x, b;
    // b = {point.x - surface.center.x, point.y - surface.center.y, point.z - surface.center.z};

    // using eigen/dense
    Eigen::MatrixXd A(3, 3);
    A << surface.l[0], surface.l[1], surface.l[2], surface.m[0], surface.m[1],
         surface.m[2], surface.normal[0], surface.normal[1], surface.normal[2];

    Eigen::VectorXd b(3);
    b << point.x - surface.center.x, point.y - surface.center.y, point.z - surface.center.z;

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    
    return {x(0), x(1), x(2)};
}

double dist(const point3D& p1, const point3D& p2) {
    return std::sqrt(
        (p1.x - p2.x) * (p1.x - p2.x) +
        (p1.y - p2.y) * (p1.y - p2.y) +
        (p1.z - p2.z) * (p1.z - p2.z)
    );
}

double area_triangle3D(const triangle3D& triangle) {
    double AC, BC, AB, s;
    
    AB = dist(triangle.a, triangle.b);
    BC = dist(triangle.b, triangle.c);
    AC = dist(triangle.a, triangle.c);
    
    s = (AB+BC+AC)/2;

    return std::sqrt(s*(s-AB)*(s-BC)*(s-AC));
}

point3D center_triangle3D(const triangle3D& triangle) {
    return point3D((triangle.a.x + triangle.b.x + triangle.c.x)/3,
                   (triangle.a.y + triangle.b.y + triangle.c.y)/3,
                   (triangle.a.z + triangle.b.z + triangle.c.z)/3);
};

std::vector<double> cross_product(const std::vector<double>& a, const std::vector<double>& b) {
    return {
        a[1] * b[2] - a[2] * b[1],  
        a[2] * b[0] - a[0] * b[2],  
        a[0] * b[1] - a[1] * b[0]   
    };
}

std::vector<std::vector<double>> merge_vectors(const std::vector<std::vector<double>>& vectors) {
    if (vectors.empty()) return {};

    size_t numRows = vectors[0].size();
    size_t numCols = vectors.size();

    for (const auto& vec : vectors) {
        if (vec.size() != numRows) {
            std::cerr << "Error: All vectors must have the same length!" << std::endl;
            return {};
        }
    }

    std::vector<std::vector<double>> matrix(numRows, std::vector<double>(numCols));

    for (size_t col = 0; col < numCols; ++col) {
        for (size_t row = 0; row < numRows; ++row) {
            matrix[row][col] = vectors[col][row];
        }
    }

    return matrix;
}