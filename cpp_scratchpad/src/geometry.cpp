#include <vector>
#include <iostream>
#include <math.h>
#include <numeric>
#include "geometry.h"
#include "utils.h"
#include <eigen/Eigen/Dense>

std::vector<bool> PIP_raycast(const std::vector<point2D>& polygon, const std::vector<point2D>& points) {

    int n = polygon.size();
    int m = points.size();
    std::vector<bool> result;
    result.reserve(m);

    for (int j=0; j<m; j++) {
        int count = 0;
        for (int i = 0; i < n; i++) {
            
            point2D p1 = polygon[i];
            point2D p2 = polygon[(i + 1) % n];
    
            if ((points[j].y > std::min(p1.y, p2.y))
                && (points[j].y <= std::max(p1.y, p2.y))
                && (points[j].x <= std::max(p1.x, p2.x))) {
        
                double xIntersect = (points[j].y - p1.y)
                                        * (p2.x - p1.x)
                                        / (p2.y - p1.y)
                                    + p1.x;
                
                if (p1.x == p2.x || points[j].x <= xIntersect) {
                    count++;
                }
            }
        }
        result[j] = (count % 2 == 1);
    }
    return result;
}

std::vector<bool> PIP_raycast_triangle(const std::vector<point2D>& triangle, const std::vector<point2D>& points) {

    int m = points.size();
    std::vector<bool> result(m, false);
    
    for (int i = 0; i < m; i++) {
        int count = 0;
        const point2D& pt = points[i];
        const point2D& p0 = triangle[0];
        const point2D& p1 = triangle[1];
        const point2D& p2 = triangle[2];
        
        // Edge from p0 to p1
        {
            double ymin = (p0.y < p1.y) ? p0.y : p1.y;
            double ymax = (p0.y < p1.y) ? p1.y : p0.y;
            double xmax = (p0.x > p1.x) ? p0.x : p1.x;
            if ((pt.y > ymin) && (pt.y <= ymax) && (pt.x <= xmax)) {
                double xIntersect = p0.x + (pt.y - p0.y) * (p1.x - p0.x) / (p1.y - p0.y);
                if (p0.x == p1.x || pt.x <= xIntersect)
                    count++;
            }
        }
        
        // Edge from p1 to p2
        {
            double ymin = (p1.y < p2.y) ? p1.y : p2.y;
            double ymax = (p1.y < p2.y) ? p2.y : p1.y;
            double xmax = (p1.x > p2.x) ? p1.x : p2.x;
            if (ymax == ymin && pt.y == ymax) {
                count++;
            }
            if ((pt.y > ymin) && (pt.y <= ymax) && (pt.x <= xmax)) {
                double xIntersect = p1.x + (pt.y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                if (p1.x == p2.x || pt.x <= xIntersect)
                    count++;
            }
        }
        
        // Edge from p2 to p0
        {
            double ymin = (p2.y < p0.y) ? p2.y : p0.y;
            double ymax = (p2.y < p0.y) ? p0.y : p2.y;
            double xmax = (p2.x > p0.x) ? p2.x : p0.x;
            if ((pt.y > ymin) && (pt.y <= ymax) && (pt.x <= xmax)) {
                double xIntersect = p2.x + (pt.y - p2.y) * (p0.x - p2.x) / (p0.y - p2.y);
                if (p2.x == p0.x || pt.x <= xIntersect)
                    count++;
            }
        }
        
        result[i] = (count & 1) != 0;
    }
    
    return result;

}

std::vector<bool> PIP_winding_triangle(const std::vector<point2D>& triangle, const std::vector<point2D>& points) {
    int m = points.size();
    std::vector<bool> result(m, false);
    for (int i = 0; i<m; i++) {
        int winding_number = 0;
        const point2D& pt = points[i];
        const point2D& p0 = triangle[0];
        const point2D& p1 = triangle[1];
        const point2D& p2 = triangle[2];

        // reference
        double ref = sign((p1.x-p0.x)*(p2.y-p0.y) - (p2.x-p0.x)*(p1.y-p0.y));
        // rest
        double dummy1 = (p0.x-pt.x)*(p1.y-pt.y) - (p1.x-pt.x)*(p0.y-pt.y);
        if (dummy1 == 0.0) {
            result[i] = 0;
            continue;
        }
        double sign1 = sign(dummy1);
        if (sign1 == -ref) {
            result[i] = 0;
            continue;
        }
        double dummy2 = (p1.x-pt.x)*(p2.y-pt.y) - (p2.x-pt.x)*(p1.y-pt.y);
        if (dummy2 == 0.0) {
            result[i] = 0;
            continue;
        }
        double sign2 = sign(dummy2);
        if (sign2 == -ref) {
            result[i] = 0;
            continue;
        }
        double dummy3 = (p2.x-pt.x)*(p0.y-pt.y) - (p0.x-pt.x)*(p2.y-pt.y);
        if (dummy3 == 0.0) {
            result[i] = 0;
            continue;
        }
        double sign3 = sign(dummy3);
        if (sign3 == -ref) {
            result[i] = 0;
            continue;
        }
        result[i] = 1;
    }
    return result;
};

std::vector<std::vector<int>> PIP_raycast_triangle_acc(const std::vector<point2D>& triangle,
     const std::vector<double>& x_coord, const std::vector<double>& y_coord) {

        std::vector<std::vector<int>> result(y_coord.size(), std::vector<int>(x_coord.size(), -1));
        // for each line on y coord compute the two intersections
        for (int i = 0; i<y_coord.size(); i++) {
            const double& y = y_coord[i];
            const point2D& p0 = triangle[0];
            const point2D& p1 = triangle[1];
            const point2D& p2 = triangle[2];
            std::vector<double> intersections;
            intersections.reserve(2);
            
            // edge 1
            double ymin_1 = (p0.y < p1.y) ? p0.y : p1.y;
            double ymax_1 = (p0.y < p1.y) ? p1.y : p0.y;
            if ((y > ymin_1) && (y <= ymax_1)) {
                double xIntersect = p0.x + (y - p0.y) * (p1.x - p0.x) / (p1.y - p0.y);
                intersections.push_back(xIntersect);
            }
            // edge 2
            double ymin_2 = (p1.y < p2.y) ? p1.y : p2.y;
            double ymax_2 = (p1.y < p2.y) ? p2.y : p1.y;
            if ((y > ymin_2) && (y <= ymax_2)) {
                double xIntersect = p1.x + (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                intersections.push_back(xIntersect);
            }
            // edge 3
            double ymin_3 = (p2.y < p0.y) ? p2.y : p0.y;
            double ymax_3 = (p2.y < p0.y) ? p0.y : p2.y;
            if ((y > ymin_3) && (y <= ymax_3)) {
                double xIntersect = p2.x + (y - p2.y) * (p0.x - p2.x) / (p0.y - p2.y);
                intersections.push_back(xIntersect);
            }
            double intersection_min = *std::min_element(intersections.begin(), intersections.end());
            double intersection_max = *std::max_element(intersections.begin(), intersections.end());
            for (int j = 0; j<x_coord.size(); j++) {
                if (x_coord[j] >= intersection_min && x_coord[j] <= intersection_max) {
                    result[i][j] = 0;
                }
            }
        }
        return result;
     };

void PIP_raycast_triangle_acc_shadowing(std::vector<std::vector<int>>& W, const std::vector<point2D>& triangle,
        const std::vector<double>& x_coord, const std::vector<double>& y_coord,
        const int& i_min, const int& i_max, const int& j_min, const int& j_max) {
            
            // for each line on y coord compute the two intersections
            for (int i = i_min; i<i_max; i++) {
                const double& y = y_coord[i];
                const point2D& p0 = triangle[0];
                const point2D& p1 = triangle[1];
                const point2D& p2 = triangle[2];
                std::vector<double> intersections;
                
                // edge 1
                double ymin_1 = (p0.y < p1.y) ? p0.y : p1.y;
                double ymax_1 = (p0.y < p1.y) ? p1.y : p0.y;
                if ((y > ymin_1) && (y <= ymax_1)) {
                    double xIntersect = p0.x + (y - p0.y) * (p1.x - p0.x) / (p1.y - p0.y);
                    intersections.push_back(xIntersect);
                }
                // edge 2
                double ymin_2 = (p1.y < p2.y) ? p1.y : p2.y;
                double ymax_2 = (p1.y < p2.y) ? p2.y : p1.y;
                if ((y > ymin_2) && (y <= ymax_2)) {
                    double xIntersect = p1.x + (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
                    intersections.push_back(xIntersect);
                }
                // edge 3
                double ymin_3 = (p2.y < p0.y) ? p2.y : p0.y;
                double ymax_3 = (p2.y < p0.y) ? p0.y : p2.y;
                if ((y > ymin_3) && (y <= ymax_3)) {
                    double xIntersect = p2.x + (y - p2.y) * (p0.x - p2.x) / (p0.y - p2.y);
                    intersections.push_back(xIntersect);
                }
                if (intersections.empty()) {
                    continue;
                }
                double intersection_min = *std::min_element(intersections.begin(), intersections.end());
                double intersection_max = *std::max_element(intersections.begin(), intersections.end());
                
                for (int j =j_min; j<j_max; j++) {
                    if (x_coord[j] >= intersection_min && x_coord[j] <= intersection_max && W[i][j]==0) {
                        W[i][j] = 1;
                    }
                }
            }    
            
    };

surface::surface(triangle3D triangle_, Eigen::Vector3d normal_, double ca_, double cd_, double cs_, std::string material_) : 
    triangle(triangle_), normal(normal_), ca(ca_), cd(cd_), cs(cs_), material(material_),
     area(area_triangle3D(triangle_)) {

        // create the plane versors: l is one of the edges while m is the cross product between n and l
        l = (triangle.b - triangle.a).normalized();
        m = normal.cross(l);
        // parallel projection to find C coordinates in local frame
        Eigen::MatrixXd A(3, 3);
        A(0, 0) =  l(0); A(1, 0) =  l(1); A(2, 0) =  l(2);
        A(0, 1) =  m(0); A(1, 1) =  m(1); A(2, 1) =  m(2);
        A(0, 2) =  normal(0); A(1, 2) =  normal(1); A(2, 2) =  normal(2);

        Eigen::VectorXd b = triangle.c - triangle.a;
        

        Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

        std::vector<double> l_coordinates = {0, (triangle.b - triangle.a).norm(), x(0)};
        std::vector<double> m_coordinates = {0, 0, x(1)};

        point2D A_shadowed(0, 0), B_shadowed((triangle.b - triangle.a).norm(), 0), C_shadowed(x(0), x(1));
        triangle_2d = {A_shadowed, B_shadowed, C_shadowed};

        l_min = *std::min_element(l_coordinates.begin(), l_coordinates.end());
        l_max = *std::max_element(l_coordinates.begin(), l_coordinates.end());
        m_min = *std::min_element(m_coordinates.begin(), m_coordinates.end());
        m_max = *std::max_element(m_coordinates.begin(), m_coordinates.end());
     };


std::vector<double> parallel_projection(const surface& surface, const Eigen::Vector3d& point, const Eigen::Vector3d& v) {
    /**
     * Function used to compute the projection of a point onto a plane, it will yield the planar coordinates 
     * and the lambda coefficient. 
     * It solves the linear system [l m v] [{c1 c2 lambda}]' = P - A,  where A is the vertex A, or M x = b
     */

    // using eigen/dense
    Eigen::Matrix3d A(3, 3);
    A(0, 0) =  surface.l[0]; A(1, 0) =  surface.l[1]; A(2, 0) =  surface.l[2];
    A(0, 1) =  surface.m[0]; A(1, 1) =  surface.m[1]; A(2, 1) =  surface.m[2];
    A(0, 2) =  v[0]; A(1, 2) =  v[1]; A(2, 2) =  v[2];

    Eigen::Vector3d b = point - surface.triangle.a;
    

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    
    return {x(0), x(1), x(2)};
}


double area_triangle3D(const triangle3D& triangle) {
    double AC, BC, AB, s;
    
    AB = (triangle.a- triangle.b).norm();
    BC = (triangle.b- triangle.c).norm();
    AC = (triangle.c- triangle.a).norm();
    
    s = (AB+BC+AC)/2;

    return std::sqrt(s*(s-AB)*(s-BC)*(s-AC));
}

Eigen::Matrix3d R_body_wind(const double& alpha, const double& beta) {

    Eigen::Matrix3d R(3, 3);
    R(0,0) = std::cos(alpha) * std::cos(beta);
    R(0,1) = -std::cos(alpha) * std::sin(beta);
    R(0,2) = -std::sin(alpha);

    R(1,0) = std::sin(beta);
    R(1,1) = std::cos(beta);
    R(1,2) = 0.0;

    R(2,0) =  std::sin(alpha) * std::cos(beta);
    R(2,1) =  -std::sin(alpha) * std::sin(beta);
    R(2,2) =  std::cos(alpha);
    return R;
};



