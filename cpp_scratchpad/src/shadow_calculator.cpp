#include "shadow_calculator.h"
#include "utils.h"
#include <algorithm>

Part::Part(const std::string& path) {
    surfaces = _load_surfaces(path);
    surfaces_number = surfaces.size();
}

std::vector<properties> Part::_load_properties(const std::string& path) {
    std::vector<properties> prop;

    // Load data
    std::vector<std::vector<double>> id = load_data(path + "id.txt");
    std::vector<std::vector<double>> pa_infrared = load_data(path + "pa_infrared.txt");
    std::vector<std::vector<double>> pd_infrared = load_data(path + "pd_infrared.txt");
    std::vector<std::vector<double>> ps_infrared = load_data(path + "ps_infrared.txt");
    std::vector<std::vector<double>> pa_optical = load_data(path + "pa_optical.txt");
    std::vector<std::vector<double>> pd_optical = load_data(path + "pd_optical.txt");
    std::vector<std::vector<double>> ps_optical = load_data(path + "ps_optical.txt");

    // Check if data was loaded correctly
    if (id.empty() || pa_infrared.empty() || pd_infrared.empty() || ps_infrared.empty() ||
        pa_optical.empty() || pd_optical.empty() || ps_optical.empty()) {
        std::cerr << "Error: One or more files could not be loaded!" << std::endl;
        return prop;  
    }

    // Ensure all matrices have the correct shape (1 row)
    if (id.size() != 1 || pa_infrared.size() != 1 || pd_infrared.size() != 1 || ps_infrared.size() != 1 ||
        pa_optical.size() != 1 || pd_optical.size() != 1 || ps_optical.size() != 1) {
        std::cerr << "Error: Matrices must have 1 row!" << std::endl;
        return prop;
    }
    size_t numSurfaces = id[0].size();  

    for (size_t i = 0; i < numSurfaces; i++) {

        prop.emplace_back(id[0][i], pa_infrared[0][i], pd_infrared[0][i], ps_infrared[0][i],
                          pa_optical[0][i], pd_optical[0][i], ps_optical[0][i]);  
    }

    return prop;  

}


std::vector<surface> Part::_load_surfaces(const std::string& path) {

    std::vector<properties> surface_properties = _load_properties(path);

    std::vector<surface> surfaces;

    // Load data
    std::vector<std::vector<double>> v_0 = load_data(path + "v0.txt");
    std::vector<std::vector<double>> v_1 = load_data(path + "v1.txt");
    std::vector<std::vector<double>> v_2 = load_data(path + "v2.txt");
    std::vector<std::vector<double>> nrl = load_data(path + "nrl.txt");

    // Check if data was loaded correctly
    if (v_0.empty() || v_1.empty() || v_2.empty() || nrl.empty()) {
        std::cerr << "Error: One or more files could not be loaded!" << std::endl;
        return surfaces;  // Return empty surfaces vector
    }

    // Ensure all matrices have the correct shape (3 rows)
    if (v_0.size() != 3 || v_1.size() != 3 || v_2.size() != 3 || nrl.size() != 3) {
        std::cerr << "Error: Matrices must have 3 rows (x, y, z coordinates)!" << std::endl;
        return surfaces;
    }

    size_t numTriangles = v_0[0].size();  // Number of triangles

    // Ensure all matrices have the same number of columns
    if (v_1[0].size() != numTriangles || v_2[0].size() != numTriangles || nrl[0].size() != numTriangles) {
        std::cerr << "Error: All matrices must have the same number of columns!" << std::endl;
        return surfaces;
    }

    // Loop through columns to construct surfaces
    for (size_t i = 0; i < numTriangles; i++) {
        // Extract triangle vertices
        point3D p0(v_0[0][i], v_0[1][i], v_0[2][i]);
        point3D p1(v_1[0][i], v_1[1][i], v_1[2][i]);
        point3D p2(v_2[0][i], v_2[1][i], v_2[2][i]);

        // Extract normal vector
        std::vector<double> normal = {nrl[0][i], nrl[1][i], nrl[2][i]};

        // Create triangle and surface
        triangle3D triangle(p0, p1, p2);
        surfaces.emplace_back(triangle, normal, surface_properties[i]);  // Add surface to vector
    }
    
    return surfaces;  // Return the vector of surfaces
}

void Part::print_info() {
    std::cout<<"-- PART INFO -----------------------------"<<std::endl;
    std::cout<<"Number of panels: "<<surfaces_number<<std::endl;
    std::cout<<"------------------------------------------"<<std::endl;
}

void Part::print_suface_data(int id_number) {
    std::cout<<"-- SURFACE INFO --------------------------"<<std::endl;
    std::cout<<"ID: "<<id_number<<std::endl;
    std::cout<<"Area: "<<surfaces[id_number-1].area<<" m^2"<<std::endl;
    std::cout<<"V1 = ["<<surfaces[id_number-1].triangle.a.x<<","<<surfaces[id_number-1].triangle.a.y
             <<","<<surfaces[id_number-1].triangle.a.z<<"] m"<<std::endl;
    std::cout<<"V2 = ["<<surfaces[id_number-1].triangle.b.x<<","<<surfaces[id_number-1].triangle.b.y
             <<","<<surfaces[id_number-1].triangle.b.z<<"] m"<<std::endl;
    std::cout<<"V3 = ["<<surfaces[id_number-1].triangle.c.x<<","<<surfaces[id_number-1].triangle.c.y
             <<","<<surfaces[id_number-1].triangle.c.z<<"] m"<<std::endl;
    std::cout<<"n = ["<<surfaces[id_number-1].normal[0]<<","<<surfaces[id_number-1].normal[1]
             <<","<<surfaces[id_number-1].normal[2]<<"] "<<std::endl;
    std::cout<<"Optical absorptivity: "<<surfaces[id_number-1].surface_properties.optical[0]<<std::endl;
    std::cout<<"Optical diffuse reflectivity: "<<surfaces[id_number-1].surface_properties.optical[1]<<std::endl;
    std::cout<<"Optical specular reflectivity: "<<surfaces[id_number-1].surface_properties.optical[2]<<std::endl;
    std::cout<<"Infrared absorptivity: "<<surfaces[id_number-1].surface_properties.infrared[0]<<std::endl;
    std::cout<<"Infrared diffuse reflectivity: "<<surfaces[id_number-1].surface_properties.infrared[1]<<std::endl;
    std::cout<<"Infrared specular reflectivity: "<<surfaces[id_number-1].surface_properties.infrared[2]<<std::endl;
    std::cout<<"------------------------------------------"<<std::endl;
    
}

bool first_discrimination(const surface& panel, const std::vector<double>& v) {
    /**
     * Function that checks if the panel is incident or not wrt a direction by
     * computing the scalar product between the vector and the norm to the panel
     */


    if (panel.normal.size() != v.size()) {
        throw std::invalid_argument("Incorrect definition of normal or direction vector!");
    }
    if (std::inner_product(panel.normal.begin(), panel.normal.end(), v.begin(), 0.0) < 0) {
        return 1;
    }
    else {
        return 0;
    }
}

bool second_discrimination(const surface& panel_shadowed, const surface& panel_shadowing, const std::vector<double>& v) {
    /**
     * Function that applies the second discrimination
     */
    std::vector<double> lambdas;

    std::vector<double> coefficients_a = parallel_projection(panel_shadowed, panel_shadowing.triangle.a, v);
    std::vector<double> coefficients_b = parallel_projection(panel_shadowed, panel_shadowing.triangle.b, v);
    std::vector<double> coefficients_c = parallel_projection(panel_shadowed, panel_shadowing.triangle.c, v);

    lambdas = {-coefficients_a[2], -coefficients_b[2], -coefficients_c[2]};

    if (*std::min_element(lambdas.begin(), lambdas.end()) >= 0.0) {
        return 0;
    }
    else {
        return 1;
    }



}

std::vector<double> compute_shadow(const Part& part, const std::vector<double>& v) {

    int n = part.surfaces_number;
    std::vector<std::vector<int>> sigma(n, std::vector<int>(n, 0));

    std::vector<double> fraction;

    // applying first discrimination
    std::vector<int> first_discrimination_indexes;
    for (int i = 0; i < n; i++) {
        if (!first_discrimination(part.surfaces[i], v)){
            first_discrimination_indexes.push_back(i);
            for (int j = 0; j < n; j++) {
                sigma[j][i] = -1;
            }
        }
    }
    for (int i = 0; i < first_discrimination_indexes.size(); i++) {
        for (int j = 0; j < n; j++) {
            sigma[first_discrimination_indexes[i]][j] = 1;
        }
    }
    std::cout<<"------------------------------------------"<<std::endl;
    std::cout<<"Plates excluded: "<<first_discrimination_indexes.size()<<std::endl;
    std::cout<<"------------------------------------------"<<std::endl;

    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++){
            // second discrimination
            if (sigma[i][j] == 0 && i != j){
                if (second_discrimination(part.surfaces[i], part.surfaces[j], v)) {
                    // apply testing
                    sigma[i][j] = 1;
                }
                else {
                    sigma[i][j] = -1;
                }
            }
        }
    }

    return fraction;

};