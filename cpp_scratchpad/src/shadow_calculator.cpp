#include "shadow_calculator.h"
#include "utils.h"
#include <algorithm>
#include <math.h> 
#include <thread>
#include <eigen/Eigen/Dense>

// default constructor
Part::Part(const std::string& path) {
    surfaces = _load_surfaces(path);
    surfaces_number = surfaces.size();
}
// copy rotation constructor
Part::Part(const Part& original, const Eigen::Matrix3d& R) : surfaces(original.surfaces),
    surfaces_number(original.surfaces_number) {

    for (auto& surf : surfaces) {
        surf.normal = (R * surf.normal).normalized();
        surf.triangle.a = R * surf.triangle.a;
        surf.triangle.b = R * surf.triangle.b;
        surf.triangle.c = R * surf.triangle.c;
        surf.l = (R * surf.l).normalized();
        surf.m = (R * surf.m).normalized();
    }
};


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
        Eigen::Vector3d p0(v_0[0][i], v_0[1][i], v_0[2][i]);
        Eigen::Vector3d p1(v_1[0][i], v_1[1][i], v_1[2][i]);
        Eigen::Vector3d p2(v_2[0][i], v_2[1][i], v_2[2][i]);

        // Extract normal vector
        Eigen::Vector3d normal(nrl[0][i], nrl[1][i], nrl[2][i]);

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
    std::cout<<"V1 = ["<<surfaces[id_number-1].triangle.a(0)<<","<<surfaces[id_number-1].triangle.a(1)
             <<","<<surfaces[id_number-1].triangle.a(2)<<"] m"<<std::endl;
    std::cout<<"V2 = ["<<surfaces[id_number-1].triangle.b(0)<<","<<surfaces[id_number-1].triangle.b(1)
             <<","<<surfaces[id_number-1].triangle.b(2)<<"] m"<<std::endl;
    std::cout<<"V3 = ["<<surfaces[id_number-1].triangle.c(0)<<","<<surfaces[id_number-1].triangle.c(1)
             <<","<<surfaces[id_number-1].triangle.c(2)<<"] m"<<std::endl;
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

bool first_discrimination(const surface& panel, const Eigen::Vector3d& v) {
    /**
     * Function that checks if the panel is incident or not wrt a direction by
     * computing the scalar product between the vector and the norm to the panel
     */
    if (panel.normal.dot(v) < 0) {
        return 1;
    }
    else {
        return 0;
    }
}

std::vector<std::vector<double>> second_discrimination(const surface& panel_shadowed, const surface& panel_shadowing, const Eigen::Vector3d& v) {
    /**
     * Function that applies the second discrimination
     */
    std::vector<double> lambdas;
    lambdas.reserve(3);

    std::vector<double> coefficients_a = parallel_projection(panel_shadowed, panel_shadowing.triangle.a, v);
    std::vector<double> coefficients_b = parallel_projection(panel_shadowed, panel_shadowing.triangle.b, v);
    std::vector<double> coefficients_c = parallel_projection(panel_shadowed, panel_shadowing.triangle.c, v);

    lambdas = {-coefficients_a[2], -coefficients_b[2], -coefficients_c[2]};

    if (*std::max_element(lambdas.begin(), lambdas.end()) <= 0.0) {
        return {{0}};
    }
    else {
        return {{1}, coefficients_a, coefficients_b, coefficients_c};
    }
}

std::vector<double> compute_shadow(const Part& part, const Eigen::Matrix3d& R, const Eigen::Vector3d& v, const int& q_max, bool multi_thread) {
    // copy part and rotate all surfaces
    Part part_rotated(part, R);

    int n = part.surfaces_number;
    std::vector<std::vector<int>> sigma(n, std::vector<int>(n, -2)); // -2 dummy value, -1 no-sh, 0 partial sh, 1 full sh

    // DISCRIMINATION
    // first discrimination
    std::vector<int> first_discrimination_indexes;
    for (int i = 0; i < n; i++) {
        if (!first_discrimination(part_rotated.surfaces[i], v)){
            first_discrimination_indexes.push_back(i);
            for (int j = 0; j < n; j++) {
                sigma[j][i] = -1; // cannot shadow other plates
            }
        }
    }
    for (int i = 0; i < first_discrimination_indexes.size(); i++) {
        for (int j = 0; j < n; j++) {
            sigma[first_discrimination_indexes[i]][j] = 1; // is completely shadowed
        }
    }
    // std::cout<<"------------------------------------------"<<std::endl;
    // std::cout<<"Plates excluded: "<<first_discrimination_indexes.size()<<std::endl;
    // std::cout<<"------------------------------------------"<<std::endl;
    
    //substitute the values on the diagonal with -1 (a plate cannot shadow itself), not if it already has a value of 1 (full shadow)
    for (int i = 0; i < n; i++) {
        if (sigma[i][i] == -2) {
            sigma[i][i] = -1;
        }
    }
    // initialize array with projections info
    std::vector<std::vector<projection>> projections(n, std::vector<projection>(n));
    for (int i = 0; i < n; i++) {
        if (sigma[i][0] == 1) {
            continue; //full shadow detected
        }
        for (int j = 0; j < n; j++){
            // second discrimination
            if (sigma[i][j] != -2) {
                continue;
            }
            std::vector<std::vector<double>> output = second_discrimination(part_rotated.surfaces[i], part_rotated.surfaces[j], v);
            if (output[0][0]) {
                // apply testing if second discrimination yiels p-sh

                // first test (max/min coordinates)
                // l,m coordinates for shadowing triangle
                std::vector<double> l_coordinates = {output[1][0], output[2][0], output[3][0]}; 
                std::vector<double> m_coordinates = {output[1][1], output[2][1], output[3][1]}; 

                projection shadowing(l_coordinates, m_coordinates);

                if ((shadowing.l_min>=part_rotated.surfaces[i].l_max)||
                    (shadowing.l_max<=part_rotated.surfaces[i].l_min)||
                    (shadowing.m_min>=part_rotated.surfaces[i].m_max)||
                    (shadowing.m_max<=part_rotated.surfaces[i].m_min)) {
                    sigma[i][j] = -1; // no-sh
                }
                else{
                    // second test (PIP testing)
                    std::vector<bool> shadowing_in_shadowed = PIP_raycast_triangle(part_rotated.surfaces[i].triangle_2d, shadowing.triangle);
                    std::vector<bool> shadowed_in_shadowing = PIP_raycast_triangle(shadowing.triangle, part_rotated.surfaces[i].triangle_2d);
                    if (std::any_of(shadowing_in_shadowed.begin(), shadowing_in_shadowed.end(), [](bool v) { return v; })){
                        sigma[i][j] = 0; // all points of the shadowing plate fall inside the shadowed one, partial sh
                        sigma[j][i] = -1; // the opposite is clearly a case of no-sh, as a shadowed plate cannot cast a shadow on the plate casting on it
                        
                    }
                    if (std::all_of(shadowed_in_shadowing.begin(), shadowed_in_shadowing.end(), [](bool v) { return v; })){
                        for (int u=0; u<n; u++){
                            sigma[i][u] = 1; // all shadowed points are inside the shadowing, full sh
                        }
                        for (int u=0; u<n; u++){
                            if (sigma[u][i] == -2) {
                                sigma[u][i] = -1; //the plate cannot cast shadow on the rest of the plates
                            }
                        }
                        
                    }
                    else {
                        sigma[i][j] = 0; // otherwise cannot decide, needs to be computed using pixelation
                        projections[i][j] = shadowing;
                    }
                }
            }
            else {
                sigma[i][j] = -1;
            }    
            
        }
    }
    // PIXELATION
    std::vector<int> to_be_pixelated;
    std::vector<double> fraction(n, 0.0);
    for (int i = 0; i<n; i++) {
        if (sigma[i][0] == 1) {
            continue;
        }
        to_be_pixelated.push_back(i);    
    }
    if (multi_thread) {
        // multithreading
        
        int n_to_be_pixeled = to_be_pixelated.size();
        int n_threads = std::thread::hardware_concurrency()-1;
        std::vector<std::thread> threads;
        std::vector<std::vector<double>> results(n_threads);
        int chunk = static_cast<int>(std::round(n_to_be_pixeled/n_threads));

        for (int i = 0; i < n_threads-1; i++) {
            // Launch thread
            std::vector<int> indexes(to_be_pixelated.begin() + i*chunk, to_be_pixelated.begin() + (i+1)*chunk);
            threads.emplace_back([&, i, indexes]() {
                results[i] = pixelation(sigma, part_rotated, q_max, projections, indexes);
            });
        }
        int i = n_threads-1;
        std::vector<int> indexes(to_be_pixelated.begin() + i*chunk, to_be_pixelated.end());
        threads.emplace_back([&, i, indexes]() {
            results[n_threads-1] = pixelation(sigma, part_rotated, q_max, projections, indexes);
        });
        
        // Join all threads
        for (auto& th : threads) {
            th.join();
        }

        // Merge results
        std::vector<double> final_result;
        for (const auto& res : results) {
            final_result.insert(final_result.end(), res.begin(), res.end());
        }
        for (int i = 0; i<n_to_be_pixeled; i++) {
            fraction[to_be_pixelated[i]] = final_result[i];
        }
    }
    else {
        // single thread
        std::vector<double> final_result = pixelation(sigma, part_rotated, q_max, projections, to_be_pixelated);
        
        int c = 0;
        for (int i = 0; i<n; i++) {
            if (c<to_be_pixelated.size()) {
                if (to_be_pixelated[c] != i) {
                    fraction[i] = 0.0;
                }
                else {
                    fraction[i] = final_result[c];
                    c++;
                } 
            }
            else {
                fraction[i] = 0.0;
            }
            
        }
        
    }
    return fraction;
    
};

std::vector<double> pixelation(const std::vector<std::vector<int>>& sigma, const Part& part, const int& q_max, 
    const std::vector<std::vector<projection>>& projections, const std::vector<int>& indexes) {
        std::vector<double> f(indexes.size(), 1.0);
        for (int ii=0; ii<indexes.size(); ii++) { 
            // cycle through all the rows in sigma
            int i = indexes[ii];
            double l_min = part.surfaces[i].l_min;
            double l_max = part.surfaces[i].l_max;
            double m_min = part.surfaces[i].m_min;
            double m_max = part.surfaces[i].m_max;
    
            double delta_l = l_max - l_min;
            double delta_m = m_max - m_min;
            // choosing number of pixels
            double n_m, n_l;
            if (delta_m>=delta_l) {
                n_m = static_cast<double>(q_max);
                n_l = static_cast<double>(std::max(q_max/10, static_cast<int>(std::round(n_m * delta_l/delta_m))));
            }
            else {
                n_l = static_cast<double>(q_max);
                n_m = static_cast<double>(std::max(q_max/10, static_cast<int>(std::round(n_l * delta_m/delta_l))));
            }
            double size_l = delta_l/n_l;
            double size_m = delta_m/n_m;
            // creating the grid
            std::vector<double> l_coord = linspace(part.surfaces[i].l_min + size_l/2, part.surfaces[i].l_max - size_l/2, n_l);
            std::vector<double> m_coord = linspace(part.surfaces[i].m_min + size_m/2, part.surfaces[i].m_max - size_m/2, n_m);
    
            std::vector<std::vector<int>> W = PIP_raycast_triangle_acc(part.surfaces[i].triangle_2d, l_coord, m_coord);
            // std::cout<<"Pixelating surface: "<<i<<"  nl: "<<n_l<<"  nm: "<<n_m<<std::endl;
            for (int j = 0; j<sigma.size(); j++) {
                if (sigma[i][j] != 0){
                    continue;
                }
                // resize the pixelation grid by computing index ranges
                int j_min = std::clamp(static_cast<int>(std::floor((projections[i][j].l_min-l_min)/size_l))-1, 0, static_cast<int>(n_l)-1);
                int j_max = std::clamp(static_cast<int>(std::ceil((projections[i][j].l_max-l_min)/size_l))+1, j_min+1, static_cast<int>(n_l));
                
                int i_min = std::clamp(static_cast<int>(std::floor((projections[i][j].m_min-m_min)/size_m))-1, 0, static_cast<int>(n_m)-1);
                int i_max = std::clamp(static_cast<int>(std::ceil((projections[i][j].m_max-m_min)/size_m))+1, i_min+1, static_cast<int>(n_m));
                
                PIP_raycast_triangle_acc_shadowing(W, projections[i][j].triangle, l_coord, m_coord, i_min, i_max, j_min, j_max);
            }
            double f_1 = 0;
            for (int p = 0; p<m_coord.size(); p++) {
                for (int q = 0; q<l_coord.size(); q++) {
                    if (W[p][q] == 1) {
                        f_1++;
                    }
                }
            }
            double f_2 = 0;
            for (int p = 0; p<m_coord.size(); p++) {
                for (int q = 0; q<l_coord.size(); q++) {
                    if (W[p][q] == 0) {
                        f_2++;
                    }
                }
            }
    
            // f[i] = (part.surfaces[i].area - f * size_l * size_m)/ part.surfaces[i].area;
            f[ii] = std::clamp(1 - f_1 / (f_1 + f_2), 1e-3, 1.0);
        }
        return f;
        
    };

std::vector<double> rp_coefficients(const Part& part, const Eigen::Matrix3d& R,
     const std::vector<double>& fraction, const Eigen::Vector3d& v) {

    Eigen::Vector3d C_op(0.0,0.0,0.0);
    Eigen::Vector3d C_ir(0.0,0.0,0.0);
    Part part_rotated(part, R);
    Eigen::Vector3d r = v;
    
    for (int i = 0; i<part_rotated.surfaces_number; i++) {
            // vectorial components
            double ca_op = part_rotated.surfaces[i].surface_properties.optical[0];
            double cd_op = part_rotated.surfaces[i].surface_properties.optical[1];
            double cs_op = part_rotated.surfaces[i].surface_properties.optical[2];
            double ca_ir = part_rotated.surfaces[i].surface_properties.infrared[0];
            double cd_ir = part_rotated.surfaces[i].surface_properties.infrared[1];
            double cs_ir = part_rotated.surfaces[i].surface_properties.infrared[2];
            Eigen::Vector3d n = part_rotated.surfaces[i].normal;
            
            double cos_theta_i = r.dot(n);
            Eigen::Vector3d dummy = -r+(2.0/3.0)*n;

            double reduced_area = fraction[i]*part_rotated.surfaces[i].area*cos_theta_i;
            // instant remission
            // C_op += -reduced_area*((ca_op + cd_op)*dummy+2*cos_theta_i*cs_op*n);
            // C_ir += -reduced_area*((ca_ir + cd_ir)*dummy+2*cos_theta_i*cs_ir*n);
            // no remission
            C_op += -reduced_area*(-r*(ca_op + cd_op) + (2.0/3.0)*n*cd_op + 2*cos_theta_i*cs_op*n);
            C_ir += -reduced_area*(-r*(ca_ir + cd_ir) + (2.0/3.0)*n*cd_op + 2*cos_theta_i*cs_ir*n);

        }
        
        
        return {C_op(0), C_op(1), C_op(2), C_ir(0), C_ir(1), C_ir(2)};
};

std::vector<double> rp_coefficients_surface(const surface& surf, const double& fraction, const Eigen::Vector3d& v) {

    double ca_op = surf.surface_properties.optical[0];
    double cd_op = surf.surface_properties.optical[1];
    double cs_op = surf.surface_properties.optical[2];
    double ca_ir = surf.surface_properties.infrared[0];
    double cd_ir = surf.surface_properties.infrared[1];
    double cs_ir = surf.surface_properties.infrared[2];
    Eigen::Vector3d n = surf.normal;
            
    double cos_theta_i = v.dot(n);
    // Eigen::Vector3d dummy = -r+(2.0/3.0)*n;

    double reduced_area = fraction*surf.area*std::max(cos_theta_i, 0.0);
    // instant remission
    // C_op += -reduced_area*((ca_op + cd_op)*dummy+2*cos_theta_i*cs_op*n);
    // C_ir += -reduced_area*((ca_ir + cd_ir)*dummy+2*cos_theta_i*cs_ir*n);
    // no remission
    Eigen::Vector3d C_op = -reduced_area*(-v*(ca_op + cd_op) + (2.0/3.0)*n*cd_op + 2*cos_theta_i*cs_op*n);
    Eigen::Vector3d C_ir = -reduced_area*(-v*(ca_ir + cd_ir) + (2.0/3.0)*n*cd_op + 2*cos_theta_i*cs_ir*n);
    return {C_op(0), C_op(1), C_op(2), C_ir(0), C_ir(1), C_ir(2)};
};


