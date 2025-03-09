#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include "utils.h"
#include "geometry.h"
#include "shadow_calculator.h"
#include <eigen/Eigen/Dense>
#include <random>
#include <chrono>
#include <thread>
#include <algorithm>


using namespace std;

int main(){
    string data_path = "../data/MRO_lowfidelity_";

    Part_v1 MRO(data_path);
    int q_max = 100;

    // MRO.print_info();
    // MRO.print_suface_data(1);
    Eigen::Vector3d v_test(0.0, 1.0, 0.0);
    Eigen::Matrix3d R_unity(3, 3);
    R_unity(0,0) = 1.0;
    R_unity(0,1) = 0.0;
    R_unity(0,2) = 0.0;

    R_unity(1,0) = 0.0;
    R_unity(1,1) = 1.0;
    R_unity(1,2) = 0.0;

    R_unity(2,0) = 0.0;
    R_unity(2,1) = 0.0;
    R_unity(2,2) = 1.0;
    
    
    Eigen::Vector3d v(-1.0, 0.0, 0.0);

    string path_test = "C:/Users/nike/Documents/ThesisProject/cpp_scratchpad/input/MRO_lowfidelity.dae";
    vector<vector<double>> materials(6);
    materials[0] = {0.9, 0.07, 0.03};
    materials[1] = {0.9, 0.07, 0.03};
    materials[2] = {0.14, 0.80, 0.06};
    materials[3] = {0.14, 0.80, 0.06};
    materials[4] = {0.14, 0.80, 0.06};
    materials[5] = {0.9, 0.05, 0.05};
    
    Part MRO_test(path_test, materials);
    Assembly MRO_test_assembly(MRO_test, MRO_test);
    Part assembled_MRO = MRO_test_assembly.combine_parts(R_unity, R_unity);
    cout<<assembled_MRO.surfaces[0].l<<endl;
    cout<<assembled_MRO.surfaces[98].l<<endl;



    // Part_v1 MRO_rotated(MRO, R);
    // MRO_rotated.print_info();
    // MRO_rotated.print_suface_data(1);


    /*
    TESTING SINGLE THREADING VS MULTI THREADING PIXELATION 
    */

    //single threading
    // auto t1_single = std::chrono::high_resolution_clock::now();
    // std::vector<double> f_single = compute_shadow(MRO, R, v, q_max, false);
    // auto t2_single = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> ms_double_single = t2_single - t1_single;
    

    // multiple threading
    // auto t1_multi = std::chrono::high_resolution_clock::now();
    // std::vector<double> f_multi = compute_shadow(MRO, v, q_max, true);
    // auto t2_multi = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> ms_double_multi = t2_multi - t1_multi;
    
    // cout<<"FRACTIONS SINGLE"<<endl;
    // for (int i = 0; i<f_single.size(); i++) {
    //     cout<<i<<"   "<<f_single[i]<<endl;
    // }
    // cout<<"FRACTIONS MULTI"<<endl;
    // for (int i = 0; i<f_multi.size(); i++) {
    //     cout<<i<<"   "<<f_multi[i]<<endl;
    // }
    // cout<<"Single: "<< ms_double_single.count()<<endl;
    // cout<<"Multi: "<< ms_double_multi.count()<<endl;

    /*
    VERIFICATION
    */
    // double pi = M_PI;
    // // vector<double> alpha_angles = {-pi/2, 0, pi/2};
    // vector<double> alpha_angles;
    // double val1 = -pi/2;
    // while (val1<= pi/2) {
    //     alpha_angles.push_back(val1);
    //     val1 += pi/36; 
    // }
    // // vector<double> beta_angles = {-pi, -pi/2, 0, pi/2, pi};
    // vector<double> beta_angles;
    // double val2 = -pi;
    // while (val2<= pi) {
    //     beta_angles.push_back(val2);
    //     val2 += pi/36; 
    // }
    
    // cout<<"size alpha: "<<alpha_angles.size()<<endl;
    // cout<<"size beta: "<<beta_angles.size()<<endl;
    // cout<<"total: "<<beta_angles.size()*alpha_angles.size()<<endl;
    // vector<Eigen::Matrix3d> Rs(alpha_angles.size()*beta_angles.size(), Eigen::Matrix3d(3, 3));
    // vector<vector<double>> angles(alpha_angles.size()*beta_angles.size(), vector<double>(2));
    // for (int i = 0; i<alpha_angles.size(); i++) {
    //     for (int j = 0; j<beta_angles.size(); j++) {
    //         angles[i*beta_angles.size()+j] = {alpha_angles[i], beta_angles[j]};
    //         Rs[i*beta_angles.size()+j] = R_body_wind(alpha_angles[i], beta_angles[j]);
    //     }
    // }
    // save_data(angles, "../output/a_2.txt");
    // cout<<"angles"<<endl;
    // cout<<angles[5292][0]<<"  "<<angles[5292][1]<<endl;
    // cout<<"rs[5292]"<<endl;
    // cout<<Rs[5292]<<endl;
    // cout<<"R(angles[5292])"<<endl;
    // cout<<R_body_wind(angles[5292][0], angles[5292][1])<<endl;
    // Eigen::Matrix3d R_buggata = R_body_wind(1e-30, 1e-30);
    // // R_buggata(0, 0) = 1.0;
    // // R_buggata(1, 1) = 1.0;
    // // R_buggata(2, 2) = 1.0;
    // R_buggata(0, 1) = 0.0;
    // cout<<"R_buggata"<<endl;
    // cout<<R_buggata<<endl;
    // vector<double> f_test_bug = compute_shadow(MRO, R_buggata, v_test, q_max, false);
    // vector<double> coeff_bug = rp_coefficients(MRO, R_buggata, f_test_bug, v_test);
    // cout<<"R_buggata"<<endl;
    // cout<<coeff_bug[0]<<endl;
    // cout<<coeff_bug[1]<<endl;
    // cout<<coeff_bug[2]<<endl;
    // vector<double> f_test = compute_shadow(MRO, Rs[5292], v_test, q_max, false);
    // vector<double> coeff = rp_coefficients(MRO, Rs[5292], f_test, v_test);
    // cout<<"Rs rotation coeffs"<<endl;
    // cout<<coeff[0]<<endl;
    // cout<<coeff[1]<<endl;
    // cout<<coeff[2]<<endl;
    // vector<double> f_test_unity = compute_shadow(MRO, R_unity, v_test, q_max, false);
    // vector<double> coeff_unity = rp_coefficients(MRO, R_unity, f_test_unity, v_test);
    // cout<<"unity rotation coeffs"<<endl;
    // cout<<coeff_unity[0]<<endl;
    // cout<<coeff_unity[1]<<endl;
    // cout<<coeff_unity[2]<<endl;
    // cout<<"unity rotation: "<<endl;
    // cout<<R_unity*v_test<<endl;
    // cout<<"bugged rotation: "<<endl;
    // cout<<R_buggata*v_test<<endl;


    // Part_v1 MRO_rotated(MRO, Rs[5293]);
    // cout<<"FRACTIONS"<<endl;
    // cout<<f_test_unity.size()<<"  "<<f_test.size()<<endl;
    // for (int i = 0; i<f_test_unity.size(); i++) {
    //     cout<<i<<" unity f: "<<f_test_unity[i]<<" rs f: "<<f_test[i]<<
    //     " Cx unity: "<<rp_coefficients_surface(MRO.surfaces[i], f_test_unity[i], v_test)[0]<<
    //     " Cx rs: "<<rp_coefficients_surface(MRO_rotated.surfaces[i], f_test[i], v_test)[0]<<endl;
    // }
    // surface surf = MRO_rotated.surfaces[6];
    // double ca_op = surf.surface_properties.optical[0];
    // double cd_op = surf.surface_properties.optical[1];
    // double cs_op = surf.surface_properties.optical[2];
    // double ca_ir = surf.surface_properties.infrared[0];
    // double cd_ir = surf.surface_properties.infrared[1];
    // double cs_ir = surf.surface_properties.infrared[2];
    // Eigen::Vector3d n = surf.normal;
            
    // double cos_theta_i = v_test.dot(n);
    // cout<<"cos theta: "<<cos_theta_i<<endl;
    // // Eigen::Vector3d dummy = -r+(2.0/3.0)*n;

    // double reduced_area = f_test[6]*surf.area;
    // cout<<"reduced area: "<<reduced_area<<endl;
    // cout<<(-v_test*(ca_op + cd_op))[0]<<endl;
    // cout<<((2.0/3.0)*cd_op*n)[0]<<endl;
    // cout<<(2*cos_theta_i*cs_op*n)[0]<<endl;
    // cout<<"normal x: "<<n[0]<<endl;
    // cout<<"normal y: "<<n[1]<<endl;
    // cout<<"normal z: "<<n[2]<<endl;
    // cout<<"normal x: "<<MRO.surfaces[6].normal[0]<<endl;
    // cout<<"normal y: "<<MRO.surfaces[6].normal[1]<<endl;
    // cout<<"normal z: "<<MRO.surfaces[6].normal[2]<<endl;

    
    // Eigen::Vector3d C_op = -reduced_area*(-v*(ca_op + cd_op) + (2.0/3.0)*n*cd_op + 2*cos_theta_i*cs_op*n);
    // Eigen::Vector3d C_ir = -reduced_area*(-v*(ca_ir + cd_ir) + (2.0/3.0)*n*cd_op + 2*cos_theta_i*cs_ir*n);




    // vector<vector<double>> coefficients(Rs.size(), vector<double>(6, 0));
    // for (int i = 0; i<Rs.size(); i++) {
    //     cout<<"Iteration: "<<i+1<<endl;
    //     vector<double> f = compute_shadow(MRO, Rs[i], v, q_max, false);
    //     coefficients[i] = rp_coefficients(MRO, Rs[i], f, v);
    // }
    // save_data(coefficients, "../output/lowfidelity/c_v1.txt");

    // multi thread
    // int n_threads = thread::hardware_concurrency()-1;
    // vector<thread> threads;
    // vector<vector<vector<double>>> results(n_threads);
    // int chunk = static_cast<int>(round(Rs.size()/n_threads));
    // cout<<"cores: "<<n_threads<<endl;

    // for (int i = 0; i < n_threads; i++) {
    //     // Launch thread
    //     threads.emplace_back([&, i]() {
            
    //         int i_max;
    //         if (i==n_threads-1) {
    //            i_max = Rs.size(); 
    //         }
    //         else {
    //             i_max = min((i+1)*chunk, static_cast<int>(Rs.size()));
    //         }
    //         cout<<"starting thread: "<<i+1<<endl;
    //         cout<<"from: "<<i*chunk<<" to: "<<i_max<<endl;
    //         vector<vector<double>> partial_results(i_max-i*chunk);
    //         for (int j = i*chunk; j<i_max; j++) {
    //             vector<double> f = compute_shadow(MRO, Rs[j], v, q_max, false);
    //             partial_results[j-i*chunk] = rp_coefficients(MRO, Rs[j], f, v);
    //         }
    //         results[i] = partial_results;
    //     });
    // }
    
    // // Join all threads
    // for (auto& th : threads) {
    //     th.join();
    // }
    // vector<vector<double>> coefficients;
    // for (auto& thread_results : results) {
    //     for (auto& line : thread_results) {
    //         coefficients.push_back(line);
    //     }
    // }
    // save_data(coefficients, "../output/c_3.txt");
    /*
    TESTING PIP ALGORITHMS
    */

    // testing PIP algo
    // point2D test(0.1, 0.2);
    // point2D A(-0.5,-0.5), B(0.5,-0.2), C(0,0.5);
    // point2D A(0,0), B(10, 0), C(6, 6);
    // vector<point2D> points = {test};
    // vector<point2D> triangle = {A, B, C};

    // // creating the grid
    // vector<double> x_coord = linspace(0.1, 9.9, 5000);
    // vector<double> y_coord = linspace(0.1, 5.9, 5000);

    // vector<point2D> points;
    // for (int i = 0; i<y_coord.size(); i++) {
    //     for (int j = 0; j<x_coord.size(); j++) {
    //         points.emplace_back(x_coord[j], y_coord[i]);
    //     }
    // }

    // // raycast normal
    // auto t1_ray = std::chrono::high_resolution_clock::now();
    // vector<bool> test_ray = PIP_raycast_triangle(triangle, points);
    // auto t2_ray = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> ms_double_raycast = t2_ray - t1_ray;
    // cout<<"normal: "<< ms_double_raycast.count()<<endl;
    
    // // ray cast acc
    // auto t1_wind = std::chrono::high_resolution_clock::now();
    // vector<vector<bool>> test_ray_acc = PIP_raycast_triangle_acc(triangle, x_coord, y_coord);
    // auto t2_wind = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> ms_double_wind = t2_wind - t1_wind;
    // cout<<"acc: "<<ms_double_wind.count()<<endl;    


    // for (int i = 0; i<30; i++) {
    //     for (int j = 0; j<100; j++) {
    //         cout<<results[i][j];
    //     }
    //     cout<<endl;
    // }

    // cout<<PIP_raycast_triangle(triangle, points)[0]<<endl;
    // cout<<PIP_winding_triangle(triangle, points)[0]<<endl;

    // vector<double> n_points = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7};
    // // benchmarking
    // for (int k = 1; k<101; k++) {


    // vector<vector<double>> times;
    // cout<<"Index: "<<k<<endl;

    // for (int j = 0; j<n_points.size(); j++) {
        
    //     random_device rd;
    //     mt19937 engine(rd());
    //     uniform_real_distribution<double> dist(-1, 1);
    //     double n = n_points[j];
    //     vector<point2D> points;
    //     points.reserve(n);

    //     for (int i = 0; i < n; ++i) {
    //         double x = dist(engine);
    //         double y = dist(engine);
    //         points.emplace_back(x, y);
    //     }
        
    //     // raycast
    //     auto t1_ray = std::chrono::high_resolution_clock::now();
    //     vector<bool> test_ray = PIP_raycast_triangle(triangle, points);
    //     auto t2_ray = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double, std::milli> ms_double_raycast = t2_ray - t1_ray;
    //     // cout<<ms_double_raycast.count()<<endl;
    
    //     //winding nuber
    //     auto t1_wind = std::chrono::high_resolution_clock::now();
    //     vector<bool> test_wind = PIP_winding_triangle(triangle, points);
    //     auto t2_wind = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double, std::milli> ms_double_wind = t2_wind - t1_wind;
    //     // cout<<ms_double_wind.count()<<endl;

    //     double t_ray = ms_double_raycast.count();
    //     double t_wind = ms_double_wind.count();
    //     times.push_back({n, t_ray, t_wind});

    // }
    
    // save_data(times, "../output4/PIP_times" + std::to_string(k) + ".txt");
    // }
    
    return 0;

}