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

    Part MRO(data_path);
    // MRO.print_info();
    // MRO.print_suface_data(1);

    Eigen::Vector3d v(-1.0, 0.0, 0.0);
    int q_max = 100;
    // Eigen::Matrix3d R;
    // R(0, 0) = 1.0; R(0, 1) = 0.0; R(0, 2) = 0.0; 
    // R(1, 0) = 0.0; R(1, 1) = 1.0; R(1, 2) = 0.0; 
    // R(2, 0) = 0.0; R(2, 1) = 0.0; R(2, 2) = 1.0; 

    // Part MRO_rotated(MRO, R);
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
    double pi = acos(-1.0);
    // vector<double> alpha_angles = {-pi/2, 0, pi/2};
    vector<double> alpha_angles;
    double val1 = -pi/2;
    while (val1<= pi/2) {
        alpha_angles.push_back(val1);
        val1 += pi/72; 
    }
    // vector<double> beta_angles = {-pi, -pi/2, 0, pi/2, pi};
    vector<double> beta_angles;
    double val2 = -pi;
    while (val2<= pi) {
        beta_angles.push_back(val2);
        val2 += pi/72; 
    }

    cout<<"size alpha: "<<alpha_angles.size()<<endl;
    cout<<"size beta: "<<beta_angles.size()<<endl;
    cout<<"total: "<<beta_angles.size()*alpha_angles.size()<<endl;
    vector<Eigen::Matrix3d> Rs(alpha_angles.size()*beta_angles.size());
    vector<vector<double>> angles(alpha_angles.size()*beta_angles.size(), vector<double>(2));
    for (int i = 0; i<alpha_angles.size(); i++) {
        for (int j = 0; j<beta_angles.size(); j++) {
            angles[i*beta_angles.size()+j] = {alpha_angles[i], beta_angles[j]};
            Rs[i*beta_angles.size()+j] = R_body_wind(alpha_angles[i], beta_angles[j]);
        }
    }
    save_data(angles, "../output/lowfidelity/a_v3.txt");
    // vector<vector<double>> coefficients(Rs.size(), vector<double>(6, 0));
    // for (int i = 0; i<Rs.size(); i++) {
    //     cout<<"Iteration: "<<i+1<<endl;
    //     vector<double> f = compute_shadow(MRO, Rs[i], v, q_max, false);
    //     coefficients[i] = rp_coefficients(MRO, Rs[i], f, v);
    // }
    // save_data(coefficients, "../output/lowfidelity/c_v1.txt");

    // multi thread
    int n_threads = thread::hardware_concurrency()-1;
    vector<thread> threads;
    vector<vector<vector<double>>> results(n_threads);
    int chunk = static_cast<int>(round(Rs.size()/n_threads));
    cout<<"cores: "<<n_threads<<endl;

    for (int i = 0; i < n_threads; i++) {
        // Launch thread
        threads.emplace_back([&, i]() {
            
            int i_max;
            if (i==n_threads-1) {
               i_max = Rs.size(); 
            }
            else {
                i_max = min((i+1)*chunk, static_cast<int>(Rs.size()));
            }
            cout<<"starting thread: "<<i+1<<endl;
            cout<<"from: "<<i*chunk<<" to: "<<i_max<<endl;
            vector<vector<double>> partial_results(i_max-i*chunk);
            for (int j = i*chunk; j<i_max; j++) {
                vector<double> f = compute_shadow(MRO, Rs[j], v, q_max, false);
                partial_results[j-i*chunk] = rp_coefficients(MRO, Rs[j], f, v);
            }
            results[i] = partial_results;
        });
    }
    
    // Join all threads
    for (auto& th : threads) {
        th.join();
    }
    vector<vector<double>> coefficients;
    for (auto& thread_results : results) {
        for (auto& line : thread_results) {
            coefficients.push_back(line);
        }
    }
    save_data(coefficients, "../output/lowfidelity/c_v3.txt");
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