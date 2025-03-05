#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <variant>
#include <fstream>
#include <iomanip>  // For std::setw

#include "geometry.h"

std::vector<std::vector<double>> load_data(const std::string& filename);



template <typename T>
void print_matrix(const std::vector<std::vector<T>>& matrix, int width = 2) {
    for (int i = matrix.size() - 1; i >= 0; --i) {  // Iterate from last row to first
        for (const auto& elem : matrix[i]) {
            std::cout << std::setw(width) << elem << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); 

  return linspaced;
}

template <typename T> int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

void save_data(const std::vector<std::vector<double>>& matrix, const std::string& filename);

#endif