#include "utils.h"

// std::vector<std::vector<double>> load_data(const std::string& filename) {
//     std::ifstream file(filename);
//     if (!file) {
//         std::cerr << "Error: Could not open file " << filename << std::endl;
//         return {};  // Return an empty matrix
//     }

//     std::vector<std::vector<double>> matrix;
//     std::string line;

//     // Read file line by line
//     while (std::getline(file, line)) {
//         std::stringstream ss(line);
//         double value;
//         std::vector<double> row;

//         while (ss >> value) {
//             row.push_back(value);
//         }

//         // Store only non-empty rows
//         if (!row.empty()) {
//             matrix.push_back(row);
//         }
//     }

//     file.close();
//     return matrix;
// }

std::vector<std::vector<double>> load_data(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {};
    }

    int rows, cols;
    std::string line;

    // Read the first line to get matrix dimensions
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        char percent;
        iss >> percent >> rows >> cols; // Read "% 1 98" for example
    } else {
        std::cerr << "Error: File is empty or missing header" << std::endl;
        return {};
    }

    // Initialize matrix
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));

    int row = 0;
    while (std::getline(file, line) && row < rows) {
        std::istringstream iss(line);
        for (int col = 0; col < cols; ++col) {
            if (!(iss >> matrix[row][col])) {
                std::cerr << "Error: Not enough data in row " << row + 1 << std::endl;
                return {};
            }
        }
        row++;
    }

    file.close();
    return matrix;
}

void save_data(const std::vector<std::vector<double>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j < row.size() - 1) file << "\t";  // Separate columns by tab
        }
        file << "\n";  // Newline after each row
    }

    file.close();
}

