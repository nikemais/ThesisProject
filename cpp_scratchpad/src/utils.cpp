#include "utils.h"

std::vector<std::vector<double>> load_data(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return {};  // Return an empty matrix
    }

    std::vector<std::vector<double>> matrix;
    std::string line;

    // Read file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double value;
        std::vector<double> row;

        while (ss >> value) {
            row.push_back(value);
        }

        // Store only non-empty rows
        if (!row.empty()) {
            matrix.push_back(row);
        }
    }

    file.close();
    return matrix;
}

