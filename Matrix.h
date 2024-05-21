#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <complex>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath> // for std::abs and std::sqrt


namespace algebra {

enum class StorageOrder { RowMajor, ColumnMajor };
enum class NormType { One, Infinity, Frobenius };

// Forward declaration of the Matrix class
template<typename T, StorageOrder Order>
class Matrix;

// operator* overload for matrix-vector multiplication
template<typename T, StorageOrder Order>
std::vector<T> operator*(const Matrix<T, Order>& A, const std::vector<T>& v) {
    // Check if the matrix and vector are compatible for multiplication
    if ((Order == StorageOrder::RowMajor && A.num_cols() != v.size()) ||
        (Order == StorageOrder::ColumnMajor && A.num_rows() != v.size()))
        throw std::invalid_argument("Matrix and vector dimensions are not compatible for multiplication.");
    
    std::vector<T> result(A.num_rows(), T{});  // Initialize the result vector with zeros as suggested

    if (A.is_compressed()) {
        if constexpr (Order == StorageOrder::RowMajor) {
            // CSR multiplication
            for (size_t i = 0; i < A.num_rows(); ++i) // row i
                for (size_t j = A.inner[i]; j < A.inner[i + 1]; ++j) // non-zero elements in row i
                    result[i] += A.values[j] * v[A.outer[j]];
        } else {
            // CSC multiplication
            auto valit = A.values.begin();
            for (size_t j = 0; j < A.num_cols(); ++j) // column j
                for (size_t i = A.inner[j]; i < A.inner[j + 1]; ++i) // non-zero elements in column j
                    result[A.outer[i]] += *(valit++) * v[j];
        }
    } else {
        // Uncompressed state multiplication: gradually add in place as suggested
        for (const auto& entry : A.data) {
            const auto& indices = entry.first;
            size_t i = indices[0];
            size_t j = indices[1];
            result[i] += entry.second * v[j];
        }
    }

    return result;
}
template<typename T, StorageOrder Order1, StorageOrder Order2>
std::vector<T> operator*(const Matrix<T, Order1>& A, const Matrix<T, Order2>& v) {
    // Ensure the second matrix is a vector-like matrix (one column)
    if (v.num_cols() != 1)
        throw std::invalid_argument("Second operand is not a vector (must have exactly one column).");

    // 'convert' the vector-like matrix to a vector, and use it for multiplication
    std::vector<T> vec(v.num_rows());
    for (size_t i = 0; i < v.num_rows(); ++i)
        vec[i] = v(i, 0);  // Copy each element from the matrix to the vector

    return A * vec;
}

template<typename T, StorageOrder Order>
void readMTX(Matrix<T, Order>& matrix, const std::string& filename) {
    // Read a matrix from a Matrix Market file (.mtx) - example:
    // %%MatrixMarket matrix coordinate real general
    // % Rows Columns Entries
    // 131 131 536
    // 1 1  1.0000000000000e+00
    // 9 1 -6.3869481600000e+00
    // 131 1  6.3554020600000e-02
    // ...

    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open file: " + filename);

    std::string line;
    bool header_passed = false;
    size_t rows, cols, entries;

    while (getline(file, line)) {
        if (line[0] == '%') // Comment line
            continue;

        if (!header_passed) {
            std::istringstream iss(line);
            if (!(iss >> rows >> cols >> entries))
                throw std::runtime_error("Invalid Matrix Market header.");
            matrix.resize(rows, cols);
            header_passed = true;
        } else {
            std::istringstream iss(line);
            size_t row, col;
            T value;
            if (!(iss >> row >> col >> value))
                throw std::runtime_error("Invalid data format.");
            matrix(row-1, col-1) = value; // adjust indices to 0-based as it all should be
        }
    }

    file.close();
}

// Matrix class definition
template<typename T, StorageOrder Order>
class Matrix {
    // Matrix properties
    size_t rows, cols; // Number of rows and columns
    bool compressed; // Flag to indicate if the matrix is compressed

    // Data storage:
    // - For non-compressed storage
    std::map<std::array<size_t, 2>, T> data; // Map of (row, col) -> value
    // - For compressed storage
    std::vector<T> values; // Non-zero values
    std::vector<size_t> outer; // Column(CSR) or row(CSC) indices
    std::vector<size_t> inner; // Row(CSR) or column(CSC) pointers

    friend std::vector<T> operator*<>(const Matrix<T, Order>& A, const std::vector<T>& v);
    friend std::vector<T> operator*<>(const Matrix<T, Order>& A, const Matrix<T, Order>& v);
    friend void readMTX<>(Matrix<T, Order>&, const std::string&);

public:
    // constructor
    Matrix(size_t r = 0, size_t c = 0);

    // getters
    size_t num_rows() const { return rows; }
    size_t num_cols() const { return cols; }

    // access/set operators
    T operator()(size_t i, size_t j) const;
    T& operator()(size_t i, size_t j);

    // methods
    void compress();
    void uncompress();
    bool is_compressed() const { return compressed; }
    void resize(size_t nrows, size_t ncols);

    template<NormType type>
    T norm() const;

    void erase(size_t i, size_t j);
    bool contains(size_t i, size_t j) const;
    void clear();

    void print() const; // Print the matrix
    void printBones() const; // Print the active data structure

private:
    T oneNorm() const;
    T infinityNorm() const;
    T frobeniusNorm() const;
};





// definitions - "Matrix.tpp"
   template<typename T, StorageOrder Order>
    Matrix<T, Order>::Matrix(std::size_t rows, std::size_t cols)
        : rows(rows), cols(cols), compressed(false) {}
    
    // Accessor and mutator methods
    template<typename T, StorageOrder Order>
    T Matrix<T,Order>::operator()(size_t i, size_t j) const {
        if (i >= rows || j >= cols)
            throw std::out_of_range("Matrix indices are out of bounds.");

        if (compressed) {
            if constexpr (Order == StorageOrder::RowMajor) {
                // Search in CSR format
                size_t start = inner[i];
                size_t end = inner[i + 1];
                for (size_t pos = start; pos < end; ++pos) {
                    if (outer[pos] == j) {
                        return values[pos];
                    }
                }
            } else {
                // Search in CSC format
                size_t start = inner[j];
                size_t end = inner[j + 1];
                for (size_t pos = start; pos < end; ++pos) {
                    if (outer[pos] == i) {
                        return values[pos];
                    }
                }
            }
            return T{}; // Return zero if element is not explicitly stored
        } else {
            // Search in uncompressed format
            auto it = data.find({i, j});
            return it != data.end() ? it->second : 0;
        }
    }
    template<typename T, StorageOrder Order>
    T& Matrix<T,Order>::operator()(size_t i, size_t j) {
        if (i >= rows || j >= cols)
            throw std::out_of_range("Matrix indices are out of bounds.");
        if (compressed){
            // Allow modification only if the element exists in compressed format
            size_t out_idx = (Order == StorageOrder::RowMajor) ? i : j;
            size_t in_idx = (Order == StorageOrder::RowMajor) ? j : i;
            size_t start = outer[out_idx];
            size_t end = outer[out_idx + 1];
            for (size_t k = start; k < end; ++k)
                if (inner[k] == in_idx)
                    return values[k]; // Return reference to the existing non-zero element

            throw std::logic_error("Adding new elements in compressed state is not allowed.");
        }

        return data[{i, j}];
    }

    // resize
    template<typename T, StorageOrder Order>
    void Matrix<T,Order>::resize(size_t r, size_t c) {
        // I decided to throw an exception if any existing elements would become out of bounds, instead of cutting them off
        // this means simply check for oob elements, if we arrive at the end simply update rows and cols (and cut inner if compressed)
        // works both for compressed and uncompressed states
        if (compressed){
            if constexpr (Order == StorageOrder::RowMajor) {
                for (size_t i = 0; i < rows; i++) {
                    for (size_t pos = inner[i]; pos < inner[i + 1]; pos++) {
                        size_t col = outer[pos];
                        if (col >= c)
                            throw std::logic_error("Resize would result in out of bounds element: (" + std::to_string(i) + ", " + std::to_string(col) + ")");
                    }
                }
                // cut inner
                if (r < rows)
                    inner.resize(r + 1);
            }
            else {
                for (size_t i = 0; i < cols; i++) {
                    for (size_t pos = inner[i]; pos < inner[i + 1]; pos++) {
                        size_t row = outer[pos];
                        if (row >= r)
                            throw std::logic_error("Resize would result in out of bounds element: (" + std::to_string(row) + ", " + std::to_string(i) + ")");
                    }
                }
                // cut inner
                if (c < cols)
                    inner.resize(c + 1);
            }
        } else {
            // uncompressed
            for (const auto& val : data) {
                size_t row = val.first[0];
                size_t col = val.first[1];
                if (row >= r || col >= c)
                    throw std::logic_error("Resize would result in out of bounds element: (" + std::to_string(row) + ", " + std::to_string(col) + ")");
            }
        }

        rows = r;
        cols = c;
    }

    // compressions
    template<typename T, StorageOrder Order>
    void Matrix<T,Order>::compress() {
        // Compress the matrix into CSR or CSC format
        if (compressed)
            throw std::logic_error("Matrix is already compressed.");

        if constexpr (Order == StorageOrder::RowMajor) {
            // Compress into CSR format
            inner.resize(rows + 1, 0);
            
            // count the non-zero entries, for each row
            for (const auto& entry : data)
                inner[entry.first[0] + 1]++;
            
            // Convert counts to indices
            for (size_t i = 1; i < inner.size(); i++)
                inner[i] += inner[i - 1];
            
            // Temporary vector to keep track of the current position in each row
            values.resize(data.size());
            outer.resize(data.size());
            std::vector<size_t> current_pos(inner.begin(), inner.end() - 1);
            for (const auto& entry : data) {
                size_t row = entry.first[0];
                size_t col = entry.first[1];
                size_t pos = current_pos[row]++;
                values[pos] = entry.second;
                outer[pos] = col;
            }
        } else {
            // Compress into CSC format
            inner.resize(cols + 1, 0);
            
            // count the non-zero entries, for each column
            for (const auto& entry : data)
                inner[entry.first[1] + 1]++;
            
            // Convert counts to indices
            for (size_t i = 1; i < inner.size(); i++)
                inner[i] += inner[i - 1];
            
            // Temporary vector to keep track of the current position in each column
            std::vector<size_t> current_pos(inner.begin(), inner.end() - 1);
            values.resize(data.size());
            outer.resize(data.size());
            for (const auto& entry : data) {
                size_t row = entry.first[0];
                size_t col = entry.first[1];
                size_t pos = current_pos[col]++;
                values[pos] = entry.second;
                outer[pos] = row;
            }
        }
        compressed = true;
        data.clear(); // Clear dynamic storage
    }
    template<typename T, StorageOrder Order>
    void Matrix<T,Order>::uncompress() {
        // Bring back the matrix to uncompressed state
        if (!compressed)
            throw std::logic_error("Matrix is already uncompressed.");

        data.clear(); // just making sure

        if constexpr (Order == StorageOrder::RowMajor) {
            // Decompress from CSR format
            for (size_t i = 0; i < rows; i++) {
                for (size_t pos = inner[i]; pos < inner[i + 1]; pos++) {
                    size_t col = outer[pos];
                    T value = values[pos];
                    data[{i, col}] = value;
                }
            }
        } else {
            // Decompress from CSC format
            for (size_t i = 0; i < cols; i++) {
                for (size_t pos = inner[i]; pos < inner[i + 1]; pos++) {
                    size_t row = outer[pos];
                    T value = values[pos];
                    data[{row, i}] = value;
                }
            }
        }

        // Clear compressed data to free memory
        values.clear();
        outer.clear();
        inner.clear();
        compressed = false;        
    }

    // erase, contains, clear
    template<typename T, StorageOrder Order>
    void Matrix<T,Order>::erase(size_t i, size_t j) {
        // remove an element
        if (i >= rows || j >= cols)
            throw std::out_of_range("Cannot erase element: indices are out of bounds.");
        if (compressed)
            throw std::logic_error("Cannot erase elements in compressed state.");
        data.erase({i, j});
    }
    template<typename T, StorageOrder Order>
    bool Matrix<T,Order>::contains(size_t i, size_t j) const {
        // check the presence of an element
        return data.find({i, j}) != data.end();
    }
    template<typename T, StorageOrder Order>
    void Matrix<T,Order>::clear() {
        data.clear();
        values.clear();
        outer.clear();
        inner.clear();
        compressed = false;
    }


    // prints
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::print() const {
        std::cout << "Matrix (" << rows << "x" << cols << "): (" 
                  << (Order == StorageOrder::RowMajor ? "RowMajor" : "ColumnMajor") << ", " 
                  << (compressed ? "compressed" : "not compressed") << ")\n";
        // just use access operator for simple printing :)
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j)
                std::cout << (*this)(i, j) << " \t";
            std::cout << std::endl;
        }
    }
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::printBones() const {
        std::cout << "Matrix (" << rows << "x" << cols << ") bones: (" 
                  << (Order == StorageOrder::RowMajor ? "RowMajor" : "ColumnMajor") << ", "
                  << (compressed ? "compressed" : "not compressed") << ")\n";
        if (compressed) {
            std::cout << " --Values: ";
            for (const auto& value : values)
                std::cout << value << " ";

            std::cout << "\n --Inner indices: (" << (Order==StorageOrder::RowMajor ? "Rows" : "Cols") << " starting indices) ";
            for (const auto& index : outer)
                std::cout << index << " ";

            std::cout << "\n --Outer indices: (" << (Order==StorageOrder::RowMajor ? "Cols" : "Rows") << " indices) ";
            for (const auto& ptr : inner)
                std::cout << ptr << " ";

            std::cout << "\n";
        } else {
            for (const auto& pair : data)
                std::cout << " - (" << pair.first[0] << ", " << pair.first[1] << "): " << pair.second << "\n";
        }
    }


    // norm
    template<typename T, StorageOrder Order>
    template<NormType type>
    T Matrix<T,Order>::norm() const {
        if constexpr (type == NormType::One)
            return oneNorm();
        else if constexpr (type == NormType::Infinity)
            return infinityNorm();
        else if constexpr (type == NormType::Frobenius)
            return frobeniusNorm();
        else
            throw std::invalid_argument("Invalid norm type.");
    }
    template<typename T, StorageOrder Order>
    T Matrix<T,Order>::oneNorm() const {
        // we store col sums, and return the max
        std::vector<T> columnSums(cols, T{});
        if (compressed) {
            if constexpr (Order == StorageOrder::RowMajor) {
                // For CSR: traverse...
                for (size_t i = 0; i < rows; ++i) {
                    size_t start = inner[i], end = inner[i + 1];
                    for (size_t pos = start; pos < end; ++pos)
                        columnSums[outer[pos]] += std::abs(values[pos]);
                }
            } else {
                // For CSC: directly calculate since it's the natural order in this case
                for (size_t j = 0; j < cols; ++j) 
                    for (size_t i = inner[j]; i < inner[j + 1]; ++i)
                        columnSums[j] += std::abs(values[i]);
            }
        } else {
            // Uncompressed: just get through the data
            for (const auto& entry : data)
                columnSums[entry.first[1]] += std::abs(entry.second);
        }
        return *std::max_element(columnSums.begin(), columnSums.end());
    }
    template<typename T, StorageOrder Order>
    T Matrix<T,Order>::infinityNorm() const {
        // we store row sums, and return the max
        std::vector<T> rowSums(rows, T{});
        if (compressed) {
            if constexpr (Order == StorageOrder::RowMajor) {
                // For CSR: directly calculate since it's the natural order for infinity norm in CSR
                for (size_t i = 0; i < rows; ++i)
                    for (size_t pos = inner[i]; pos < inner[i + 1]; ++pos)
                        rowSums[i] += std::abs(values[pos]);
            } else {
                // For CSC: traverse...
                for (size_t j = 0; j < cols; ++j) {
                    size_t start = inner[j], end = inner[j + 1];
                    for (size_t pos = start; pos < end; ++pos)
                        rowSums[outer[pos]] += std::abs(values[pos]);
                }
            }
        } else {
            for (const auto& entry : data)
                rowSums[entry.first[0]] += std::abs(entry.second);
        }
        return *std::max_element(rowSums.begin(), rowSums.end());
    }
    template<typename T, StorageOrder Order>
    T Matrix<T,Order>::frobeniusNorm() const {
        // sum of squares of all elements
        // std::norm returns the squared magnitude of the complex number
        T sum = T{};
        if (compressed) {
            for (const T& value : values)
                sum += std::norm(value);
        } else {
            for (const auto& entry : data)
                sum += std::norm(entry.second);
        }
        return std::sqrt(sum);
    }
} // namespace algebra

#endif // MATRIX_H
