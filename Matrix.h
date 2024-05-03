#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <complex>
#include <algorithm>


namespace algebra {

enum class StorageOrder { RowMajor, ColumnMajor };

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
        if (Order == StorageOrder::RowMajor) {
            // Compressed Sparse Row (CSR) multiplication
            for (size_t i = 0; i < A.num_rows(); ++i) {
                for (size_t j = A.inner[i]; j < A.inner[i + 1]; ++j) {
                    size_t col = A.outer[j];
                    result[i] += A.values[j] * v[col];
                }
            }
        } else {
            // Compressed Sparse Column (CSC) multiplication
            for (size_t col = 0; col < A.num_cols(); ++col) {
                T vj = v[col];
                if (vj != T{})
                    for (size_t idx = A.inner[col]; idx < A.inner[col + 1]; ++idx)
                        result[A.outer[idx]] += A.values[idx] * vj;
            }
        }
    } else {
        // Uncompressed state multiplication
        for (const auto& entry : A.data) {
            const auto& indices = entry.first;
            size_t i = indices[0];
            size_t j = indices[1];
            result[i] += entry.second * v[j];
        }
    }

    return result;
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
    std::vector<size_t> inner; // Row or column pointers
    std::vector<size_t> outer; // Column or row indices

    friend std::vector<T> operator*<>(const Matrix<T, Order>& A, const std::vector<T>& v);

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
    void resize(size_t nrows, size_t ncols);
    void compress();
    void uncompress();
    bool is_compressed() const { return compressed; }
    void print() const; // Print the matrix
    void printBones() const; // Print the active data structure

private:
    void clear();
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
            if (Order == StorageOrder::RowMajor) {
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
        if (compressed) //TODO allow edit for existing non-zero elements
            throw std::logic_error("Modification not allowed in compressed state.");

        // TODO If the element to set is zero, remove it from the map instead... how

        return data[{i, j}];
    }

    // resize
    template<typename T, StorageOrder Order>
    void Matrix<T,Order>::resize(size_t r, size_t c) {
        // I decided to only allow resizing for uncompressed matrices (not very clear directions - maybe I should decmpress first?)
        // Also I decided to throw an exception if any existing elements would become out of bounds, instead of cutting them off
        // TODO improve
        if (compressed)
            throw std::logic_error("Resize not allowed in compressed state.");
        
        // Check if any existing elements would become out of bounds
        for (const auto& val : data) {
            size_t row = val.first[0];
            size_t col = val.first[1];
            if (row >= r || col >= c)
                throw std::logic_error("Resize would result in out of bounds element: (" + std::to_string(row) + ", " + std::to_string(col) + ")");
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

        if (Order == StorageOrder::RowMajor) {
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

        if (Order == StorageOrder::RowMajor) {
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

    // clear
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
        std::cout << "Matrix (" << rows << "x" << cols << "):\n";
        // just use access operator :)
        for (std::size_t i = 0; i < rows; ++i) {
            for (std::size_t j = 0; j < cols; ++j)
                std::cout << (*this)(i, j) << " \t";
            std::cout << std::endl;
        }
    }
    template<typename T, StorageOrder Order>
    void Matrix<T, Order>::printBones() const {
        std::cout << "Matrix bones: sz=(" << rows << ", " << cols 
                  << ") - The matrix is " << (compressed ? "compressed" : "not compressed") << "\n";
        if (compressed) {
            std::cout << " - Values: ";
            for (const auto& value : values)
                std::cout << value << " ";

            std::cout << "\n - Inner indices: (" << (Order==StorageOrder::RowMajor ? "Rows" : "Cols") << " starting indices) ";
            for (const auto& index : outer)
                std::cout << index << " ";

            std::cout << "\n - Outer indices: (" << (Order==StorageOrder::RowMajor ? "Cols" : "Rows") << " indices) ";
            for (const auto& ptr : inner)
                std::cout << ptr << " ";

            std::cout << "\n";
        } else {
            for (const auto& pair : data)
                std::cout << " - (" << pair.first[0] << ", " << pair.first[1] << "): " << pair.second << "\n";
        }
    }


} // namespace algebra

#endif // MATRIX_H