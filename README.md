# PACS Challenge 2 - A Sparse Matrix

This repo implements a simple Sparse Matrix class in C++. The challenge required the implementation to be capable of handling both compressed and uncompressed matrix formats (CSR/CSC), supporting operations with vectors, and reading matrices in the [matrix market format](https://math.nist.gov/MatrixMarket/formats.html#mtx).

## Matrix Implementation

The `Matrix` class handles sparse matrices efficiently with both compressed and uncompressed storage. The class is templated with a data type `T` and a storage order (`RowMajor` or `ColumnMajor`). The implementation is designed to optimize operations depending on the matrix's storage state (compressed or uncompressed) and storage order.

The `Matrix.h` file contains the declaration and definition of the Matrix class. The definitions are not in a separated file due to the use of templates, this is one of the suggested methods to handle this in the lectures. The long methods definitions are put at the end of file for clarity, except for the friends methods that are put on top.

### Data structure

The class maintains data either in a `std::map` for uncompressed storage or using three `std::vector`s for compressed storage, reflecting CSR or CSC formats.

### Methods

Among the methods:

- **operator()(size_t i, size_t j)**:
    - (const) Accesses the element at row `i` and column `j`. Throws `std::out_of_range` if the indices are out of bounds.
    - (non-const) Edit the matrix element. Note that, in the current implementation, if the element to set is 0 it is still added to the storage, since there is no way to check what values is to be set (we are only passing the reference). This could possibly be solved with a method to check and remove zeros if present, called by other methods, or maybe with a proxy class, but I decided not to complicate things. Instead of setting an element to 0, the user should use the `erase(i,j)` method.
- **compress()**: Converts the matrix into its compressed form. In `RowMajor`, this means using Compressed Sparse Row (CSR) format; in `ColumnMajor`, Compressed Sparse Column (CSC) format. Free the uncompressed data storage.
- **uncompress()**: Converts the matrix back to its uncompressed state, freeing the data mambers used for compressed memory.
- **is_compressed()**: Checks if the matrix is currently in a compressed state.
- **resize(size_t nrows, size_t ncols)**: Resizes the matrix to have `nrows` rows and `ncols` columns, preserving existing data if within new bounds, throws error if resizing would cause data loss.
- **erase(size_t i, size_t j)**: Removes the element at row `i` and column `j` if the matrix is uncompressed. Throws an error if the matrix is compressed or if indices are out of bounds.
- **contains(size_t i, size_t j) const**: Checks if an element exists at row `i` and column `j`.
- **clear()**: Clears all data from the matrix, resetting it to an empty state.
- **print() const**: Prints the entire matrix to the console in a readable format.
- **printBones() const**: Provides a debug view of the matrix, showing the underlying active data structure (uncompressed or compressed).

Additionally, there are two friends:

- **readMTX(Matrix<T, Order>& matrix, const std::string& filename)**: Reads a matrix from a Matrix Market file (.mtx), populating it according to the file contents.

- **operator\*(const Matrix<T, Order>& A, const std::vector<T>& v)**: Multiplies the matrix `A` with the vector `v` and returns the resulting vector. This function checks for dimension compatibility and uses the appropriate algorithm based on whether `A` is compressed.

- **operator\*(const Matrix<T, Order>& A, const Matrix<T, Order>& v)**: overload for `operator*` to accept a single-column Matrix as a vector. Converts the Matrix `v` to a `std::vector<T>` so it can call the main method. It uses a doubled StorageOrder template to allow any Order combination of `A` and `v` inputs



## Testing with main.cpp
`main.cpp` serves as the test driver for the Matrix class.

First it initializes a custom matrix with some values, to test various functionalities implemented:

- Basic matrix creation and operations in both `RowMajor` and `ColumnMajor` storage orders.
- Compression and decompression of matrices.
- Exception handling in edge cases like modifications in a compressed state.
- Matrix resizing and element manipulation.
- Performance measurements of matrix-vector multiplication under different conditions.

Then, it imports a MTX file to test the multiplication. A template function `timeMultiplication` is defined to measure the execution time of matrix-vector multiplications. It uses `std::chrono` for high-resolution timing, providing a clear output of the time taken for specified multiplication tasks in all combinations of storage order and storage state.


## Build and Run
A `Makefile` is provided to compile and run the project. It includes flags for optimization to ensure accurate performance testing. To build and test the project, use the following commands:
```bash
make
./main
