#include "Matrix.h"

int main() {
    // basic tests
    algebra::Matrix<double, algebra::StorageOrder::RowMajor> mat(4, 4);
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.14;
    mat(3, 2) = 3.;
    mat(3, 3) = 4.0;
    mat.print();
    mat.printBones();

    // compress/uncompress tests
    std::cout << "\n\n ** COMPRESSION TESTS **" << std::endl;
    std::cout << "Compressing matrix...";
    mat.compress();
    std::cout << "done.\n";

    try {
        std::cout << "Trying to modify matrix value at (0,0) while compressed..." << std::endl;
        mat(0, 0) = 10; // error since matrix is compressed
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    mat.print();
    mat.printBones();

    std::cout << "\nUncompressing matrix...";
    mat.uncompress();
    std::cout << "done.\n";
    std::cout << "Trying to modify matrix value at (0,0) while uncompressed..." << std::endl;
    mat(0, 0) = 10;
    std::cout << "Matrix value at (0,0): " << mat(0, 0) << std::endl;
    mat.print();
    mat.printBones();

    std::cout << "\n\n ** RESIZE TEST **" << std::endl;
    try {
        std::cout << "Resizing matrix to 2x3" << std::endl;
        mat.resize(2, 3);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    mat.printBones();


    // A*v test
    std::cout << "\n\n ** MATRIX-VECTOR MULTIPLICATION TEST **" << std::endl;
    std::vector<double> vec = {1, 2, 0, 50};
    mat.print();
    std::cout << "Vector: ";
    for (double val : vec) std::cout << val << " ";
    std::cout << std::endl;
    // Perform multiplication in uncompressed state
    try {
        auto result_uncompressed = mat * vec;
        std::cout << "Uncompressed matrix-vector multiplication result:\n";
        for (double val : result_uncompressed) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    // Perform multiplication in compressed state
    mat.compress();
    try {
        auto result_compressed = mat * vec;
        std::cout << "Compressed matrix-vector multiplication result:\n";
        for (double val : result_compressed) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }


    return 0;
}