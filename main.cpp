#include "Matrix.h"

int main() {
    std::cout << "  ** PACS Challenge 2 - A Sparse Matrix **  " << std::endl << std::endl;
    // basic tests
    algebra::Matrix<double, algebra::StorageOrder::ColumnMajor> mat(4, 4);
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.14;
    mat(3, 2) = 3.;
    mat(3, 3) = 4.0;
    mat.print();
    mat.printBones();

    // compress/uncompress tests
    std::cout << "\n\n ** COMPRESSION TESTS **" << std::endl;
    std::cout << " - Compressing matrix...";
    mat.compress();
    std::cout << " done.\n";
    mat.print();
    mat.printBones();

    try {
        // edit existing values
        std::cout << "\n - Trying to modify matrix value at (0,0) while compressed..." << std::endl;
        mat(0, 0) = 10; // ok
        std::cout << "Matrix value at (0,0): " << mat(0, 0) << std::endl;
        mat.print();
        // add new values
        std::cout << " - Trying to add NEW matrix value at (1,1) while compressed..." << std::endl;
        mat(1, 1) = 10; // error
        std::cout << "Matrix value at (1,1): " << mat(0, 0) << std::endl;
        mat.print();
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    std::cout << "\n - Uncompressing matrix...";
    mat.uncompress();
    std::cout << " done.\n";
    std::cout << " - Trying to modify matrix value at (0,0) while uncompressed..." << std::endl;
    mat(0, 0) = 11;
    std::cout << "Matrix value at (0,0): " << mat(0, 0) << std::endl;
    mat.print();




    std::cout << "\n\n ** RESIZE TEST **" << std::endl;
    try {
        std::cout << "Resizing matrix to 2x3" << std::endl;
        mat.resize(2, 3);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    // remove oob first
    // mat.compress();
    try {
        std::cout << "Removing oob elements...";
        for(size_t i = 0; i < mat.num_rows(); i++)
            for(size_t j = 0; j < mat.num_cols(); j++)
                if (i >= 2 || j >= 3)
                    mat.erase(i, j);
        std::cout << " done.\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    mat.print();
    try {
        std::cout << "Resizing matrix to 2x3" << std::endl;
        mat.resize(2, 3);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    mat.printBones();


    // A*v test
    std::cout << "\n\n ** MATRIX-VECTOR MULTIPLICATION TEST **" << std::endl;
    std::vector<double> vec = {1, 2, 0, 5};
    try {
        mat.resize(4, 4);
        mat(3,3) = .1;
        mat.print();
        mat.printBones();
        std::cout << "\nVector: ";
        for (double val : vec) std::cout << val << " ";
        std::cout << std::endl;
        
        // Perform multiplication in uncompressed state
        auto result_uncompressed = mat * vec;
        std::cout << "Uncompressed matrix-vector multiplication result:\n";
        for (double val : result_uncompressed) {
            std::cout << val << " ";
        }
        std::cout << "\n";
        // Perform multiplication in compressed state
        mat.compress();
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