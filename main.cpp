#include "Matrix.h"
#include <chrono>

// time helper
template <typename T, algebra::StorageOrder Order>
void timeMultiplication(algebra::Matrix<T, Order>& mat, std::vector<T>& vec, const std::string& description) {
    auto start = std::chrono::high_resolution_clock::now();
    auto result = mat * vec;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<T, std::milli> elapsed = end - start;
    std::cout << "Time taken for " << description << " multiplication: " << elapsed.count() << " ms\n";
}

int main() {

    std::cout << "test with complex variables" << std::endl;
    algebra::Matrix<std::complex<double>, algebra::StorageOrder::RowMajor> matComplex(2,2);
    matComplex(0,0) = {1,0};
    matComplex(1,1) = {0,1};

    matComplex.print();
    matComplex.printBones();
    
    std::vector<std::complex<double>> vComplex = {{1,0},{0,1}};
    auto xComplex = matComplex*vComplex;
    for (std::complex<double> val : xComplex) std::cout << val << " "; std::cout << "\n";
   

    std::cout << "  ** PACS Challenge 2 - A Sparse Matrix **  " << std::endl << std::endl;
    // basic tests
    algebra::Matrix<double, algebra::StorageOrder::RowMajor> matR(4, 4);
    algebra::Matrix<double, algebra::StorageOrder::ColumnMajor> matC(4, 4);
    matR(0, 0) = 1.;    matC(0, 0) = 1.;
    matR(0, 1) = 2.;    matC(0, 1) = 2.;
    matR(1, 0) = 3.14;  matC(1, 0) = 3.14;
    matR(3, 2) = 3.;    matC(3, 2) = 3.;
    matR(3, 3) = 4.;    matC(3, 3) = 4.;

    matR.print();
    matR.printBones();
    matC.print();
    matC.printBones();

    // compress/uncompress tests
    std::cout << "\n\n ** COMPRESSION TESTS **" << std::endl;
    std::cout << " - Compressing matrix...";
    matC.compress(); matR.compress();
    std::cout << " done.\n";
    matC.print(); matR.print();
    matC.printBones(); matR.printBones();

    try {
        // edit existing values
        std::cout << "\n - Trying to modify matrix value at (0,0) while compressed... ";
        matC(0, 0) = 10; matR(0, 0) = 10; // ok
        std::cout << "OK, matrix value at (0,0): " << matC(0, 0) << std::endl;
        matC.print(); matR.print();
        // add new values
        std::cout << " - Trying to add NEW matrix value at (1,1) while compressed... ";
        matC(1, 1) = 10; matR(1, 1) = 10; // error
        std::cout << "OK, matrix value at (1,1): " << matC(0, 0) << std::endl;
        matC.print(); matR.print();
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        matC.print(); matR.print();
    }

    std::cout << "\n - Uncompressing matrix...";
    matC.uncompress(); matR.uncompress();
    std::cout << " done.\n";
    std::cout << " - Trying to modify matrix value at (0,0) while uncompressed..." << std::endl;
    matC(0, 0) = 11; matR(0, 0) = 11;
    std::cout << "Matrix value at (0,0): " << matC(0, 0) << std::endl;
    matC.print(); matR.print();




    std::cout << "\n\n ** RESIZE TEST **" << std::endl;
    try {
        std::cout << "Resizing matrix to 2x3" << std::endl;
        matC.resize(2, 3); matR.resize(2, 3);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    // remove oob first
    // mat.compress();
    try {
        std::cout << "Removing oob elements...";
        for(size_t i = 0; i < matC.num_rows(); i++)
            for(size_t j = 0; j < matC.num_cols(); j++)
                if (i >= 2 || j >= 3)
                    matC.erase(i, j);
        for(size_t i = 0; i < matR.num_rows(); i++)
            for(size_t j = 0; j < matR.num_cols(); j++)
                if (i >= 2 || j >= 3)
                    matR.erase(i, j);
        std::cout << " done.\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    matC.compress();
    matC.print(); matR.print();
    try {
        std::cout << "Resizing matrix to 2x3" << std::endl;
        matC.resize(2, 3);
        matC.print();
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    matC.printBones();
    matC.uncompress();



    // A*v test
    std::cout << "\n\n ** MATRIX-VECTOR MULTIPLICATION TEST **" << std::endl;
    // std::vector<double> vec = {1, 2, 0, 5};
    algebra::Matrix<double, algebra::StorageOrder::ColumnMajor> vec(4, 1);
    vec(0, 0)=1; vec(1, 0)=2; vec(2, 0)=0; vec(3, 0)=5;
    vec.print();
    try {
        matR.resize(4, 4);  matC.resize(4, 4);
        matR(3,3) = .1;     matC(3,3) = .1;
        matR.print();
        
        // Perform multiplication in uncompressed state
        auto result_uncompressedR = matR * vec;
        auto result_uncompressedC = matC * vec;
        std::cout << "Matrix-vector multiplication result (RowMajor Uncompressed):\n";
        for (double val : result_uncompressedR) std::cout << val << " "; std::cout << "\n";
        std::cout << "Matrix-vector multiplication result (ColumnMajor Uncompressed):\n";
        for (double val : result_uncompressedC) std::cout << val << " "; std::cout << "\n";
        // Perform multiplication in compressed state
        matR.compress();
        auto result_compressedR = matR * vec;
        auto result_compressedC = matC * vec;
        std::cout << "Matrix-vector multiplication result (RowMajor Compressed):\n";
        for (double val : result_compressedR) std::cout << val << " "; std::cout << "\n";
        std::cout << "Matrix-vector multiplication result (ColumnMajor Compressed):\n";
        for (double val : result_compressedC) std::cout << val << " "; std::cout << "\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }



    // norm tests
    std::cout << "\n\n ** NORM TESTS **" << std::endl;
    matC.print();
    std::cout << " - 1-norm: " << matC.norm<algebra::NormType::One>() << std::endl;
    std::cout << " - Inf norm: " << matC.norm<algebra::NormType::Infinity>() << std::endl;
    std::cout << " - Frobenius norm: " << matC.norm<algebra::NormType::Frobenius>() << std::endl;
    matR.print();
    std::cout << " - 1-norm: " << matR.norm<algebra::NormType::One>() << std::endl;
    std::cout << " - Inf norm: " << matR.norm<algebra::NormType::Infinity>() << std::endl;
    std::cout << " - Frobenius norm: " << matR.norm<algebra::NormType::Frobenius>() << std::endl;


    // read matrix from file (mtx format)(RowMajor)
    std::cout << "\n\n ** MTX FILE TESTS **" << std::endl;
    std::string filename = "lnsp_131.mtx";
    std::cout << "Reading matrix from file: " << filename << std::endl;
    algebra::Matrix<double, algebra::StorageOrder::RowMajor> mtxR;
    readMTX(mtxR, filename);
    std::cout << "Matrix read from file: size (" << mtxR.num_rows() << ", " << mtxR.num_cols() << ")" << std::endl;
    std::cout << "Multiplying matrix by a random vector..." << std::endl;
    // generate a vector of random values
    std::vector<double> vecr(mtxR.num_cols());
    std::generate(vecr.begin(), vecr.end(), []() { return (double)rand() / RAND_MAX; });
    // >> time
    timeMultiplication(mtxR, vecr, "RowMajor uncompressed"); // uncompressed
    mtxR.compress();
    timeMultiplication(mtxR, vecr, "RowMajor compressed"); // compressed

    // read matrix from file (mtx format)(ColumnMajor)
    algebra::Matrix<double, algebra::StorageOrder::ColumnMajor> mtxC;
    readMTX(mtxC, filename);
    // >> time
    timeMultiplication(mtxC, vecr, "ColumnMajor uncompressed"); // uncompressed
    mtxC.compress();
    timeMultiplication(mtxC, vecr, "ColumnMajor compressed"); // compressed


    return 0;
}
