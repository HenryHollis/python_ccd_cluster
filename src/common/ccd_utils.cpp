//
// Created by Henry Hollis on 1/28/24.
//

#include "ccd_utils.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
double ccd_utils::calcCCDsimple(const std::vector<double> &ref, int num_ref_rows,
                                const std::vector<double> &emat, size_t emat_row, size_t emat_col,
                                bool scale) {

    std::vector<double> cormat = calcCorMat(emat, emat_row, emat_col); //calculate cormat of expression matrix
    if (num_ref_rows != emat_row || cormat.empty() || ref.empty() ) {
        throw std::invalid_argument("Matrices must be of the same size for calcCCDsimple");
    }
    double upperTriDiff = 0.0;
    //loop through triangular matrix and accumulate the difference between entries of ref and cormat
    for (size_t i = 0; i < num_ref_rows; ++i) {
        for (size_t j = i; j < num_ref_rows; ++j) {
            upperTriDiff += pow(ref[i * num_ref_rows + j] - cormat[i * num_ref_rows + j], 2);
        }
    }
    double ccd = sqrt(upperTriDiff);

    if (!ref.empty() ) {
        // Number of columns is the size of any row (assuming all rows have the same size)
        if (scale) {
            size_t nPairs = choose(num_ref_rows, 2);
            ccd /= static_cast<double>(nPairs);
        }
    } else
        std::cerr << "Matrix ref is empty or has empty rows." << std::endl;

    return ccd;
}

std::vector<double> ccd_utils::calcCorMat(const std::vector<double> &rect, size_t numRows, size_t numCols) {
    /* Takes rectangular matrix and calculated the gene x gene correlation matrix
     *
     */

    // Initialize the flat correlation matrix with zeros
    std::vector<double> correlationMatrix(numRows*numRows ,0.0);

    // Calculate pairwise correlations for each row
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = i; j < numRows; ++j) {
            if (i == j) {
                // Diagonal elements (correlation with itself) are always 1
                correlationMatrix[i * numRows + j] = 1.0;
            } else {
                // Pass a single row to the function
                std::vector<double> rowToProcessi(rect.begin() + i * numCols, rect.begin() + (i + 1) * numCols);
                std::vector<double> rowToProcessj(rect.begin() + j * numCols, rect.begin() + (j + 1) * numCols);
//                correlationMatrix[i * numRows + j] = correlationMatrix[i * numRows + j] = cor(rowToProcessi, rowToProcessj);

                std::vector<double> rank_x = rankVector(rowToProcessi);
                std::vector<double> rank_y = rankVector(rowToProcessj);
                correlationMatrix[i * numRows + j] = correlationMatrix[i * numRows + j] = cor(rank_x, rank_y);
            }
        }
    }

    return correlationMatrix;
}


long ccd_utils::choose(size_t n, int k) {
    if (0 == k)
        return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

/**
 * Spearman Correlation Code
 * Inspired from https://www.geeksforgeeks.org/program-spearmans-rank-correlation/
 */
// Function returns the rank vector
// of the set of observations
std::vector<double> ccd_utils::rankVector(const std::vector<double> &X) {

    int N = X.size();

    // Rank Vector
    std::vector<double> Rank_X(N);

    for(int i = 0; i < N; i++)
    {
        int r = 1, s = 1;

        // Count no of smaller elements
        // in 0 to i-1
        for(int j = 0; j < i; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Count no of smaller elements
        // in i+1 to N-1
        for (int j = i+1; j < N; j++) {
            if (X[j] < X[i] ) r++;
            if (X[j] == X[i] ) s++;
        }

        // Use Fractional Rank formula
        // fractional_rank = r + (n-1)/2
        Rank_X[i] = r + (s-1) * 0.5;
    }

    // Return Rank Vector
    return Rank_X;
}

// function that returns
// Pearson correlation coefficient.
double ccd_utils::cor(const std::vector<double> &X, const std::vector<double> &Y) {
    int n = X.size();
    double sum_X = 0, sum_Y = 0,
            sum_XY = 0;
    double squareSum_X = 0,
            squareSum_Y = 0;

    for (int i = 0; i < n; i++)
    {
        // sum of elements of array X.
        sum_X = sum_X + X[i];

        // sum of elements of array Y.
        sum_Y = sum_Y + Y[i];

        // sum of X[i] * Y[i].
        sum_XY = sum_XY + X[i] * Y[i];

        // sum of square of array elements.
        squareSum_X = squareSum_X +
                      X[i] * X[i];
        squareSum_Y = squareSum_Y +
                      Y[i] * Y[i];
    }

    // use formula for calculating
    // correlation coefficient.
    double corr = (double)(n * sum_XY -
                           sum_X * sum_Y) /
                  sqrt((n * squareSum_X -
                        sum_X * sum_X) *
                       (n * squareSum_Y -
                        sum_Y * sum_Y));

    return corr;
}




std::vector<double> ccd_utils::sliceColumns(const std::vector<double> &matrix, const std::vector<size_t> &columnsToAccess,
                                            size_t nrow, size_t ncol)  {
    if (matrix.empty() || columnsToAccess.empty()) {
        std::cerr << "Matrix is empty or columns to slice are empty." << std::endl;
        return {}; // Return an empty matrix if either the matrix or columns are empty
    }

    // Create a vector to store the elements of the specified columns
    std::vector<double> columnsToProcess;

    // Copy the elements of the specified columns into the processing vector
    for (int row = 0; row < nrow; ++row) {
        for (int col : columnsToAccess) {
            if (col >= ncol) {
                std::cerr << "Column index out of range" <<std::endl;
                throw std::out_of_range("Column index out of range");
            }
            columnsToProcess.push_back(matrix[row * ncol + col]);
        }
    }


    return columnsToProcess;
}