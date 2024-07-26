//
// Created by Henry Hollis on 1/28/24.
//

#include "ccd_utils.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

double ccd_utils::calcCCDsimple(const std::vector<double> &ref, int num_ref_rows,
                                const std::vector<double> &emat, size_t emat_row, size_t emat_col,
                                bool scale) {

    std::vector<double> cormat = calcCorMat(emat, emat_row, emat_col); //calculate cormat of expression matrix
    if (num_ref_rows != emat_row || cormat.empty() || ref.empty() ) {
        throw std::invalid_argument("Matrices must be of the same size for calcCCDsimple");
    }
    double ccdUpperTriDiff = 0.0;
    //loop through triangular matrix and accumulate the difference between entries of ref and cormat
    for (size_t i = 0; i < num_ref_rows; ++i) {
        for (size_t j = i; j < num_ref_rows; ++j) {
            ccdUpperTriDiff += pow(ref[i * num_ref_rows + j] - cormat[i * num_ref_rows + j], 2);
        }
    }
    double ccd = sqrt(ccdUpperTriDiff);
    
 
    // Number of columns is the size of any row (assuming all rows have the same size)
    if (scale) {
        size_t nPairs = choose(num_ref_rows, 2);
        ccd /= static_cast<double>(nPairs);
    }


    return ccd;
}

double ccd_utils::calcCCS(const std::vector<double> &ref, int num_ref_rows,
                                const std::vector<double> &emat, size_t emat_row, size_t emat_col) {

    std::vector<double> cormat = calcCorMat(emat, emat_row, emat_col); //calculate cormat of expression matrix
    if (num_ref_rows != emat_row || cormat.empty() || ref.empty() ) {
        throw std::invalid_argument("Matrices must be of the same size for calcCCDsimple");
    }
    double ccdUpperTriDiff = 0.0;
    double nullUpperTriDiff = 0.0;
    //loop through triangular matrix and accumulate the difference between entries of ref and cormat
    //also accumulate the difference between cormat and the identity matrix
    for (size_t i = 0; i < num_ref_rows; ++i) {
        for (size_t j = i; j < num_ref_rows; ++j) {
            ccdUpperTriDiff += pow(ref[i * num_ref_rows + j] - cormat[i * num_ref_rows + j], 2);
            if(i != j) //avoids the main diagonal 
                nullUpperTriDiff += pow(cormat[i * num_ref_rows + j], 2);
            
        }
    }
    double ccd = sqrt(ccdUpperTriDiff);
    double ncd = sqrt(nullUpperTriDiff);


    return ncd - ccd;
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

// Function to sum columns based on group membership
int ccd_utils::sumColumnsByGroup(const std::vector<double>& matrix, size_t rows, size_t cols, const std::vector<int>& membership,
                      std::vector<double>& result) {
    // Determine the number of unique groups
    std::unordered_map<int, int> groupIndex;
    int groupCount = 0;
    for (int group : membership) {
        if (groupIndex.find(group) == groupIndex.end()) {
            groupIndex[group] = groupCount++;
        }
    }

    // Initialize the result vector (flat vector)
    result.assign(rows * groupCount, 0.0);

    // Sum columns based on group membership
    for (int col = 0; col < cols; ++col) {
        int group = membership[col];
        int resultCol = groupIndex[group];
        for (int row = 0; row < rows; ++row) {
            result[row * groupCount + resultCol] += matrix[row * cols + col];
        }
    }

    return(groupCount);
}

int ccd_utils::sliceAndSumByGroup(const std::vector<double>& matrix, const std::vector<size_t>& columnsToAccess,
                       size_t nrow, size_t ncol, const std::vector<int>& membership, std::vector<double>& result) {
    if (matrix.empty() || columnsToAccess.empty()) {
        std::cerr << "Matrix or columns to slice are empty." << std::endl;
        return -1; // Return an error code if either the matrix or columns are empty
    }

    // Determine the number of unique groups for columns to access and initialize the result vector
    std::unordered_map<int, int> groupIndex;
    int groupCount = 0;
    for (size_t col : columnsToAccess) {
        if (col >= membership.size()) {
            std::cerr << "Column index out of range in membership vector" << std::endl;
            throw std::out_of_range("Column index out of range in membership vector");
        }
        int group = membership[col];
        if (groupIndex.find(group) == groupIndex.end()) {
            groupIndex[group] = groupCount++;
        }
    }

    // Initialize the result vector (flat vector)
    result.assign(nrow * groupCount, 0.0);

    // Sum columns based on group membership
    size_t colsToProcess = columnsToAccess.size();
    for (size_t row = 0; row < nrow; ++row) {
        for (size_t colIdx = 0; colIdx < colsToProcess; ++colIdx) {
            size_t col = columnsToAccess[colIdx];
            if (col >= ncol) {
                std::cerr << "Column index out of range in matrix" << std::endl;
                throw std::out_of_range("Column index out of range in matrix");
            }
            int group = membership[col];
            int resultCol = groupIndex[group];
            result[row * groupCount + resultCol] += matrix[row * ncol + col];
        }
    }

    return groupCount;
}