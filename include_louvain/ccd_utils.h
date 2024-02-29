//
// Created by Henry Hollis on 1/28/24.
//

#ifndef LOUVAIN_CCD_CCD_UTILS_H
#define LOUVAIN_CCD_CCD_UTILS_H
#include <vector>
class ccd_utils {
public:
    static std::vector<double> sliceColumns(const std::vector<double> &matrix, const std::vector<size_t> &columnsToAccess,
                                            size_t nrow, size_t ncol);

    static double calcCCDsimple(const std::vector<double> &ref, int num_ref_rows,
                                const std::vector<double> &emat, size_t emat_row, size_t emat_col,
                                bool scale);

    static std::vector<double> calcCorMat(const std::vector<double> &rect, size_t numRows, size_t numCols);
    static long choose(size_t n, int k);

    static std::vector<double> rankVector(const std::vector<double> &X);

    static double cor(const std::vector<double> &X, const std::vector<double> &Y);

};
#endif //LOUVAIN_CCD_CCD_UTILS_H
