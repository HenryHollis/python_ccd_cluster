//
// Created by Henry Hollis on 1/27/24.
//

#ifndef LOUVAIN_CCD_CCDMODULARITYVERTEXPARTITION_H
#define LOUVAIN_CCD_CCDMODULARITYVERTEXPARTITION_H
#define SAMPLES_IN_COMMUNITY 3 //number of samples that contribute to correlation matrix.
#include <unordered_map>
#include <numeric>
#include <functional>


#include "TreeNode.h"
#include "MutableVertexPartition.h"

class ccdModularityVertexPartition : public MutableVertexPartition
{
public:
    //constructors
    // Constructor with matrix initialization
    ccdModularityVertexPartition(Graph* graph,
                                 vector<size_t> const& membership,
                                 const std::vector<double>& geneSampleMatrix);
    ccdModularityVertexPartition(Graph* graph,
                                 vector<size_t> const& membership);
    ccdModularityVertexPartition(Graph* graph);

    virtual ~ccdModularityVertexPartition();
    virtual ccdModularityVertexPartition* create(Graph* graph);
    virtual ccdModularityVertexPartition* create(Graph* graph, vector<size_t> const& membership);
    virtual ccdModularityVertexPartition* create(Graph* graph, vector<size_t> const& membership, const std::vector<double>& geneSampleMatrix);

    virtual double diff_move(size_t v, size_t new_comm);
    virtual double quality();
    //  Setter method for the matrix
    void setGeneSampleMatrix(const std::vector<double>& geneSampleMatrix, size_t rows, size_t cols);
    void setRefMatrix(const std::vector<double>& refMat, size_t rows, size_t cols);
    void setSubjectGroup(const std::vector<int> &subject_group);
    void setCCSweight(float weight);
    void setCellsInComm(size_t cellsInComm);
    void setCellsPerSamp(size_t cellsPerSamp);
    void setSamplesInComm(size_t samplesInComm);

    // Getter for geneSampleMatrix
    const std::vector<double>& getGeneMatrix();
    const std::vector<double> & getRefMatrix();

    void move_node(size_t v,size_t new_comm) override;
    void relabel_communities(vector<size_t> const& new_comm_id) override;
protected:
private:
    // Matrix representing genes and samples
    std::vector<double> geneSampleMatrix;
    size_t geneMatRows;
    size_t geneMatCols;    
    //Reference matrix
    std::vector<double> refMatrix;
    size_t refMatRows;
    size_t refMatCols;
    float ccsWeight;
    size_t cellsInComm;
    size_t cellsPerSamp;
    size_t samplesInComm;
    struct vecHash {
        size_t operator()(const std::vector<size_t>& v) const;
    };
    std::unordered_map<std::vector<size_t>, double, vecHash> ccdCache;
    std::vector<TreeNode*>tree;

};

#endif //LOUVAIN_CCD_CCDMODULARITYVERTEXPARTITION_H
