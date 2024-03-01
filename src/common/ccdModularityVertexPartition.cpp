#include "ccdModularityVertexPartition.h"
#include "ccd_utils.h"

#ifdef DEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif
ccdModularityVertexPartition::ccdModularityVertexPartition(Graph *graph, const vector<size_t> &membership,
                                                           const vector<double> &geneSampleMatrix):
        MutableVertexPartition(graph, membership),
        geneSampleMatrix(geneSampleMatrix){

}

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph,
                                                           vector<size_t> const& membership) :
        MutableVertexPartition(graph,
                               membership)
{ }

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph) :
        MutableVertexPartition(graph)
{ }

ccdModularityVertexPartition::~ccdModularityVertexPartition()
{ }

ccdModularityVertexPartition* ccdModularityVertexPartition::create(Graph* graph)
{
    std::vector<double> GeneMatrix = this->geneSampleMatrix;
    size_t geneMatRows = this-> geneMatRows;
    size_t geneMatCols = this-> geneMatCols;

    std::vector<double> RefMat = this->refMatrix;
    size_t refMatRows = this-> refMatRows;
    size_t refMatCols = this-> refMatCols;

    auto* tmp = new ccdModularityVertexPartition(graph);
    tmp->geneSampleMatrix = GeneMatrix;
    tmp->refMatrix = RefMat;
    tmp->refMatRows = refMatRows;
    tmp->refMatCols = refMatCols;
    tmp->geneMatRows = geneMatRows;
    tmp->geneMatCols = geneMatCols;
    return tmp;

}

ccdModularityVertexPartition* ccdModularityVertexPartition::create(Graph* graph, vector<size_t> const& membership)
{
    std::vector<double> GeneMatrix = this->geneSampleMatrix;
    size_t geneMatRows = this-> geneMatRows;
    size_t geneMatCols = this-> geneMatCols;

    std::vector<double> refMat = this->refMatrix;
    size_t refMatRows = this-> refMatRows;
    size_t refMatCols = this-> refMatCols;

    auto* tmp = new  ccdModularityVertexPartition(graph, membership);
    tmp->geneSampleMatrix = GeneMatrix;
    tmp->refMatrix = refMat;
    tmp->refMatRows = refMatRows;
    tmp->refMatCols = refMatCols;
    tmp->geneMatRows = geneMatRows;
    tmp->geneMatCols = geneMatCols;

    return tmp;
}

ccdModularityVertexPartition *ccdModularityVertexPartition::create(Graph *graph, const vector<size_t> &membership,
                                                                   const vector<double> &geneSampleMatrix) {
    return new ccdModularityVertexPartition(graph, membership, geneSampleMatrix);
}

void ccdModularityVertexPartition::setGeneSampleMatrix(const vector<double> &geneSampleMatrix, size_t rows, size_t cols) {
    if( !geneSampleMatrix.empty()){
        this->geneSampleMatrix = geneSampleMatrix;
        this->geneMatRows = rows;
        this->geneMatCols = cols;
    } else
        throw std::invalid_argument("Gene expression must not be empty");
}

void ccdModularityVertexPartition::setRefMatrix(const vector<double> &refMat, size_t rows, size_t cols) {
    if( !refMat.empty()){
        this->refMatrix = refMat;
        this->refMatRows = rows;
        this->refMatCols = cols;
    } else
        throw std::invalid_argument("Reference Matrix must not be empty");
}

const std::vector<double> &ccdModularityVertexPartition::getGeneMatrix() {
    return geneSampleMatrix;
}
const std::vector<double> &ccdModularityVertexPartition::getRefMatrix() {
    return refMatrix;
}

/*****************************************************************************
  Returns the difference in modularity if we move a node to a new community
*****************************************************************************/
double ccdModularityVertexPartition::diff_move(size_t v, size_t new_comm)
{
#ifdef DEBUG
    cerr << "double ccdModularityVertexPartition::diff_move(" << v << ", " << new_comm << ")" << endl;
#endif
    size_t old_comm = this->_membership[v]; //what community is v in?
    double diff = 0.0;
    double ccd_diff = NAN;
    double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());

    std::vector<size_t> Nodes_in_old_comm = this->get_community(old_comm);
    std::vector<size_t> Nodes_in_new_comm = this->get_community(new_comm);
    Nodes_in_new_comm.push_back(v); // the nodes in the new community union v
    double old_ccd = INFINITY;
    double new_ccd = INFINITY;

    if (total_weight == 0.0)
        return 0.0;
    if (new_comm != old_comm)
    {
        // ********CALC CCD*************
        std::vector<double> emat = this->getGeneMatrix(); //Get the expression matrix associated with the partition object
        std::vector<double> refmat = this->getRefMatrix();
        // calculate ccd in old community if enough nodes are aggregated into c's community:
        if (CCD_COMM_SIZE < Nodes_in_old_comm.size()) {
            std::vector<double> comm_emat = ccd_utils::sliceColumns(emat, Nodes_in_old_comm, this->geneMatRows,Nodes_in_old_comm.size() );
            old_ccd = ccd_utils::calcCCDsimple(refmat, this->refMatRows, comm_emat, this->geneMatRows,Nodes_in_old_comm.size(), false);
        }
        //calc ccd of adding v into new community
        if (CCD_COMM_SIZE < Nodes_in_new_comm.size()) {
            std::vector<double> comm_emat = ccd_utils::sliceColumns(emat,  Nodes_in_new_comm, this->geneMatRows, Nodes_in_new_comm.size());
            new_ccd = ccd_utils::calcCCDsimple(refmat, this->refMatRows, comm_emat, this->geneMatRows,Nodes_in_new_comm.size(), false);

        }
        //****************************
#ifdef DEBUG
        cerr << "\t" << "old_comm: " << old_comm << endl;
#endif
        double w_to_old = this->weight_to_comm(v, old_comm); //sum of edges going to old community
#ifdef DEBUG
        cerr << "\t" << "w_to_old: " << w_to_old << endl;
#endif
        double w_from_old = this->weight_from_comm(v, old_comm); //sum of edges coming from old community
#ifdef DEBUG
        cerr << "\t" << "w_from_old: " << w_from_old << endl;
#endif
        double w_to_new = this->weight_to_comm(v, new_comm); //sum of edges going to new community
#ifdef DEBUG
        cerr << "\t" << "w_to_new: " << w_to_new << endl;
#endif
        double w_from_new = this->weight_from_comm(v, new_comm); //sum of edges coming from new community
#ifdef DEBUG
        cerr << "\t" << "w_from_new: " << w_from_new << endl;
#endif
        double k_out = this->graph->strength(v, IGRAPH_OUT); //sum of all edges leaving node v
#ifdef DEBUG
        cerr << "\t" << "k_out: " << k_out << endl;
#endif
        double k_in = this->graph->strength(v, IGRAPH_IN); //sum of all edges coming into v
#ifdef DEBUG
        cerr << "\t" << "k_in: " << k_in << endl;
#endif
        double self_weight = this->graph->node_self_weight(v);
#ifdef DEBUG
        cerr << "\t" << "self_weight: " << self_weight << endl;
#endif
        double K_out_old = this->total_weight_from_comm(old_comm); //total weights of edges leaving OLD community
#ifdef DEBUG
        cerr << "\t" << "K_out_old: " << K_out_old << endl;
#endif
        double K_in_old = this->total_weight_to_comm(old_comm);  //total weights of edges Entering OLD community
#ifdef DEBUG
        cerr << "\t" << "K_in_old: " << K_in_old << endl;
#endif
        double K_out_new = this->total_weight_from_comm(new_comm) + k_out;
#ifdef DEBUG
        cerr << "\t" << "K_out_new: " << K_out_new << endl;
#endif
        double K_in_new = this->total_weight_to_comm(new_comm) + k_in;
#ifdef DEBUG
        cerr << "\t" << "K_in_new: " << K_in_new << endl;
      cerr << "\t" << "total_weight: " << total_weight << endl;
#endif
        double diff_old = (w_to_old - k_out*K_in_old/total_weight) + \
               (w_from_old - k_in*K_out_old/total_weight);
#ifdef DEBUG
        cerr << "\t" << "diff_old: " << diff_old << endl;
#endif
        double diff_new = (w_to_new + self_weight - k_out*K_in_new/total_weight) + \
               (w_from_new + self_weight - k_in*K_out_new/total_weight);
#ifdef DEBUG
        cerr << "\t" << "diff_new: " << diff_new << endl;
#endif
        diff = diff_new - diff_old;
#ifdef DEBUG
        cerr << "\t" << "diff: " << diff << endl;
#endif
        ccd_diff = old_ccd - new_ccd; //negative number returns smaller score
    }
#ifdef DEBUG
    cerr << "exit double ccdModularityVertexPartition::diff_move((" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
#endif
    double m;
    if (this->graph->is_directed())
        m = this->graph->total_weight();
    else
        m = 2*this->graph->total_weight();
    ccd_diff = isfinite(ccd_diff) ? ccd_diff : 0.0;
    return diff/m  + 0.1 * ccd_diff;
}


/*****************************************************************************
  Give the modularity of the partition.

  We here use the unscaled version of modularity, in other words, we don"t
  normalise by the number of edges.
******************************************************************************/
double ccdModularityVertexPartition::quality()
{
#ifdef DEBUG
    cerr << "double ccdModularityVertexPartition::quality()" << endl;
#endif
    double mod = 0.0;

    double m;
    if (this->graph->is_directed())
        m = this->graph->total_weight();
    else
        m = 2*this->graph->total_weight();

    if (m == 0)
        return 0.0;

    for (size_t c = 0; c < this->n_communities(); c++)
    {
        double w = this->total_weight_in_comm(c);
        double w_out = this->total_weight_from_comm(c);
        double w_in = this->total_weight_to_comm(c);
#ifdef DEBUG
        double csize = this->csize(c);
      cerr << "\t" << "Comm: " << c << ", size=" << csize << ", w=" << w << ", w_out=" << w_out << ", w_in=" << w_in << "." << endl;
#endif
        mod += w - w_out*w_in/((this->graph->is_directed() ? 1.0 : 4.0)*this->graph->total_weight());
    }
    double q = (2.0 - this->graph->is_directed())*mod;
#ifdef DEBUG
    cerr << "exit double ccdModularityVertexPartition::quality()" << endl;
    cerr << "return " << q/m << endl << endl;
#endif
    return q/m;
}
