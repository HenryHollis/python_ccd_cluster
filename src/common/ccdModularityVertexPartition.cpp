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
        MutableVertexPartition(graph, membership)
{}

ccdModularityVertexPartition::ccdModularityVertexPartition(Graph* graph) :
        MutableVertexPartition(graph)
{
    // Create initial singleton communities.
    for (size_t nodeID : _membership){
        TreeNode* comm = new TreeNode(nodeID);
        TreeNode* vert = new TreeNode(nodeID);
        comm->addChild(vert);
        this->tree.push_back(comm);

    }
}

ccdModularityVertexPartition::~ccdModularityVertexPartition()
{
}

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
    vector<TreeNode*> newTree; //the new tree we build for the collapse partition
    for(int i = 0; i < this->tree.size(); i++){  //for each element of old tree vector
        size_t comm_id = this->tree[i]->id;      //grab the ID associated with element i
        TreeNode* comm =  new TreeNode(comm_id); //create a new TreeNode with same ID
        comm->addChild(this->tree[i]);        //Add children of this->tree to new community
        newTree.push_back(comm);                   //Add the whole thing to my new tree vector
    }
    tmp->tree = newTree;
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

    //TODO: Will need to add code for tree member when using Leiden clustering.
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

void ccdModularityVertexPartition::setSubjectGroup(const std::vector<int> &subject_group) {
    //Marches through tree object and sets group of leaves according to subject_group
    if(!subject_group.empty()){
        for (int i = 0; i < subject_group.size(); i++)
            this->tree[i]->children[0]->group = subject_group[i];
    }else
        throw std::invalid_argument("Subject Grouping vector must not be empty");
}

const std::vector<double> &ccdModularityVertexPartition::getGeneMatrix() {
    return geneSampleMatrix;
}
const std::vector<double> &ccdModularityVertexPartition::getRefMatrix() {
    return refMatrix;
}

size_t ccdModularityVertexPartition::vecHash::operator()(const std::vector<size_t>& v) const {
    return std::accumulate(v.begin(), v.end(), 0, [](size_t a, size_t b) { return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2)); });
}

void ccdModularityVertexPartition::move_node(size_t v, size_t new_comm) {
    size_t old_comm = _membership[v];
    this->tree = move_node_tree(this->tree,old_comm, new_comm, v);
    MutableVertexPartition::move_node(v, new_comm);

}

void ccdModularityVertexPartition::relabel_communities(const vector<size_t> &new_comm_id) {
    MutableVertexPartition::relabel_communities(new_comm_id);
    //Relabel elements of this->tree:
    for (size_t i = 0; i < this->tree.size(); i++)
        this->tree[i]->id = new_comm_id[this->tree[i]->id];
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
    double ccd_diff;
    double total_weight = this->graph->total_weight()*(2.0 - this->graph->is_directed());

    vector<TreeNode*> verts = searchTreeVec(this->tree, old_comm)->getChildren(); //all verts in old community
    vector<TreeNode*>vert_leaves = searchTreeVec(verts, v)->getLeaves();  //get nodes under vertex v
    vector<size_t> nodes_in_v = get_ids_from_tree(vert_leaves);
    vector<int> Groups_in_v = get_group_from_tree(vert_leaves);

    vector<TreeNode*> TreeNodes_in_old_comm_v = searchTreeVec(this->tree, old_comm)->getLeaves();
    vector<size_t> Nodes_in_old_comm_v = get_ids_from_tree(TreeNodes_in_old_comm_v);
    vector<int> Groups_in_old_comm_v = get_group_from_tree(TreeNodes_in_old_comm_v);

    TreeNode* new_comm_TreeNode = searchTreeVec(this->tree, new_comm);
    vector<size_t> Nodes_in_new_comm_no_v;
    vector<int> Groups_in_new_comm_no_v;

    if(new_comm_TreeNode){  //If the new node is not empty...
        //Get all leaves of the nodes in the new comm
        vector<TreeNode*> TreeNodes_in_new_comm_no_v = searchTreeVec(this->tree, new_comm)->getLeaves();
        Nodes_in_new_comm_no_v = get_ids_from_tree(TreeNodes_in_new_comm_no_v);
        Groups_in_new_comm_no_v = get_group_from_tree(TreeNodes_in_new_comm_no_v);

    }

    std::vector<size_t> Nodes_in_new_comm_v;
    std::vector<int> Groups_in_new_comm_v;

    Nodes_in_new_comm_v.assign(Nodes_in_new_comm_no_v.begin(), Nodes_in_new_comm_no_v.end()); //deep copy
    Nodes_in_new_comm_v.insert(Nodes_in_new_comm_v.end(), std::begin(nodes_in_v), std::end(nodes_in_v));
    Groups_in_new_comm_v.assign(Groups_in_new_comm_no_v.begin(), Groups_in_new_comm_no_v.end()); //deep copy of grp list
    Groups_in_new_comm_v.insert(Groups_in_new_comm_v.end(), std::begin(Groups_in_v), std::end(Groups_in_v));

    std::vector<size_t> Nodes_in_old_comm_no_v;
    std::vector<int> Groups_in_old_comm_no_v;
    Nodes_in_old_comm_no_v.assign(Nodes_in_old_comm_v.begin(), Nodes_in_old_comm_v.end()); //deep copy
    Groups_in_old_comm_no_v.assign(Groups_in_old_comm_v.begin(), Groups_in_old_comm_v.end());

    // Define a lambda function to check if elements in 'nodes_in_v' are in another array
    auto is_in_array_to_delete = [&](int val) {
        return std::find(std::begin(nodes_in_v), std::end(nodes_in_v), val) != std::end(nodes_in_v);
    };

    // Remove elements from both vectors based on the indices removed from Nodes_in_old_comm_no_v.
    auto it1 = Nodes_in_old_comm_no_v.begin();
    auto it2 = Groups_in_old_comm_no_v.begin();

    while (it1 != Nodes_in_old_comm_no_v.end() && it2 != Groups_in_old_comm_no_v.end()) {
        if (is_in_array_to_delete(*it1)) { //If elmnt of nodes_old_comm_no_v needs deletion, rm from grps_old_comm_no_v too
            it1 = Nodes_in_old_comm_no_v.erase(it1);
            it2 = Groups_in_old_comm_no_v.erase(it2);
        } else {
            ++it1;
            ++it2;
        }
    }

    // Use std::remove_if with the lambda function to remove elements from vec
    // Nodes_in_old_comm_no_v.erase(std::remove_if(Nodes_in_old_comm_no_v.begin(), Nodes_in_old_comm_no_v.end(), is_in_array_to_delete), Nodes_in_old_comm_no_v.end());

    //Change in ccd should be [ccd(new+v) + ccd(old - v)] - [ccd(old + v) + ccd(new - v)]
    double old_ccd_v = 0.;
    double new_ccd_no_v = 0.;  //WHATS the ccd of a random matrix?
    double old_ccd_no_v = 0.;
    double new_ccd_w_v = 0.;
    if (total_weight == 0.0)
        return 0.0;
    if (new_comm != old_comm)
    {
        // ********CALC CCD*************
        std::vector<double> emat = this->getGeneMatrix(); //Get the expression matrix associated with the partition object
        std::vector<double> refmat = this->getRefMatrix();
        // calculate ccd in old community if enough nodes are aggregated into c's community:
        if (CCD_COMM_SIZE < Nodes_in_old_comm_v.size()) {
            auto it = this->ccdCache.find(Nodes_in_old_comm_v);
            // if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                // old_ccd_v =  it->second;
            // } else{
                //calculate the result and store it
                try{
                    std::vector<double> comm_emat_old_v = ccd_utils::sliceColumns(emat, Nodes_in_old_comm_v, this->geneMatRows, this->geneMatCols);
                    std::vector<double> comm_emat_old_v_grp_sumd;
                    int num_groups_old_v = ccd_utils::sumColumnsByGroup(comm_emat_old_v, this->geneMatRows, Nodes_in_old_comm_v.size(), Groups_in_old_comm_v, comm_emat_old_v_grp_sumd);
                    if (CCD_COMM_SIZE < num_groups_old_v)
                        old_ccd_v = ccd_utils::calcCCS(refmat, this->refMatRows, comm_emat_old_v_grp_sumd, this->geneMatRows, num_groups_old_v);
                    // this->ccdCache[Nodes_in_old_comm_v] = old_ccd_v;
                }catch (const std::out_of_range& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }

            // }
        }
        if (CCD_COMM_SIZE < Nodes_in_old_comm_no_v.size()) {
            auto it = this->ccdCache.find(Nodes_in_old_comm_no_v);
            // if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                // old_ccd_no_v =  it->second;
            // } else{
                //calculate the result and store it
                try{
                    std::vector<double> comm_emat_old_no_v = ccd_utils::sliceColumns(emat, Nodes_in_old_comm_no_v, this->geneMatRows, this->geneMatCols);
                    std::vector<double> comm_emat_old_no_v_grp_sumd;
                    int num_groups_old_no_v = ccd_utils::sumColumnsByGroup(comm_emat_old_no_v, this->geneMatRows, Nodes_in_old_comm_no_v.size(), Groups_in_old_comm_no_v, comm_emat_old_no_v_grp_sumd);
                    if (CCD_COMM_SIZE < num_groups_old_no_v )
                        old_ccd_no_v = ccd_utils::calcCCS(refmat, this->refMatRows, comm_emat_old_no_v_grp_sumd, this->geneMatRows, num_groups_old_no_v);
                    // this->ccdCache[Nodes_in_old_comm_no_v] = old_ccd_no_v;
                }catch (const std::out_of_range& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }

            // }
        }
        //calc ccd of adding v into new community
        if (CCD_COMM_SIZE < Nodes_in_new_comm_v.size()) {
            auto it = this->ccdCache.find(Nodes_in_new_comm_v);
            // if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                // new_ccd_w_v = it->second;
            // }else{
                //    calculate the result and store it
                try{
                    std::vector<double> comm_emat_new_v = ccd_utils::sliceColumns(emat,  Nodes_in_new_comm_v, this->geneMatRows, this->geneMatCols);
                    std::vector<double> comm_emat_new_v_grp_sumd;
                    int num_groups_new_v = ccd_utils::sumColumnsByGroup(comm_emat_new_v, this->geneMatRows, Nodes_in_new_comm_v.size(), Groups_in_new_comm_v, comm_emat_new_v_grp_sumd);
                       cout<<"comm_emat_new_v_grp_sumd:"<<endl;
                       for (int i = 0; i < 12; ++i) {
                            for (int j = 0; j < num_groups_new_v; ++j) {
                                std::cout << comm_emat_new_v_grp_sumd[i * num_groups_new_v + j] << " ";
                            }
                            std::cout << std::endl;
                        }
                        cout<<endl;
                    if (CCD_COMM_SIZE < num_groups_new_v )
                             new_ccd_w_v = ccd_utils::calcCCS(refmat, this->refMatRows, comm_emat_new_v_grp_sumd, this->geneMatRows, num_groups_new_v);
                    
                    // this->ccdCache[Nodes_in_new_comm_v] = new_ccd_w_v;
                }catch (const std::out_of_range& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }


            // }
        }
        //calc ccd of adding v into new community
        if (CCD_COMM_SIZE < Nodes_in_new_comm_no_v.size()) {
            auto it = this->ccdCache.find(Nodes_in_new_comm_no_v);
            // if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                // new_ccd_no_v = it->second;
            // }else{
                //    calculate the result and store it
                try{
                    std::vector<double> comm_emat_new_no_v = ccd_utils::sliceColumns(emat,  Nodes_in_new_comm_no_v, this->geneMatRows, this->geneMatCols);
                    vector<double> comm_emat_new_no_v_grp_sumd;
                    int num_groups_new_no_v = ccd_utils::sumColumnsByGroup(comm_emat_new_no_v, this->geneMatRows, Nodes_in_new_comm_no_v.size(), Groups_in_new_comm_no_v, comm_emat_new_no_v_grp_sumd);
                    if (CCD_COMM_SIZE < num_groups_new_no_v )
                             new_ccd_no_v = ccd_utils::calcCCS(refmat, this->refMatRows, comm_emat_new_no_v_grp_sumd, this->geneMatRows, num_groups_new_no_v);
                    // this->ccdCache[Nodes_in_new_comm_no_v] = new_ccd_no_v;
                }catch (const std::out_of_range& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }


            // }
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

    }
    new_ccd_no_v *= (Nodes_in_new_comm_no_v.size()+1);
    // new_ccd_no_v = 1. / (1. + new_ccd_no_v);
    // new_ccd_no_v = (Nodes_in_new_comm_no_v.size()+1) / (1. + new_ccd_no_v);
    new_ccd_w_v *= (Nodes_in_new_comm_v.size()+1);
    // new_ccd_w_v = 1. / (1. + new_ccd_w_v);
    // new_ccd_w_v = (Nodes_in_new_comm_v.size()+1) / (1. + new_ccd_w_v);
    old_ccd_no_v *= (Nodes_in_old_comm_no_v.size()+1);
    // old_ccd_no_v = 1. / (1. + old_ccd_no_v);
    // old_ccd_no_v = (Nodes_in_old_comm_no_v.size()+1) / (1. + old_ccd_no_v);
    old_ccd_v *= (Nodes_in_old_comm_v.size()+1);
    // old_ccd_v = 1./ (1. + old_ccd_v);
    // old_ccd_v = (Nodes_in_old_comm_v.size()+1) / (1. + old_ccd_v);

//     //   ccd_diff = (old_ccd_v + new_ccd_no_v) - (new_ccd_w_v + old_ccd_no_v); //negative number returns smaller score
    ccd_diff = (new_ccd_w_v + old_ccd_no_v) - (old_ccd_v + new_ccd_no_v) ; //negative number returns smaller score
    cout<<"v: "<<v<<endl;
    cout<<"old comm:"<<old_comm <<" --> new comm: "<<new_comm<<endl;
    cout<<"Nodes in v: " << nodes_in_v.size();
    // for(size_t node : nodes_in_v){cout<<node<<" ";}
    cout<<"\nNodes in old comm v: " << Nodes_in_old_comm_v.size();
    // for(size_t node : Nodes_in_old_comm_no_v){cout<<node<<" ";}
    cout<<" ccd(): " << old_ccd_v; 
    cout<<"\nNodes in new comm NO v: " << Nodes_in_new_comm_no_v.size();
    // for(size_t node : Nodes_in_old_comm_v){cout<<node<<" ";}
    cout<<" ccd(): "<< new_ccd_no_v; 
    cout<<"\nNodes in old comm NO v: " << Nodes_in_old_comm_no_v.size();
    // for(size_t node : Nodes_in_new_comm_v){cout<<node<<" ";}
    cout<<" ccd(): "<< old_ccd_no_v; 
    cout<<"\nNodes in new comm v: " << Nodes_in_new_comm_v.size();
    // for(size_t node : Nodes_in_new_comm_no_v){cout<<node<<" ";}
    cout<<" ccd(): "<< new_ccd_w_v <<endl; 
    // std::cout <<"v: " << v<< "; new comm: " << new_comm <<"; old_com:" << old_comm <<"; old ccd w v:" << old_ccd_v <<"; old ccd no v:" << old_ccd_no_v  <<"; new_ccd_w_v:" <<  new_ccd_w_v << "; new_ccd_no_v:" << new_ccd_no_v << "; ccd_diff:" <<ccd_diff << endl;

#ifdef DEBUG
    cerr << "exit double ccdModularityVertexPartition::diff_move((" << v << ", " << new_comm << ")" << endl;
    cerr << "return " << diff << endl << endl;
#endif
    double m;
    if (this->graph->is_directed())
        m = this->graph->total_weight();
    else
        m = 2*this->graph->total_weight();

//        int min_comm_involved = std::min(Nodes_in_old_comm_no_v.size(),Nodes_in_new_comm_v.size());
//        int total_nodes = this->graph->vcount();
//        double frac = min_comm_involved/total_nodes;
// double result = diff/m  + frac * ccd_diff;
    double result = diff/m  + .1 * ccd_diff;

    std::cout << "ccd_diff: " << ccd_diff << " mod: " << diff/m <<" res: " << result << endl;
    if (std::isnan(new_ccd_w_v)) {
                            std::cout << "new_ccd_w_v is NaN.\n Nodes in new_comm_v:" << std::endl;
                            for(size_t node : Nodes_in_new_comm_v){cout<<node<<" ";}
                            std::cout <<std::endl;
                            std::cout << "Samples in new_comm_v:" << std::endl;
                            for(size_t node : Groups_in_new_comm_v){cout<<node<<" ";} 
                            std::cout<<std::endl;
                            std::cout << "nodes in v:" << std::endl;
                            for(size_t node : nodes_in_v){cout<<node<<" ";}
                            std::cout<<std::endl;
                            std::cout << "groups in v:" << std::endl;
                            for(size_t node : Groups_in_v){cout<<node<<" ";}
                        std::cout <<std::endl;
                    } 
    return result;
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
    double ccd_score = 0.0;
    double m;
    if (this->graph->is_directed())
        m = this->graph->total_weight();
    else
        m = 2*this->graph->total_weight();

    if (m == 0)
        return 0.0;

    for (size_t c = 0; c < this->n_communities(); c++)
    {
        double c_ccd = 0.;
        TreeNode* treeComm = this->tree[c];
        vector<TreeNode*> leaves = treeComm->getLeaves();
        vector<size_t> nodes_in_c = get_ids_from_tree(leaves);
        if (CCD_COMM_SIZE < nodes_in_c.size()) {
            auto it = this->ccdCache.find(nodes_in_c);
            if (it != this->ccdCache.end()) {
                // Result is already in the cache, return it
                c_ccd =  it->second;
            } else{
                //calculate the result and store it
                try{
                    std::vector<double> sliced_emat = ccd_utils::sliceColumns(this->getGeneMatrix(), nodes_in_c, this->geneMatRows, this->geneMatCols);
                    c_ccd = ccd_utils::calcCCS(this->getRefMatrix(), this->refMatRows, sliced_emat, this->geneMatRows, nodes_in_c.size());
                    this->ccdCache[nodes_in_c] = c_ccd;
                }catch (const std::out_of_range& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }

            }
        }
        ccd_score += nodes_in_c.size() / (1. + c_ccd);


        //ccd_score +=
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
    ccd_score = ccd_score/n_communities();
    cout<< "Global CCD Metric: " << ccd_score <<endl;
    return q/m;
}