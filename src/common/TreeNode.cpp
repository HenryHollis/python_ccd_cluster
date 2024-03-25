//
// Created by Henry Hollis on 3/20/24.
//
#include <TreeNode.h>
void TreeNode::addChild(TreeNode* child) {
    children.push_back(child);
    child->parent = this;
}
void TreeNode::removeChild(TreeNode* child) {
    children.erase(std::remove(children.begin(), children.end(), child), children.end());
    child->parent = nullptr;
}

TreeNode* TreeNode::findChildById(size_t id) {
    for (auto* child : children) {
        if (child->id == id) {
            return child;
        }
    }
    return nullptr;
}

vector<TreeNode*> TreeNode::getChildren() {
        return children;
}

vector<TreeNode*> TreeNode::getLeaves() {
    vector<TreeNode*> communities;
    getLeavesHelper(this, communities);
    return communities;
}

void TreeNode::getLeavesHelper(const TreeNode* node, vector<TreeNode*>& communities)const {
    if (node->children.empty()) {
        communities.push_back(const_cast<TreeNode*>(node));
    } else {
        for (const auto& child : node->children) {
            getLeavesHelper(child, communities);
        }
    }
}

TreeNode* searchTreeVec(const vector<TreeNode*>& communities, size_t id) {
    for (TreeNode* leaf : communities) {
        if (leaf->id == id) {
            return leaf;
        }
    }
    return nullptr;
}

vector<TreeNode*> move_node_tree(vector<TreeNode*>& communities, size_t from_node_id, size_t  to_node_id, size_t childID) {
    TreeNode* from_node = searchTreeVec(communities, from_node_id);
    TreeNode* to_node = searchTreeVec(communities, to_node_id);
    if (!to_node){  //In case a new community is needed, create it
        to_node = new TreeNode(to_node_id);
        communities.push_back(to_node);
    }
    if (from_node){
        TreeNode* childToMove = from_node->findChildById(childID);
        if (childToMove) {
            // Remove the child branch from its current parent
            childToMove->parent->removeChild(childToMove);

            // Add the child branch as a child of toNode
            to_node->addChild(childToMove);

            // Check if from_node is now childless and remove it from communities if needed
            if (from_node->children.size() == 0) {
                auto it = std::find(communities.begin(), communities.end(), from_node);
                if (it != communities.end()) {
                    communities.erase(it);
                }
            }

        }else
            cerr<<"Child with ID: "<< childID << " not found under node with ID: "<< from_node_id << "."<<endl;

    }else
        cerr<<"from_node nodeID is not found in vector provided."<<endl;
    return communities;

}

vector<size_t> get_ids_from_tree(vector<TreeNode*> &communities){
    vector<size_t> res;
    res.reserve(communities.size());
    for(TreeNode* comm : communities){
        res.push_back(comm->id);
    }
    return res;
}

void printTree(const vector<TreeNode*>& communities, int depth) {
    for (const auto& node : communities) {
        for (int i = 0; i < depth; ++i) {
            cout << "  ";
        }
        if (node->parent) {
            cout << "|-- " << node->id << " (Parent: " << node->parent->id << ")" << endl;
        } else {
            cout << "|-- " << node->id << " (Root)" << endl;
        }
        printTree(node->getChildren(), depth + 1);
    }
}

vector<TreeNode*> mergeNodes(vector<TreeNode*>& communities, size_t id1, size_t  id2, size_t parentID) {
    TreeNode* node1 = searchTreeVec(communities, id1);
    TreeNode* node2 = searchTreeVec(communities, id2);
    if (node1 && node2) {
        TreeNode *parent = new TreeNode(parentID);  // Placeholder ID for internal nodes
        parent->addChild(node1);
        if (node1 != node2)
            parent->addChild(node2);

        // Remove node1 and node2 from communities vector
        communities.erase(std::remove(communities.begin(), communities.end(), node1), communities.end());
        if (node1 != node2)
            communities.erase(std::remove(communities.begin(), communities.end(), node2), communities.end());

        // Add the new parent node to the communities vector
        communities.push_back(parent);

    }else
        cerr<<"One or both nodeIDs provided were not found in vector provided."<<endl;

    return communities;
}