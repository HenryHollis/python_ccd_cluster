#include <iostream>
#include <vector>

using namespace std;
class TreeNode {
    public:
        size_t id;
        TreeNode *parent;
        vector<TreeNode*> children;
        int numLeaves; // Number of leaves under this node
        int group;  // Group identifier

        TreeNode(size_t id) : id(id), numLeaves(0), parent(nullptr) {}

        void addChild(TreeNode *child);

        void removeChild(TreeNode *child);

        TreeNode *findChildById(size_t id);

        vector<TreeNode *> getChildren();

        vector<TreeNode *> getLeaves();
    private:
        void getLeavesHelper(const TreeNode *node, vector<TreeNode *> &communities) const;
        void updateNumLeaves();
        int countLeaves(TreeNode* node) ;
};

TreeNode *searchTreeVec(const vector<TreeNode *> &communities, size_t id);

vector<TreeNode *> move_node_tree(vector<TreeNode *> &communities, size_t from_node_id, size_t to_node_id, size_t childID);

void printTree(const vector<TreeNode *> &communities, int depth = 0);

vector<TreeNode *> mergeNodes(vector<TreeNode *> &communities, size_t id1, size_t id2, size_t parentID);

vector<size_t> get_ids_from_tree(vector<TreeNode*> &communities);
vector<int> get_group_from_tree(vector<TreeNode*> &communities);


bool checkCommNodeCount(const vector<TreeNode*>& leaves, size_t cutoff);