#ifndef GENETREE_H
#define GENETREE_H

#include "iqtree.h"
#include "genenode.h"

class GeneTree: public IQTree { 
public:
    GeneTree();

    GeneTree(Alignment *aln);

    GeneTree(const string &treeString);

    ~GeneTree();

    GeneNode *newNode(int node_id = -1, const char* node_name = NULL) override;

    GeneNode *newNode(int node_id, int node_name) override;

    void setOriginalNodeIndex(const map<string, int> &seqNameToIndex);

    void addSplitOfNode(GeneNode *node);

    void removeSplitOfNode(GeneNode *node);

    Split* getLeavesSplit();

    void dfsEncodeSplits(GeneNode *node, GeneNode *parent);

    void encodeSplits();

    void resetAfterMergeSCM();

    void makeFirstChildOfRoot(GeneNode *node);

    void reroot(GeneNode *node);

    void deroot();

    void collapseEdge(GeneNode *head, GeneNode *tail);

    void collapseSubtree(GeneNode *node, bool do_collapse);
    
    void readTreeString(const string &treeString);

    void printResultTree(string fileName, bool isAppend);
    
    string getResultTreeString();

    void reInitializeTree(GeneNode *node = NULL, GeneNode* parent = NULL);

    GeneNode* getAnyOtherLeaf(GeneNode *node = NULL, GeneNode *parent = NULL);

    void setRootLeaf(char *my_root);

    void dfsFixTree(GeneNode *node, GeneNode *parent);

    void getPolytomies(vector<GeneNode*> &polytomies, GeneNode *node = NULL, GeneNode *parent = NULL);

    void getRelabelMap(map<string, string> &relabel, map<string, GeneNode*> &delabel, int &label, GeneNode *node, GeneNode *parent = NULL);

    void relabelAndCollapse(const map<string, string> &relabel);

    void relabelTree(const map<string, string> &relabel, GeneNode *node = NULL, GeneNode *parent = NULL);

    void delabelTree(map<string, GeneNode*> &delabel, GeneNode *node = NULL, GeneNode *parent = NULL);

    void deleteDuplicateNode(GeneNode *node = NULL, GeneNode *parent = NULL);

    void doubleCheckUniqueName(map<string, int> &checkUniqueLabel, GeneNode *node = NULL, GeneNode *parent = NULL);

    void setNodeIdByMapName(const map<string, int> &seqNameToIndex);

    string getBootstrapTree(int index);

    void computeGeneConcordance(vector<GeneTree*> trees, map<string,string> &meanings);

    void computeSiteConcordance(map<string,string> &meanings);

    void computeSiteConcordance(pairNode &branch, int nquartets);
   
    void extractQuadSubtrees(vector<Split*> &subtrees, BranchVector &branches, Node *node = NULL, Node *dad = NULL);

    map<int, string> bootstrapTrees;

    Params treeParams;

    int totalLeafNum = 0;

    map<Split, GeneNode*> splitToNode;

    // pseudo root for unrooted tree
    GeneNode *seedNode = NULL;
};

#endif