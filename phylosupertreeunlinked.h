#ifndef PHYLOSUPERTREEUNLINKED_H
#define PHYLOSUPERTREEUNLINKED_H 

#include "phylosupertree.h"
#include "genetree.h"
#include "strictconsensusmerge.h"

class PhyloSuperTreeUnlinked : public PhyloSuperTree {
public:
    PhyloSuperTreeUnlinked(Params &params);

    PhyloSuperTreeUnlinked(Params &params, const StrVector &sourceTrees);

    ~PhyloSuperTreeUnlinked();

    bool isSuperTreeUnlinked() override {
        return true;
    }

    StrVector getAllSeqNames();

    void runGeneTreesReconstruction();

    void printGeneTrees();

    void dfsMRP(Node* u, Node* pa, int &time, vector<pair<int, int>> &eulerInternalBranch, map<string, int> &leafIndex);

    void buildMRPMatrix();

    void doMRP();

    void printResultWithMRPTree();

    void doSCM();

    void printResultWithSCMTree();

    void printScoreWithConAln(GeneTree *tree, string treeType);

    StrVector createBootstrapGeneTrees();

    Alignment *conAln = NULL;

    Alignment *mrpAln = NULL;

    GeneTree *mrpTree = NULL;

    GeneTree *scmTree = NULL;

    int index;

    string firstSCMTree;

    StrVector allSeqNames;

    map<string, int> seqNameToIndex;

    vector<int> gene_tree_assigned;
};

#endif
