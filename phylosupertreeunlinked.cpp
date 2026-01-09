#include "phylosupertreeunlinked.h"

PhyloSuperTreeUnlinked::PhyloSuperTreeUnlinked(Params &params): PhyloSuperTree(params, true) {
    this->params = &params;

    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        tree->treeParams = *(this->params);
    }

    if (params.aln_file) {
        VerboseMode saved_mode;
        saved_mode = verbose_mode;
        verbose_mode = VB_QUIET;

        conAln = new Alignment(params.aln_file, params.sequence_type, params.intype);
        conAln->checkGappySeq();

        verbose_mode = saved_mode;
    }

    // if (params.gene_trees_file) {
    //     // Print gene trees here because we need bifurcating trees
    //     printGeneTrees();
    // }
}

PhyloSuperTreeUnlinked::PhyloSuperTreeUnlinked(Params &params, const StrVector &sourceTrees): PhyloSuperTree() {
    this->params = &params;

    for (auto treeString: sourceTrees) {
        GeneTree* tree = new GeneTree(treeString);
        tree->treeParams = *(this->params);
        push_back(tree);
    }
}

PhyloSuperTreeUnlinked::~PhyloSuperTreeUnlinked() {
    if (conAln) {
        delete conAln;
    }
    if (mrpAln) {
        delete mrpAln;
    }
    if (mrpTree) {
        delete mrpTree;
    }
    // if (scmTree) {
    //     delete scmTree;
    // }
}

StrVector PhyloSuperTreeUnlinked::getAllSeqNames() {
    if (allSeqNames.size() > 0) {
        return allSeqNames;
    }

    if (conAln) {
        for (int i = 0; i < conAln->getNSeq(); ++i) {
            string seqName = conAln->getSeqName(i);
            allSeqNames.push_back(seqName);
            seqNameToIndex[seqName] = i;
        }
    } else {
        for (auto it = begin(); it != end(); it++) {
            GeneTree* tree = (GeneTree*)(*it);
            NodeVector taxa;
            tree->getTaxa(taxa);
            for (auto taxon: taxa) {
                string seqName = taxon->name;
                if (seqNameToIndex.find(seqName) == seqNameToIndex.end()) {
                    allSeqNames.push_back(seqName);
                    seqNameToIndex[seqName] = allSeqNames.size() - 1;
                }
            }
        }
    }

    return allSeqNames;
}

void PhyloSuperTreeUnlinked::runGeneTreesReconstruction() {
    for (auto it = begin(); it != end(); it++) {
        cout << "----------     Reconstructing gene tree " << (it - begin()) << "     ----------\n";
        VerboseMode saved_mode;
        saved_mode = verbose_mode;
        verbose_mode = VB_QUIET;

        GeneTree* tree = (GeneTree*)(*it);
        runOptimizeAndReconstruction(tree->treeParams, tree);
        
        verbose_mode = saved_mode;
        cout << "\n---------- Reconstruction of gene tree " << (it - begin()) << " done ----------\n\n";
    }
}

void PhyloSuperTreeUnlinked::printGeneTrees() {
    string treeFile(this->params->out_prefix);
    treeFile += ".gene_trees";
    // open treeFile and remove all
    ofstream outFile(treeFile.c_str());
    outFile.close();
    
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        tree->setRootLeaf(NULL);
        tree->printResultTree(treeFile, true);
    }
}

void PhyloSuperTreeUnlinked::dfsMRP(Node* u, Node* pa, int &time, vector<pair<int, int>> &eulerInternalBranch, map<string, int> &leafIndex) {
    int in = ++time;
    FOR_NEIGHBOR_IT(u, pa, it){
        Node* v = (*it)->node;
        dfsMRP(v, u, time, eulerInternalBranch, leafIndex);
    }
    int out = time;

    if (u->isLeaf()) {
        assert(in == out);
        leafIndex[u->name] = in;
    } else if (!pa->isLeaf()) {
        eulerInternalBranch.push_back(make_pair(in, out));
    }
}

void PhyloSuperTreeUnlinked::buildMRPMatrix() {
    StrVector seqNames = getAllSeqNames();
    StrVector sequences(seqNames.size());
    
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        int time = 0;
        vector<pair<int, int>> eulerInternalBranch;
        map<string, int> leafIndex;

        tree->setRootLeaf(NULL);
        assert(tree->root->isLeaf());

        leafIndex[tree->root->name] = 0;
        dfsMRP(tree->root->neighbors[0]->node, tree->root, time, eulerInternalBranch, leafIndex);

        for (int i = 0; i < seqNames.size(); ++i) {
            if (leafIndex.find(seqNames[i]) != leafIndex.end()) {

                int leafInd = leafIndex[seqNames[i]];
                for (auto [in, out] : eulerInternalBranch) {
                    sequences[i] += (leafInd >= in && leafInd <= out) ? '1' : '0';
                }

            } else {
                sequences[i] += string(eulerInternalBranch.size(), '?');
            }
        }
    }

    for (int i = 0; i < seqNames.size(); ++i) {
        assert(seqNameToIndex[seqNames[i]] == i);
        seqNames[i] = to_string(i);
    }

    mrpAln = new Alignment(seqNames, sequences, params->sequence_type);

    cout << "\nMRP: MRP matrix built with " << mrpAln->getNSeq() << " sequences and " << mrpAln->getNSite() << " characters\n";
    // mrpAln->printPhylip(cout);
}

void PhyloSuperTreeUnlinked::doMRP() {
    buildMRPMatrix();
    mrpTree = new GeneTree(mrpAln);
    mrpTree->treeParams = *(this->params);
    mrpTree->treeParams.gbo_replicates = 0;
    runOptimizeAndReconstruction(mrpTree->treeParams, mrpTree);
    
    switch (params->mrp_type) {
        case MRPType::MRP_GREEDY: {
            StrVector bestTrees = mrpTree->candidateTrees.getHighestScoringTrees(params->popSize);
            StringIntMap treels;
            
            for (int i = 0; i < bestTrees.size(); ++i) {
                treels[bestTrees[i]] = i;
            }
            
            IntVector weight(treels.size(), 1);

            string greedyTree = computeConsensusTreeNoFileIO(treels, weight, params->tree_max_count, 
                params->split_threshold, params->split_weight_threshold, params);

            mrpTree->readTreeString(greedyTree);
            break;
        }
        case MRPType::MRP_RANDOM: {
            mrpTree->readTreeString(mrpTree->candidateTrees.getRandCandTree());
            break;
        }
        case MRPType::MRP_BEST: {
            // Do nothing
            break;
        }
    }
    NodeVector taxa;
    mrpTree->getTaxa(taxa);
    for (auto taxon: taxa) {
        taxon->name = allSeqNames[stoi(taxon->name)];
    }
}

void PhyloSuperTreeUnlinked::printResultWithMRPTree() {
    assert(mrpTree);

    string treeFile(this->params->out_prefix);
    treeFile += ".treefile";
    mrpTree->printResultTree(treeFile, false);

    string drawFile(this->params->out_prefix);
    drawFile += ".draw";

    ofstream out;
    out.exceptions(ios::failbit | ios::badbit);
    out.open(drawFile.c_str());
    
    out << "MRP TREE\n--------------------------------------------------------\n\n";
    mrpTree->drawTree(out, WT_BR_SCALE | WT_SORT_TAXA);
    out << "\n\n";
    
    out.close();

    printScoreWithConAln(mrpTree, "MRP");
}

void PhyloSuperTreeUnlinked::doSCM() {
    StrVector sourcesTree;
    int scaffoldDensity = 0;
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        scaffoldDensity = max(scaffoldDensity, tree->getNumTaxa());
        sourcesTree.push_back(tree->getTreeString());
    }
    getAllSeqNames();

    cout << fixed << setprecision(2) << "SCM: Scaffold density: " << 1.0 * scaffoldDensity / allSeqNames.size() << "\n";

    StrictConsensusMerge scm(sourcesTree, seqNameToIndex);
    scmTree = scm.getSCMTree();
    scmTree->reInitializeTree();

    cout << fixed << setprecision(2) << "SCM: Resolution of SCM Tree: " << 1.0 * (scmTree->nodeNum - scmTree->leafNum - 1) / (scmTree->leafNum - 3) << "\n"; 

    firstSCMTree = scmTree->getTreeString();

    if (params->mrp_type == MRPType::MRP_NONE) {
        delete scmTree;
        return;
    }

    vector<GeneNode*> polytomies;
    scmTree->getPolytomies(polytomies);

    cout << "SCM: Refining SCM Tree\n";
    
    VerboseMode saved_mode;
    saved_mode = verbose_mode;
    verbose_mode = VB_QUIET;

    int maxDegree = 0;

    for (auto polytomy: polytomies) {
        maxDegree = max(maxDegree, polytomy->degree());

        map<string, string> relabel;
        map<string, GeneNode*> delabel;
        int label = -1;
        scmTree->getRelabelMap(relabel, delabel, label, polytomy);

        Params params = *(this->params);
        PhyloSuperTreeUnlinked *newTree = new PhyloSuperTreeUnlinked(params, sourcesTree);

        for (int i = 0; i < newTree->size(); ++i) {
            GeneTree* tree = (GeneTree*)(*newTree)[i];
            tree->relabelAndCollapse(relabel);

            if (tree->leafNum < 4 || tree->branchNum <= tree->leafNum) {
                iter_swap(newTree->begin() + i, newTree->end() - 1);
                delete newTree->back();
                newTree->pop_back();
                --i;
            }
        }
        
        if (!newTree->empty()) {
            newTree->doMRP();
            GeneTree *mrpTree = newTree->mrpTree;
            mrpTree->delabelTree(delabel);
            mrpTree->root = NULL;
            assert(polytomy->degree() == 0);
            delete polytomy;
        }
        
        delete newTree;
    }

    verbose_mode = saved_mode;

    scmTree->reInitializeTree();

    cout << "SCM: Refining SCM Tree successfully with " 
    + to_string(polytomies.size()) + " polytomies"
    + " and max degree " + to_string(maxDegree) << "\n";

    cout << fixed << setprecision(2) << "SCM: Resolution of refined SCM Tree: " << 1.0 * (scmTree->nodeNum - scmTree->leafNum - 1) / (scmTree->leafNum - 3) << "\n";
}

void PhyloSuperTreeUnlinked::printResultWithSCMTree() {
    GeneTree *tmpTree = new GeneTree(firstSCMTree);

    string treeFile(this->params->out_prefix);

    tmpTree->printResultTree(treeFile + ".scm", false);

    string drawFile(this->params->out_prefix);
    drawFile += ".draw";

    ofstream out;
    out.exceptions(ios::failbit | ios::badbit);
    out.open(drawFile.c_str());
    
    out << "STRICT CONSENSUS MERGER TREE\n--------------------------------------------------------\n\n";
    tmpTree->drawTree(out, WT_SORT_TAXA);
    out << "\n\n";
    
    if (params->mrp_type != MRPType::MRP_NONE) {
        scmTree->printResultTree(treeFile + ".treefile", false);

        out << "SCM + MRP TREE\n--------------------------------------------------------\n\n";
        scmTree->drawTree(out, WT_BR_SCALE | WT_SORT_TAXA);
        out << "\n\n";
    
        printScoreWithConAln(scmTree, "SCM");
    }
    
    out.close();
    delete tmpTree;
}

void PhyloSuperTreeUnlinked::printScoreWithConAln(GeneTree *tree, string treeType) {
    if (conAln == NULL) {
        return;
    }
    
    IQTree *tmpTree = new IQTree(conAln);
    tmpTree->copyTree(tree);

    Params newParams = *(this->params);
    newParams.gbo_replicates = 0;
    tmpTree->setParams(newParams);

    cout << "\nSCORE OF " + treeType + " TREE: " << tmpTree->computeParsimony() << "\n";
    
    delete tmpTree;
}

StrVector PhyloSuperTreeUnlinked::createBootstrapGeneTrees() {
    StrVector bootstrapGeneTrees;
    for (auto it = begin(); it != end(); it++) {
        GeneTree* tree = (GeneTree*)(*it);
        int randIdx = random_int(tree->boot_trees.size());
        bootstrapGeneTrees.push_back(tree->getBootstrapTree(randIdx));
    }
    return bootstrapGeneTrees;
}