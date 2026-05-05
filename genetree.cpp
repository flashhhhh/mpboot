#include "genetree.h"

#define PUT_MEANING(value, description) meanings.insert({#value, description})

GeneTree::GeneTree() : IQTree() {
    // Constructor implementation
}

GeneTree::GeneTree(Alignment *aln) : IQTree(aln) {
    // Constructor implementation
}

GeneTree::GeneTree(const string &treeString) : IQTree() {
    readTreeString(treeString);
    
    seedNode = (GeneNode*) root->neighbors[0]->node;
    assert(!seedNode->isLeaf());

    // Collapse all node has degree 2, except the seedNode
    dfsFixTree(seedNode, NULL);
}

void GeneTree::dfsFixTree(GeneNode *node, GeneNode *parent) {
    node->parent = parent;
    
    if (node->isLeaf()) {
        return;
    }

    for (auto child: node->getChildren()) {
        dfsFixTree(child, node);
    }

    if (node->parent != NULL && node->degree() == 2) {
        collapseEdge(node->parent, node);
    }
}

GeneTree::~GeneTree() {
    // Destructor implementation
}

GeneNode *GeneTree::newNode(int node_id, const char* node_name) {
    return new GeneNode(node_id, node_name);
}
GeneNode *GeneTree::newNode(int node_id, int node_name) {
    return new GeneNode(node_id, node_name);
}

void GeneTree::setOriginalNodeIndex(const map<string, int> &seqNameToIndex) {
    totalLeafNum = seqNameToIndex.size();
    
    NodeVector taxa;
    getTaxa(taxa);

    for (auto taxon : taxa) {
        assert(taxon->isLeaf() && seqNameToIndex.find(taxon->name) != seqNameToIndex.end());
        ((GeneNode*) (taxon))->setOriginalIndex(seqNameToIndex.at(taxon->name));
    }
}

void GeneTree::addSplitOfNode(GeneNode *node) {
    assert(node->split);
    assert(splitToNode.find(*node->split) == splitToNode.end());
    splitToNode[*node->split] = node;
}

void GeneTree::removeSplitOfNode(GeneNode *node) {
    assert(node->split);
    assert(splitToNode[*node->split] == node);
    assert(splitToNode.erase(*node->split) == 1);
}

Split* GeneTree::getLeavesSplit() {
    assert(seedNode != NULL);
    if (splitToNode.empty() || seedNode->split == NULL) {
        assert(splitToNode.empty() && seedNode->split == NULL);
        encodeSplits();
    }
    return seedNode->split;
}

void GeneTree::dfsEncodeSplits(GeneNode *node, GeneNode *parent) {
    node->parent = parent;
    assert(node->getNumChildren() != 1);
    
    if (node->split) {
        delete node->split;
    }
    node->split = new Split(totalLeafNum, 0);

    if (node->isLeaf()) {
        node->split->addTaxon(((GeneNode*) node)->getOriginalIndex());
    }

    for (auto child: node->getChildren()) {
        dfsEncodeSplits(child, node);
        *node->split += *child->split;
    }

    addSplitOfNode(node);
}

void GeneTree::encodeSplits() {
    assert(seedNode != NULL && totalLeafNum > 0);

    if (seedNode->split != NULL) {
        return;
    }
    
    splitToNode.clear();

    if (seedNode->getNumChildren() == 2) {
        deroot();
    }

    dfsEncodeSplits((GeneNode*) seedNode, NULL);
}

void GeneTree::resetAfterMergeSCM() {
    delete seedNode->split;
    seedNode->split = NULL;
    splitToNode.clear();
}

void GeneTree::makeFirstChildOfRoot(GeneNode *node) {
    assert(node->parent != NULL);
    reroot(node->parent);
    assert(seedNode == node->parent);

    for (int i = 0 ; i < seedNode->neighbors.size(); ++i) {
        if (seedNode->neighbors[i]->node == node) {
            swap(seedNode->neighbors[i], seedNode->neighbors[0]);
            return;
        }
    }
    assert(0);
}

void GeneTree::reroot(GeneNode *node) {
    GeneNode *oldPar = node->parent;
    if (oldPar == NULL) {
        assert(seedNode == node);
        return;
    }

    if (oldPar != seedNode) {
        reroot(oldPar);
    }

    assert(oldPar == seedNode && oldPar->getNumChildren() > 2);

    removeSplitOfNode(oldPar);
    removeSplitOfNode(node);

    oldPar->removeChild(node);
    *oldPar->split -= *node->split;

    node->addChild(oldPar);
    *node->split += *oldPar->split;
        
    addSplitOfNode(oldPar);
    addSplitOfNode(node);

    seedNode = node;
}

void GeneTree::deroot() {
    if (!seedNode || seedNode->degree() != 2) {
        return;
    }

    // cout << "Derooting tree" << endl;

    GeneNode *child0 = (GeneNode*) seedNode->neighbors[0]->node;
    GeneNode *child1 = (GeneNode*) seedNode->neighbors[1]->node;

    if (child0->degree() < 3 && child1->degree() < 3) {
        cout << "Error tree's topology, both children have degree < 3" << endl;
        assert(0);
    }
    if (child0->degree() < 3) {
        swap(child0, child1);
    }

    collapseEdge(seedNode, child0);
}

void GeneTree::collapseEdge(GeneNode *head, GeneNode *tail) {
    assert(tail->isLeaf() == false);
    assert(head->parent == tail || tail->parent == head);
    
    bool down = (tail->parent == head);
    
    // Collapse should be from head = parent to tail = child in this version
    assert(down);

    if (down) {
        head->removeChild(tail);
    
        for (auto neighbor: tail->getChildren()) {
            tail->removeChild(neighbor);
            head->addChild(neighbor);
        }
    } else {
        tail->removeChild(head);
        GeneNode *parent = tail->parent;
        if (parent) {
            parent->removeChild(tail);
            parent->addChild(head);
        }
    }

    delete tail;
}

void GeneTree::collapseSubtree(GeneNode *node, bool do_collapse) {
    if (node->isLeaf()) {
        return;
    }

    for (auto child: node->getChildren()) {
        collapseSubtree(child, true);
    }

    if (do_collapse) {
        assert(node->parent != NULL);
        collapseEdge(node->parent, node);
    }
}

void GeneTree::readTreeString(const string &treeString) {
    stringstream str;
	str << treeString;
	str.seekg(0, ios::beg);
	freeNode();
	readTree(str, rooted);
}

void GeneTree::printResultTree(string fileName, bool isAppend) {
    printTree(fileName.c_str(), WT_SORT_TAXA | WT_NEWLINE | (isAppend ? WT_APPEND : 0));
}

string GeneTree::getResultTreeString() {
    std::ostringstream oss;
    printTree(oss, WT_SORT_TAXA | WT_NEWLINE);
    return oss.str();
}

void GeneTree::reInitializeTree(GeneNode *node, GeneNode* parent) {
    if (!node) {
        setRootLeaf(NULL);
        node = (GeneNode*) root;
        nodeNum = getNumTaxa();
        leafNum = 0;
        branchNum = 0;
    }

    node->id = (node->isLeaf() ? leafNum++ : nodeNum++);
    if (node->isLeaf() == false) {
        node->name = "";
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        (*it)->id = branchNum;
        (*it)->node->findNeighbor(node)->id = branchNum;
        branchNum++;
        reInitializeTree((GeneNode*)(*it)->node, node);
    }
}

GeneNode* GeneTree::getAnyOtherLeaf(GeneNode *node, GeneNode *parent) {
    if (!node) {
        assert(root);
        node = (GeneNode*) root;
    } else if (node->isLeaf()) {
        return (node->name != root->name) ? node : NULL;
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        GeneNode *leaf = getAnyOtherLeaf(child, node);
        if (leaf) {
            return leaf;
        }
    }
    return NULL;
}

void GeneTree::setRootLeaf(char *my_root) {
    string root_name;
    if (my_root) {
        root_name = my_root;
    } else {
        if (root && root->isLeaf()) {
            return;
        }
        if (aln) {
            root_name = aln->getSeqName(0);
        } else {
            root = getAnyOtherLeaf();
            assert(root && root->isLeaf());
            return;
        }
    }
    root = findLeafName(root_name);
    assert(root);
}

void GeneTree::getPolytomies(vector<GeneNode*> &polytomies, GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    }

    if (node->degree() > 3) {
        polytomies.push_back(node);
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        getPolytomies(polytomies, child, node);
    }
}

void GeneTree::getRelabelMap(map<string, string> &relabel, map<string, GeneNode*> &delabel, int &label, GeneNode *node, GeneNode *parent) {
    node->parent = parent;

    if (node->isLeaf()) {
        relabel[node->name] = to_string(label);
        assert(parent != NULL);
    }

    for (auto child: node->getChildren()) {
        if (parent == NULL) {
            delabel[to_string(++label)] = child;
        }
        getRelabelMap(relabel, delabel, label, child, node);
    }
}

void GeneTree::relabelAndCollapse(const map<string, string> &relabel) {
    setRootLeaf(NULL);
    relabelTree(relabel);
    deleteDuplicateNode();
    GeneNode *newRoot = getAnyOtherLeaf();
    if (newRoot != NULL) {
        assert(root->name != newRoot->name);
        root = newRoot;
        deleteDuplicateNode();
    } else {
        assert(getNumTaxa() <= 2);
    }
    reInitializeTree();
    if (leafNum > 2) {
        map<string, int> checkUniqueLabel;
        doubleCheckUniqueName(checkUniqueLabel);
    }
}

void GeneTree::doubleCheckUniqueName(map<string, int> &checkUniqueLabel, GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    }

    if (node->isLeaf()) {
        assert(checkUniqueLabel.find(node->name) == checkUniqueLabel.end());
        checkUniqueLabel[node->name] = 1;
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        doubleCheckUniqueName(checkUniqueLabel, child, node);
    }
}

void GeneTree::relabelTree(const map<string, string> &relabel, GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    }

    if (node->isLeaf()) {
        assert(relabel.find(node->name) != relabel.end());
        node->name = relabel.at(node->name);
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        GeneNode *child = (GeneNode*) (*it)->node;
        relabelTree(relabel, child, node);
    }
}

void GeneTree::delabelTree(map<string, GeneNode*> &delabel, GeneNode *node, GeneNode *parent) {
    if (!node) {
        assert(root && root->isLeaf());
        node = (GeneNode*) root->neighbors[0]->node;
    }

    node->parent = parent;

    if (node->isLeaf()) {
        assert(delabel.find(node->name) != delabel.end());
        GeneNode *subtree = delabel[node->name];
        delabel.erase(node->name);

        assert(subtree->parent != NULL && parent != NULL);
        
        subtree->parent->removeChild(subtree);
        parent->removeChild(node);
        parent->addChild(subtree);
        
        delete node;
        return;
    }
    
    for (auto child: node->getChildren()) {
        delabelTree(delabel, child, node);
    }

    if (parent == NULL && delabel.empty() == false) {
        for (auto [_, subtree]: delabel) {
            assert(subtree->parent != NULL);
            subtree->parent->removeChild(subtree);
            node->addChild(subtree);
        }
        delabel.clear();
    }
}

void GeneTree::deleteDuplicateNode(GeneNode *node, GeneNode *parent) {
    if (!node) {
        node = (GeneNode*) root;
    } else if (node->isLeaf()) {
        node->parent = parent;
        return;
    }
    
    node->parent = parent;

    for (auto child: node->getChildren()) {
        deleteDuplicateNode(child, node);
    }

    bool doDelete = true;

    GeneNode *prevChild = NULL;

    for (auto child: node->getChildren()) {
        if (child->isLeaf() == false) {
            doDelete = false;
            break;
        }

        if (prevChild == NULL) {
            prevChild = child;
        } else if (child->name != prevChild->name) {
            doDelete = false;
            break;
        }
    }

    if (doDelete && parent != NULL) {
        parent->removeChild(node);
        node->removeChild(prevChild);
        parent->addChild(prevChild);

        for (auto child: node->getChildren()) {
            delete child;
        }
        delete node;
    }
}

void GeneTree::setNodeIdByMapName(const map<string, int> &seqNameToIndex) {
    NodeVector taxa;
    getTaxa(taxa);
    
    assert(taxa.size() == seqNameToIndex.size());

    for (auto taxon: taxa) {
        string seqName = taxon->name;
        assert(seqNameToIndex.find(seqName) != seqNameToIndex.end());
        taxon->id = seqNameToIndex.at(seqName);
    }
}

string GeneTree::getBootstrapTree(int index) {
    if (bootstrapTrees.empty()) {
        for (const auto &[tree, id]: treels) {
            assert(aln);
            
            GeneTree *tmpTree = new GeneTree(tree);
            
            NodeVector taxa;
            tmpTree->getTaxa(taxa);
            
            for (auto taxon: taxa) {
                taxon->name = aln->getSeqName(stoi(taxon->name));
            }

            bootstrapTrees[id] = tmpTree->getTreeString();
            
            delete tmpTree;
        }
    }

    assert(index >= 0 && index < boot_trees.size());
    assert(bootstrapTrees.find(boot_trees[index]) != bootstrapTrees.end());

    return bootstrapTrees[boot_trees[index]];
}

void GeneTree::extractQuadSubtrees(vector<Split*> &subtrees, BranchVector &branches, Node *node, Node *dad) {
	if (!node) node = root;
	FOR_NEIGHBOR_IT(node, dad, it) {
		extractQuadSubtrees(subtrees, branches, (*it)->node, node);
		if ((*it)->node->isLeaf()) continue;
		// internal branch
		assert(node->degree() == 3 && (*it)->node->degree() == 3);
		int cnt = 0;
		Node *child = (*it)->node;
        string treestrings[4];
        int nodeid = 0;
		FOR_NEIGHBOR_DECLARE(child, node, it2) {
			Split *sp = new Split(leafNum);
			getTaxa(*sp, (*it2)->node, child);
			subtrees.push_back(sp);
			cnt += sp->countTaxa();
		}
		FOR_NEIGHBOR(node, child, it2) {
			Split *sp = new Split(leafNum);
			getTaxa(*sp, (*it2)->node, node);
			subtrees.push_back(sp);
			cnt += sp->countTaxa();
		}
		assert(cnt == leafNum);
        branches.push_back({node, child});
	}
}

void GeneTree::computeGeneConcordance(vector<GeneTree*> trees, map<string,string> &meanings) {
    StrVector names;
    getTaxaName(names);

    StringIntMap name_map;
    for (auto stri = names.begin(); stri != names.end(); stri++)
        name_map[*stri] = stri - names.begin();

    BranchVector branches;
    vector<Split*> subtrees;
    extractQuadSubtrees(subtrees, branches, root->neighbors[0]->node);

    IntVector decisive_counts; // number of decisive trees
    decisive_counts.resize(branches.size(), 0);

    IntVector supports[3]; // number of trees supporting 3 alternative splits
    supports[0].resize(branches.size(), 0);
    supports[1].resize(branches.size(), 0);
    supports[2].resize(branches.size(), 0);

    string prefix[3] = {"gC", "gD1", "gD2"};

    for (int treeid = 0; treeid < trees.size(); treeid++) {
        GeneTree *tree = trees[treeid];

        StrVector taxname;
        tree->getTaxaName(taxname);

        // create the map from taxa between 2 trees
        Split taxa_mask(leafNum);
        for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
            if (name_map.find(*it) == name_map.end()) {
                if (*it == "__root__") {
                    cout << "WARNING : By default, trees without a multifurcation at the root are treated as rooted." << endl;
                    cout << "          You may need to change your tree structure." << endl;
                }
                outError("Taxon not found in full tree: ", *it);
            }
            taxa_mask.addTaxon(name_map[*it]);
        }

        // make the taxa ordering right before converting to split system
        taxname.clear();
        int smallid = 0;
        for (int taxid = 0; taxid < leafNum; taxid++) {
            if (taxa_mask.containTaxon(taxid)) {
                taxname.push_back(names[taxid]);
                tree->findLeafName(names[taxid])->id = smallid++;
            }
        }
        assert(taxname.size() == tree->leafNum);

        SplitGraph sg;
        //NodeVector nodes;
        tree->convertSplits(sg);

        SplitIntMap hash_ss;
        for (auto sit = sg.begin(); sit != sg.end(); sit++) {
            hash_ss.insertSplit((*sit), 1);
        }

        // now scan through all splits in current tree
        for (int id = 0, qid = 0; qid < subtrees.size(); id++, qid += 4) {
            Neighbor *nei = branches[id].second->findNeighbor(branches[id].first);

            bool decisive = true;
            for (int i = 0; i < 4; i++) {
                if (!taxa_mask.overlap(*subtrees[qid+i])) {
                    decisive = false;
                    break;
                }
            }

            if (!decisive) continue;
            
            decisive_counts[id]++;
            for (int i = 0; i < 3; i++) {
                Split this_split = *subtrees[qid]; // current split
                this_split += *subtrees[qid+i+1];
                Split *subsp = this_split.extractSubSplit(taxa_mask);
                if (subsp->shouldInvert())
                    subsp->invert();
                if (hash_ss.findSplit(subsp)) {
                    supports[i][id]++;
                }
                delete subsp;
            }
        }
    }

    for (int i = 0; i < branches.size(); i++) {
        if (decisive_counts[i] == 0)
            continue;
        Neighbor *nei = branches[i].second->findNeighbor(branches[i].first);
        int gN = decisive_counts[i];
        int gCF_N = supports[0][i];
        int gDF1_N = supports[1][i];
        int gDF2_N = supports[2][i];
        int gDFP_N = gN - gCF_N - gDF1_N - gDF2_N;
        double gCF = round((double)gCF_N/gN * 10000)/100;
        double gDF1 = round((double)gDF1_N/gN * 10000)/100;
        double gDF2 = round((double)gDF2_N/gN * 10000)/100;
        double gDFP = round((double)gDFP_N/gN * 10000)/100;
        PUT_ATTR(nei, gCF);
        PUT_ATTR(nei, gDF1);
        PUT_ATTR(nei, gDF2);
        PUT_ATTR(nei, gDFP);
        PUT_ATTR(nei, gN);
        PUT_ATTR(nei, gCF_N);
        PUT_ATTR(nei, gDF1_N);
        PUT_ATTR(nei, gDF2_N);
        PUT_ATTR(nei, gDFP_N);
        stringstream g_factors;
        g_factors << gCF << "/" << gDF1 << "/" << gDF2 << "/" << gDFP;
        nei->putAttr("gCF/gDF1/gDF2/gDFP", g_factors.str());
        stringstream g_factors_N;
        g_factors_N << gCF_N << "/" << gDF1_N << "/" << gDF2_N << "/" << gDFP_N;
        nei->putAttr("gCF_N/gDF1_N/gDF2_N/gDFP_N", g_factors_N.str());

        stringstream tmp;
        tmp.precision(3);
        tmp << (double)supports[0][i]/decisive_counts[i]*100;
        if (verbose_mode >= VB_MED)
            tmp << "%" << decisive_counts[i];
        
        Node *node = branches[i].second;
        if (params->newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']')
                node->name += ",[gCF=" + tmp.str() + "]";
            else
                node->name = node->name.substr(0, node->name.length()-1) + ",gCF=" + tmp.str() + "]";
        } else {
            if (!node->name.empty())
                node->name.append("/");
            node->name.append(tmp.str());
        }
    }
    for (vector<Split*>::reverse_iterator it = subtrees.rbegin(); it != subtrees.rend(); it++)
        delete (*it);

    PUT_MEANING(gCF, "Gene concordance factor (=gCF_N/gN %)");
    PUT_MEANING(gDF1, "Gene discordance factor for NNI-1 branch (=gDF1_N/gN %)");
    PUT_MEANING(gDF2, "Gene discordance factor for NNI-2 branch (=gDF2_N/gN %)");
    PUT_MEANING(gDFP, "Gene discordance factor due to polyphyly (=gDFP_N/gN %)");
    PUT_MEANING(gN, "Number of trees decisive for the branch");
    PUT_MEANING(gCF_N, "Number of trees concordant with the branch");
    PUT_MEANING(gDF1_N, "Number of trees concordant with NNI-1 branch");
    PUT_MEANING(gDF2_N, "Number of trees concordant with NNI-2 branch");
    PUT_MEANING(gDFP_N, "Number of trees decisive but discordant due to polyphyly");
    meanings.insert({"*NOTE*", "(gCF+gDF1+gDF2+gDFP) = 100% and (gCF_N+gDF1_N+gDF2_N+gDFP_N) = gN"});
}

void GeneTree::computeSiteConcordance(map<string,string> &meanings) {
    BranchVector branches;
    getInnerBranches(branches);

    for (auto ii = 0; ii < branches.size(); ii++) {
        BranchVector::iterator it = branches.begin()+ii;
        computeSiteConcordance((*it), params->site_concordance);
        Neighbor *nei = it->second->findNeighbor(it->first);
        double sCF = 0.0;
        if (!GET_ATTR(nei, sCF))
            continue;

        stringstream tmp;
        tmp.precision(3);
        tmp << sCF;
        string sup_str = tmp.str();
        Node *node = it->second;
        if (params->newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']') {
                node->name += ",[sCF=" + sup_str + "]";
            } else
                node->name = node->name.substr(0, node->name.length()-1) + ",sCF=" + sup_str + "]";
        } else {
            if (!node->name.empty())
                node->name += "/";
            node->name += sup_str;
        }
    }

    PUT_MEANING(sCF, "Site concordance factor averaged over " + convertIntToString(params->site_concordance) +  " quartets (=sCF_N/sN %)");
    PUT_MEANING(sN, "Number of informative sites averaged over " + convertIntToString(params->site_concordance) +  " quartets");
    PUT_MEANING(sDF1, "Site discordance factor for alternative quartet 1 (=sDF1_N/sN %)");
    PUT_MEANING(sDF2, "Site discordance factor for alternative quartet 2 (=sDF2_N/sN %)");
    PUT_MEANING(sCF_N, "sCF in absolute number of sites");
    PUT_MEANING(sDF1_N, "sDF1 in absolute number of sites");
    PUT_MEANING(sDF2_N, "sDF2 in absolute number of sites");
}

void GeneTree::computeSiteConcordance(pairNode &branch, int nquartets) {
    vector<IntVector> left_taxa, right_taxa;

    // extract the taxa from the two left subtrees
    left_taxa.resize(branch.first->neighbors.size()-1);
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, it) {
        if (rooted && (*it)->node == root)
            return;
        getTaxaID(left_taxa[id], (*it)->node, branch.first);
        id++;
    }
    assert(id == left_taxa.size());

    // extract the taxa from the two right subtrees
    right_taxa.resize(branch.second->neighbors.size()-1);
    id = 0;
    FOR_NEIGHBOR(branch.second, branch.first, it) {
        if (rooted && (*it)->node == root)
            return;
        getTaxaID(right_taxa[id], (*it)->node, branch.second);
        id++;
    }
    assert(id == right_taxa.size());

    if (rooted) {
        vector<IntVector>::iterator it;
        for (it = left_taxa.begin(); it != left_taxa.end(); it++)
            for (auto it2 = it->begin(); it2 != it->end(); it2++)
                if (*it2 == leafNum-1) {
                    it->erase(it2);
                    break;
                }
        for (it = right_taxa.begin(); it != right_taxa.end(); it++)
            for (auto it2 = it->begin(); it2 != it->end(); it2++)
                if (*it2 == leafNum-1) {
                    it->erase(it2);
                    break;
                }
    }

    double sCF = 0.0; // concordance factor
    double sDF1 = 0.0;
    double sDF2 = 0.0;
    double sN = 0.0;
    size_t sum_sites = 0;
    double sCF_N = 0, sDF1_N = 0, sDF2_N = 0;
    vector<int64_t> support;
    support.resize(3, 0);

    Neighbor *nei = branch.second->findNeighbor(branch.first);
    for (size_t i = 0; i < nquartets; ++i) {
        // get a random quartet
        IntVector quartet;
        quartet.resize(4);
        int left_id0 = 0, left_id1 = 1, right_id0 = 0, right_id1 = 1;
        if (left_taxa.size() > 2) {
            left_id0 = random_int(left_taxa.size());
            do {
                left_id1 = random_int(left_taxa.size());
            } while (left_id0 == left_id1);
        }
        if (right_taxa.size() > 2) {
            right_id0 = random_int(right_taxa.size());
            do {
                right_id1 = random_int(right_taxa.size());
            } while (right_id0 == right_id1);
        }
        quartet[0] = left_taxa[left_id0][random_int(left_taxa[left_id0].size())];
        quartet[1] = left_taxa[left_id1][random_int(left_taxa[left_id1].size())];
        quartet[2] = right_taxa[right_id0][random_int(right_taxa[right_id0].size())];
        quartet[3] = right_taxa[right_id1][random_int(right_taxa[right_id1].size())];

        support[0] = support[1] = support[2] = 0;
        aln->computeQuartetSupports(quartet, support);
        size_t sum = support[0] + support[1] + support[2];
        sum_sites += sum;
        if (sum > 0) {
            sCF += ((double)support[0]) / sum;
            sDF1 += ((double)support[1]) / sum;
            sDF2 += ((double)support[2]) / sum;
            sCF_N += support[0];
            sDF1_N += support[1];
            sDF2_N += support[2];
        }
    }

    sN = (double)sum_sites / nquartets;
    // rounding
    sCF = round(sCF / nquartets * 10000)/100;
    sDF1 = round(sDF1 / nquartets * 10000)/100;
    sDF2 = round(sDF2 / nquartets * 10000)/100;
    sCF_N = round(sCF_N / nquartets * 100)/100;
    sDF1_N = round(sDF1_N / nquartets * 100)/100;
    sDF2_N = round(sDF2_N / nquartets * 100)/100;
    PUT_ATTR(nei, sCF);
    PUT_ATTR(nei, sN);
    PUT_ATTR(nei, sDF1);
    PUT_ATTR(nei, sDF2);
    PUT_ATTR(nei, sCF_N);
    PUT_ATTR(nei, sDF1_N);
    PUT_ATTR(nei, sDF2_N);
    stringstream s_factors;
    s_factors << sCF << "/" << sDF1 << "/" << sDF2;
    nei->putAttr("sCF/sDF1/sDF2", s_factors.str());
    stringstream s_factors_N;
    s_factors_N << sCF_N << "/" << sDF1_N << "/" << sDF2_N;
    nei->putAttr("sCF_N/sDF1_N/sDF2_N", s_factors_N.str());
    // insert key-value for partition-wise con/discordant sites
    string keys[] = {"sC", "sD1", "sD2"};
    for (size_t i = 3; i < support.size(); ++i) {
        if (support[i] >= 0)
            nei->putAttr(keys[i%3] + convertIntToString(i/3), (double)support[i]/nquartets);
        else
            nei->putAttr(keys[i%3] + convertIntToString(i/3), "NA");
    }
}