#include "Tree.h"
#include "lca_preprocessing.h"
#include <sdsl/bit_vectors.hpp>

#include "RMQ.h"

// Build depths array
void RMQ::preprocess_depths(Tree* tree) {
    depths = std::vector<int>(tree->get_nodes_num());
    depths[0] = 0;
    for (int i = 1; i < tree->get_nodes_num(); i++) {
        depths[i] = depths[tree->get_node(i)->parent->id] + 1;
    }
}

// Build and preprocess centroid path related information
void RMQ::preprocess_centroid_paths(Tree* tree) {
    cp_roots = std::vector<Tree::Node*>(tree->get_nodes_num());
    cp_rmqs = std::vector<gen_rmq_t*>();
    cp_indices = std::vector<int>(tree->get_nodes_num());
    cp_depths = std::vector<int>(tree->get_nodes_num());
    
    // Build and preprocess centroid path related information
    cp_roots[0] = tree->get_root();
    cp_rmqs.push_back(new gen_rmq_t);
    cp_rmqs[0]->v.push_back(-tree->get_root()->weight);
    cp_indices[0] = 0;
    cp_depths[0] = 0;
    
    // Assumes that nodes are iterated top-down
    for (int i = 1; i < tree->get_nodes_num(); i++) {
        Tree::Node* node = tree->get_node(i);
        if (node->pos_in_parent == 0) {
            cp_roots[i] = cp_roots[node->parent->id];
            cp_indices[i] = cp_indices[node->parent->id];
            cp_rmqs[cp_indices[i]]->v.push_back(-node->weight);
            cp_depths[i] = cp_depths[node->parent->id];
        } else {
            // Create new centroid path
            cp_roots[i] = node;
            
            gen_rmq_t* cp_rmq = new gen_rmq_t;
            cp_rmq->v.push_back(-node->weight);
            cp_rmqs.push_back(cp_rmq);
            
            cp_indices[i] = cp_rmqs.size() - 1;
            cp_depths[i] = cp_depths[node->parent->id] + 1;
        }
    }
    
    for (gen_rmq_t* cp_rmq : cp_rmqs) {
        if (cp_rmq->v.size() > 1) {
            general_rmq_preprocess(cp_rmq);
        }
    }
}

// Build centroid subpath rmq for each leaf
void RMQ::preprocess_leaf_rmqs(Tree* tree) {
    preprocess_centroid_paths(tree);
    preprocess_depths(tree);
    
    // Assumes all taxa are represented in tree
    std::vector<gen_rmq_t*> leaf_rmqs = std::vector<gen_rmq_t*>(Tree::get_taxas_num());
    
    for (int taxon = 0; taxon < Tree::get_taxas_num(); taxon++) {
        gen_rmq_t* leaf_rmq = new gen_rmq_t;

        Tree::Node* curr_node = tree->get_leaf(taxon);
        
        // Get the max weight of each centroid subpath on path from leaf to root
        while (true) {
            Tree::Node* cp_root = cp_roots[curr_node->id];

            int depth_in_cp = depths[curr_node->id] - depths[cp_root->id];
            gen_rmq_t* curr_rmq = cp_rmqs[cp_indices[curr_node->id]];
            
            // If node is already the root of centroid path, then no need to query the rmq
            int max_weight = depth_in_cp == 0 ? curr_rmq->v[0] : curr_rmq->v[general_rmq(curr_rmq, 0, depth_in_cp)];
            
            leaf_rmq->v.push_back(max_weight);
            
            if (cp_root == tree->get_root()) {
                break;
            } else {
                curr_node = cp_root->parent;
            }
        }
        
        if (leaf_rmq->v.size() > 1) {
            general_rmq_preprocess(leaf_rmq);
        }
        
        leaf_rmqs[taxon] = leaf_rmq;
    }
    
    // Point each node to a leaf rmq (for any leaf in its leafset)
    leaf_rmq_for_nodes = std::vector<gen_rmq_t*>(tree->get_nodes_num());

    // Assumes nodes are iterated bottom-up
    for (int i = tree->get_nodes_num() - 1; i >= 0; i--) {
        Tree::Node* curr_node = tree->get_node(i);
        if (curr_node->is_leaf()) {
            leaf_rmq_for_nodes[i] = leaf_rmqs[curr_node->taxa];
        } else {
            leaf_rmq_for_nodes[i] = leaf_rmq_for_nodes[curr_node->children[0]->id];
        }
    }
}

// Numbers leaves from left to right, starting at first_available index and populates min_leaf_idx_in_side_trees_cps_int
// Returns the new first_available_index
int RMQ::preprocess_min_leaf_idx_in_side_trees_helper(Tree::Node* node, int first_available_index) {
    if (node->is_leaf()) {
        min_leaf_idx_in_side_trees_cps_int[cp_indices[node->id]].push_back(first_available_index);
        return first_available_index + 1;
    }
    
    // This gives index of first leaf in side trees
    first_available_index = preprocess_min_leaf_idx_in_side_trees_helper(node->children[0], first_available_index);
    min_leaf_idx_in_side_trees_cps_int[cp_indices[node->id]].push_back(first_available_index);
    
    for (int i = 1; i < node->children.size(); i++) {
        first_available_index = preprocess_min_leaf_idx_in_side_trees_helper(node->children[i], first_available_index);
    }
    
    return first_available_index;
}

void RMQ::preprocess_min_leaf_idx_in_side_trees(Tree* tree) {
    min_leaf_idx_in_side_trees_cps_int = std::vector<std::vector<int>>(tree->get_leaves_num());
    preprocess_min_leaf_idx_in_side_trees_helper(tree->get_root(), 0);
    
    // Centroid paths were populated bottom to top, so need to be reversed
    for (int i = 0; i < min_leaf_idx_in_side_trees_cps_int.size(); i++) {
        std::reverse(min_leaf_idx_in_side_trees_cps_int[i].begin(), min_leaf_idx_in_side_trees_cps_int[i].end());
    }
}

void RMQ::populate_min_leaf_idx_in_side_trees_cps(Tree* tree) {
    min_leaf_idx_in_side_trees_cp_ranks = std::vector<sdsl::rank_support_v<1>>();

    for (int i = 0; i < min_leaf_idx_in_side_trees_cps_int.size(); i++) {
        // Each entry is offset by the smallest one to reduce space taken by the bit vector
        int bit_vector_size = min_leaf_idx_in_side_trees_cps_int[i][0] - min_leaf_idx_in_side_trees_cps_int[i].back() + 1;
        sdsl::bit_vector* min_leaf_idx_in_side_trees_cp = new sdsl::bit_vector(bit_vector_size, 0);
        
        for (int min_leaf_idx_in_side_trees : min_leaf_idx_in_side_trees_cps_int[i]) {
            // We also subtract 1 from each entry. The reason is illustrated below:
            // Say the entries are 0, 1, 4
            // Then the leafs in the side trees of the second node in the centroid path are 1, 2, 3
            // All of these should give us the same answer when a rank query is made; however 1 gives a different result
            // To account for this we subtract 1 from the entries
            if (min_leaf_idx_in_side_trees - min_leaf_idx_in_side_trees_cps_int[i].back() != 0) {
                (*min_leaf_idx_in_side_trees_cp)[min_leaf_idx_in_side_trees - min_leaf_idx_in_side_trees_cps_int[i].back() - 1] = 1;
            }
        }
        
        // Create a rank support structure for the bit vector
        min_leaf_idx_in_side_trees_cp_ranks.push_back(*new sdsl::rank_support_v<1>(min_leaf_idx_in_side_trees_cp));
    }
}

void RMQ::preprocess_leaf_rank_structures(Tree* tree) {
    // Number of centroid paths is the same as number of leaves
    preprocess_min_leaf_idx_in_side_trees(tree);
    populate_min_leaf_idx_in_side_trees_cps(tree);
}

RMQ::RMQ(Tree* tree) {
    std::cout << tree->to_string() << std::endl;
    preprocess_leaf_rmqs(tree);
    preprocess_leaf_rank_structures(tree);
}

// Max weight from v to the root of its centroid path
int RMQ::rmq_q1(Tree::Node* v, Tree::Node* w) {
    Tree::Node* cp_root = cp_roots[v->id];
    int depth_in_cp_v = depths[v->id] - depths[cp_root->id];
    gen_rmq_t* curr_rmq = cp_rmqs[cp_indices[v->id]];
    
    if (depth_in_cp_v == 0) {
        // If node is already the root of centroid path, then no need to query the rmq
        // Negate answer since values in rmqs are negative
        return -curr_rmq->v[0];
    } else if (cp_indices[v->id] == cp_indices[w->id]) {
        // If v and w are in the same centroid path, then query for the path between them
        int depth_in_cp_w = depths[w->id] - depths[cp_root->id];
        return -curr_rmq->v[general_rmq(curr_rmq, depth_in_cp_w, depth_in_cp_v)];
    } else {
        // The path from v to root of centroid path must be between v and w
        return -curr_rmq->v[general_rmq(curr_rmq, 0, depth_in_cp_v)];
    }
}

// Max weight from parent of cp root of v to shallowest cp root that is a proper descendant of w
int RMQ::rmq_q2_qg1(Tree::Node* v, Tree::Node* w) {
    if (cp_depths[v->id] - cp_depths[w->id] < 2) {
        // No centroid subpaths between v and w
        return 0;
    }
    
    gen_rmq_t* curr_rmq = leaf_rmq_for_nodes[v->id];
    // The leaf rmq must have at least 3 entries, so no need to check for size > 1
    return -curr_rmq->v[general_rmq(curr_rmq, cp_depths[w->id] + 1, cp_depths[v->id] - 1)];
}

// Max weight from w to deepest node in centroid path of w that is an ancestor of v
int RMQ::rmq_qg(Tree::Node* v, Tree::Node* w) {
    if (cp_depths[v->id] - cp_depths[w->id] == 0) {
        // v and w are on the same centroid path
        return 0;
    }
    
    // Get index of leftmost leaf in subtree at v
    int idx_of_leaf = min_leaf_idx_in_side_trees_cps_int[cp_indices[v->id]].back();
    int w_cp_idx = cp_indices[w->id];

    // Rank query gives index from bottom of centroid path, subtract from length of centroid path - 1 to get real index
    int deepest_node = min_leaf_idx_in_side_trees_cps_int[w_cp_idx].size() - 1 -
                       min_leaf_idx_in_side_trees_cp_ranks[w_cp_idx](idx_of_leaf - min_leaf_idx_in_side_trees_cps_int[w_cp_idx].back());
    gen_rmq_t* curr_rmq = cp_rmqs[w_cp_idx];
    
    return -curr_rmq->v[general_rmq(curr_rmq, deepest_node, depths[w->id] - depths[cp_roots[w->id]->id])];
}

int RMQ::rmq(Tree::Node* v, Tree::Node* w) {
    return std::max({rmq_q1(v, w), rmq_q2_qg1(v, w), rmq_qg(v, w)});
}
