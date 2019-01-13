/*
 * freqdiff.h
 *
 *  Created on: 21 Oct 2015
 *      Author: Mesh
 */

#ifndef FREQDIFF_H_
#define FREQDIFF_H_

#include <iostream>
#include <queue>
#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include "taxas_ranges.h"
#include "lca_preprocessing.h"
#include "utils.h"

#include "Tree.h"
#include "RMQ.h"

struct node_bitvec_t {
	Tree::Node* node;
	boost::dynamic_bitset<>* bitvector;

	node_bitvec_t() : node(NULL), bitvector(NULL) {}
	node_bitvec_t(Tree::Node* node, boost::dynamic_bitset<>* bitvector) : node(node), bitvector(bitvector) {}
	//~node_bitvec_t() { delete bitvector; } FIXME
};
bool operator < (const node_bitvec_t& bv1, const node_bitvec_t& bv2) {
	return *(bv1.bitvector) < *(bv2.bitvector);
}

struct subpath_query_info_t {
	size_t nodes_num;
	Tree::Node** cp_roots;
	gen_rmq_t** cp_rmqs;
	int* depths;

	subpath_query_info_t(size_t nodes_num) : nodes_num(nodes_num), cp_roots(new Tree::Node*[nodes_num]),
			cp_rmqs(new gen_rmq_t*[nodes_num]), depths(new int[nodes_num]) {}
	~subpath_query_info_t() {
		delete[] cp_roots;
		for (size_t i = 0; i < nodes_num; i++) {
			delete cp_rmqs[i];
		}
		delete[] cp_rmqs;
		delete[] depths;
	}
};

int* start,* stop;
int* e,* m;
std::vector<Tree::Node*>* rsort_lists;

bool* marked;

Tree::Node** left,** right;
size_t* orig_pos_in_parent;

size_t* counter;
bool* BT;
int* leaf_p_index;

int* vleft,* vright;
int* pointer;
int* levels,* ids;
int* parent;
bool* exists;
Tree::Node** tree_nodes;
Tree::Node empty_node = Tree::Node(-1, -1, 0);


// calculate cluster weights using kn^2 method
void calc_w_kn2(std::vector<Tree*>& trees) {

	size_t n = Tree::get_taxas_num();

	int tot_int_nodes = 0;
	for (Tree* tree : trees) {
		tot_int_nodes += tree->get_nodes_num() - n - 1;
	}

	// Generate bit vectors
	node_bitvec_t* node_bitvecs = new node_bitvec_t[tot_int_nodes];
	int bitvc = 0;
	for (Tree* tree : trees) {
		taxas_ranges_t* tr = build_taxas_ranges(tree);

		size_t nodes_num = tree->get_nodes_num();
		for (size_t i = 1; i < nodes_num; i++) {
			Tree::Node* node = tree->get_node(i);
			if (!node->is_leaf()) {
				// fill bit vector for current node
				node_bitvec_t node_bitvec(node, new boost::dynamic_bitset<>(Tree::get_taxas_num()));
				for (int j = tr->intervals[node->id].start; j <= tr->intervals[node->id].end; j++) {
					node_bitvec.bitvector->set(tr->taxas[j]);
				}
				node_bitvecs[bitvc] = node_bitvec;
				bitvc++;
			}
		}

		delete tr;
	}

	// Sort bit vectors
	std::sort(node_bitvecs, node_bitvecs+tot_int_nodes);

	// count adjacent equal bitvectors and set weights
	for (int i = 0; i < bitvc; ) {
		int w = 1;
		while (i+w < bitvc && *(node_bitvecs[i].bitvector) == *(node_bitvecs[i+w].bitvector)) w++;
		for (int j = i; j < i+w; j++) {
			node_bitvecs[j].node->weight = w;
		}
		i += w;
	}

	for (int i = 0; i < bitvc; i++) {
		delete node_bitvecs[i].bitvector;
	}
	delete[] node_bitvecs;
}


void calc_w_k2n(std::vector<Tree*>& trees) {

	size_t n = Tree::get_taxas_num();

	int tot_int_nodes = 0;
	for (Tree* tree : trees) {
		tot_int_nodes += tree->get_nodes_num();
	}

	int* tree_start = new int[trees.size()];
	int* cache_taxa = new int[tot_int_nodes];
	int* cache_parent_ids = new int[tot_int_nodes];
	int* cache_size = new int[tot_int_nodes];
	int* cache_weight = new int[tot_int_nodes];
	std::fill(cache_weight, cache_weight+tot_int_nodes, 0);

	int pos = 0;
	for (size_t i = 0; i < trees.size(); i++) {
		Tree* tree = trees[i];
		tree_start[i] = pos;
		for (size_t j = 0; j < tree->get_nodes_num(); j++) {
			cache_taxa[pos] = tree->get_node(j)->taxa;
			if (j > 0) {
				cache_parent_ids[pos] = tree->get_node(j)->parent->id;
			}
			cache_size[pos] = tree->get_node(j)->size;
			pos++;
		}
	}

	int* e = new int[n];
	std::pair<int,int>* clusters = new std::pair<int,int>[n*2];
	std::pair<int,int>* H = new std::pair<int,int>[n*2];
	int* H2node = new int[n*2];
	std::pair<int,int>* minmax = new std::pair<int,int>[n*2];
	for (size_t i = 0; i < trees.size(); i++) {
		std::fill(H, H+Tree::get_taxas_num()*2, std::pair<int,int>(INT32_MAX, 0));

		Tree* tree1 = trees[i];

		// calculate e
		int count = 0;
		for (size_t j = 0; j < tree1->get_nodes_num(); j++) {
			if (tree1->get_node(j)->is_leaf()) {
				e[tree1->get_node(j)->taxa] = count++;
			}
			tree1->get_node(j)->weight++;
		}

		// calculate L and R for each cluster
		for (int j = tree1->get_nodes_num()-1; j > 0; j--) {
			Tree::Node* node = tree1->get_node(j);
			if (node->is_leaf()) {
				clusters[j].first = clusters[j].second = e[node->taxa];
			} else {
				clusters[j].first = clusters[node->children[0]->id].first;
				clusters[j].second = clusters[(*node->children.rbegin())->id].second;

				int pos;
				pos = node->pos_in_parent == node->parent->get_children_num()-1 ? clusters[j].first : clusters[j].second;
				H[pos] = clusters[j];
				H2node[pos] = tree_start[i]+j;
			}
		}

		// for each tree, see if cluster is in our ref tree, and if so increase its count by 1
		for (size_t j = i+1; j < trees.size(); j++) {
			Tree* tree2 = trees[j];

			std::fill(minmax, minmax+tree2->get_nodes_num(), std::pair<int,int>(INT32_MAX, 0));
			for (int k = tree2->get_nodes_num()-1; k > 0; k--) {
				int taxa = cache_taxa[tree_start[j]+k];
				int parent_id = cache_parent_ids[tree_start[j]+k];
				int size = cache_size[tree_start[j]+k];

				if (taxa != Tree::Node::NONE) {
					minmax[k].first = minmax[k].second = e[taxa];
				}

				if (minmax[parent_id].first > minmax[k].first) {
					minmax[parent_id].first = minmax[k].first;
				}
				if (minmax[parent_id].second < minmax[k].second) {
					minmax[parent_id].second = minmax[k].second;
				}

				if (taxa == Tree::Node::NONE && minmax[k].second-minmax[k].first+1 == size) {
					if (H[minmax[k].first] == minmax[k]) {
						cache_weight[H2node[minmax[k].first]]++;
						cache_weight[tree_start[j]+k]++;
					} else if (H[minmax[k].second] == minmax[k]) {
						cache_weight[H2node[minmax[k].second]]++;
						cache_weight[tree_start[j]+k]++;
					}
				}
			}
		}
	}

	for (size_t i = 0; i < trees.size(); i++) {
		Tree* tree = trees[i];
		for (size_t j = 0; j < tree->get_nodes_num(); j++) {
			tree->get_node(j)->weight += cache_weight[tree_start[i]+j];
		}
	}

	delete[] tree_start;
	delete[] cache_taxa;
	delete[] cache_parent_ids;
	delete[] cache_size;
	delete[] cache_weight;

	delete[] e;
	delete[] clusters;
	delete[] H;
	delete[] H2node;
	delete[] minmax;
}


// for node_id in T1:
// start[node_id] = smallest rank (i.e. leftmost leaf) in T2 among \Lambda(T1[node_id])
// stop[node_id] = biggest rank (i.e. rightmost leaf) in T2 among \Lambda(T1[node_id])
void compute_start_stop(Tree* tree1, Tree* tree2, int* t2_leaves_ranks) {
	// for each cluster in t1, find its leftmost and rightmost leaves in t2
	for (int i = tree1->get_nodes_num()-1; i >= 0; i--) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf()) {
			start[i] = stop[i] = t2_leaves_ranks[node->taxa];
		} else {
			start[i] = INT32_MAX;
			stop[i] = 0;
			for (Tree::Node* child : node->children) {
				if (start[i] > start[child->id])
					start[i] = start[child->id];
				if (stop[i] < stop[child->id])
					stop[i] = stop[child->id];
			}
		}
	}
}


void filter_clusters_n2(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, taxas_ranges_t* t2_tr, lca_t* t2_lcas,
		bool* to_del) {

	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);
	delete t2_ranks;

	// mark clusters in t1 to be deleted if a heavier incompatible cluster is in t2
	for (size_t i = 1; i < tree1->get_nodes_num(); i++) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf())
			continue;

		Tree::Node* leftmost_leaf = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* rightmost_leaf = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		int lcaX = lca(t2_lcas, leftmost_leaf->id, rightmost_leaf->id);
		marked[lcaX] = true;
		for (int j = t1_tr->intervals[i].start; j <= t1_tr->intervals[i].end; j++) {
			//mark all ancestors
			Tree::Node* curr = tree2->get_leaf(t1_tr->taxas[j]);
			while (!marked[curr->id]) {
				marked[curr->id] = true;
				curr = curr->parent;
			}
		}

		for (int j = t2_tr->intervals[lcaX].start; j <= t2_tr->intervals[lcaX].end; j++) {
			Tree::Node* leaf = tree2->get_leaf(t2_tr->taxas[j]);
			if (!marked[leaf->id]) { // leaf not in X, mark all ancestors up to lcaX
				Tree::Node* curr = leaf->parent;
				while (curr->id != lcaX) {
					if (marked[curr->id] && node->weight <= curr->weight) {
						// weight of node i in t1 is lower than an incompatible cluster in t2
						// mark for deletion
						to_del[i] = true;
						break;
					}
					curr = curr->parent;
				}
			}
			if (to_del[i])
				break;
		}

		marked[lcaX] = false;
		for (int j = t1_tr->intervals[i].start; j <= t1_tr->intervals[i].end; j++) {
			//mark all ancestors
			Tree::Node* curr = tree2->get_leaf(t1_tr->taxas[j]);
			while (marked[curr->id]) {
				marked[curr->id] = false;
				curr = curr->parent;
			}
		}
	}
}


subpath_query_info_t* preprocess_subpaths_queries(Tree* tree) {
	subpath_query_info_t* subpath_query_info = new subpath_query_info_t(Tree::get_taxas_num()*2);

	// pointer from each node to the root of its cp
	Tree::Node** cp_roots = subpath_query_info->cp_roots;
	cp_roots[0] = tree->get_root();
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		Tree::Node* node = tree->get_node(i);
		if (node->pos_in_parent == 0) {
			cp_roots[i] = cp_roots[node->parent->id];
		} else {
			cp_roots[i] = node;
		}
	}

	gen_rmq_t** cp_rmqs = subpath_query_info->cp_rmqs;
	gen_rmq_t* nullp = NULL;

    // This initialisation looks to be taking O(n) instead of O(m) time - that probably violates the efficiency guarantee
	std::fill(cp_rmqs, cp_rmqs+Tree::get_taxas_num()*2, nullp);
	for (size_t i = 0; i < tree->get_nodes_num(); i++) {
		if (cp_rmqs[cp_roots[i]->id] == NULL) {
			cp_rmqs[cp_roots[i]->id] = new gen_rmq_t;
		}
		cp_rmqs[cp_roots[i]->id]->v.push_back(-tree->get_node(i)->weight);
	}
	for (size_t i = 0; i < tree->get_nodes_num(); i++) {
		if (cp_rmqs[i] != NULL && cp_rmqs[i]->v.size() > 1) {
			general_rmq_preprocess(cp_rmqs[i]); // TODO: better to remove v from struct
		}
	}

	int* depths = subpath_query_info->depths;
	depths[0] = 0; // root depth
	for (size_t i = 1; i < tree->get_nodes_num(); i++) {
		depths[i] = 1 + depths[tree->get_node(i)->parent->id];
	}

	return subpath_query_info;
}

int max_subpath_query(subpath_query_info_t* subpq_info, Tree::Node* ancestor, Tree::Node* descendant) { // TODO: we might need only ids
	Tree::Node* curr = descendant;
	int res = 0;
	while (subpq_info->cp_roots[curr->secondary_id]->id != subpq_info->cp_roots[ancestor->secondary_id]->id) {
		gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->secondary_id]->id];
		int query_endp = subpq_info->depths[curr->secondary_id] - subpq_info->depths[subpq_info->cp_roots[curr->secondary_id]->id];
		if (query_endp == 0) {
			res = std::min(res, currpath_rmq->v[0]);
		} else {
			res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, 0, query_endp)]);
		}
		curr = subpq_info->cp_roots[curr->secondary_id]->parent;
	}
	gen_rmq_t* currpath_rmq = subpq_info->cp_rmqs[subpq_info->cp_roots[curr->secondary_id]->id];
	int query_startp = subpq_info->depths[ancestor->secondary_id] - subpq_info->depths[subpq_info->cp_roots[curr->secondary_id]->id];
	int query_endp = subpq_info->depths[curr->secondary_id] - subpq_info->depths[subpq_info->cp_roots[curr->secondary_id]->id];
	if (query_startp < query_endp) {
		res = std::min(res, currpath_rmq->v[general_rmq(currpath_rmq, query_startp+1, query_endp)]);
	}
	return -res;
}


subpath_query_info_t* subpq_info; // TODO: temporary
// leaves in marked must be sorted in left-to-right order
Tree* contract_tree_fast(Tree* tree, lca_t* lca_prep, std::vector<int>& marked) {	//TODO: think about externalizing lca_t
	if (marked.empty()) return NULL;

	int count = 0;

	levels[count] = tree->get_leaf(marked[0])->depth;
	ids[count] = tree->get_leaf(marked[0])->secondary_id;
	count++;
	for (size_t i = 1; i < marked.size(); i++) {
		int lca_id = lca(lca_prep, tree->get_leaf(marked[i-1])->id, tree->get_leaf(marked[i])->id);
		levels[count] = tree->get_node(lca_id)->depth;
		ids[count] = tree->get_node(lca_id)->secondary_id;
		count++;
		levels[count] = tree->get_leaf(marked[i])->depth;
		ids[count] = tree->get_leaf(marked[i])->secondary_id;
		count++;
	}

	//std::vector<int> pointer(levels.size(), -1);
	for (int i = 0; i < count; i++) pointer[i] = i;

	std::stack<int> Sl;
	for (int i = 0; i < count; i++) {
		while (!Sl.empty() && levels[i] <= levels[Sl.top()]) {
			Sl.pop();
		}

		if (!Sl.empty()) {
			vleft[i] = Sl.top();
		} else {
			vleft[i] = -1;
		}
		Sl.push(i);
	}

	std::stack<int> Sr;
	for (int i = count-1; i >= 0; i--) {
		while (!Sr.empty() && levels[i] <= levels[Sr.top()]) {
			Sr.pop();
		}

		if (!Sr.empty()) {
			vright[i] = Sr.top();
		} else {
			vright[i] = -1;
		}
		Sr.push(i);
	}

	int root_pos = -1;
	std::fill(exists, exists+count, true);
	for (int i = 0; i < count; i++) {
		parent[i] = -1;
		if (vleft[i] == -1 && vright[i] == -1) { //root
			if (root_pos == -1) {
				root_pos = i;
			}
		} else if (vleft[i] == -1) {
			parent[i] = vright[i];
		} else if (vright[i] == -1) {
			parent[i] = pointer[vleft[i]];
		} else {
			if (levels[vleft[i]] >= levels[vright[i]]) {
				parent[i] = pointer[vleft[i]];
				if (levels[vleft[i]] == levels[vright[i]]) { // deals with non-binary nodes
					pointer[vright[i]] = pointer[vleft[i]];
					exists[vright[i]] = false;
				}
			} else {
				parent[i] = pointer[vright[i]];
			}
		}
	}

	Tree* new_tree = new Tree(count*2);
	Tree::Node* nullp = NULL;
	std::fill(tree_nodes, tree_nodes+count, nullp);
	for (int i = 0; i < count; i++) {
		// create node
		if (i%2 == 0 || exists[i]) {
			if (tree_nodes[i] == NULL) {
				if (i%2 == 0) { //a leaf
					tree_nodes[i] = new_tree->add_node(marked[i/2]);
				} else {
					tree_nodes[i] = new_tree->add_node();
				}
			}
			tree_nodes[i]->weight = tree->get_node(ids[i])->weight;
			tree_nodes[i]->secondary_id = ids[i];
		}

		// connect to parent
		if (tree_nodes[i] != NULL && parent[i] != -1) {
			if (tree_nodes[parent[i]] == NULL) {
				tree_nodes[parent[i]] = new_tree->add_node();
			}
			tree_nodes[parent[i]]->add_child(tree_nodes[i]);
		}
	}

	// insert special nodes
	size_t newtree_nodes = new_tree->get_nodes_num();
	for (size_t i = 0; i < newtree_nodes; i++) {
		Tree::Node* curr_node = new_tree->get_node(i);

		if (curr_node->is_root()) continue;

		Tree::Node* desc_par = tree->get_node(curr_node->secondary_id)->parent;
		int msq = max_subpath_query(subpq_info, curr_node->parent, desc_par);
		if (msq > 0) {
			Tree::Node* sp_node = new_tree->add_node();
			sp_node->weight = msq;
			sp_node->secondary_id = curr_node->parent->secondary_id;
			curr_node->parent->null_child(curr_node->pos_in_parent);
			curr_node->parent->add_child(sp_node);
			sp_node->add_child(curr_node);
		}
	}
	new_tree->fix_tree(tree_nodes[root_pos]);

	return new_tree;
}


Tree* orig_t2 = NULL; // TODO: temporary
bool is_special(Tree::Node* node) {
	return node->size < orig_t2->get_node(node->secondary_id)->size;
}
void filter_clusters_nlog2n(Tree::Node* t1_root, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas, lca_t* orig_t2_lcas,
		bool* to_del) {

	std::priority_queue<std::pair<int, int> > BTw;

	Tree::Node* curr = t1_root;
	while (!curr->is_leaf()) {
		curr = curr->children[0];
	}

	int beta = 0;

	Tree::Node* rim1 = tree2->get_leaf(curr->taxa);
	curr = curr->parent;
	Tree::Node* p_2 = curr;

	/* Handle side trees */
	int p_index = 0;

	// build sorted sets of leaves
	leaf_p_index[t1_tr->taxas[t1_tr->intervals[t1_root->id].start]] = -1; // p_1 must be ignored

	std::vector<Tree::Node*> st_roots;
	while (curr != t1_root->parent) {
		for (size_t i = 1; i < curr->get_children_num(); i++) {
			Tree::Node* child = curr->children[i];
			for (int j = t1_tr->intervals[child->id].start; j <= t1_tr->intervals[child->id].end; j++) {
				leaf_p_index[t1_tr->taxas[j]] = p_index;
			}
			st_roots.push_back(child);

			p_index++;
		}
		curr = curr->parent;
	}
	taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
	std::vector<int>* leaves_sets = new std::vector<int>[p_index];
	for (size_t i = 0; i < t2_tr->taxas_num; i++) {
		if (leaf_p_index[t2_tr->taxas[i]] >= 0) {
			leaves_sets[leaf_p_index[t2_tr->taxas[i]]].push_back(t2_tr->taxas[i]);
		}
	}

	// build restrictions and recurse
	for (int i = 0; i < p_index; i++) {
		if (leaves_sets[i].size() > 1) {
			Tree* contracted = contract_tree_fast(orig_t2, orig_t2_lcas, leaves_sets[i]);
			lca_t* lca_ct2 = lca_preprocess(contracted);
			filter_clusters_nlog2n(st_roots[i], contracted, t1_tr, lca_ct2, orig_t2_lcas, to_del);
			delete lca_ct2;
			delete contracted;
		}
	}
	delete t2_tr;
	delete[] leaves_sets;

	std::fill(BT, BT+tree2->get_nodes_num(), false);
	std::fill(counter, counter+tree2->get_nodes_num(), 0);

	// T_2[p_1] should is "initialized" separately
	counter[rim1->parent->id] = 1;
	rim1 = rim1->parent;

	/* Handle centroid path */
	int curr_t1_start = t1_tr->intervals[t1_root->id].start + 1;
	curr = p_2;
	while (curr != t1_root->parent) {
		int curr_t1_end = t1_tr->intervals[curr->id].end;

		Tree::Node* ri = rim1;
		for (int i = curr_t1_start; i <= curr_t1_end; i++) {
			ri = tree2->get_node(lca(t2_lcas, ri->id, tree2->get_leaf(t1_tr->taxas[i])->id));
		}

		if (counter[rim1->id] == rim1->size && ri != rim1 && rim1->get_children_num() != 1) {
			// if children of rim1 are all "marked", should not go into BT
			counter[rim1->parent->id] += rim1->size;
			if (rim1->weight > beta && is_special(rim1)) beta = rim1->weight;
			rim1 = rim1->parent;
		}
		while (rim1 != ri) {
			BT[rim1->id] = true;
			BTw.push(std::make_pair(rim1->weight, rim1->id));
			if (rim1->weight > beta && is_special(rim1)) beta = rim1->weight;
			rim1 = rim1->parent;
		}

		for (int i = curr_t1_start; i <= curr_t1_end; i++) {
			Tree::Node* x = tree2->get_leaf(t1_tr->taxas[i]);
			while (!BT[x->id] && x != ri) {
				BT[x->id] = true;
				BTw.push(std::make_pair(x->weight, x->id));
				if (x->weight > beta && is_special(x)) beta = x->weight; // TODO: and is special node
				x = x->parent;
			}
		}
		for (int i = curr_t1_start; i <= curr_t1_end; i++) {
			Tree::Node* x = tree2->get_leaf(t1_tr->taxas[i]);
			counter[x->id]++;
			while (x != ri && counter[x->id] == x->size && !is_special(x)) { // TODO: check
				counter[x->parent->id] += x->size;
				BT[x->id] = false;
				x = x->parent;
			}
		}

		std::pair<int, int> top;
		while (!BTw.empty()) {

			top = BTw.top();
			if (BT[top.second]) break;
			BTw.pop();
		}
		int M = BTw.empty() ? 0 : top.first;
		if (curr->weight <= std::max(M, beta)) {
			to_del[curr->id] = true;
		}

		curr = curr->parent;
		curr_t1_start = curr_t1_end + 1;
	}
}


void compute_m(Tree::Node* node, int* e, int* m, std::vector<Tree::Node*>* rsort_lists) {
	if (node->is_leaf()) {
		m[node->id] = e[node->taxa];
	}
	for (Tree::Node* child : node->children) {
		compute_m(child, e, m, rsort_lists);
		if (m[node->id] > m[child->id]) {
			m[node->id] = m[child->id];
		}
	}
	if (!node->is_root()) {
		rsort_lists[m[node->id]].push_back(node);
	}
}


// See Section 2.4 of paper [XXX]
void merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas) {
	//compute e
	for (size_t j = 0; j < Tree::get_taxas_num(); j++) {
		e[t1_tr->taxas[j]] = j;
	}
	std::fill(m, m+tree2->get_nodes_num(), INT32_MAX);
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) rsort_lists[i].clear();

	compute_m(tree2->get_root(), e, m, rsort_lists);

	// sort tree2
	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		tree2->get_node(i)->clear_children();
	}
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		for (auto it = rsort_lists[i].begin(); it != rsort_lists[i].end(); it++) {
			(*it)->parent->add_child(*it);
		}
	}

	taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);
	delete[] t2_ranks;

	// calc x_left and x_right
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		Tree::Node* curr = tree2->get_leaf(i);
		Tree::Node* parent = curr->parent;
		while (parent != NULL && *(parent->children.begin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		left[i] = curr;

		curr = tree2->get_leaf(i);
		parent = curr->parent;
		while (parent != NULL && *(parent->children.rbegin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		right[i] = curr;
	}

	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		orig_pos_in_parent[i] = tree2->get_node(i)->pos_in_parent;
	}

	for (int i = tree1->get_nodes_num()-1; i >= 1; i--) {
		Tree::Node* a = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* b = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		if (a == b) continue; // tree1->node i is a leaf

		Tree::Node* ru = tree2->get_node(lca(t2_lcas, a->id, b->id));

		Tree::Node* a_left = left[a->taxa];
		Tree::Node* b_right = right[b->taxa];

		size_t du_pos = (a_left->depth > ru->depth) ? orig_pos_in_parent[a_left->id] : 0;
		size_t eu_pos = (b_right->depth > ru->depth) ? orig_pos_in_parent[b_right->id] : ru->get_children_num()-1;
		if (du_pos == 0 && eu_pos == ru->get_children_num()-1) continue;

		Tree::Node* newnode = tree2->add_node();
		newnode->weight = tree1->get_node(i)->weight;
		newnode->size = tree1->get_node(i)->size;
		for (size_t j = du_pos; j <= eu_pos; j++) {
			if (ru->children[j] != NULL) {
				newnode->add_child(ru->children[j]);
				// lazy child deletion
				ru->null_child(j);
			}
		}
		ru->set_child(newnode, du_pos);
	}
	tree2->fix_tree();

	delete t2_tr;
}

// Returns true if node has a marked leaf in subtree
// Populates is_lca, true if node is lca of marked leaves
bool mark_lcas(Tree::Node* node, std::vector<bool>* is_marked_leaf, std::vector<bool>* is_lca) {
	if (is_marked_leaf->at(node->id)) {
		is_lca->at(node->id) = true;
		return true;
	}

	int num_children_with_marked_leaves = 0;

	for (Tree::Node* child : node->children) {
		if (mark_lcas(child, is_marked_leaf, is_lca)) {
			num_children_with_marked_leaves++;
		}
	}

	if (num_children_with_marked_leaves >= 2) {
		is_lca->at(node->id) = true;
	}

	return num_children_with_marked_leaves > 0;
}

// Builds the subtree of tree restricted to some leaves given marked lcas
Tree::Node* build_restricted_subtree(Tree* tree, Tree::Node* root, std::vector<bool>* is_lca, std::vector<Tree::Node*>* assoc_nodes) {
    bool root_is_lca = is_lca->at(root->id);
    
    if (root->is_leaf()) {
        if (root_is_lca) {
            // If root is a marked leaf then create new node
            // Leaves are labelled with 1, this is only useful for trees with a single leaf
            return assoc_nodes->at(root->id) = tree->add_node(root->taxa, 1);
        } else {
            // Else no subtree
            return NULL;
        }
    }

    // If there is only one child with non-NULL subtree, then we can just return that subtree
    // If more than one child have non-NULL subtrees, then root is an lca and we need to connect each of these subtrees to the node
    Tree::Node* node = root_is_lca ? tree->add_node() : NULL;
    for (Tree::Node* child : root->children) {
        Tree::Node* child_subtree = build_restricted_subtree(tree, child, is_lca, assoc_nodes);
        
        if (child_subtree != NULL) {
            if (root_is_lca) {
                node->add_child(child_subtree);
            } else {
                assoc_nodes->at(root->id) = child_subtree;
                return child_subtree;
            }
        }
    }
    
    // At this point, either root is an lca or has no subtree
    if (root_is_lca) {
        return assoc_nodes->at(root->id) = node;
    } else {
        return NULL;
    }
}

// Builds the subtree of tree restricted to marked_leaves
// Also populates the associated node for each node in tree, defined as shallowest descendant of the node in the restricted subtree
Tree* restricted_subtree(Tree* tree, std::vector<int>* marked_taxa, std::vector<Tree::Node*>* assoc_nodes) {
	if (marked_taxa->empty()) return NULL;

	std::vector<bool> is_lca(tree->get_nodes_num(), false);
	std::vector<bool> is_marked_leaf(tree->get_nodes_num(), false);

    for (int taxon : *marked_taxa) {
		is_marked_leaf[tree->get_leaf(taxon)->id] = true;
	}

	mark_lcas(tree->get_root(), &is_marked_leaf, &is_lca);

	// Build new tree with only lcas
    Tree* new_tree = new Tree(std::count(is_lca.begin(), is_lca.end(), true));
    Tree::Node* root = build_restricted_subtree(new_tree, tree->get_root(), &is_lca, assoc_nodes);
    new_tree->fix_tree(root);
    
    return new_tree;
}

// Stably counting sorts elements of a vector along axis given by index
std::vector<std::vector<int>> counting_sort(std::vector<std::vector<int>>& vector, int index) {
    // Asummed to be non-negative
    int max_label = 0;
    for (int i = 0; i < vector.size(); i++) {
        max_label = max_label > vector[i][index] ? max_label : vector[i][index];
    }
    
    // Each bucket will hold a list of elements
    std::vector<std::vector<std::vector<int>>> buckets(max_label + 1);
    
    // Bucketize
    for (int i = 0; i < vector.size(); i++) {
        buckets[vector[i][index]].push_back(vector[i]);
    }
    
    std::vector<std::vector<int>> sorted_vector;
    
    // Reassemble
    for (int i = 0; i < buckets.size(); i++) {
        for (int j = 0; j < buckets[i].size(); j++) {
            sorted_vector.push_back(buckets[i][j]);
        }
    }

    return sorted_vector;
}

// Labels tree using the labels of the associated nodes in the restricted subtrees
void label_trees_using_assoc_nodes(std::vector<Tree*> trees, std::vector<std::vector<Tree::Node*>>& assoc_nodes_vector1, std::vector<std::vector<Tree::Node*>>& assoc_nodes_vector2) {
    // A label pair looks like (left_label, right_label, tree_index, node_id)
    std::vector<std::vector<int>> label_pairs;
    
    for (int tree_index = 0; tree_index < trees.size(); tree_index++) {
        for (int node_id = 0; node_id < trees[tree_index]->get_nodes_num(); node_id++) {
            label_pairs.push_back({assoc_nodes_vector1[tree_index][node_id]->label, assoc_nodes_vector2[tree_index][node_id]->label, tree_index, node_id});
        }
    }

    // Stable counting sort along each of the labels
    label_pairs = counting_sort(label_pairs, 1);
    label_pairs = counting_sort(label_pairs, 0);
    
    // Identical nodes stores (tree_index, node_id) of identical nodes
    std::vector<std::vector<int>> identical_nodes;
    int label = 1;
    
    for (int i = 0; i <= label_pairs.size(); i++) {
        // If identical nodes is not empty and this pair is different from the previous one then label the identical nodes
        // Also do so if i == label_pairs.size() since there are no more pairs to be checked
        if (i == label_pairs.size() ||
            (!identical_nodes.empty() &&
             (label_pairs[i - 1][0] != label_pairs[i][0] || label_pairs[i - 1][1] != label_pairs[i][1])
             )
            ) {

            for (int j = 0; j < identical_nodes.size(); j++) {
                trees[identical_nodes[j][0]]->get_node(identical_nodes[j][1])->label = label;
            }
            
            label++;
            identical_nodes.clear();
        }
        
        if (i < label_pairs.size()) {
            identical_nodes.push_back({label_pairs[i][2], label_pairs[i][3]});
        }
    }
}

// Divide marked_leaves into 2 parts and recursively label
void label_trees_helper(std::vector<Tree*>& trees, std::vector<int>& marked_taxa) {
    if (marked_taxa.size() == 1) {
        return;
    }
    
    std::vector<int> marked_taxa1, marked_taxa2;
    
    int i = 0;
    for (; i < marked_taxa.size() / 2; i++) {
        marked_taxa1.push_back(marked_taxa[i]);
    }
    for (; i < marked_taxa.size(); i++) {
        marked_taxa2.push_back(marked_taxa[i]);
    }
    
    std::vector<Tree*> trees1, trees2;
    std::vector<std::vector<Tree::Node*>> assoc_nodes_vector1, assoc_nodes_vector2;
    
    for (Tree* tree : trees) {
        // Default assoc node is the empty node, representing no marked taxa in leafset
        std::vector<Tree::Node*> assoc_nodes1(tree->get_nodes_num(), &empty_node);
        std::vector<Tree::Node*> assoc_nodes2(tree->get_nodes_num(), &empty_node);

        trees1.push_back(restricted_subtree(tree, &marked_taxa1, &assoc_nodes1));
        trees2.push_back(restricted_subtree(tree, &marked_taxa2, &assoc_nodes2));
        
        assoc_nodes_vector1.push_back(assoc_nodes1);
        assoc_nodes_vector2.push_back(assoc_nodes2);
    }
    
    label_trees_helper(trees1, marked_taxa1);
    label_trees_helper(trees2, marked_taxa2);
    
    label_trees_using_assoc_nodes(trees, assoc_nodes_vector1, assoc_nodes_vector2);
}

// Labels nodes in the trees, giving identical labels to nodes with identical leafsets
void label_trees(std::vector<Tree*>& trees) {
    std::vector<int> marked_taxa;
    for (int i = 0; i < Tree::get_taxas_num(); i++) {
        marked_taxa.push_back(i);
    }
    
    label_trees_helper(trees, marked_taxa);
}

void calc_w_knlogn(std::vector<Tree*>& trees) {
    label_trees(trees);
    
    // A label looks like (label, tree_index, node_id)
    std::vector<std::vector<int>> labels;
    
    for (int tree_index = 0; tree_index < trees.size(); tree_index++) {
        for (int node_id = 0; node_id < trees[tree_index]->get_nodes_num(); node_id++) {
            labels.push_back({trees[tree_index]->get_node(node_id)->label, tree_index, node_id});
        }
    }
    
    labels = counting_sort(labels, 0);
    
    // Identical nodes stores (tree_index, node_id) of identical nodes
    std::vector<std::vector<int>> identical_nodes;
    
    for (int i = 0; i <= labels.size(); i++) {
        // If identical nodes is not empty and this label is different from the previous one then compute freq for the identical nodes
        // Also do so if i == labels.size() since there are no more labels to be checked
        if (i == labels.size() ||
            (!identical_nodes.empty() && labels[i - 1][0] != labels[i][0])
            ) {

            for (int j = 0; j < identical_nodes.size(); j++) {
                trees[identical_nodes[j][0]]->get_node(identical_nodes[j][1])->weight = identical_nodes.size();
            }
            
            identical_nodes.clear();
        }
        
        if (i < labels.size()) {
            identical_nodes.push_back({labels[i][1], labels[i][2]});
        }
    }
}

// Returns rootsOfSubtrees(\leafset(T_A))
boost::unordered_set<Tree::Node*>* filter_clusters_nlogn_helper(Tree::Node* root_T_A, Tree* T_B, RMQ* rmq_T_B, bool* to_del_T_A) {
    Tree::Node* curr = root_T_A;
    std::vector<boost::unordered_set<Tree::Node*>> lower_boundaries = std::vector<boost::unordered_set<Tree::Node*>>();
    
    // recursively call filter_clusters on side trees and build up lower boundaries
    while (!curr->is_leaf()) {
        boost::unordered_set<Tree::Node*> lower_boundary = boost::unordered_set<Tree::Node*>();

        for (int i = 1; i < curr->get_children_num(); i++) {
            lower_boundary.merge(*filter_clusters_nlogn_helper(curr->children[i], T_B, rmq_T_B, to_del_T_A));
        }
        
        lower_boundaries.push_back(lower_boundary);
        curr = curr->children[0];
    }

    // Was built top to bottom, but more convenient in the other direction
    std::reverse(lower_boundaries.begin(), lower_boundaries.end());

    // First lca, leaf in T_B labelled by leaf in T_A
    Tree::Node* l_i = T_B->get_leaf(curr->taxa);
    l_i->counter = 0;
    l_i->parent->counter = 1;
    int lower_boundary_counter = -1; // Start at -1 since no lower boundary for leaf
    
    boost::unordered_set<Tree::Node*>* roots = new boost::unordered_set<Tree::Node*>({l_i});
    boost::heap::fibonacci_heap<int> incompatible = boost::heap::fibonacci_heap<int>();
    std::unordered_map<int, boost::heap::fibonacci_heap<int>::handle_type>* heap_handles = new std::unordered_map<int, boost::heap::fibonacci_heap<int>::handle_type>;

    
    while (curr != root_T_A) {
        curr = curr->parent;
        lower_boundary_counter++;
        
        Tree::Node* prev_l_i = l_i;
        for (Tree::Node* v : lower_boundaries[lower_boundary_counter]) {
            l_i = T_B->get_node(lca(rmq_T_B->lca_prep, l_i->id, v->id));
        }
        
        // Add nodes from prev_l_i to l_i
        // Include prev_l_i if incompatible with previous centroid path node
        if (prev_l_i->counter == prev_l_i->get_children_num()) {
            heap_handles->insert({prev_l_i->parent->id, incompatible.push(rmq_T_B->max_weight_path(prev_l_i->parent, l_i))});
        } else {
            heap_handles->insert({prev_l_i->id, incompatible.push(rmq_T_B->max_weight_path(prev_l_i, l_i))});
        }
        
        for (Tree::Node* v : lower_boundaries[lower_boundary_counter]) {
            v->parent->counter++;
            roots->insert(v);
            
            v = v->parent;
            
            while (v->counter == v->get_children_num()) {
                v->parent->counter++;

                roots->insert(v);
                for (Tree::Node* child : v->children) {
                    roots->erase(child);
                }

                // If v in incompatible then remove it
                if (heap_handles->find(v->id) != heap_handles->end()) {
                    // Increase key then delete-max
                    incompatible.increase(heap_handles->at(v->id), incompatible.top() + 1);
                    incompatible.pop();
                    heap_handles->erase(v->id);
                }
                
                v = v->parent;
            }
            
            // Add max weight from v until l_i to incompatible
            heap_handles->insert({v->id, incompatible.push(rmq_T_B->max_weight_path(v, l_i))});
        }
        
        int max_weight = heap_handles->size() == 0 ? 0 : incompatible.top();
        if (curr->weight <= max_weight) {
            to_del_T_A[curr->id] = true;
        }
    }
    
    for (Tree::Node* root : *roots) {
        root->parent->counter = 0;
    }
    
    return roots;
}

void filter_clusters_nlogn(Tree* T_A, Tree* T_B, lca_t* lca_T_B, bool* to_del_T_A) {
    RMQ* rmq_T_B = new RMQ(T_B, lca_T_B);

    filter_clusters_nlogn_helper(T_A->get_root(), T_B, rmq_T_B, to_del_T_A);
}

Tree* freqdiff(std::vector<Tree*>& trees, bool centroid_paths) {
    start = new int[Tree::get_taxas_num()*2];
	stop = new int[Tree::get_taxas_num()*2];
	e = new int[Tree::get_taxas_num()];
	m = new int[Tree::get_taxas_num()*2];
	rsort_lists = new std::vector<Tree::Node*>[Tree::get_taxas_num()];

	marked = new bool[Tree::get_taxas_num()*2];
    std::fill(marked, marked+Tree::get_taxas_num()*2, false);
	bool* to_del_t = new bool[Tree::get_taxas_num()*2];
	bool* to_del_ti = new bool[Tree::get_taxas_num()*2];

	left = new Tree::Node*[Tree::get_taxas_num()];
	right = new Tree::Node*[Tree::get_taxas_num()];
	orig_pos_in_parent = new size_t[Tree::get_taxas_num()*2];

	BT = new bool[Tree::get_taxas_num()*2];
	counter = new size_t[Tree::get_taxas_num()*2];
	leaf_p_index = new int[Tree::get_taxas_num()];

	// tree contract
	vleft = new int[Tree::get_taxas_num()*2];
	vright = new int[Tree::get_taxas_num()*2];
	pointer = new int[Tree::get_taxas_num()*2];
	levels = new int[Tree::get_taxas_num()*2];
	ids = new int[Tree::get_taxas_num()*2];
	parent = new int[Tree::get_taxas_num()*2];
	exists = new bool[Tree::get_taxas_num()*2];
	tree_nodes = new Tree::Node*[Tree::get_taxas_num()*2];

	for (size_t i = 0; i < trees.size(); i++) {
		trees[i]->reorder();
	}

    calc_w_knlogn(trees);

    // TODO: Move inside filter clusters maybe?
    lca_t** lca_preps = new lca_t*[trees.size()];
	for (size_t i = 0; i < trees.size(); i++) {
		lca_preps[i] = lca_preprocess(trees[i]);
	}

	subpath_query_info_t** subpath_query_ti = NULL;
	if (centroid_paths) {
		subpath_query_ti = new subpath_query_info_t*[trees.size()];
		for (size_t i = 0; i < trees.size(); i++) {
			subpath_query_ti[i] = preprocess_subpaths_queries(trees[i]);
		}
	}
    
	Tree* T = new Tree(trees[0]);
	for (size_t i = 1; i < trees.size(); i++) {
		Tree* Ti = new Tree(trees[i]);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(Ti);

		taxas_ranges_t* tr_T = NULL;
		tr_T = build_taxas_ranges(T);

		lca_t* lca_T = lca_preprocess(T);

		// filter clusters
		std::fill(to_del_ti, to_del_ti+Ti->get_nodes_num(), false);
		if (centroid_paths) {
			orig_t2 = T; // TODO: temporary
			subpq_info = preprocess_subpaths_queries(T);
			filter_clusters_nlog2n(Ti->get_root(), T, tr_Ti, lca_T, lca_T, to_del_ti);
			delete subpq_info;
		} else {
			filter_clusters_n2(Ti, T, tr_Ti, tr_T, lca_T, to_del_ti);
		}

		std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
		if (centroid_paths) {
			orig_t2 = Ti; // TODO: temporary
			subpq_info = subpath_query_ti[i];
			filter_clusters_nlog2n(T->get_root(), Ti, tr_T, lca_preps[i], lca_preps[i], to_del_t);
		} else {
			filter_clusters_n2(T, Ti, tr_T, tr_Ti, lca_preps[i], to_del_t);
		}

		Ti->delete_nodes(to_del_ti);
		T->delete_nodes(to_del_t);

		delete lca_T;
		delete tr_T;
		delete tr_Ti;

		lca_T = lca_preprocess(T);
		tr_Ti = build_taxas_ranges(Ti);

		merge_trees(Ti, T, tr_Ti, lca_T);

		delete lca_T;
		delete tr_Ti;

		// reset secondary ids
		for (size_t j = 0; j < T->get_nodes_num(); j++) {
			T->get_node(j)->secondary_id = T->get_node(j)->id;
		}

		delete Ti;
	}

	std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
	taxas_ranges_t* tr_T = build_taxas_ranges(T);
	for (size_t i = 0; i < trees.size(); i++) {
		if (centroid_paths) {
			orig_t2 = trees[i];
			subpq_info = subpath_query_ti[i];
			filter_clusters_nlog2n(T->get_root(), trees[i], tr_T, lca_preps[i], lca_preps[i], to_del_t);
		} else {
			taxas_ranges_t* tr_Ti = build_taxas_ranges(trees[i]);
			filter_clusters_n2(T, trees[i], tr_T, tr_Ti, lca_preps[i], to_del_t);
			delete tr_Ti;
		}

		delete trees[i];
	}
	T->delete_nodes(to_del_t);
	delete tr_T;

	delete[] start;
	delete[] stop;
	delete[] e;
	delete[] m;
	delete[] rsort_lists;

	delete[] marked;
	delete[] to_del_ti;
	delete[] to_del_t;

	delete[] orig_pos_in_parent;
	delete[] left;
	delete[] right;

	for (size_t i = 0; i < trees.size(); i++) {
		delete lca_preps[i];
	}
	delete[] lca_preps;

	if (centroid_paths) {
		for (size_t i = 0; i < trees.size(); i++) {
			delete subpath_query_ti[i];
		}
		delete[] subpath_query_ti;
	}

	return T;
}


#endif /* FREQDIFF_H_ */
