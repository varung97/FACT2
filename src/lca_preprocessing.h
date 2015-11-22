#ifndef LCA_PREPROC_H_
#define LCA_PREPROC_H_

#include <vector>
#include <cmath>
#include <iostream>

#include "utils.h"


struct rmq_t {
	std::vector<std::vector<int> > M;
	std::vector<int> v;
	int block_size;
	std::vector<int> addresses;
	std::vector<int**> prep_blocks;
};

struct lca_t {
	std::vector<int> E, R;
	rmq_t* rmq_prep;

	lca_t() : rmq_prep(new rmq_t) {}
};


void rmq_preprocess(rmq_t* rmq_prep, std::vector<int>& v) {
	size_t size = v.size();
	int block_size = int_log2(size);
	int blocks = size / block_size + (size % block_size != 0);

	std::vector<int> Ap(blocks, INT32_MAX);
	std::vector<int> B(blocks, 0);
	for (int i = 0; i < blocks; i++) {
		for (int j = 0; j < block_size; j++) {
			if (Ap[i] > v[i*block_size+j]) {
				Ap[i] = v[i*block_size+j];
				B[i] = i*block_size + j;
			}
		}
	}

	rmq_prep->M.push_back(std::vector<int>(Ap.size(), 0));
	for (unsigned int j = 0; j < Ap.size(); j++) {
		rmq_prep->M[0][j] = B[j];
	}
	for (int i = 1; (1 << i) <= (int) Ap.size(); i++) {
		rmq_prep->M.push_back(std::vector<int>(Ap.size()-(1 << i)+1, 0));
		for (unsigned int j = 0; j < Ap.size()-(1 << i)+1; j++) {
			if (v[rmq_prep->M[i-1][j]] <= v[rmq_prep->M[i-1][j+(1 << (i-1))]]) {
				rmq_prep->M[i][j] = rmq_prep->M[i-1][j];
			} else {
				rmq_prep->M[i][j] = rmq_prep->M[i-1][j+(1 << (i-1))];
			}
		}
	}

	rmq_prep->prep_blocks.resize(size, NULL);
	for (int i = 0; i < blocks; i++) {
		int address = 0;
		for (int j = 1; j < block_size; j++) {
			address <<= 1;
			address |= (v[i*block_size + j] > v[i*block_size + j-1]);
		}

		rmq_prep->addresses.push_back(address);

		if (rmq_prep->prep_blocks[address] == NULL) {
			rmq_prep->prep_blocks[address] = alloc_int_matrix(block_size);
			for (int j = 0; j < block_size; j++) {
				rmq_prep->prep_blocks[address][j][j] = j;
				for (int k = j+1; k < block_size; k++) {
					rmq_prep->prep_blocks[address][j][k] = rmq_prep->prep_blocks[address][j][k-1];
					if (v[i*block_size + rmq_prep->prep_blocks[address][j][k-1]] > v[i*block_size + k]) {
						rmq_prep->prep_blocks[address][j][k] = k;
					}
				}
			}
		}

		for (int j = 0; j < block_size; j++) {
			for (int k = j; k < block_size; k++) {
				int min_pos = j;
				for (int h = j; h <= k; h++) {
					if (v[i*block_size + min_pos] > v[i*block_size + h]) {
						min_pos = h;
					}
				}
			}
		}
	}

	//rmq_prep->v = v;
	rmq_prep->block_size = block_size;
}

int rmq2(rmq_t* rmq_prep, int a, int b) {

	int ba = (a/rmq_prep->block_size) + (a % rmq_prep->block_size != 0);
	int bb = (b/rmq_prep->block_size) - (b % rmq_prep->block_size != rmq_prep->block_size-1);

	int range = bb-ba+1;

	// There is no block to be checked, so it does not matter what I return until it is a valid index in [a-b]
	if (range == 0) return a;

	int k = 0;
	while (range >>= 1) ++k;

	if (rmq_prep->v[rmq_prep->M[k][ba]] <= rmq_prep->v[rmq_prep->M[k][bb-(1<<k)+1]]) {
		return rmq_prep->M[k][ba];
	} else {
		return rmq_prep->M[k][bb-(1<<k)+1];
	}
}

int rmq(rmq_t* rmq_prep, int a, int b) {
	if (a/rmq_prep->block_size == b/rmq_prep->block_size) { //same block
		int block_idx = a/rmq_prep->block_size;
		return block_idx*rmq_prep->block_size + rmq_prep->prep_blocks[rmq_prep->addresses[a/rmq_prep->block_size]]
		                             [a%rmq_prep->block_size][b%rmq_prep->block_size];
	}

	int min_pos = rmq2(rmq_prep, a, b);

	if (a % rmq_prep->block_size != 0) {
		int temp_min = rmq_prep->prep_blocks[rmq_prep->addresses[a/rmq_prep->block_size]]
		                                     [a%rmq_prep->block_size][rmq_prep->block_size-1];
		if (rmq_prep->v[min_pos] > rmq_prep->v[(a/rmq_prep->block_size)*rmq_prep->block_size + temp_min]) {
			min_pos = (a/rmq_prep->block_size)*rmq_prep->block_size + temp_min;
		}
	}
	if (b % rmq_prep->block_size != rmq_prep->block_size-1) {
		int temp_min = rmq_prep->prep_blocks[rmq_prep->addresses[b/rmq_prep->block_size]]
		                                     [0][b%rmq_prep->block_size];
		if (rmq_prep->v[min_pos] > rmq_prep->v[(b/rmq_prep->block_size)*rmq_prep->block_size + temp_min]) {
			min_pos = (b/rmq_prep->block_size)*rmq_prep->block_size + temp_min;
		}
	}

	return min_pos;
}


int lca(lca_t* lca_prep, int u, int v) {
	if (lca_prep->R[u] < lca_prep->R[v]) {
		return lca_prep->E[rmq(lca_prep->rmq_prep, lca_prep->R[u], lca_prep->R[v])];
	} else {
		return lca_prep->E[rmq(lca_prep->rmq_prep, lca_prep->R[v], lca_prep->R[u])];
	}
}

void eulerian_walk(Tree::Node* node, std::vector<int>& E, std::vector<int>& L, std::vector<int>& R, int depth) {
	E.push_back(node->id);
	L.push_back(depth);
	if (R[node->id] == -1) R[node->id] = E.size()-1;
	for (Tree::Node* child : node->children) {
		eulerian_walk(child, E, L, R, depth+1);
		E.push_back(node->id);
		L.push_back(depth);
	}
}

void resize_to_logmul(std::vector<int>& v) {
	int log_size = int_log2(v.size());
	if (v.size() % log_size != 0)
		v.resize((v.size()/log_size + 1)*log_size, v[v.size()-1]+1);
}

lca_t* lca_preprocess(Tree* t) {
	lca_t* lca_prep = new lca_t;

	lca_prep->R.resize(t->get_nodes_num(), -1);
	eulerian_walk(t->get_root(), lca_prep->E, lca_prep->rmq_prep->v, lca_prep->R, 0);

	resize_to_logmul(lca_prep->E);
	resize_to_logmul(lca_prep->rmq_prep->v);

	rmq_preprocess(lca_prep->rmq_prep, lca_prep->rmq_prep->v);
	return lca_prep;
}


#endif
