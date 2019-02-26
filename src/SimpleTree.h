#ifndef SIMPLETREE_H_
#define SIMPLETREE_H_

#include <vector>

class SimpleTree {
public:
	SimpleTree(size_t nodes_num_hint = 0);
	virtual ~SimpleTree();

	class SimpleNode;

	SimpleNode* get_node(int i);
	SimpleNode* get_root();
    void set_root(SimpleNode* root);
	size_t get_nodes_num();

	SimpleTree::SimpleNode* add_node(int taxa = -1, int label = 1);

private:
	std::vector<SimpleNode*> nodes;
    SimpleNode* root;
};

class SimpleTree::SimpleNode {
public:
	static const int NONE = -1;

	std::vector<SimpleNode*> children;

	SimpleNode* parent;
	int id;
	int taxa;
    int label;

	SimpleNode(int id, int taxa, int label);

	void add_child(SimpleNode* child);

	bool is_leaf();
	bool is_root();
};

#endif
