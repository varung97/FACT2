import sys


class Node:
    def __init__(self):
        self.value = None
        self.taxa = None
        self.children = []

    def is_leaf(self):
        return self.taxa is not None

    def to_newick(self):
        if self.is_leaf():
            return str(self.taxa)

        return '(' + ','.join([str(child) for child in self.children]) + ')'

    def __repr__(self):
        return self.to_newick()

def from_newick(newick):
    root = Node()

    if newick.isdigit():
        root.value = int(newick)
        root.taxa = root.value
        return root

    newick = newick[1:-1] + ','

    num_brackets = 0
    children = []
    child = ''

    for char in newick:
        if num_brackets == 0 and char == ',':
                children.append(child)
                child = ''
        else:
            child += char
            if char == '(':
                num_brackets += 1
            elif char == ')':
                num_brackets -= 1

    for child in children:
        root.children.append(from_newick(child))

    root.children.sort(key=lambda node: node.value)
    root.value = root.children[0].value

    return root

if __name__ == '__main__':
    print(str(from_newick(sys.argv[1])) == str(from_newick(sys.argv[2])))
