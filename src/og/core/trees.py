

num = len


class BaseNode(object):
    def __init__(self, id=None):
        self.id = id
        self.descendants = []

    def __hash__(self):
        return id(self)


    
class TreeNode(BaseNode):
    def __init__(self, id=None, depth=-1):
        BaseNode.__init__(self, id)
        self.depth = depth
        
    def __gt__(self, other):
        return self.depth < other.depth

    def __eq__(self, other):
        return self.depth == other.depth

    def __ne__(self, other):
        return self.depth != other.depth
    
    def __ge__(self, other):
        return self < other or self == other

    def __lt__(self, other):
        return self.depth > other.depth

    def __le__(self, other):
        return self > other or self == other

    def __str__(self):
        return str(self.id)

    def __repr__(self):
        return '%s(%s)' % (
            self.__class__.__name__, self.id,
        )
        
    

class BaseTree(object):
    def __init__(self):
        self.nodes = []

        

class AdjacencyTree(BaseTree):
    def __init__(self):
        BaseTree.__init__(self)
        self.internal_indices = []
        self.terminal_indices = []

        
    def __iter__(self):
        nodes_stack = [0]        
        nodes_seen = [0] * num(self.nodes)
        while nodes_stack:
            node = nodes_stack.pop()
            if nodes_seen[node]:
                continue
            if self.nodes[node].descendants:
                nodes_stack.extend(reversed(self.nodes[node].descendants))
            nodes_seen[node] = 1
            
            yield self.nodes[node]

            
    def __reversed__(self):
        nodes_list = []
        nodes_stack = [0]
        nodes_seen = [0] * num(self.nodes)
        while nodes_stack:
            node = nodes_stack.pop()
            if nodes_seen[node]:
                continue
            if self.nodes[node].descendants:
                nodes_stack.extend(reversed(self.nodes[node].descendants))
            nodes_seen[node] = 1

            nodes_list.append(self.nodes[node])
            
        return reversed(nodes_list)
                    

    def get_edges(self, indices=False, reverse=False):
        edges_list = []
        nodes_stack = [(None, 0)]        
        nodes_seen = [0] * num(self.nodes)
        while nodes_stack:
            node_anc, node_dsc = nodes_stack.pop()
            if nodes_seen[node_dsc]:
                continue
            if self.nodes[node_dsc].descendants:
                nodes_stack.extend(
                    map(lambda node_n: (node_dsc, node_n),
                        reversed(self.nodes[node_dsc].descendants)
                    )
                )
            nodes_seen[node_dsc] = 1

            if indices:
                edges_list.append(
                    (node_anc, node_dsc)
                )
            else:
                if node_anc is None:
                    edges_list.append(
                        (node_anc, self.nodes[node_dsc])
                    )
                else:
                    edges_list.append(
                        (self.nodes[node_anc], self.nodes[node_dsc])
                    )

        if reverse:
            return reversed(edges_list)
        else:
            return iter(edges_list)
    
    @property
    def internal_nodes(self):
        return map(lambda i: self.nodes[i], self.internal_indices)

    
    @property
    def terminal_nodes(self):
        return map(lambda i: self.nodes[i], self.terminal_indices)



