#!/usr/bin/env python

# TODO: add proper parsing of quoted strings
# TODO: add to_string() method

import re

from og.constants import (
    _COLON,
    _COMMA,
    _DASH,
    _EMPTY,
    _LPAREN,
    _RPAREN,
    _SEMICOLON,
    _SPACE,
    _2QUOTE,
    _1QUOTE
)
from og.core.trees import AdjacencyTree, TreeNode
from math import isinf

_INAPPROPRIATE_WHITESPACE = 'Input newick string contains inappropriate use of white space'

_PARENTHESES = re.compile(r'''^\(.+\)(?:[0-9a-zA-Z!"#$%&'()*+,\-./:;<=>?@\[\]^_`{\|}~]+)?(?:\:\d+(?:\-\d+)?)?;?$''')
_SUBCLADE_NEW = _LPAREN
_SUBCLADE_END = _RPAREN
_SUBCLADE_NODE_SEP = _COMMA
_SUBCLADE_TERMINATORS = set((_SUBCLADE_NEW, _SUBCLADE_END, _SUBCLADE_NODE_SEP)) 
_SUBCLADE_NODE_ANNOT_SEP = _COLON
_BRANCH_LENGTH_INTERVAL_SEP = _DASH
_NODE = 'n{:d}'.format

num = len


def _format_id(id):
    return _EMPTY if id is None else str(id)
        

class NewickTreeNode(TreeNode):
    def __init__(self, id=None, length=None, depth=-1):
        TreeNode.__init__(self, id, depth)
        self.length = length

        
    def __repr__(self):
        return '%s(%s, length=%s, depth=%s)' % (
            self.__class__.__name__,
            str(self.id),
            str(self.length),
            str(self.depth)
        )
        


class NewickTree(AdjacencyTree):
    def __init__(self, string=None, default_length=None):
        AdjacencyTree.__init__(self)
        self._default_length = default_length
        self.from_string(string)

        
    def __str__(self):
        return self.to_string()

    
    def _parse_branch_length(self, string):
        length = string.lstrip(_SUBCLADE_NODE_ANNOT_SEP)
        if len(length) > 0:
            length = float(length)
        return length


    def _format_branch_length(self, length):
        if length is None:
            return _EMPTY
        return str(length)

    
    def _parse_node_string(self, string, node=None):
        if node is None:
            node = NewickTreeNode()

        string = string.strip()
        for i in range(len(string)):
            if string[~i] == _SUBCLADE_NODE_ANNOT_SEP:
                node.length = self._parse_branch_length(string[~i:])
                string = string[:~i].strip()
                break
        if len(string) > 0:
            node.id = string.strip(_1QUOTE + _2QUOTE)
        if not node.length:
            node.length = self._default_length
        return node
 

    def _format_node_string(self, node, fmt_id=_format_id):
        if node is None:
            return _EMPTY
        
        node_id = fmt_id(node.id)
        length = self._format_branch_length(node.length)
        if length:
            return _SUBCLADE_NODE_ANNOT_SEP.join((node_id,length))
        else:
            return node_id
            
    
    def from_string(self, tree_string):
        if tree_string is None:
            return
        
        tree_string = tree_string.strip().rstrip(_SPACE + _SEMICOLON)
        if tree_string.count(_SUBCLADE_NEW) != tree_string.count(_SUBCLADE_END):
            raise ValueError('Input newick string contains mismatched number of opening and closing parentheses')
        if not _PARENTHESES.match(tree_string):
            raise ValueError('Input newick string lacks matching outer parentheses')

        N = 0     # tree node index
        I = 0     # for consecutive internal node IDs
        i = 0     # tree string index
        P = None  # previous operation
        L = len(tree_string)
        depth = 0
        node_stack = list()
        while i < L:
            if tree_string[i] == _SPACE:
                if P != _SUBCLADE_NODE_SEP:
                    raise ValueError(_INAPPROPRIATE_WHITESPACE)

            elif tree_string[i] == _SUBCLADE_NEW:
                # new internal node
                if N > 0:
                    self.nodes[node_stack[-1]].descendants.append(N)
                self.internal_indices.append(N)
                self.nodes.append(NewickTreeNode(I, self._default_length, depth))
                node_stack.append(N)
                depth += 1
                P = _SUBCLADE_NEW
                b = i + 1
                I += 1
                N += 1

            elif tree_string[i] == _SUBCLADE_END:
                # end of internal node
                n = node_stack.pop()
                if P != _SUBCLADE_END:
                    # terminal node
                    self.nodes.append(self._parse_node_string(tree_string[b:i]))
                    self.nodes[n].descendants.append(N)
                    self.nodes[-1].depth = depth
                    self.terminal_indices.append(N)
                    N += 1
                depth -= 1
                    
                # must scan ahead for branch lengths
                b = i + 1
                e = i + 1
                while e < L and tree_string[e] not in _SUBCLADE_TERMINATORS:
                    if tree_string[e] == _SPACE:
                        raise ValueError(_INAPPROPRIATE_WHITESPACE)
                    e += 1
                
                if b < e:
                    self._parse_node_string(tree_string[b:e], node=self.nodes[n])
                
                P = _SUBCLADE_END
                b = e

            elif tree_string[i] == _SUBCLADE_NODE_SEP:
                if P != _SUBCLADE_END:
                    # terminal node
                    self.nodes.append(self._parse_node_string(tree_string[b:i]))
                    self.nodes[node_stack[-1]].descendants.append(N)
                    self.nodes[-1].depth = depth
                    self.terminal_indices.append(N)
                    N += 1
                
                P = _SUBCLADE_NODE_SEP
                b = i + 1

            i += 1

            
    def to_string(self, internal_ids=False, lengths=True):
        def _formatter(id):
            if id is None:
                return _EMPTY
            if isinstance(id, int):
                return _NODE(id) if internal_ids else _EMPTY
            else:
                return _format_id(id)
            
        tree_strings = [''] * num(self.nodes)
        for dsc_index in self.terminal_indices:
            tree_strings[dsc_index] = self._format_node_string(self.nodes[dsc_index])
            
        for anc_index in reversed(self.internal_indices):
            anc_node = self.nodes[anc_index]
            tree_strings[anc_index] = '%s%s%s%s' % (
                _SUBCLADE_NEW,
                _SUBCLADE_NODE_SEP.join(
                    map(tree_strings.__getitem__, anc_node.descendants)
                ),
                _SUBCLADE_END,
                self._format_node_string(self.nodes[anc_index], _formatter)
            )
        return tree_strings[0]
            
            

class _BranchLengthInterval(object):
    def __init__(self, minimum, maximum): 
        self.minimum = minimum
        self.maximum = maximum

    def __repr__(self):
        return '%s(minimum=%s, maximum=%s)' % (
            self.__class__.__name__,
            str(self.minimum),
            str(self.maximum)
        )
        
    def __str__(self):
        return '[%s, %s]' % (str(self.minimum), str(self.maximum))
    

    
class IntervalNewickTree(NewickTree):
    def __init__(self, string=None, default_length=(None, None)):
        if num(default_length) != 2:
            raise ValueError('Expected a 2-element list or 2-tuple for default_length')
        # Must bootstrap self._default_length to set NewickTree._default_length
        self._default_length = _BranchLengthInterval(default_length[0], default_length[1])  
        NewickTree.__init__(self, string, self._default_length)


    def _get_branch_length_interval(self, lengths):
        interval = _BranchLengthInterval(
            self._default_length.minimum,
            self._default_length.maximum
        )
        if num(lengths) > 0:
            if lengths[0] is not None and \
               lengths[0] is not _EMPTY:
                interval.minimum = float(lengths[0])
            if num(lengths) == 1:
                interval.maximum = float(lengths[0])            
        if num(lengths) > 1:
            if lengths[1] is not None and \
               lengths[1] is not _EMPTY:
                interval.maximum = float(lengths[1])
        return interval


    def _parse_branch_length(self, string):
        return self._get_branch_length_interval(
            string.strip().lstrip(_SUBCLADE_NODE_ANNOT_SEP).split(_BRANCH_LENGTH_INTERVAL_SEP)
        )

    
    def _format_branch_length(self, length):
        _min = str(length.minimum)
        _max = str(length.maximum)
        if isinf(length.minimum):
            _min = _EMPTY
        if isinf(length.maximum):
            _max = _EMPTY
        if _min is _EMPTY and _max is _EMPTY:
            return _EMPTY
        else:
            return _BRANCH_LENGTH_INTERVAL_SEP.join((_min,_max))

    
            
if __name__ == '__main__':
    import sys
    from math import inf as _POS_INF
    _NEG_INF = -1.0 * _POS_INF
    tree = IntervalNewickTree(sys.argv[1], default_length=(_NEG_INF, _POS_INF))
    for node in tree:
        print(repr(node), node.descendants)

    print(tree.to_string(internal_ids=True))
    print(tree.to_string(internal_ids=False))
    print(tree)
