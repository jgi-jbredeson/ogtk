
# Invertebrate chromosomes have the following possible letter codes:
# A1, A2, B1, B2, B3, C1, C2, D, E, F, G, H, I, J1, J2, K, L, M, N, O1, O2, P,
#   Q, R
# in the files you have sometimes these are separated by "x" because there are
# some invertebrate chromosomes that are combinations of these codes. For
# example, I think there is a Sca chromosome that is "B1xR".  This should be
# read as the genes on this chromosome having ancestry B1 | R.

# Vertebrates/chordates have the following possible codes: 
# A1, A2, B (=B1|B2|B3), C(=C1|C2), D, E, F, G, H, I, J=(J1|J2), K, L, M, N,
#   O(=O1|O2), P, Q and possibly R.

import os
import re
import sys

__authors__ = 'Jessen V. Bredeson, Daniel S. Rokhsar'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Assign ancestral Chordate LG identity'


from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.newick import IntervalNewickTree
from og.constants import (
    _RPAREN,
    _LPAREN,
    _COMMA,
    _EMPTY,
    _PIPE
)

_ambiguous = set(('?',))
_specific_codes = re.compile('((?:A[12]?)|(?:B[123]?)|(?:C[12]?)|D|E|[Ff]|G|H|I|(?:J[12]?)|K|L|M|N|(?:O[12]?)|P|Q|R)')
_alpha_numeric  = re.compile('([A-Za-z])(\d{0,1})')
_alpha_only = re.compile('[A-Za-z]')
_homoeolog_notation = re.compile('\[[12][ab]\]')
_translocation_chars = re.compile('[x\*\+]')
_strip_digits = _alpha_only.findall

num = len

_CODE_EXPANSION_MAP = {
    'A' : set(('A1','A2')),
    'A1': set(('A1',)),
    'A2': set(('A2',)),
    'B' : set(('B1','B2','B3')),
    'B1': set(('B1',)),
    'B2': set(('B2',)),
    'B3': set(('B3',)),
    'C' : set(('C1','C2')),
    'C1': set(('C1',)),
    'C2': set(('C2',)),
    'D' : set(('D',)),
    'E' : set(('E',)),
    'F' : set(('F',)),
    'f' : set(('F',)),
    'G' : set(('G',)),
    'H' : set(('H',)),
    'I' : set(('I',)),
    'J' : set(('J1','J2')),
    'J1': set(('J1',)),
    'J2': set(('J2',)),
    'K' : set(('K',)),
    'L' : set(('L',)),
    'M' : set(('M',)),
    'N' : set(('N',)),
    'O' : set(('O1','O2')),
    'O1': set(('O1',)),
    'O2': set(('O2',)),
    'P' : set(('P',)),
    'Q' : set(('Q',)),
    'R' : set(('R',)),
}

_SEQUENCE_CODE_MAP = {
    'Bfl': {
        'Bfl01':'A',
        'Bfl10':'B1',
        'Bfl16':'B2',
        'Bfl18':'B3',
        'Bfl06':'D',
        'Bfl05':'E',
        'Bfl07':'F',
        'Bfl11':'G',
        'Bfl13':'H',
        'Bfl04':'I*O1',
        'Bfl17':'J1',
        'Bfl02':'J2*C1',
        'Bfl09':'K',
        'Bfl15':'L',
        'Bfl08':'M',
        'Bfl12':'N',
        'Bfl19':'O2',
        'Bfl14':'P',
        'Bfl03':'Q*C2',
        'BflSca1189':'(G)',
        'BflSca168':'(G)',
        'BflSca1511':'(K)',
    },
    'Bla': {
        'Bla1':'G',  # and K
        'Bla2':'C1', # and J2
        'Bla3':'O1', # and I
        'Bla4':'A2', # and A1
        'Bla5':'E',
        'Bla6':'F',
        'Bla7':'D',
        'Bla8':'H',
        'Bla9':'M',
        'Bla10':'P',
        'Bla11':'L',
        'Bla12':'B2',
        'Bla13':'B1',
        'Bla14':'Q',
        'Bla15':'B3',
        'Bla16':'O2',
        'Bla17':'N',
        'Bla18':'J1',
        'Bla19':'C2',
    },
    'Llo': {
        'Llo3':'A1',
        'Llo18':'A2',
        'Llo16':'B1',
        'Llo13':'B2',
        'Llo15':'B3',
        'Llo2':'C1*G',
        'Llo14':'C2',
        'Llo11':'D',
        'Llo7':'E',
        'Llo9':'F',
        # '':'G',  # Llo2
        'Llo1':'H*Q',
        'Llo6':'I',
        'Llo19':'J1',
        'Llo4':'J2*L',
        'Llo5':'K*O2',
        # '':'L',  # Llo4
        'Llo17':'M',
        'Llo10':'N',
        'Llo8':'O1*R',
        # '':'O2',  # Llo5
        'Llo12':'P',
        # '':'Q',  # Llo1
        # '':'R'  # Llo8
    },
    'Lva': {
        'Lva7':'A1',
        'Lva8':'A2',
        'Lva13':'B1',
        'Lva3':'B2xExC2',
        'Lva5':'C1',
        'Lva2':'D*G',
        'Lva1':'F*(B3xJ1)',
        'Lva10':'H',
        'Lva12':'I',
        'Lva17':'J2',
        'Lva4':'K',
        'Lva11':'L',
        'Lva15':'M',
        'Lva16':'N',
        'Lva9':'O1',
        'Lva19':'O2',
        'Lva14':'P',
        'Lva18':'Q',
        'Lva6':'R',
    },
    'Ofu': {
        'Ofu8':'A1',
        'Ofu10':'A2*C1',
        'Ofu7':'B1*B3*C2',
        'Ofu1':'B2*J1*J2*L*P',
        # '':'B3',  # Ofu7
        # '':'C1',
        # '':'C2',
        'Ofu3':'D',
        'Ofu5':'E*I',
        'Ofu2':'F',
        'Ofu12':'G',
        'Ofu9':'H*Q',
        # '':'I',
        # '':'J1',
        # '':'J2',
        'Ofu11':'K*O2',
        # '':'L',
        'Ofu6':'M*O1*R',
        'Ofu4':'N',
        # '':'O1',
        # '':'O2',
        # '':'P',
        # '':'Q',
        # '':'R',
    },
    'Pma': {
        'Pma3':'A1',
        'Pma13':'A2',
        'Pma15':'B1',
        'Pma2':'B2*M',
        'Pma19':'B3',
        'Pma9':'C1',
        'Pma17':'C2',
        'Pma11':'D',
        'Pma7':'E',
        'Pma6':'F',
        'Pma8':'G',
        'Pma1':'H*Q',
        'Pma12':'I',
        'Pma18':'J1',
        'Pma4':'J2*L',
        'Pma5':'K*O2',
        # 'Pma4':'L',
        # 'Pma2':'M',
        'Pma14':'N',
        'Pma10':'O1*R',
        # 'Pma5':'O2',
        'Pma16':'P',
        # 'Pma1':'Q',
        # 'Pma10':'R',
    },
    'Sca': {
        'Sca01':'A1',
        'Sca14':'A2',
        'Sca05':'B1xR',
        'Sca11':'(B2xC2)1',  # and 'D'
        'Sca23':'(B2xC2)2',
        'Sca21':'B3',
        'Sca10':'C1',
        'Sca13':'D',
        'Sca03':'E',
        'Sca09':'F',
        'Sca06':'G',
        'Sca02':'H',
        'Sca07':'I',
        'Sca15':'J1',
        'Sca12':'J2',
        'Sca04':'K',
        'Sca08':'L',
        'Sca17':'M',
        'Sca20':'N',
        'Sca19':'O1',
        'Sca16':'O2',
        'Sca22':'P',
        'Sca18':'Q',
    }
}


def parseLGs(LGs, specific=False):
    _codes = _specific_codes if specific else _alpha_only
    _LGs = set()
    for LG in LGs:
        LG = _homoeolog_notation.sub(_EMPTY, LG)
        LG = _translocation_chars.sub(_PIPE, LG)
        for code in _codes.findall(LG):
            try:
                _LGs.update(_CODE_EXPANSION_MAP[code])
            except KeyError:
                raise KeyError("Invalid code `%s` in `%s`" % (code, LG)) from None
    return _LGs



def inferLGlogically(sample_list, group, specific=False):
    sample_indices = {s.id:s.index for s in sample_list}
    
    bla_code = parseLGs(group[sample_indices['Bla_LG']], specific)
    lva_code = parseLGs(group[sample_indices['Lva_LG']], specific)
    sca_code = parseLGs(group[sample_indices['Sca_LG']], specific)

    llo_code = parseLGs(group[sample_indices['Llo_LG']], specific)
    pma_code = parseLGs(group[sample_indices['Pma_LG']], specific)
    ofu_code = parseLGs(group[sample_indices['Ofu_LG']], specific)

    bla_intersects_lva = bla_code & lva_code
    bla_intersects_sca = bla_code & sca_code
    bla_intersects_llo = bla_code & llo_code
    bla_intersects_pma = bla_code & pma_code
    bla_intersects_ofu = bla_code & ofu_code

    lva_intersects_llo = lva_code & llo_code
    lva_intersects_pma = lva_code & pma_code
    lva_intersects_ofu = lva_code & ofu_code

    sca_intersects_llo = sca_code & llo_code
    sca_intersects_pma = sca_code & pma_code
    sca_intersects_ofu = sca_code & ofu_code


    comments = []    
    consensus = set()
    if not consensus and bla_intersects_lva:
        consensus = bla_intersects_lva
    if not consensus and bla_intersects_sca:
        consensus = bla_intersects_sca
    if not consensus and bla_intersects_llo:
        consensus = bla_intersects_llo
    if not consensus and bla_intersects_pma:
        consensus = bla_intersects_pma
    if not consensus and bla_intersects_ofu:
        consensus = bla_intersects_ofu
        
    if not consensus and lva_intersects_llo:
        consensus = lva_intersects_llo
    if not consensus and lva_intersects_pma:
        consensus = lva_intersects_pma
    if not consensus and lva_intersects_ofu:
        consensus = lva_intersects_ofu

    if not consensus and sca_intersects_llo:
        consensus = sca_intersects_llo
    if not consensus and sca_intersects_pma:
        consensus = sca_intersects_pma
    if not consensus and sca_intersects_ofu:
        consensus = sca_intersects_ofu

    if num(consensus) == 0:
        comments.append('no-consensus')
    else:
        if num(consensus) > 1:
            comments.append('multi-consensus')
        num_intersects = 0
        num_present = 0
        if bla_code:
            if (bla_code & consensus) == consensus:
                num_intersects += 1
            else:
                comments.append('Bla-translocation')
            num_present += 1            
        if lva_code:
            if (lva_code & consensus) == consensus:
                num_intersects += 1
            else:
                comments.append('Lva-translocation')
            num_present += 1            
        if sca_code:
            if (sca_code & consensus) == consensus:
                num_intersects += 1
            else:
                comments.append('Sca-translocation')
            num_present += 1
        if llo_code:
            if (llo_code & consensus) == consensus:
                num_intersects += 1
            else:
                comments.append('Llo-translocation')
            num_present += 1
        if pma_code:
            if (pma_code & consensus) == consensus:
                num_intersects += 1
            else:
                comments.append('Pma-translocation')
            num_present += 1
        if ofu_code:
            if (ofu_code & consensus) == consensus:
                num_intersects += 1
            else:
                comments.append('Ofu-translocation')
            num_present += 1
            
        if num_intersects == num_present:
            comments = ['complete-consensus'] + comments
        else:
            comments = ['partial-consensus'] + comments
        
    return consensus, comments



def collapseLGs(codes):
    _codes = set()
    if codes is None or \
       codes is _ambiguous:
        return _codes
    for code in codes:
        match = _alpha_numeric.match(code)
        if match:
            alpha = match.group(1)
            if ((num(_CODE_EXPANSION_MAP[alpha]) > 1) and \
                (num(_CODE_EXPANSION_MAP[alpha] & codes) > 1)):
                _codes.add(alpha)
            else:
                _codes.add(code)
        else:
            raise AssertionError(code)
    return _codes
                


def append_samples_LG_to_header(ortho_samples, sample_list)
    for sample in sample_list:
        ortho_samples.append(sample.id + '_LG')


def append_samples_LG_to_groups(ortho_samples, sample_list, group):
    sample_indices = {s.id:s.index for s in ortho_samples}
    while len(group) < len(ortho_samples):
        group.append(set())
        
    for sample in sample_list:
        for member in group[sample.index]:
            try:
                group[sample_indices[sample.id + '_LG']].add(
                    _SEQUENCE_CODE_MAP[sample.id][member]
                )
            except KeyError:
                pass
    

def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <in.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    
        
def RULE_BASED(argv):
    if len(argv) != 1:
        usage('Unexpected number of arguments')
    
    ortho = OrthoFinderOrthogroups(argv[0])

    sample_list = list(ortho.samples)
    append_samples_LG_to_header(ortho.samples, sample_list)
    ortho.samples.append('Cns_LG')
    ortho.samples.append('Comments')
    
    for group in ortho.groups:
        append_samples_LG_to_groups(ortho.samples, sample_list, group)
        
        codes, comments = inferLGlogically(ortho.samples, group, specific=True)
        # if codes is _ambiguous:
        #   codes, comments = inferLGlogically(ortho.groups[g], samples_index, specific=False)
           
        codes = sorted(collapseLGs(codes))
        codes = codes or _ambiguous
            
        group.append((_PIPE.join(codes),))
        group.append((_COMMA.join(comments),))

    ortho.to_file(sys.stdout)



def inferLGphylogenetically(sample_indices, group, tree, include=set(), specific=False):
    if num(include) == 0:
        include = set(map(lambda node: node.id, tree.terminal_nodes))
    comments = []
    codes = [None] * num(tree.nodes)
    for dsc_index in tree.terminal_indices:
        dsc_node = tree.nodes[dsc_index]
        codes[dsc_index] = parseLGs(group[sample_indices[dsc_node.id]], specific)
        
    for anc_index, dsc_index in tree.get_edges(indices=True, reverse=True):
        if anc_index is None:
            break

        anc_node = tree.nodes[anc_index]
        dsc_node = tree.nodes[dsc_index]

        if dsc_node.id in include:
            include.add(anc_node.id)
        else:
            continue
        if codes[anc_index] is None:
            codes[anc_index] = set()
        if codes[anc_index].intersection(codes[dsc_index]):
            codes[anc_index].intersection_update(codes[dsc_index])
        else:
            codes[anc_index].update(codes[dsc_index])
                            
    cns_index = 0
    consensus = collapseLGs(codes[cns_index])
    if num(consensus) == 1:
        # cns_index = root_index 
        # iter down tree and note any translocations
        for dsc_index in tree.terminal_indices:
            if num(codes[dsc_index]) > 0 and \
               num(codes[dsc_index].intersection(codes[cns_index])) < 1:
                comments.append('%s-translocation' % str(tree.nodes[dsc_index].id))
    # elif num(codes[crown_index]) == 1:
    #     cns_index = crown_index
    #     comments.append('ingroup-crown-consensus')
    else:
        consensus = _ambiguous

    return consensus, comments


    
def PHYLOGENY_BASED(argv):
    print('[WARNING]: Running experimental, phylogenetic-based consensus procedure', file=sys.stderr)
    ortho = OrthoFinderOrthogroups(argv[0])
    tree = IntervalNewickTree(argv[1])

    sample_list = list(ortho.samples)
    
    outgroups = set()
    for anc_index, dsc_index in tree.get_edges(indices=True, reverse=True):
        if anc_index is None:
            break
        
        anc_node = tree.nodes[anc_index]
        dsc_node = tree.nodes[dsc_index]
        if str(dsc_node.id).endswith('*'):
            dsc_node.id = dsc_node.id.rstrip('*')
            outgroups.add(dsc_node.id)
        if dsc_node.id in outgroups:
            outgroups.add(anc_node.id)
        else:
            crown_index = anc_index
    
    append_samples_LG_to_header(ortho.samples, sample_list)
    ortho.samples.append('Cns_LG')
    ortho.samples.append('Comments')
    
    for group in ortho.groups:
        append_samples_LG_to_groups(ortho.samples, sample_list, group)

        consensus, comments = inferLGphylogenetically(sample_indices, group, tree, specific=True)
        if consensus is _ambiguous:
            consensus, comments = inferLGphylogenetically(sample_indices, group, tree, specific=False)
            
        group.append(('|'.join(sorted(consensus)),))
        group.append((';'.join(comments),))

    ortho.to_file(sys.stdout)



main = RULE_BASED  # PHYLOGENY_BASED  #RULE_BASED


main(sys.argv[1:])
