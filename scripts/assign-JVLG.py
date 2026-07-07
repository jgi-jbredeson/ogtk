#!/usr/bin/env python

import os
import sys
from og.core.orthogroups.parsers import OrthoFinderOrthogroups

# with edits from Dan to correct some errors, and to add 'Group 1or2' and 'Group 13or14' for two ambiguous cases.

__authors__ = 'Jessen V. Bredeson, Daniel S. Rokhsar'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Assign ancestral jawed-vertebrate LG identity'

# Vertebrates/chordates have the following possible codes: 
# A1, A2, B (=B1|B2|B3), C(=C1|C2), D, E, F, G, H, I, J=(J1|J2), K, L, M, N, O(=O1|O2), P, Q and possibly R.
# drop \[[12][ab]\] in verts, but not inverts


num = len
AncLG = dict()

# Group 1 (FIQ2a)
# (('Cmi01' OR ('Cpl32' OR 'Cpl01' OR 'Cpl02') OR ('Ler01 OR Ler03')) AND ('Gga04'))
AncLG[1] = (
    'FIQ[2a]',
    (
        set(('Cmi01','Cpl32', 'Cpl01', 'Cpl02', 'Ler01','Ler03')),
        set(('Gga04',)),
        set(('GgaZ',))
    )
)

# Group 2 (CL2a)
# (('Cmi01' OR ('Cpl32' OR 'Cpl01' OR 'Cpl02') OR ('Ler01 OR Ler03')) AND ('GgaZ'))
AncLG[2] = (
    'CL[2a]',
    (
        set(('Cmi01','Cpl32','Cpl01','Cpl02','Ler01','Ler03')),
        set(('GgaZ',)),
        set(('Gga04',))
    )
)

# Group 100 (accounts for ambiguity of genes in Group 1 or Group 2)
AncLG[100] = (
    'FIQ[2a]|CL[2a]',
    (
        set(('Cmi01','Cpl32','Cpl01','Cpl02','Ler01','Ler03')),
        set(('Loc02','Loc04','Xtr01')),
        set(('GgaZ','Gga04'))
    )
)

# Group 3 (DBJ1a)
# (('Cmi02' OR ('Cpl04' OR 'Cpl05' OR 'Cpl29') OR ('Ler02' OR 'Ler04')) AND (('Loc09' OR 'Loc11') OR 'Gga02' OR 'Xtr06')
AncLG[3] = (
    'DBJ[1a]',
    (
        set(('Cmi02','Cpl04','Cpl05','Cpl29','Ler02','Ler04')),
        set(('Loc09','Loc11','Gga02','Xtr06')),
    )
)

# Group 4 (A1KJ2a)
# ('Cmi03' OR ('Cpl03' OR 'Cpl09') OR ('Ler05' OR 'Ler08')) AND (('Loc01' OR 'Loc16') OR 'Gga03' OR 'Xtr05')
AncLG[4] = (
    'A1KJ[2a]',
    (
        set(('Cmi03','Cpl03','Cpl09','Ler05','Ler08')),
        set(('Loc01','Loc16','Gga03','Xtr05')),
    )
)
    
# Group 5 (FKN1a plus A1A2?)
# ('Cmi04' OR ('Cpl06' OR 'Cpl12') OR ('Ler06' OR 'Ler45' OR 'Ler13')) AND (('Loc03' OR 'Loc17' OR 'Loc14') OR 'Gga01' OR 'Xtr02')
AncLG[5] = (
    'FKN[1a]+A2',
    (
        set(('Cmi04','Cpl06','Cpl12','Ler06','Ler45','Ler13')),
        set(('Loc03','Loc17','Loc14','Gga01','Xtr02')),
    )
)
    
# Group 6 LM1a
# ('Cmi05' OR 'Cpl11' OR 'Ler10') AND ('Loc10' OR 'Gga08' OR 'Xtr04')
AncLG[6] = (
    'LM[1a]',
    (
        set(('Cmi05','Cpl11','Ler10')),
        set(('Loc10','Gga08','Xtr04')),
    )
)
    
# Group 7 (EO2a)
# ('Cmi05' OR ('Cpl19' OR 'Cpl23') OR ('Ler19' OR 'Ler22')) AND ('Loc08' OR 'Gga01' OR 'Xtr03')
AncLG[7] = (
    'EO[2a]',
    (
        set(('Cmi05','Cpl19','Cpl23','Ler19','Ler22')),
        set(('Loc08','Gga01','Xtr03')),
    )
)
    
# Group 8 (B2a)
# ('Cmi06' OR 'Cpl07' OR 'Ler07') AND ('Loc12' OR 'Gga07' OR 'Xtr09')
AncLG[8] = (
    'B[2a]',
    (
        set(('Cmi06','Cpl07','Ler07')),
        set(('Loc12','Gga07','Xtr09')),
    )
)
    
# Group 9 (IQ1a)
# ('Cmi07' OR ('Cpl22' OR 'Cpl38') OR ('Ler15' OR 'Ler34')) AND ('Loc05' OR 'Gga06' OR 'Xtr07')
AncLG[9] = (
    'IQ[1a]',
    (
        set(('Cmi07','Cpl22','Cpl38','Ler15','Ler34')),
        set(('Loc05','Gga06','Xtr07')),
    )
)
    
# Group 10 (G1a)
# ('Cmi08' OR 'Cpl25' OR 'Ler25') AND ('Loc20' OR 'Gga15' OR 'Xtr01')
AncLG[10] = (
    'G[1a]',
    (
        set(('Cmi08','Cpl25','Ler25')),
        set(('Loc20','Gga15','Xtr01')),
    )
)
    
# Group 11 (M2a)
# ('Cmi08' OR 'Cpl30' OR 'Ler31') AND ('Loc21' OR 'Gga17' OR 'Xtr08')
AncLG[11] = (
    'M[2a]',
    (
        set(('Cmi08','Cpl30','Ler31')),
        set(('Loc21','Gga17','Xtr08')),
    )
)
    
# Group 12 (A1J2b)
# ('Cmi09' OR 'Cpl10' OR 'Cpl12' OR 'Ler09') AND ('Loc07' OR 'Gga05' OR 'Xtr08')
AncLG[12] = (
    'A1',
    (
        set(('Cmi09','Cpl10','Cpl2','Ler09')),
        set(('Loc07','Gga05','Xtr08')),
    )
)
    
# Group 13 (JK2b)
# ('Cmi10' OR 'Cpl27' OR 'Ler26') AND ('Loc06' OR 'Gga23')
AncLG[13] = (
    'JK[2b]',
    (
        set(('Cmi10','Cpl27','Ler26')),
        set(('Loc06','Gga23')),
    ),
    (
        set(('Cpl27','Ler26')),
        set(('Xtr02',)),
    )
)

# Group 13.1 (JK2b)
# ('Cpl27' OR 'Ler26') AND ('Xtr02')
# AncLG[13.1] = (
#     'JK[2b]',
#     (
#         set(('Cpl27','Ler26')),
#         set(('Xtr02',)),
#     )
# )
    
# Group 14 (G2a)
# ('Cmi10' OR 'Cpl28' OR 'Ler28') AND ('Loc22' OR 'Gga19')
AncLG[14] = (
    'G[2a]',
    (
        set(('Cmi10','Cpl28','Ler28')),
        set(('Loc22','Gga19')),
    ),
    (
        set(('Cpl28','Ler28')),
        set(('Xtr02',)),
    )
)

# Group 14.1 (G2a)
# ('Cpl28' OR 'Ler28') AND ('Loc22' OR 'Gga19')
# AncLG[14.1] = (
#     'G[2a]',
#     (
#         set(('Cpl28','Ler28')),
#         set(('Xtr02',)),
#     )
# )
    
# Group 15 (NP2a)
# ('Cmi11' OR 'Cpl13' OR 'Ler14') AND ('Loc14' OR 'Gga09' OR 'Xtr05')
AncLG[15] = (
    'NP[2a]',
    (
        set(('Cmi11','Cpl13','Ler14')),
        set(('Loc14','Gga09','Xtr05')),
    )
)
    
# Group 16 (E1a)
# ('Cmi12' OR 'Cpl18' OR 'Ler16') AND ('Loc05' OR 'Gga12' OR 'Xtr04')
AncLG[16] = (
    'E[1a]',
    (
        set(('Cmi12','Cpl18','Ler16')),
        set(('Loc05','Gga12','Xtr04')),
    )
)
    
# Group 17 (FIQ2b +?)
# ('Cmi13' OR 'Cpl14' OR 'Ler11') AND ('Loc06' OR 'Gga13' OR 'Xtr03')
AncLG[17] = (
    'FIQ[2b]+?',
    (
        set(('Cmi13','Cpl14','Ler11')),
        set(('Loc06','Gga13','Xtr03')),
    )
)
    
# Group 18 (D2a)
# ('Cmi14' OR 'Cpl17' OR 'Ler17') AND ('Loc23' OR 'Gga11' OR 'Xtr04')
AncLG[18] = (
    'D[2a]',
    (
        set(('Cmi14','Cpl17','Ler17')),
        set(('Loc23','Gga11','Xtr04')),
    )
)
    
# Group 19 (FKN1b)
# ('Cmi15' OR 'Cpl15' OR 'Ler12') AND ('Loc07' OR 'Gga04' OR 'Xtr08')
AncLG[19] = (
    'FKN[1b]',
    (
        set(('Cmi15','Cpl15','Ler12')),
        set(('Loc07','Gga04','Xtr08')),
    )
)
    
# Group 20 (Hf1a)
# ('Cmi16' OR 'Cpl21' OR 'Ler20') AND ('Loc13' OR 'Gga14' OR 'Xtr09')
AncLG[20] = (
    'Hf[1a]',
    (
        set(('Cmi16','Cpl21','Ler20')),
        set(('Loc13','Gga14','Xtr09')),
    )
)
    
# Group 21 (O1a)
# ('Cmi17' OR 'Cpl16' OR 'Ler18') AND ('Loc27' OR 'Gga05' OR 'Xtr04')
AncLG[21] = (
    'O[1a]',
    (
        set(('Cmi17','Cpl16','Ler18')),
        set(('Loc27','Gga05','Xtr04')),
    )
)
    
# Group 22 (DJ1b)
# ('Cmi18' OR 'Cpl20' OR 'Ler21') AND ('Loc18' OR 'Gga20' OR 'Xtr10')
AncLG[22] = (
    'DJ[1b]',
    (
        set(('Cmi18','Cpl20','Ler21')),
        set(('Loc18','Gga20','Xtr10')),
    )
)
    
# Group 23 (EO2b)
# ('Cmi19' OR 'Cpl26' OR 'Ler24') AND ('Loc03' OR 'Gga26' OR 'Xtr02')
AncLG[23] = (
    'EO[2b]',
    (
        set(('Cmi19','Cpl26','Ler24')),
        set(('Loc03','Gga26','Xtr02')),
    )
)
    
# Group 24 (Hf2a)
# ('Cmi20' OR 'Cpl24' OR 'Ler23') AND ('Loc10' OR 'Gga18' OR 'Xtr10')
AncLG[24] = (
    'Hf[2a]',
    (
        set(('Cmi20','Cpl24','Ler23')),
        set(('Loc10','Gga18','Xtr10')),
    )
)
    
# Group 25 (CL2b)
# ('Cmi21' OR 'Cpl31' OR 'Ler29') AND ('Loc19' OR 'Gga28' OR 'Xtr01')
AncLG[25] = (
    'CL[2b]',
    (
        set(('Cmi21','Cpl31','Ler29')),
        set(('Loc19','Gga28','Xtr01')),
    )
)
    
# Group 26 (C1a)
# (('Cmi22' OR 'Cmi26') OR ('Cpl36' OR 'Cpl40' OR 'Cpl32') OR ('Ler33' OR 'Ler36')) AND ('Loc03' OR 'Gga10' OR 'Xtr03')
AncLG[26] = (
    'C[1a]',
    (
        set(('Cmi22','Cmi26','Cpl36','Cpl40','Cpl32','Ler33','Ler36')),
        set(('Loc03','Gga10','Xtr03')),
    )
)
    
# Group 27 (B2b)
# ('Cmi23' OR 'Cpl33' OR 'Ler27') AND ('Loc15' OR 'Gga27' OR 'Xtr10')
AncLG[27] = (
    'B[2b]',
    (
        set(('Cmi23','Cpl33','Ler27')),
        set(('Loc15','Gga27','Xtr10')),
    )
)
    
# Group 28 (P1a)
# ('Cmi24' OR 'Cpl34' OR 'Ler30') AND ('Loc25' OR 'Gga21' OR 'Xtr07')
AncLG[28] = (
    'P[1a]',
    (
        set(('Cmi24','Cpl34','Ler30')),
        set(('Loc25','Gga21','Xtr07')),
    )
)
    
# Group 29 (A2R)
# ('Cmi25' OR ('Cpl35' OR 'Cpl36') OR 'Ler32') AND ('Loc26' OR 'Gga24' OR 'Xtr07')
AncLG[29] = (
    'A2R',
    (
        set(('Cmi25','Cpl35','Cpl36','Ler32')),
        set(('Loc26','Gga24','Xtr07')),
    )
)

# Group 30 (B1b)
# ('Cmi27' OR 'Cpl43' OR 'Ler40') AND ('Loc04' OR 'Gga34' OR 'Xtr02')
AncLG[30] = (
    'B[1b]',
    (
        set(('Cmi27','Cpl43','Ler40')),
        set(('Loc04','Gga34','Xtr02')),
    )
)
    
# Group 31 (IQ1b)
# ('Cmi28' OR 'Cpl42' OR 'Ler35') AND ('Loc01' OR 'Gga22' OR 'Xtr03')
AncLG[31] = (
    'IQ[1b]',
    (
        set(('Cmi28','Cpl42','Ler35')),
        set(('Loc01','Gga22','Xtr03')),
    )
)
    
# Group 32 (LM1b)
# ('Cmi29' OR 'Cpl08' OR 'Ler39') AND ('Loc06' OR 'Gga30' OR 'Xtr03')
AncLG[32] = (
    'LM[1b]',
    (
        set(('Cmi29','Cpl08','Ler39')),
        set(('Loc06','Gga30','Xtr03')),
    )
)
    
# Group 33 (HF2b)
# ('Cmi31' OR 'Cpl39' OR 'Ler37') AND ('Loc12' OR 'Gga01' OR 'Xtr04')
AncLG[33] = (
    'HF[2b]',
    (
        set(('Cmi31','Cpl39','Ler37')),
        set(('Loc12','Gga01','Xtr04')),
    )
)
    
# Group 34 (O1b)
# ('Cmi31' OR 'Cpl39' OR 'Ler37') AND ('Gga31' OR 'Xtr07')
AncLG[34] = (
    'O[1b]',
    (
        set(('Cmi31','Cpl39','Ler37')),
        set(('Gga31','Xtr07')),
    )
)
    
# Group 35 (C1b)
# ('Cpl51' OR ('LerSca69' OR 'Ler24')) AND ('Loc24' OR 'Gga25' OR 'Xtr08')
AncLG[35] = (
    'C[1b]',
    (
        set(('Cpl51','Ler24','LerSca69')),
        set(('Loc24','Gga25','Xtr08')),
    )
)
    
# Group 36 (A2?)
# ('Cpl47') AND ('Loc24' OR  'Xtr07')
AncLG[36] = (
    'A2?',
    (
        set(('Cpl47',)),
        set(('Loc24','Xtr07')),
    )
)
    
# Group 37 (M2b)
# ('Cmi30' OR 'Cpl37' OR ('LerSca57' OR LerSca65')) AND ('Gga16' OR 'Xtr08')
AncLG[37] = (
    'M[2b]',
    (
        set(('Cmi30','Cpl37','LerSca57','LerSca65')),
        set(('Gga16','Xtr08')),
    )
)

# Group 38 (A12b)
# (('CmiSca45' OR 'CmiSca68') OR 'Cpl41' OR ('Ler38')) AND ('Loc02' OR 'Gga38' OR 'Xtr08'))
AncLG[38] = (
    'A1[2b]',
    (
        set(('CmiSca45','CmiSca68','Cpl41','Ler38')),
        set(('Loc02','Gga38','Xtr08')),
    )
)
    
# Group 39 (NP2b)
# (('CmiSca41' OR 'CmiSca60') OR ('Cpl48' OR 'Cpl43') OR ('LerSca58')) AND ('Loc02' OR 'Xtr03')
AncLG[39] = (
    'NP[2b]',
    (
        set(('CmiSca41','CmiSca60','Cpl48','Cpl43','LerSca58')),
        set(('Loc02','Xtr03')),
    )
)
    
# Group 40 (A1)
# ('CmiSca59' OR 'Cpl45' OR 'Ler44') AND ('Loc28' OR 'Gga39' OR 'Xtr04')
AncLG[40] = (
    'A1',
    (
        set(('CmiSca59','Cpl45','Ler44')),
        set(('Loc28','Gga39','Xtr04')),
    )
)

# Group 41 E(1beta) **** NEW ************
# This may require some tweaking, as most cartilaginous fish assemblies are a collection of
# scaffolds "SpeSca" when collapsed, and not just one scaffold...
AncLG[41] = (
    'E[1b]',
    (
        set(('Cpl39','LerSca530','LerSca1269','LerSca264')),
        set(('Loc01','Xtr08'))
    )
)

# AncLG[42] = (
#     'G[b]',
#     (
#         set(()),
#         set(('Gga05','Xtr09'))
#     )


def _pre(n):
    return 'pre%d' % n


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write('\n')
    stream.write('Program: %s (%s)\n' % (__program__, __purpose__))
    stream.write('Version: %s %s\n' % (__pkgname__, __version__))
    stream.write('Contact: %s\n' % __contact__)
    stream.write('\n')
    stream.write('Usage:   %s [options] <in.tsv> <out-prefix>\n' % __program__)
    stream.write('\n')
    stream.write('\n%s' % message)
    sys.exit(exitcode)
    

def main(argv):
    if len(argv) != 2:
        usage('Unexpected number of arguments')
        
    ortho = OrthoFinderOrthogroups(argv[0])
    prefix = argv[1]

    samples = [s.id for s in ortho.samples]
    
    lab = (0, 1, 'M')
    OGs = (
        OrthoFinderOrthogroups(samples=(samples + ['JV_N','JV_LG'])),
        OrthoFinderOrthogroups(samples=(samples + ['JV_N','JV_LG'])),
        OrthoFinderOrthogroups(samples=(samples + ['JV_N','JV_LG']))
    )
    
    for group in ortho.groups:
        group_set = set()
        for sample in ortho.samples:
            group_set.update(group[sample.index])

        membership_letters = []
        membership_numbers = []
        for lg in AncLG:
            for conditions in AncLG[lg][1:]:
                if ((group_set & conditions[0]) and (group_set & conditions[1])):
                    if (num(conditions) > 2 and (group_set & conditions[2])):
                        continue
                    membership_letters.append(AncLG[lg][0])
                    membership_numbers.append(lg)
                    break

        og = OGs[min(2, num(membership_numbers))]
        
        if og == 0:
            membership_letters.append('NONE')
            membership_numbers.append(0)

        assigned = og.new_group(append=True)
        assigned.id = group.id
        for sample in ortho.samples:
            assigned[sample.index] = group[sample.index]
        assigned[-2] = ('|'.join(map(str, membership_numbers)),)
        assigned[-1] = ('|'.join(membership_letters),)

    for n, og in enumerate(OGs):
        og.to_file('%s.%s.tsv' % (prefix, lab[n]))

        
if __name__ == '__main__':
    main(sys.argv[1:])
