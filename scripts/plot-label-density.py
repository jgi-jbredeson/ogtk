#!/usr/bin/env python

import os
import sys
import math
import numpy
import bisect
import pandas
import getopt
import matplotlib.pyplot as plotter

from og.constants import dict
from og.core.utils import index_list
from matplotlib.colors import to_hex
from matplotlib import colormaps as cmaps

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Plot label densities along the genome'

num = len

_bed_field_names = ['chr','beg','end','name','score','strand']
_bed_field_types = [str, int, int, str, float, str]

_map_field_names = ['name','label']
_map_field_types = [str, str]

_valid_output_types = ('pdf','png','ps','eps','svg')
_valid_coordinate_systems = ('genomic','ordinal')

_palette_priority = ['tab10','Set3','Dark2','Pastel1','tab20b','tab20c',
                     'Accent','Set2','Set1','Pastel2','tab20','Paired']
_default_colors = []
for palette in _palette_priority:
    _default_colors.extend(reversed(cmaps[palette].colors))


    
def get_num_bins(num_data, bin_width=1, bin_shift=1):
    return max(1, math.ceil(float(num_data) / bin_shift))



def get_freqs(data, labels, bin_width=1, bin_shift=1):
    """
    The data object is a pandas.Series and the labels object is a dictionary
    containing labels as keys and indices as values
    """
    num_bins = get_num_bins(num(data), bin_width, bin_shift)
    label_freq = numpy.zeros((num_bins, num(labels)), dtype=numpy.float64)
    label_max = numpy.zeros(num_bins)
    idx_begs = []
    idx_ends = []
    for i in range(num_bins):
        bin_beg = bin_shift * i
        bin_end = bin_shift * i + bin_width

        bin_loci = data.iloc[bin_beg:bin_end]
            
        for label in labels:
            label_freq[i, labels[label]] = num(bin_loci[bin_loci == label])
            
        label_freq[i,:] = label_freq[i,:] / max(1, label_freq[i,:].sum())

        idx_begs.append(bin_loci.index[0])
        idx_ends.append(bin_loci.index[-1])

    return label_freq, numpy.array(idx_begs), numpy.array(idx_ends)



def get_blocks(data, labels, bin_width=1, bin_shift=1):
    num_bins = get_num_bins(num(data), bin_width, bin_shift)
    bin_freq, _, _ = get_freqs(data, labels, bin_width, bin_shift)
    bin_max = bin_freq.argmax(1)

    blocks = []
    block_beg = 0
    block_end = 0
    min_delta = 1.0
    for j in range(1, num_bins):
        i = max(0, j - bin_width // bin_shift)
        if bin_max[i] != bin_max[j]:
            delta = abs(bin_freq[i, bin_max[i]] - bin_freq[j, bin_max[j]])
            if delta < min_delta:
                min_delta = delta
                block_end = j
        elif block_end:
            blocks.append((data.index[block_beg], data.index[block_end]))
            block_beg = block_end
            block_end = 0
            min_delta = 1.0

    blocks.append((data.index[block_beg], data.index[-1] + 1))
        
    return blocks



def infer_output_type(output_file, output_type=None):
    if output_type is None:
        output_type = output_file.rsplit('.', maxsplit=1)
        if len(output_type) > 1:
            output_type = output_type[-1].lower()
        else:
            output_type = _valid_output_types[0]
    if output_type not in _valid_output_types:
        usage('Unsupported image type: `%s`' % output_type)
    return output_type



def plot_color_legend(lg_colors, file=None):
    if file:
        legend_type = infer_output_type(file)
    else:
        legend_type = _valid_output_types[0]
        
    num_colors = num(lg_colors)
    num_side = math.ceil(math.sqrt(num_colors))
    fig, ax = plotter.subplots(num_side, num_side, figsize=(num_side, num_side))
    fig.tight_layout()
    i = 0
    j = 0
    k = 0
    for lg in lg_colors:
        if num_side <= i:
            j += 1
            i = 0
        ax[i][j].bar([0], [1], 1, color=lg_colors[lg], align='edge')
        ax[i][j].set_xbound(0, 1)
        ax[i][j].set_ybound(0, 1)
        ax[i][j].xaxis.set_tick_params(labelbottom=False)
        ax[i][j].yaxis.set_tick_params(labelbottom=False)
        ax[i][j].set_xticks([])
        ax[i][j].set_yticks([])
        ax[i][j].set_title(lg, loc='left')
        i += 1
        k += 1

    for n in range(k, num_side**2):
        if num_side <= i:
            j += 1
            i = 0
        ax[i][j].set_axis_off()
        i += 1

    if file:
        plotter.savefig(file, format=legend_type)
    else:
        plotter.show()



def write_color_legend(lg_colors, file):
    with open(file, 'wt') as colors_file:
        for lg in lg_colors:
            colors_file.write('%s\t"%s"\n' % (lg, to_hex(lg_colors[lg])))


            
def get_chrom_sizes(locus_bed, system='genomic'):
    chr_size = {}
    for chr_name in locus_bed['chr'].unique():
        chr_loci = locus_bed[locus_bed['chr'] == chr_name]
        chr_size[chr_name] = chr_loci['end'].max()

    return chr_size



def get_chrom_offsets(locus_bed, system='genomic'):
    sum_size = 0
    chr_offset = {}
    for chr_name in locus_bed['chr'].unique():
        chr_loci = locus_bed[locus_bed['chr'] == chr_name]
        chr_offset[chr_name] = sum_size
        sum_size += chr_loci['end'].max()
    return chr_offset



def read_color_legend(color_legend_name):
    colordict = {}
    color_legend = open(color_legend_name, 'r')
    for line in color_legend:
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        fields = line.split('\t')
    
        colordict[fields[0]] = fields[1].strip(""""'""")

    color_legend.close()
    
    return colordict



def read_table(file_name, header=None, types=None):
    input_table = pandas.read_table(file_name, header=None, sep="\t")

    if header:
        if types:
            assert num(types) == num(header), \
                'Unequal number of column names and types'
        
        num_columns = min(num(input_table.columns), num(header))
        if num(input_table.columns) > num_columns:
            # more columns than header, cut down table:
            input_table = input_table[list(range(num_columns))]
        else:
            # fewer columns than header, cut down header:
            header = header[:num_columns]
            types = types[:num_columns]
            
        # input_table = input_table[list(range(num(header)))]
        input_table.columns = header
        input_table = input_table.astype(dict(zip(header, types)))
        
    return input_table



def index_bed(input_bed, reference_bed=None):
    """Assumes reference_bed is NOT already indexed"""
    indexed_bed = []
    if reference_bed is None:
        for chr_name in input_bed['chr'].unique():
            chr_loci = input_bed[input_bed['chr'] == chr_name].copy()
            chr_loci.loc[:,'beg'] = list(range(0, num(chr_loci)))
            chr_loci.loc[:,'end'] = list(range(1, num(chr_loci)+1))
            indexed_bed.append(chr_loci)
    else:
        # find nearest indices for loci in input_bed that correspond to reference_bed
        for chr_name in reference_bed['chr'].unique():
            chr_loci = input_bed[input_bed['chr'] == chr_name].copy()
            ref_loci = reference_bed[reference_bed['chr'] == chr_name]
            ref_beg  = list(ref_loci.beg)
            ref_end  = list(ref_loci.end)
            for i in chr_loci.index:
                chr_loci.loc[i,'beg'] = bisect.bisect_right(ref_end, chr_loci.loc[i,'beg'])
                chr_loci.loc[i,'end'] = bisect.bisect_left(ref_beg, chr_loci.loc[i,'end'])
            indexed_bed.append(chr_loci)
    return pandas.concat(indexed_bed, ignore_index=True)
        


def sort_bed(input_bed, order):
    sorted_bed = []
    for chr_name in order:
        sorted_bed.append(
            input_bed[input_bed['chr'] == chr_name].copy()
        )

    return pandas.concat(sorted_bed, ignore_index=True)
        


def plot_axes(axes, max_size, label=None, xpad=0.04, ypad=0.04, xlab=True):
    axes.set_xbound(-1*xpad*max_size, (1+xpad)*max_size)
    axes.set_ybound(-1*ypad, 1.00)
    axes.spines['bottom'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.yaxis.set_ticks([])
    if not xlab:
        # axes.xaxis.set_ticks([])
        axes.xaxis.set_ticklabels([])
    axes.set_ylabel(label)


    
def plot_chr_glyph(axes, max_size, chr_size, centromere=None):
    if centromere is None:
        _plot_chr_without_centromere(axes, max_size, chr_size, ypos=-0.1)
    else:
        _plot_chr_with_centromere(axes, max_size, chr_size, centromere, ypos=-0.1)

        
        
def _plot_chr_with_centromere(axes, max_size, chr_size, centromere, ypos=0):
    #TODO: handling multiple centromeres
    padding = 0.01 * max_size
    Larm = centromere - padding
    Rarm = (chr_size - centromere) - padding
    
    axes.hlines(ypos, 0, Larm, linewidth=8, color='black', capstyle='round')
    axes.hlines(ypos, centromere + padding, centromere + padding + Rarm, linewidth=8, color='black', capstyle='round')
    axes.plot(centromere, ypos, marker='.', markersize=12.5, color="lightgrey")


    
def _plot_chr_without_centromere(axes, max_size, chr_size, ypos=0):
    axes.hlines(ypos, 0, chr_size, linewidth=8, color='black', capstyle='round')



def get_centromeres(centromeres, chr_name):
    if centromeres is not None:
        cen_loci = centromeres[centromeres['chr'] == chr_name]
        if num(cen_loci):
            return 0.5 * (cen_loci.loc[:,'beg'] + cen_loci.loc[:,'end'])

    return None



def get_colors(locus_map, file=None):
    label_colors = {}
    label_index = 0
    if file is not None:
        label_colors = read_color_legend(file)
        
    for label in locus_map['label'].unique():
        if label not in label_colors:
            label_colors[label] = _default_colors[label_index]
            label_index += 1

    return label_colors



def get_bin_freqs(chr_loci, bin_width=10, bin_shift=10,
                  label_colors=None, label_index=None, adaptive_binning=False):
    if label_index is None:
        label_index = index_list(label_colors)

    if adaptive_binning:
        blocks = get_blocks(chr_loci['label'], label_index, bin_width, bin_shift=1)
    else:
        blocks = [(chr_loci.index[0], chr_loci.index[-1] + 1)]
    
    bin_freq = None
    bin_begs = []
    bin_ends = []
    for block_beg, block_end in blocks:
        block_loci = chr_loci.loc[block_beg:block_end]
        
        _bin_freq, _bin_begs, _bin_ends = \
            get_freqs(
                block_loci['label'],
                label_index,
                bin_width=bin_width,
                bin_shift=bin_width
            )
        
        if bin_freq is None:
            bin_freq = _bin_freq
        else:
            bin_freq = numpy.vstack((bin_freq, _bin_freq))
            
        bin_begs.extend(_bin_begs)
        bin_ends.extend(_bin_ends)
        
    assert num(bin_freq) == num(bin_begs) == num(bin_ends), \
        'Mismatched number of frequencies and bin positions'
    
    chr_freq = bin_freq.transpose()
    chr_begs = numpy.zeros(num(bin_begs), dtype=numpy.int64)
    chr_ends = numpy.zeros(num(bin_ends), dtype=numpy.int64)
    for i in range(num(bin_ends)):
        chr_begs[i] = chr_ends[i-1]
        chr_ends[i] = chr_loci['end'][bin_ends[i]]
        
    return chr_freq, chr_begs, chr_ends



def plot_chr_freqs(axes, chr_loci, bin_width=10, bin_shift=10,
                   label_colors=None, label_index=None, adaptive_binning=False):
    if label_index is None:
        label_index = index_list(label_colors)

    chr_freq, chr_begs, chr_ends = \
        get_bin_freqs(
            chr_loci,
            bin_width=bin_width,
            bin_shift=bin_width,
            label_colors=label_colors,
            adaptive_binning=adaptive_binning
        )

    bottom = numpy.zeros(num(chr_ends))
    for label in label_index:
        axes.bar(
            x=chr_begs, 
            height=chr_freq[label_index[label],:], 
            width=chr_ends-chr_begs,
            bottom=bottom,
            color=label_colors[label],
            align='edge'
        )
        bottom += chr_freq[label_index[label],:]


    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = '' if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <locus.bed> <locus-to-lg.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0         10         20       30        40        50        60        70        80
    stream.write("  -A,--adaptive-binning\n")
    stream.write("     Optimizes binning boundaries where the majority label changes\n")
    stream.write("     class. This is achieved by modulating the number of loci in each\n")
    stream.write("     window (i.e., window size becomes dynamic).\n")
    stream.write("\n")
    stream.write("  -C,--input-color-legend <file>\n")
    stream.write("     Input label-to-color map file name. A two-column tab-separated\n")
    stream.write("     table with label in the first column and color in second. Colors\n")
    stream.write("     can be matplotlib.colormap names or quoted RGB hex values.\n")
    stream.write("     (default: tab10)\n")
    stream.write("\n")
    stream.write("  -c,--coordinate-system <str>\n")
    stream.write("     Output plot in coordinate system. Enumerative: {genomic,ordinal}\n")
    stream.write("     (default: genomic)\n")
    stream.write("\n")
    stream.write("  -l,--output-color-legend <str>\n")
    stream.write("     Output the label-to-color map to two files with the specified\n")
    stream.write("     file prefix. This option outputs both a tab-separated text file\n")
    stream.write("     and an image file of the legend in the format designated.\n")
    stream.write("\n")
    stream.write("  -O,--output-type <str>\n")
    stream.write("     Output plot image to file in the specified file format.\n")
    stream.write("     Enumerative: {%s}\n" % ','.join(_valid_output_types))
    stream.write("     (default: %s)\n" % _valid_output_types[0])
    stream.write("\n")
    stream.write("  -o,--output-file <str>\n")
    stream.write("     Output plot image to file with the specified file name.\n")
    stream.write("     (default: input file prefix)\n")
    stream.write("\n")
    stream.write("  -w,--bin-width <uint>\n")
    stream.write("     Plot label class densities in bins of the specified number of\n")
    stream.write("     loci.\n")
    stream.write("     (default: 10)\n")
    stream.write("\n")
    stream.write("  -x,--centromere-bed <file>\n")
    stream.write("     Input BED-formatted file of centromere positions. A circle is\n")
    stream.write("     added to each chromosome glyph to represent the centromere's\n")
    stream.write("     position.\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit.\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0         10         20       30        40        50        60        70        80
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)



def main(argv):
    short_options = 'hAC:c:O:o:w:l:x:'
    long_options = (
        'help',
        'adaptive-binning',
        'bin-width=',
        'coordinate-system=',
        'input-color-legend=','color-map-file=',
        'centromere-bed=',
        'output-file=',
        'output-type=',
        'output-color-legend=',
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    ypadding = 0.4
    bin_width = 10
    output_type = None
    output_file = None
    adaptive_binning = False
    input_color_legend_name = None
    output_color_legend_name = None
    centromere_bed = None
    centromere_bed_name = None
    coordinate_system = _valid_coordinate_systems[0]
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-w','--bin-width'):
            bin_width = int(value)
        elif flag in ('-A','--adaptive-binning'):
            adaptive_binning = True
        elif flag in ('-C','--input-color-legend','--color-map-file'):
            input_color_legend_name = value
        elif flag in ('-c','--coordinate-system'):
            coordinate_system = value
        elif flag in ('-l','--output-color-legend'):
            output_color_legend_name = value
        elif flag in ('-O','--output-type'):
            output_type = value
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-x','--centromere-bed'):
            centromere_bed_name = value

    if coordinate_system not in _valid_coordinate_systems:
        usage('Unsupported coordinate system: `%s`' % coordinate_system)
    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    if output_file is None:
        output_type = infer_output_type('', output_type)
        if arguments[0].lower().endswith('.bed'):
            output_file = arguments[0][:-4]
        output_file = '%s.%s' % (
            output_file,
            output_type
        )
    else:
        output_type = infer_output_type(output_file, output_type)
        
    locus_bed_name = arguments[0]
    locus_map_name = arguments[1]

    locus_bed = read_table(locus_bed_name, header=_bed_field_names, types=_bed_field_types)
    locus_map = read_table(locus_map_name, header=_map_field_names, types=_map_field_types)

    if centromere_bed_name is not None:
        centromere_bed = read_table(centromere_bed_name, header=_bed_field_names, types=_bed_field_types)
    
    if coordinate_system == 'ordinal':
        if centromere_bed is not None:
            centromere_bed = index_bed(centromere_bed, locus_bed)
        locus_bed = index_bed(locus_bed)


    chr_size = get_chrom_sizes(locus_bed)
    chr_offset = get_chrom_offsets(locus_bed)
    
    labeled_loci = pandas.merge(locus_bed, locus_map, on='name', how='inner')
    labeled_loci = sort_bed(labeled_loci, chr_offset)

    label_colors = get_colors(locus_map, file=input_color_legend_name)
        
    #TODO: need to check that there is a color for every label and assign if not.
        
    max_size = max(chr_size.values())
    num_chr = num(locus_bed['chr'].unique())
    max_index = num_chr - 1
    
    fig, ax = plotter.subplots(num_chr, figsize=(8, num_chr))
    plotter.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)

    for chr_index, chr_name in enumerate(locus_bed['chr'].unique()):
        chr_loci = labeled_loci[labeled_loci['chr'] == chr_name]
        cen_loci = get_centromeres(centromere_bed, chr_name)

        if num(chr_loci):    
            plot_chr_freqs(
                ax[chr_index],
                chr_loci,
                bin_width=bin_width,
                bin_shift=bin_width,
                label_colors=label_colors,
                adaptive_binning=adaptive_binning
            )
            
        plot_chr_glyph(
            ax[chr_index],
            max_size,
            chr_size[chr_name],
            centromere=cen_loci
        )
        plot_axes(
            ax[chr_index],
            max_size,
            label=chr_name,
            ypad=ypadding,
            xlab=(chr_index == max_index)
        )

    fig.tight_layout()     
    plotter.savefig(output_file, format=output_type)

    if output_color_legend_name:
        plot_color_legend(
            label_colors,
            file='%s.%s' % (
                output_color_legend_name,
                output_type
            )
        )
        write_color_legend(
            label_colors,
            file='%s.tsv' % (
                output_color_legend_name,
            )
        )


        
if __name__ == '__main__':
    main(sys.argv[1:])

