"""
Benchling2SBOLv

"""

# Versions should comply with PEP440.  For a discussion on single-sourcing
# the version across setup.py and the project code, see
# https://packaging.python.org/en/latest/single_source_version.html
__version__ = '0.1.0'

import copy

import benchlingclient
import dnaplotlib
from matplotlib import pyplot

# ``ANN_PARTS_MAPPING`` maps benchling annotations to dnaplotlib's part types.
# Each element of this list is a dictionary, where the value given by the
# ``annotation`` key is a dictionary that specifies all features that a
# benchling annotation should comply with to be represented by the part type
# given by the key ``part``.
# Elements of this list are checked in order. If an annotation matches two
# elements in this list, it will be represented by the part type specified by
# the first element.
ANN_PARTS_MAPPING = [
    {'annotation': {'type': 'Promoter'},
     'part': 'Promoter'},
    {'annotation': {'type': 'ncRNA', 'name': 'RiboJ'},
     'part': 'Ribonuclease'},
    {'annotation': {'type': 'RBS'},
     'part': 'RBS'},
    {'annotation': {'type': 'CDS'},
     'part': 'CDS'},
    {'annotation': {'type': 'Terminator'},
     'part': 'Terminator'},
     ]

# ``RENDER_OPT`` specifies the rendering options for each part type. With the
# exception of ``backbone_linewidth``, every key is the name of a dnaplotlib
# part type. Each value given by a key is a dictionary with options to be passed
# to dnaplotlib, via a part's ``opts`` parameter.
# A few caveats:
#   - a value with the key ``label`` is added by ``plot_sequence()``, either
#     from the part's name, or as specified in ``plot_sequence()``'s argument
#     ``labels`` if present.
#   - options ``label_x_offset`` and ``label_y_offset`` are multiplied by -1 for
#     a part oriented in reverse.
#   - options specified in ``RENDER_OPT[key]`` as functions are evaluated with
#     the calculated part's label width as an argument before being passed.
#   - CDS' ``color`` and ``label_color`` will be replaced as specified in
#     ``plot_sequence()``'s arguments ``cds_colors`` and ``cds_label_colors``.
#   - Options under the keys ``CDS`` and ``CDSFragment`` are used for CDS
#     fragments except for the last one when a CDS is split via
#     ``plot_sequence()``'s ``cds_split_char`` argument. In this case, when
#     the same option is specified under ``CDS`` and ``CDSFragment``, the one
#     under ``CDSFragment`` takes precedence.
RENDER_OPT = {
    'backbone_linewidth': 1.5,
    'Promoter': {'start_pad': lambda lw: max((lw - 10)/2, 0) + 2,
                 'end_pad': lambda lw: max((lw - 10)/2, 0) + 2,
                 'label_y_offset': -4,
                 },
    'Ribonuclease': {'start_pad': lambda lw: max((lw - 5)/2, 0) + 2,
                     'end_pad': lambda lw: max((lw - 5)/2, 0) + 2,
                     'label_y_offset': -4,
            },
    'RBS': {'start_pad': lambda lw: max((lw - 10)/2, 0) + 2,
            'end_pad': lambda lw: max((lw - 10)/2, 0) + 2,
            'label_y_offset': -4,
            },
    'CDS': {'start_pad': 2,
            'end_pad': 2,
            'x_extent': lambda lw: max(lw, 10) + 14,
            'label_x_offset': -1.5,
            'label_style': 'italic',
            },
    'CDSFragment': {'x_extent': lambda lw: max(lw, 10) + 12,
                    'label_x_offset': 0,
                    },
    'Terminator': {'start_pad': lambda lw: max((lw - 10)/2, 0) + 2,
                   'end_pad': lambda lw: max((lw - 10)/2, 0) + 2,
                   'label_y_offset': -4,
                   'label_size': 7,
                   },
    '5ChromosomalLocus': {'dashed_end': True,
                          'label_y_offset': 4,
                          'label_style': 'italic',
                          },
    '3ChromosomalLocus': {'dashed_end': True,
                          'label_y_offset': 4,
                          'label_style': 'italic',
                          },
}

def check_annotation_features(annotation, features):
    """
    Check whether an annotation's features match the specified features.

    Parameters
    ----------
    annotation : benchlingclient.Annotation
        Annotation whose features will be checked
    features : dict
        Dictionary with ``feature_name: feature_value`` pairs to be
        checked.

    Returns
    -------
    match : bool
        True if the annotation's features match, False otherwise.

    """
    match = True
    for feature_name, feature_value in features.items():
        if annotation.__getattribute__(feature_name) != feature_value:
            match = False
            break

    return match

def plot_sequence(seq=None,
                  seq_name=None,
                  start_position=None,
                  end_position=None,
                  labels={},
                  ignore_names=[],
                  cds_split_char='',
                  cds_colors={},
                  cds_label_colors={},
                  chromosomal_locus=None,
                  ax=None,
                  ax_x_extent=250,
                  ax_x_alignment='center',
                  ax_ylim=(-15, 15),
                  savefig=None):
    """
    Plot a specified benchling sequence as SBOL visual

    Parameters
    ----------
    seq : benchlingclient.DNASequence
        Sequence to be plotted. Can be omitted if `seq_name` is provided.
    seq_name : str, optional
        Name of the sequence to load and plot. Ignored if `seq` is
        specified.
    start_position, end_position : int, optional
        Only annotations in the sequence completely contained by the range
        given by these two parameters will be plotted.
    labels : dict, optional
        Dictionary with ``name: label`` pairs that specify a glyph's label,
        when the label is different than the name. If a part's name is not
        a key of `labels`, the name will be used as the glyph's label.
    ignore_names : list, optional
        Names of parts that should not be plotted.
    cds_split_char : str, optionsl
        If specified, a CDS whose name contains `cds_split_char` as a
        substring will be split along this substring and shown as a
        multipart CDS. A multipart CDS looks like a single CDS divided into
        many fragments, each with its own label and color that can be
        specified in `labels`, `cds_colors`, and `cds_label_colors`.
    cds_colors : dict, optional
        Dictionary with ``name: color`` pairs that specify a CDS' face
        color, when the color is different than the one specified by
        ``RENDER_OPT['CDS']``.
    cds_label_colors : dict, optional
        Dictionary with ``name: label_color`` pairs that specify a CDS'
        label color, when the color is different than the one specified by
        ``RENDER_OPT['CDS']``.
    chromosomal_locus : str, optional
        If specified, a pair of "ChromosomalLocus" glyphs will be added to
        both sides of the rendered design, and a label with the text given
        by `chromosomal_locus` on the left glyph.
    ax : matplotlib.axes, optional
        Axes to draw into.
    ax_x_extent : float, optinal
        Range covered by the x axis limits.
    ax_x_alignment : {'left', 'center', 'right'}
        Alignment of the rendered diagram in the x axis.
    ax_ylim : tuple-like, optional
        Y axis limits.
    savefig : str, optional
        If specified, save figure with a filename given by `savefig`.

    Raises
    ------
    ValueError
        If more than one sequence with the name given by `seq_name` was
        found.

    """
    # If seq is provided, the following is not executed.
    # If not, load from seq_name
    if seq is None:
        if seq_name is not None:
            # Load sequence from benchling
            seq_list = benchlingclient.DNASequence.list_all(name=seq_name)
            # Test that only one sequence has been found
            if len(seq_list) > 1:
                raise ValueError("more than one sequence found with name {}".\
                    format(seq_name))
            elif len(seq_list) < 1:
                raise ValueError("no sequence with name {} found".\
                    format(seq_name))
            seq = seq_list[0]
        else:
            # No sequence or sequence name provided, raise exception
            raise ValueError("seq or seq_name should be provided")

    # Get annotations
    annotations = seq.annotations.copy()

    # Remove annotations outside of range specified
    if start_position is None:
        start_position = 0
    if end_position is None:
        end_position = seq.length - 1
    annotations_filtered = []
    for annotation in annotations:
        if (annotation.start >= start_position) and \
                (annotation.end <= end_position):
            annotations_filtered.append(annotation)
    annotations = annotations_filtered

    # Sort by start position
    annotations = sorted(annotations, key=lambda x: x.start)

    # Initialize plot
    if ax is None:
        fig, ax = pyplot.subplots()
    # Set axis limits and aspect right away
    # This allows for a more precise estimation of label widths later on
    ax.set_xlim((0, ax_x_extent))
    ax.set_ylim(ax_ylim)
    ax.set_aspect('equal')
    
    # Iterate and extract parts to be plotted
    parts = []
    for annotation in annotations:
        # Iterate through mapping and select appropriate part type
        match = False
        for mapping in ANN_PARTS_MAPPING:
            if check_annotation_features(annotation, mapping['annotation']):
                match = True
                break
        # Only add part if matching part has been found
        if match:
            part_type = mapping['part']
            part_name = annotation.name
            # Part will be split if type is CDS and a cds_split_char has been
            # specified.
            # All resulting parts, but the last one, will be given the temporary
            # part type "CDSFragment". This is so arrowheads and padding are
            # modified later to make a continuous multipart CDS glyph.
            if cds_split_char and mapping['part']=='CDS':
                part_names = part_name.split(cds_split_char)
                # Remove parts flagged to be ignored
                part_names = [p for p in part_names if p not in ignore_names]
                if len(part_names)==0:
                    continue
                # Reverse if orientation is reverse
                if annotation.strand==-1:
                    part_names = part_names[::-1]
                for part_index, part_name in enumerate(part_names):
                    # Define new part
                    part = {}
                    if annotation.strand==-1:
                        if part_index <= 0:
                            part['type'] = 'CDS'
                        else:
                            part['type'] = 'CDSFragment'
                    else:
                        if part_index < (len(part_names) - 1):
                            part['type'] = 'CDSFragment'
                        else:
                            part['type'] = 'CDS'
                    part['name'] = part_name
                    # Orientation ('fwd') will be False if the reverse strand is
                    # explicitly specified in the annotation, True if any other
                    # value, not specified if None.
                    if annotation.strand is not None:
                        part['fwd'] = not (annotation.strand==-1)
                    # Save part
                    parts.append(part)
            else:
                # Check if part name is flagged to be ignored
                if part_name in ignore_names:
                    continue
                # Define new part
                part  = {}
                part['type'] = part_type
                part['name'] = part_name
                # Orientation ('fwd') will be False if the reverse strand is
                # explicitly specified in the annotation, True if any other
                # value, not specified if None.
                if annotation.strand is not None:
                    part['fwd'] = not (annotation.strand==-1)
                # Save part
                parts.append(part)

    # Construct plotting options for each part
    for part_index, part in enumerate(parts):
        # Initialize dictionary with rendering options
        opts = {}

        # GENERAL TYPE-DEPENDENT RENDERING OPTIONS
        # ========================================

        # Get opts from RENDER_OPT if available
        if part['type']=='CDSFragment':
            # CDSFragment combines options from 'CDS' and 'CDSFragment', with
            # the latter taking precedence.
            part_type_opts = RENDER_OPT.get('CDSFragment', {})
            part_type_opts_cds = RENDER_OPT.get('CDS', {})
            for k, v in part_type_opts_cds.items():
                if k not in part_type_opts:
                    part_type_opts[k] = v
        else:
            part_type_opts = RENDER_OPT.get(part['type'], {})

        # PART-SPECIFIC RENDERING OPTIONS
        # ===============================

        # The following deals with issues specific to the type "CDSFragment"
        # 'CDSFragment' does not have an arrowhead
        if part['type']=='CDSFragment':
            opts['arrowhead_length'] = 0
            opts['arrowhead_height'] = 0
        # Correction for label_x_offset
        # Removing padding will affect how the CDS looks but not the label,
        # therefore the label offset has to be corrected.
        label_x_offset_c = 0
        if ('fwd' in part) and (part['fwd']==False):
            # Start padding is zero if there is a 'CDS' or CDSFragment' to the
            # left
            if part['type']=='CDSFragment':
                if (part_index > 0) and \
                        (parts[part_index-1]['type'] in ['CDS', 'CDSFragment']):
                    opts['end_pad'] = 0
                    label_x_offset_c += part_type_opts.get('end_pad', 0)/2
            # End padding is zero if there is a 'CDSFragment' to the right
            if part['type'] in ['CDS', 'CDSFragment']:
                if (part_index < (len(parts)-1)) and \
                        (parts[part_index+1]['type']=='CDSFragment'):
                    opts['start_pad'] = 0
                    label_x_offset_c -= part_type_opts.get('start_pad', 0)/2
        else:
            # Start padding is zero if there is a 'CDSFragment' to the left
            if part['type'] in ['CDS', 'CDSFragment']:
                if (part_index > 0) and \
                        (parts[part_index-1]['type']=='CDSFragment'):
                    opts['start_pad'] = 0
                    label_x_offset_c -= part_type_opts.get('start_pad', 0)/2
            # End padding is zero if there is a 'CDS' or 'CDSFragment' to the
            # right
            if part['type']=='CDSFragment':
                if (part_index < (len(parts)-1)) and \
                        (parts[part_index+1]['type'] in ['CDS', 'CDSFragment']):
                    opts['end_pad'] = 0
                    label_x_offset_c += part_type_opts.get('end_pad', 0)/2

        # Colors for CDS glyph and label
        if part['type'] in ['CDS', 'CDSFragment']:
            if part['name'] in cds_colors:
                opts['color'] = cds_colors[part['name']]
            if part['name'] in cds_label_colors:
                opts['label_color'] = cds_label_colors[part['name']]

        # Some elements in RENDER_OPT can be functions, in which case the
        # actual value is computed by calling that function with the label width
        # as an argument.
        # Label is taken from the labels dictionary, or from the name.
        label = labels.get(part['name'], part['name'])
        opts['label'] = label
        # Method to calculate label_width from "https://stackoverflow.com/\
        # questions/24581194/matplotlib-text-bounding-box-dimensions"
        t = ax.text(0,
                    0,
                    label,
                    fontsize=part_type_opts.get('label_size', 7),
                    fontstyle=part_type_opts.get('label_style', 'normal'))
        bb = t.get_window_extent(renderer=ax.figure.canvas.get_renderer())
        bb_datacoords = bb.transformed(ax.transData.inverted())
        label_width = bb_datacoords.width
        t.remove()

        # Iterate over options
        for k, v in part_type_opts.items():
            if k in opts:
                # Don't override specific options established above
                continue
            elif hasattr(v, '__call__'):
                # Evaluate option based on label width
                opts[k] = v(label_width)
            else:
                # Special case: label_y_offset is modified depending on the
                # orientation of the part
                if k=='label_y_offset':
                    if ('fwd' in part) and (not part['fwd']):
                        opts[k] = -v
                    else:
                        opts[k] = v
                # Special case: label_x_offset is modified depending on the
                # orientation of the part and an offset correction 
                elif k=='label_x_offset':
                    if ('fwd' in part) and (not part['fwd']):
                        opts[k] = -(v + label_x_offset_c)
                    else:
                        opts[k] = v + label_x_offset_c
                else:
                    opts[k] = v

        # Add options
        part['opts'] = opts

    # 'CDSFragment' is not a real dnaplotlib part type. Change it to 'CDS'
    for part in parts:
        if part['type']=='CDSFragment':
            part['type'] = 'CDS'

    # Add chromosomal locus parts if specified
    if chromosomal_locus is not None:
        # 5' glyph
        cl5 = {}
        cl5['type'] = '5ChromosomalLocus'
        cl5['name'] = 'cl5_{}'.format(chromosomal_locus)
        cl5['fwd'] = True
        opts = RENDER_OPT.get('5ChromosomalLocus', {}).copy()
        opts['label'] = chromosomal_locus
        opts['linewidth'] = RENDER_OPT.get('backbone_linewidth', 1)
        cl5['opts'] = opts
        # 3' glyph
        cl3 = {}
        cl3['type'] = '3ChromosomalLocus'
        cl3['name'] = '3cl_{}'.format(chromosomal_locus)
        opts = RENDER_OPT.get('3ChromosomalLocus', {}).copy()
        opts['linewidth'] = RENDER_OPT.get('backbone_linewidth', 1)
        cl3['opts'] = opts
        # Save glyphs
        parts = [cl5] + parts + [cl3]

    # Create the DNAplotlib renderer
    dr = dnaplotlib.DNARenderer(
        linewidth=RENDER_OPT.get('backbone_linewidth', 1))

    # Redend the DNA to axis
    start, end = dr.renderDNA(ax, parts, dr.SBOL_part_renderers())
    # Set x axis limits depending on alignment
    # Note that the extent is always ax_x_extent
    if ax_x_alignment=='left':
        ax.set_xlim((start, start + ax_x_extent))
    if ax_x_alignment=='right':
        ax.set_xlim((end - ax_x_extent, end))
    if ax_x_alignment=='center':
        ax.set_xlim(((start + end - ax_x_extent)/2,
                     (start + end + ax_x_extent)/2))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    if savefig is not None:
        pyplot.savefig(savefig, bbox_inches='tight', dpi=300)

def plot_sequences(seqs=None,
                   seq_names=None,
                   start_position=None,
                   end_position=None,
                   labels={},
                   ignore_names=[],
                   cds_split_char='',
                   cds_colors={},
                   cds_label_colors={},
                   chromosomal_locus=None,
                   ax_x_extent=250,
                   ax_x_alignment='center',
                   ax_ylim=(-15, 15),
                   hspace=0,
                   figsize=None,
                   savefig=None):
    """
    Plot several benchling sequences as SBOL visual.

    Parameters
    ----------
    seqs : list of benchlingclient.DNASequence
        List of sequences to be plotted. Can be omitted if `seq_names` is
        provided.
    seq_names : list of str, optional
        Names of the sequences to load and plot. Ignored if `seqs` is
        specified.
    hspace : float, optional
        Vertical space to be kept between sequences. The default is zero,
        which, if `figsize` has the same aspect ratio than all axes stacked
        vertically, perfectly stacks the sequences' underlying axes.
        `hspace` can be positive or negative, which increases or reduces
        this distance, respectively.
    figsize : tuple, optional
        Size of the figure to be created. If not specified, get the width
        from the current defaults, and calculate the height to match the
        aspect ratio of all axes stacked vertically.
    savefig : str, optional
        If specified, save figure with a filename given by `savefig`.

    Other parameters
    ----------------
    All parameters in ``plot_sequence()``, with the exception of `seq`,
    `seq_name`, and `savefig` can be passed to this function. These will
    then be directly passed to ``plot_sequence()`` when it is called to
    plot each diagram.

    """
    # If seqs is provided, the following is not executed.
    # If not, load from seq_names
    if seqs is None:
        if seq_names is not None:
            seqs = []
            for seq_name in seq_names:
                # Load sequence from benchling
                seq_list = benchlingclient.DNASequence.list_all(name=seq_name)
                # Test that only one sequence has been found
                if len(seq_list) > 1:
                    raise ValueError("more than one sequence found with name "
                        "{}".format(seq_name))
                elif len(seq_list) < 1:
                    raise ValueError("no sequence with name {} found".\
                        format(seq_name))
                seqs.append(seq_list[0])
        else:
            # No sequence or sequence name provided, raise exception
            raise ValueError("seqs or seq_names should be provided")

    # Initialize figure
    if figsize is None:
        fig_width = pyplot.rcParams.get('figure.figsize')[0]
        fig_height = fig_width*(ax_ylim[1] - ax_ylim[0])/ax_x_extent*len(seqs)
        figsize = (fig_width, fig_height)
    fig = pyplot.figure(figsize=figsize)
    # Plot each sequence in a separate axes
    for seq_index, seq in enumerate(seqs):
        ax = fig.add_subplot(len(seqs), 1, seq_index + 1)
        plot_sequence(seq=seq,
                      ax=ax,
                      start_position=start_position,
                      end_position=end_position,
                      labels=labels,
                      ignore_names=ignore_names,
                      cds_split_char=cds_split_char,
                      cds_colors=cds_colors,
                      cds_label_colors=cds_label_colors,
                      chromosomal_locus=chromosomal_locus,
                      ax_x_extent=ax_x_extent,
                      ax_x_alignment=ax_x_alignment,
                      ax_ylim=ax_ylim)
    # Adjust vertical space between subplots
    fig.subplots_adjust(hspace=hspace)

    # Save figure if specified
    if savefig is not None:
        pyplot.savefig(savefig, bbox_inches='tight', dpi=300)
