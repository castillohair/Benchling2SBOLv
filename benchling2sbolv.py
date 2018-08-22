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
     'part': 'RNACleavageSite'},
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
# Values that are not passed directly are the following:
#   - a value with the key ``label`` is added by ``plot_sequence()``, either
#     from the part's name, or as specified in ``plot_sequence()``'s argument
#     ``labels`` if present.
#   - option ``label_y_offset`` is multiplied by -1 for a part oriented in
#     reverse.
#   - options specified in ``RENDER_OPT[key]`` as functions are evaluated with
#     the calculated part's label width as an argument before being passed.
#   - CDS' ``color`` and ``label_color`` will be replaced as specified in
#     ``plot_sequence()``'s arguments ``cds_colors`` and ``cds_label_colors``.
RENDER_OPT = {
    'backbone_linewidth': 1.5,
    'Promoter': {'start_pad': lambda lw: max((lw - 10)/2, 0) + 2,
                 'end_pad': lambda lw: max((lw - 10)/2, 0) + 2,
                 'label_y_offset': -4,
                 },
    'RNACleavageSite': {'start_pad': lambda lw: max((lw - 5)/2, 0) + 2,
                        'end_pad': lambda lw: max((lw - 5)/2, 0) + 2,
                        'label_y_offset': -4,
            },
    'RBS': {'start_pad': lambda lw: max((lw - 10)/2, 0) + 2,
            'end_pad': lambda lw: max((lw - 10)/2, 0) + 2,
            'label_y_offset': -4,
            },
    'CDS': {'start_pad': 2,
            'end_pad': 2,
            'x_extent': lambda lw: max(lw, 10) + 15,
            'label_x_offset': -1,
            'label_style': 'italic',
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
                  cds_colors={},
                  cds_label_colors={},
                  chromosomal_locus=None,
                  ax=None,
                  ax_xlim=(0, 250),
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
    ax_xlim, ax_ylim : tuple-like, optional
        Axes limits.
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
    ax.set_xlim(ax_xlim)
    ax.set_ylim(ax_ylim)
    ax.set_aspect('equal')
    
    # Iterate and select parts to be plotted
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
            # Define new part
            part  = {}
            part['type'] = mapping['part']
            part['name'] = annotation.name
            # Orientation ('fwd') will be False if the reverse strand is
            # explicitly specified in the annotation, True if any other value,
            # not specified if None.
            if annotation.strand is not None:
                part['fwd'] = not (annotation.strand==-1)
            # Label is taken from the labels dictionary, or from the name.
            label = labels.get(part['name'], part['name'])
            # Get rendering options
            opts = {}
            opts['label'] = label
            # Colors for CDS glyph and label
            if part['type']=='CDS':
                if part['name'] in cds_colors:
                    opts['color'] = cds_colors[part['name']]
                if part['name'] in cds_label_colors:
                    opts['label_color'] = cds_label_colors[part['name']]
            # Get opts from RENDER_OPT if available
            part_type_opts = RENDER_OPT.get(part['type'], {})
            # Some elements in RENDER_OPT can be functions, in which case the
            # actual value is computed by calling that function with label_width
            # as an argument.
            # First, calculate label_width
            # Method from "https://stackoverflow.com/questions/24581194/\
            # matplotlib-text-bounding-box-dimensions"
            label_size = part_type_opts.get('label_size', 7)
            label_style = part_type_opts.get('label_style', 'normal')
            t = ax.text(0,
                        0,
                        label,
                        fontsize=label_size,
                        fontstyle=label_style)
            bb = t.get_window_extent(renderer=ax.figure.canvas.get_renderer())
            bb_datacoords = bb.transformed(ax.transData.inverted())
            label_width = bb_datacoords.width
            t.remove()
            # Iterate over options
            for k, v in part_type_opts.items():
                if hasattr(v, '__call__'):
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
                    # Special case: if type is CDS, and color not already set
                    # from cds_colors and cds_label_colors, set them here.
                    elif (k=='color') and \
                            (part['type']=='CDS') and \
                            ('color' not in opts):
                        opts['color'] = v
                    elif (k=='label_color') and \
                            (part['type']=='CDS') and \
                            ('label_color' not in opts):
                        opts['label_color'] = v
                    else:
                        opts[k] = v
            part['opts'] = opts
            # Save part
            parts.append(part)

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
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    if savefig is not None:
        pyplot.savefig(savefig, bbox_inches='tight', dpi=300)
