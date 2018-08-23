import benchlingclient
import benchling2sbolv

# Specify benchling key
benchlingclient.LOGIN_KEY = 'write_your_key_here!'

# Plot many sequences
# Sequence names are 'pSR43_2', 'pSR43_3', 'pSR43_4'.
# From bases 700 to 5200
# Colors and label colors for CDSs ho1, pcyA, and sfgfp are custom.
# Save as "example_many_seqs.png"
seq = benchling2sbolv.plot_sequences(
    seq_names=['pSR43_2', 'pSR43_3', 'pSR43_4'],
    start_position=700,
    end_position=5200,
    cds_colors={'ho1': '#08519c',
                'pcyA': '#08519c',
                'ccaS': '#008800',
                },
    cds_label_colors={'ho1': '#ffffff',
                       'pcyA': '#ffffff',
                       'ccaS': '#ffffff',
                       },
    savefig='example_many_seqs.png')
