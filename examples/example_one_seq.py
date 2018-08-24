import matplotlib
import benchlingclient
import benchling2sbolv

# The following is necessary so that latex-rendered text is not italicized by
# default.
matplotlib.rcParams['mathtext.default'] = 'regular'

# Specify benchling key
benchlingclient.LOGIN_KEY = 'write_your_key_here!'

# Plot one sequence
# Name: pSR58_6, from base 700 to 3000.
# Promoter 'PcpcG2-172' has a special latex-formatted label so that it
# looks prettier.
# Color of CDSs 'ccaR' and 'sfgfp' are changed to different shades of
# green.
# Save figure as "example_one_seq.png"
seq = benchling2sbolv.plot_sequence(
    seq_name='pSR58_6',
    start_position=700,
    end_position=3000,
    labels={'PcpcG2-172': r'$P_{\mathit{cpcG2}-172}$'},
    cds_colors={'ccaR': '#aaffaa',
                'sfgfp': '#00FF00',
                },
    savefig='example_one_seq.png')
