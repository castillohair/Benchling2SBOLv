import benchlingclient
import benchling2sbolv

# Specify benchling key
benchlingclient.LOGIN_KEY = 'write_your_key_here!'

seq = benchling2sbolv.plot_sequence(seq_name='lSC0328',
                                    start_position=3000,
                                    end_position=5600,
                                    chromosomal_locus='amyE',
                                    cds_colors={'sfgfp': (0, 0.4, 0.1),
                                                },
                                    cds_label_colors={'sfgfp': (1, 1, 1),
                                                       },
                                    savefig='example.png')

