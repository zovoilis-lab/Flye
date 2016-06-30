#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
File with configurations
"""

vals = {
        "circular_window" : 50000,
        "err_rate_threshold" : 0.24,
        "profiles" : {
            "pacbio" : {
                "subs_matrix" : "pacbio_substitutions.mat",
                "hopo_matrix" : "pacbio_homopolymers.mat"
            },
            "nano" : {
                "subs_matrix" : "nano_substitutions.mat",
                "hopo_matrix" : "nano_homopolymers.mat"
            }
        }
    }
