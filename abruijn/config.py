#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
File with configurations
"""

vals = {
        "circular_window" : 50000,
        "err_rate_threshold" : 0.24,
        "min_alignment_length" : 5000,
        "simple_kmer_length" : 4,
        "solid_kmer_length" : 10,
        "max_bubble_length" : 1000,
        "big_genome" : 50 * 1024 * 1024,

        "err_modes" : {
            "pacbio" : {
                "subs_matrix" : "pacbio_substitutions.mat",
                "hopo_matrix" : "pacbio_homopolymers.mat",
                "solid_missmatch" : 0.2,
                "solid_indel" : 0.2,
            },
            "nano" : {
                "subs_matrix" : "nano_substitutions.mat",
                "hopo_matrix" : "nano_homopolymers.mat",
                "solid_missmatch" : 0.3,
                "solid_indel" : 0.3,
            },
            "pacbio_hi_err" : {
                "subs_matrix" : "p6c4_substitutions.mat",
                "hopo_matrix" : "p6c4_homopolymers.mat",
                "solid_missmatch" : 0.25,
                "solid_indel" : 0.25,
            }
        }
    }
