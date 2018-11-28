#(c) 2016-2018 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Modifies repeat graph using the Tresle output
"""

import logging

logger = logging.getLogger()


class Connection:
    __slots__ = ("id", "path", "sequence")
    def __init__(self, id=None, path=[], sequence=""):
        self.id = id
        self.path = path
        self.sequence = sequence


def apply_changes(repeat_graph, trestle_results, resolved_repeats_file):
    #repeat_graph.output_dot("before.dot")
    connections = _get_connections(trestle_results)
    _add_sequences(connections)
    for conn in connections:
        repeat_graph.separate_path(conn.path, conn.id, conn.sequence)

    #repeat_graph.output_dot("after.dot")


def _get_connections(trestle_results):
    connections = []
    resolved_repeats = set()
    with open(trestle_results, "r") as f:
        for line in f:
            if line.startswith("Repeat"): continue

            tokens = line.strip().split()
            repeat_id, bridged = int(tokens[0]), tokens[5]
            if bridged == "True" and abs(repeat_id) not in resolved_repeats:
                resolved_repeats.add(abs(repeat_id))

                path_1, path_2 = tokens[9].split("+")
                in_1, out_1 = path_1.split("|")
                in_2, out_2 = path_2.split("|")
                connection_1 = [int(in_1[2:]), repeat_id,
                                int(out_1[3:])]    #stripping in/out prefixes
                connection_2 = [int(in_2[2:]), repeat_id,
                                int(out_2[3:])]    #stripping in/out prefixes

                logger.info("Repeat {0}: {1}, {2}"
                    .format(repeat_id, connection_1, connection_2))

                new_seq_id = "trestle_resolved_" + str(repeat_id) + "_copy_"
                connections.extend([Connection(new_seq_id + "1", connection_1, ""),
                                    Connection(new_seq_id + "2", connection_2, "")])

    return connections


def _add_sequences(connections):
    for conn in connections:
        conn.sequence = "ACGT"
