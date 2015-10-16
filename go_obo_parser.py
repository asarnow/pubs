#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
A constant-space parser for the GeneOntology OBO v1.2 format

Version 1.1
Based on Uli Koehler's original 1.0 Version.

(C) 2013 Uli Koehler
(C) 2015 Daniel Asarnow
"""
from collections import defaultdict


__author__ = "Uli Koehler"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__ = "Apache v2.0"


def process_go_term(go_term):
    """
    In an object representing a GO term, replace single-element lists with
    their only member.
    Returns the modified object as a dictionary.
    """
    ret = dict(go_term)  # Input is a defaultdict, might express unexpected behaviour
    for key, value in ret.iteritems():
        if len(value) == 1:
            ret[key] = value[0]
    return ret


def parse_go_obo(filename):
    """
    Parses a Gene Ontology dump in OBO v1.2 format.
    Yields each
    Keyword arguments:
        filename: The filename to read
    """
    with open(filename, "r") as infile:
        current_go_term = None
        for line in infile:
            line = line.strip()
            if not line:
                continue  # Skip empty
            if line == "[Term]":
                if current_go_term:
                    yield process_go_term(current_go_term)
                current_go_term = defaultdict(list)
            elif line == "[Typedef]":
                # Skip [Typedef sections]
                current_go_term = None
            else:  # Not [Term]
                # Only process if we're inside a [Term] environment
                if current_go_term is None:
                    continue
                key, sep, val = line.partition(":")
                current_go_term[key].append(val.strip())
        # Add last term
        if current_go_term is not None:
            yield process_go_term(current_go_term)


if __name__ == "__main__":
    """Print out the number of GO objects in the given GO OBO file"""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='The input file in GO OBO v1.2 format.')
    args = parser.parse_args()
    # Iterate over GO terms
    term_counter = 0
    for goTerm in parse_go_obo(args.infile):
        term_counter += 1
    print "Found %d GO terms" % term_counter
