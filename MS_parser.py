#!/usr/bin/env python2.7
import os
import os.path
import sys


def main(args):
    fname = args[1]
    with open(fname, 'rU') as f:
        # with closes as well as opening
        head = f.readline().rstrip(os.linesep).split('\t')
        # creating list of headers
        data = {h: [] for h in head}
        for l in f:
            tok = l.rstrip(os.linesep).split('\t')
            # os.linesep represents line seperator on running os
            if len(tok) == 1:
                continue
            # continue tells forloop to go to next iteration so that if something is
            # not there there's no error when you call 5 and its not there
            assert len(tok) == len(head)  # aborts script if this isn't true, makes sure that blank cells are accounted for
            for i in xrange(0, len(tok)):
                data[head[i]].append(tok[i])
    print data
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

    # it's checking if the top level interpreter is running, if so run our thing
