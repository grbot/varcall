#!/usr/bin/env python2.7

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--inBED", default="${inBED}", help="")
parser.add_argument("--chrm", default="${chrm}", help="")
parser.add_argument("--outBED", default="${outBED}", help="")
args = parser.parse_args()


def split_bed(inBED, chrm, outBED):
    """
    :param inBED:
    :param chrm:
    :param outBED:
    :return:
    """
    output = open(outBED, 'w')
    for line in open(inBED):
        data = line.strip().split()
        if str(data[0]) == str(chrm):
            output.writelines(line)
    output.close()


if __name__ == '__main__':
    split_bed(args.inBED, args.chrm, args.outBED)
