# -*- coding: utf-8 -*-
# @Time : 2020/11/22 9:37
# @Author : Zhongyi Hua
# @FileName: core_func.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import argparse
import os
from ssr_utils import fasta_ssrs, build_rep_set
from plugin_utils import PrimerDesign, ScreenSSR, combine2


def getArgs():
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input sequence file. fasta format')
    required.add_argument('-o', '--output', required=True, metavar='<FILE>', help='Output table file path.')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-c', '--circular', action='store_true', default=False,
                          help='The sequences are circular sequence')
    optional.add_argument('-m', '--min-motif-size', type=int, nargs='+', metavar='<INT>', default=[9, 6, 6, 6, 5, 5],
                          help='Minimum size of a repeat motif in bp. Default 9, 6, 6, 6, 5, 5 for mono-, di-, tri-, '
                               'tetra-, penta-, hexa- nucleotide repeat, respectively')
    optional.add_argument('-p', '--primer', action='store_true', default=False,
                          help='Design primer at the same time')
    optional.add_argument('-s', '--screen', action='store_true', default=False,
                          help='Screen out the SSR could used in identification at the same time')
    optional.add_argument('-a', '--assembly', nargs='+',
                          help='Confirm the SSR could used in Identification at the same time')
    # Multiprocessing threads
    optional.add_argument('-t', '--threads', type=int, metavar='<INT>', default=1, help='Number of threads to run the process on. Default is 1.')

    args = parser.parse_args()

    if args.screen and (not args.assembly):
        parser.error("-a must with -s/--screen")

    return args


def main():
    args = getArgs()

    location = os.path.abspath(os.path.join(__file__, '..'))
    motif_size_dict = {_[0]: _[1] for _ in zip([1, 2, 3, 4, 5, 6], args.min_motif_size)}
    repeats_info = build_rep_set(open(os.path.join(location, 'all_repeats_1-6nt.txt'), 'r'), motif_size_dict)
    fasta_ssrs(args, repeats_info)
    if args.screen:
        screen = ScreenSSR(args)
        screen.blast_pair()
        screen.blast_parse()
        screen.remove_tmp_file()
    if args.primer:
        primer_ = PrimerDesign(args)
        primer_.prepare_cfg(circular=args.circular)
        primer_.design_primer()
        primer_.parse_output()
        primer_.remove_tmp_file()
    if args.screen and args.primer:
        # filter primer for screened SSR
        combine2(args)


if __name__ == '__main__':
    main()
