# -*- coding: utf-8 -*-
# @Time : 2020/11/22 8:50
# @Author : Zhongyi Hua
# @FileName: ssr_utils.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

from Bio import SeqIO
from os import remove as del_file
import multiprocessing as multi
from tqdm import tqdm


def build_rep_set(repeat_file, unit_cutoff=None, motif_length_cutoff=None):
    """
        Outputs the repeats info dictionary used by the get_ssrs function.
        Takes list of repeat motifs from repeats file(output by generate_repeats function) as input.
        Creates a dictionary with expanded repeat as the key and (class, motif_length, strand) as values.
        Works either by "length_cutoff=" or by "unit_cutoff=" arguments.
    """
    repeats_out = dict()
    cutoffs = set()
    if unit_cutoff is None:
        unit_cutoff = {1: 9, 2: 6, 3: 6, 4: 6, 5: 5, 6: 4}
    for line in repeat_file:
        motif_dict = dict()
        L = line.strip().split('\t')
        motif_length = int(L[2])
        if motif_length < motif_length_cutoff:
            continue
        motif = L[0]
        motif = motif*unit_cutoff[motif_length]
        cutoffs.add(len(motif))
        motif_dict['class'] = L[1]
        motif_dict['motif_length'] = motif_length
        motif_dict['strand'] = L[3]
        repeats_out[motif] = motif_dict
    repeats_out['cutoff'] = sorted(list(cutoffs))

    return repeats_out


def get_ssrs(seq_record, repeats_info, out):
    """Native function that identifies repeats in fasta files."""
    if type(out) == str:
        out_file = open(out, 'w')
    else:
        out_file = out
    length_cutoffs = repeats_info['cutoff']
    input_seq = str(seq_record.seq).upper()
    input_seq_length = len(input_seq)
    for length_cutoff in length_cutoffs:
        fallback = length_cutoff - 1
        sub_start = 0  # substring start
        sub_stop = sub_start + length_cutoff  # substring stop
        while sub_stop <= input_seq_length:
            sub_stop = sub_start + length_cutoff
            sub_seq = input_seq[sub_start:sub_stop]
            if sub_seq in repeats_info:
                match = True
                repeat_data = repeats_info[sub_seq]
                motif_length = repeat_data['motif_length']
                rep_class = repeat_data['class']
                strand = repeat_data['strand']
                offset = length_cutoff % motif_length
                repeat_seq = input_seq[sub_start+offset:sub_start+offset+motif_length]
                i = 0
                while match:
                    j = sub_stop
                    if sub_stop >= input_seq_length:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        print(seq_record.id, sub_start, sub_stop, rep_class, match_length, strand, num_units, sub_seq[:motif_length], sep="\t", file=out_file)
                        sub_start = sub_stop - fallback
                    elif input_seq[j] == repeat_seq[i]:
                        sub_stop += 1
                        i += 1
                        if i >= motif_length:
                            i = 0
                    else:
                        match = False
                        match_length = sub_stop - sub_start
                        num_units = int(match_length/motif_length)
                        print(seq_record.id, sub_start, sub_stop, rep_class, match_length, strand, num_units, sub_seq[:motif_length], sep="\t", file=out_file)
                        sub_start = sub_stop - fallback
            else:
                sub_start += 1
    if type(out) == str:
        out_file.close()


def fasta_ssrs(args, repeats_info):
    """

    :param args:
    :param repeats_info:  return from build_rep_set function
    :return:
    """
    handle = open(args.input, 'r')

    records = [ _ for _ in SeqIO.parse(handle, 'fasta')]
    num_records = len(records)

    file_output = open(args.output, 'w')
    print('\t'.join(['seqid', 'start', 'end', 'class', 'length', 'strand', 'units', 'motif']), file=file_output)

    if args.threads > 1:
        i = 0
        pool = multi.Pool(processes=args.threads)
        for record in records:
            out_name = './temp_%s.tsv' % (i)
            i += 1
            pool.apply_async(get_ssrs, (record, repeats_info, out_name,))

        pool.close()
        pool.join()

        # Concat all the output files into one.
        temp_outs = tqdm(range(num_records), total=num_records)

        for o in temp_outs:
            name = './temp_%s.tsv' % (o)
            temp_outs.set_description("Concatenating file: %d " % (o))
            with open(name, 'r') as fh:
                for line in fh:
                    print(line.strip(), file=file_output)
            del_file(name)

    elif args.threads == 1:
        records = tqdm(records, total=num_records)
        for record in records:
            records.set_description("Processing %s" % (record.id))
            get_ssrs(record, repeats_info, file_output)

