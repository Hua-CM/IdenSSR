# -*- coding: utf-8 -*-
# @Time : 2020/11/23 9:40
# @Author : Zhongyi Hua
# @FileName: plugin_utils.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import pandas as pd
import os
from tempfile import mktemp
from Bio import SeqIO
from collections import deque
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


def get_seq(_seq, _start, _end):
    if _start <= 200:
        return _seq.seq[len(_seq.seq)-201+_start:] + _seq.seq[0: _end+200]
    elif _end >= len(_seq.seq)-200:
        return _seq.seq[_start-201:] + _seq.seq[:200-(len(_seq.seq)-_end)]
    else:
        return _seq.seq[_start-201: _end+200]


def del_directory(_dir):
    for r_, _d_, f_ in os.walk(_dir):
        for files in f_:
            os.remove(os.path.join(r_, files))
        os.removedirs(r_)


def combine2(args):
    tb_primer = pd.read_table('_primer'.join(os.path.splitext(args.out)))
    tb_screen_out = pd.read_table('_screen'.join(os.path.splitext(args.out)))
    tb_primer['tmp_query'] = tb_primer.apply(lambda x: '_'.join([str(x['seqid']), str(x['start'], str(x['stop']))]),
                                             axis=1)
    tb_screen_out['tmp_query'] = tb_screen_out.apply(lambda x: '_'.join([str(x['seqid']), str(x['start'],str(x['stop']))]),
                                                     axis=1)
    tb_primer[tb_primer['tmp_query'].isin(tb_screen_out['tmp_query'].to_list())].to_csv(
        '_filter_primer'.join(os.path.splitext(args.out))
    )


class PrimerDesign:
    location = os.path.abspath(os.path.join(__file__, '..'))

    def __init__(self, args):
        self.sequences = {_.id: _.seq for _ in SeqIO.parse(args.input, 'fasta')}
        self.lengths = {_id: len(_seq.seq) for _id, _seq in self.sequences.items()}
        self.ssr_info = pd.read_table(args.output)
        self.output = '_primer'.join(os.path.splitext(args.out))
        self._config = mktemp()
        self._p3out = mktemp()

    def prepare_cfg(self, circular=False):
        if not circular:
            meta_tb = self.ssr_info[(self.ssr_info['start'] > 200) & (self.ssr_info['end'] < self.lengths[self.ssr_info[
                'seqid']]-200)]
        else:
            meta_tb = self.ssr_info
        _write_list = []
        for _idx, _item in meta_tb.iterrows():
            genome_seq = self.sequences[_item['seqid']]
            _sequence_id = "SEQUENCE_ID=" + str(genome_seq.id) + "_" + str(_item['start']) + "_" + str(_item['end'])
            _sequence_template = "SEQUENCE_TEMPLATE=" + str(get_seq(genome_seq, _item['start'], _item['end']))
            _write_list.append(_sequence_id)
            _write_list.append(_sequence_template)
            _write_list.append('SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=1,180,' +
                               str(181 + int(_item['end']) - int(_item['start'])) +
                               ',200')
            _write_list.append('=')
        with open(self._config, 'w') as f2:
            f2.write('=\n')
            f2.write('\n'.join(_write_list))
            f2.write("\n")

    def design_primer(self):
        _cmd = ' '.join(
                       ['primer3_core',
                        '--p3_settings_file='+os.path.join(PrimerDesign.location, 'SSR_primer3_settings.txt'),
                        '--output='+self._p3out,
                        self._config]
                       )
        os.system(_cmd)

    def parse_output(self):
        _a = open(self._p3out, 'r')
        _file = deque(_a.read().splitlines())
        _module = []
        _record_list = []
        for _line in _file:
            _module.append(_line)
            if _line == '=':
                _model_dict = dict(_.split('=') for _ in _module)
                if 'SEQUENCE_ID' in _model_dict:
                    if not _model_dict['PRIMER_PAIR_NUM_RETURNED'] == '0':
                        _record = {'seqid': '_'.join(_model_dict['SEQUENCE_ID'].split('_')[:-2]),
                                   'start': _model_dict['SEQUENCE_ID'].split('_')[-2],
                                   'end': _model_dict['SEQUENCE_ID'].split('_')[-1],
                                   'Product Size': _model_dict['PRIMER_PAIR_0_PRODUCT_SIZE'],
                                   'Forward':  _model_dict['PRIMER_LEFT_0_SEQUENCE'],
                                   'Forward_position': _model_dict['PRIMER_LEFT_0'].split(',')[0],
                                   'Forward_length': _model_dict['PRIMER_LEFT_0'].split(',')[1],
                                   'Forward_TM': _model_dict['PRIMER_LEFT_0_TM'],
                                   'Reverse': _model_dict['PRIMER_RIGHT_0_SEQUENCE'],
                                   'Reverse_position': _model_dict['PRIMER_RIGHT_0'].split(',')[0],
                                   'Reverse_length': _model_dict['PRIMER_RIGHT_0'].split(',')[1],
                                   'Reverse_TM': _model_dict['PRIMER_RIGHT_0_TM']
                                   }
                        _record_list.append(_record)
                    else:
                        _record = {'Assembly': '_'.join(_model_dict['SEQUENCE_ID'].split('_')[:-2]),
                                   'GenomePosition': '-'.join(_model_dict['SEQUENCE_ID'].split('_')[-2:]),
                                   'Seq': _model_dict['SEQUENCE_TEMPLATE']
                                   }
                        _record_list.append(_record)
                    _module = []
                else:
                    continue
        pd.DataFrame(_record_list).to_csv(self.output, sep='\t', index=False)

    def remove_tmp_file(self):
        os.remove(self._config)
        os.remove(self._p3out)


class ScreenSSR:
    def __init__(self, args):
        self._seq_path = args.input
        self._ssr_info = args.output
        self._assembly = args.assembly
        self._circular = args.circular
        self._threads = args.threads
        self._tmpdir = mktemp()

    def blast_pair(self):
        # make a temporary directory for BLAST
        os.mkdir(self._tmpdir)
        sequences = {_.id: _.seq for _ in SeqIO.parse(self._seq_path, 'fasta')}
        lengths = {_id: len(_seq.seq) for _id, _seq in sequences.items()}
        ssr_info = pd.read_table(self._ssr_info)
        ssr_info = ssr_info[ssr_info['seqid'].isin(self._assembly)]
        if not self._circular:
            ssr_info = ssr_info[(ssr_info['start'] > 200) &
                                (ssr_info['end'] < lengths[ssr_info['seqid']] - 200)]
        # prepare sequences
        seqlist = []
        for _idx, _item in ssr_info.iterrows():
            genome_seq = sequences[_item['seqid']]
            _sequence_id = str(genome_seq.id) + "_" + str(_item['start']) + "_" + str(_item['end'])
            _sequence_template = str(get_seq(genome_seq, _item['start'], _item['end']))
            seqlist.append('>' + _sequence_id)
            seqlist.append(_sequence_template)
        with open(os.path.join(self._tmpdir, 'query.fasta'), 'w') as f:
            f.write('\n'.join(seqlist))
            f.write('\n')
        # BLAST
        database_cmd = NcbimakeblastdbCommandline(
            dbtype='nucl',
            input_file=os.path.join(self._tmpdir, 'query.fasta'))
        blastn_cmd = NcbiblastnCommandline(
            query=os.path.join(self._tmpdir, 'query.fasta'),
            db=os.path.join(self._tmpdir, 'query.fasta'),
            dust='no',
            outfmt='\"6 qacc sacc length pident evalue\"',
            num_threads=self._threads,
            evalue='1e-3',
            out=os.path.join(self._tmpdir, 'blast_result.txt'),
            max_hsps=1)
        database_cmd()
        blastn_cmd()

    def blast_parse(self):
        blast_result = pd.read_table(os.path.join(self._tmpdir, 'blast_result.txt'),
                                     names=['query', 'subject', 'length', 'identity', 'evalue'])
        ssr_info = pd.read_table(self._ssr_info)
        # filter blast result
        query_asm = self._assembly[0]
        blast_result = blast_result[(blast_result['query'] == query_asm) & ~(blast_result['subject'] == query_asm)]
        #
        query_list = list(set(blast_result['query'].to_list()))
        keep_list = []
        # get keep list
        for _query in query_list:
            _query_table = blast_result[blast_result['query'] == _query]
            q_start, q_end = _query.split('_')[-2:]
            query_motif, query_unit = ssr_info[(ssr_info['seqid'] == _query) &
                                               (ssr_info['start'] == q_start) &
                                               (ssr_info['end'] == q_end),
                                               ['motif', 'units']]
            subject_list = [query_asm]
            subject_id_list = []
            for _idx, _row in _query_table.iterrows():
                s_id = '_'.join(_row['subject'].split('_')[:-2])
                s_start, s_end = _row['subject'].split('_')[-2:]
                subject_motif, subject_unit = ssr_info[(ssr_info['seqid'] == s_id) &
                                                       (ssr_info['start'] == s_start) &
                                                       (ssr_info['end'] == s_end),
                                                       ['motif', 'units']]
                if (subject_motif == query_motif) and (not query_unit == subject_unit):
                    subject_list.append(s_id)
                    subject_id_list.append(_row['subject'])
            if set(subject_list) == set(self._assembly):
                keep_list.append(_query)
                keep_list += subject_id_list
        # screen out keep list in ssr_info
        ssr_info['tmp_query'] = ssr_info.apply(lambda x: '_'.join([str(x['seqid']), str(x['start'], str(x['stop']))]),
                                               axis=1)
        ssr_info = ssr_info[ssr_info['tmp_query'].isin(keep_list)]
        ssr_info.drop(columns='tmp_query', inplace=True)
        # output
        ssr_info.to_csv('_screen'.join(os.path.splitext(self._ssr_info)), sep='\t', index=False)

    def remove_tmp_file(self):
        del_directory(self._tmpdir)
