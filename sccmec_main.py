import sys
sys.path.append("/Volume/biocode2/")

import os
import re
import ast
import time
import json
import random
import csv
import shutil
import pickle
import math
import pandas as pd
import numpy as np
import operator
import Levenshtein

from typing import Any, List
from matplotlib.ticker import AutoMinorLocator
from itertools import repeat, groupby
from functools import reduce
from tqdm import tqdm
from statistics import median, mean
from scipy.stats import binom, chi2, chi2_contingency

from collections import Counter, OrderedDict, defaultdict

from multiprocessing import Process, Pool, Manager
from main_classes.gbkGenome import gbkGenome, LocationTuple
from main_classes.main_module import read_file, write_file, download_genome, ensure_directory, mutate_sequence, random_seq, read_json, write_json, get_blast_record
from main_classes.motifs_library import BackgroundDist, BackgroundModel

from mdm_df import defense_finder
from mdm_df import defense_finder_posttreat

import subprocess

from matplotlib.lines import Line2D

import matplotlib.colors
import matplotlib.pyplot as plt
from matplotlib.container import BarContainer
from matplotlib.ticker import FuncFormatter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqRecord import SeqRecord
from Bio import Entrez

np.set_printoptions(linewidth=np.inf, threshold=np.inf, suppress=True)

pd.set_option('display.expand_frame_repr', False)

pd.set_option("display.max_rows", 6000, 
              "display.max_columns", 6000, 
              "display.max_colwidth", 30)

path_pfamdb = "/Volume/biodata2/PfamDB/PfamDB_35.0/Pfam-A.hmm"

pm = "/Volume/biodata2/"
pr = pm + "sccmec_snail_analysis/"
pp = pm + "sccmec_snail_analysis/ANALYSIS_RESULTS_FOR_PAPER/"
pg = pm + "genbank_collection/"

Entrez.email = '---'

def load_dvars(pr):

    with open(pr + 'dvars.json', 'r') as fp:

        dvars = json.load(fp)

    return dvars

def calculate_nucleotide_distance(lt_1: LocationTuple, 
                                  lt_2: LocationTuple, 
                                  sequence: Seq):

    a = lt_1[0]
    b = lt_1[1]
    c = lt_2[0]
    d = lt_2[1]
    M = len(sequence)-1

    if lt_1[2] == 1: orig_dist = d - a
    else: orig_dist = b-c

    if orig_dist < 0: span = orig_dist + M
    else: span = orig_dist

    return span

def multiprocess_with_pbar(funct, main_list, n_processes = 64, description = 'None'):

    def listener(q_tempory, q_results, total_tasks):

        pbar = tqdm(total=total_tasks, desc=description, position=0, leave=True)
        for item in iter(q_tempory.get, None):
            pbar.update()
            q_results.put(item)

    processes = []
    q_tempory = Manager().Queue()
    q_results = Manager().Queue()

    total_tasks = len(main_list)
    sub_lists = np.array_split(main_list, n_processes)  

    proc = Process(target=listener, args=(q_tempory, q_results, total_tasks))
    proc.start()

    for a_list in sub_lists:
        P = Process(target=funct, args=(a_list, q_tempory))
        processes.append(P)

    for p in processes: p.start()

    for p in processes: p.join()

    q_tempory.put(None) 

    proc.join()

    results = [q_results.get() for _ in range(q_results.qsize())]

    return results

def levenshtein_tester():

    m = 1000
    n = int(0.05 * m)

    s1 = random_seq(m)
    s2 = mutate_sequence(n, s1, 'mutate')

    r = Levenshtein.ratio(s1, s2)
    print(r)

def lv3(query, SecRec):

    
    r = Levenshtein.ratio(query, SecRec.seq)

    if r > 0.95:
        return SecRec.id, 1
    else:
        return SecRec.id, 0

def levenshtein_searcher(query, dynamic_acc_pos_list):

    p = Pool(processes=64) 
    processing_list = zip(repeat(query), dynamic_acc_pos_list)

    st = time.perf_counter()
    res = p.starmap(lv3, processing_list)
    ed = time.perf_counter()

    hits = [t[0] for t in res if t[1] == 1]

    print(len(res), len(dynamic_acc_pos_list), len(hits), end=" - ")
    

    return hits

def message(*args, **kwargs):

    global message_enabled

    if message_enabled: print(*args, **kwargs)

def run_defense_finder(path_ngbrs, path_defin, outdir, 
                       dbtype='ordered_replicon', workers=16, 
                       preserve_raw=False):

    def build_systems(genes):

        system_groups = [list(it) for k, it in groupby(genes, lambda val: val['sys_id'])]
        systems = []

        for system_group in system_groups:

            item = {}
            first_item = system_group[0]
            last_item = system_group[-1]
            item['sys_id'] = first_item['sys_id']
            item['sys_beg'] = first_item['hit_id']
            item['sys_end'] = last_item['hit_id']
            item['type'] = first_item['type']
            item['subtype'] = first_item['subtype']
            item['protein_in_syst'] = reduce(lambda acc, s: acc + ',' + s['hit_id'], system_group, '')[1:]
            item['genes_count'] = len(system_group)
            item['name_of_profiles_in_sys'] = reduce(lambda acc, s: acc + ',' + s['gene_name'], system_group, '')[1:]
            systems.append(item)

        return systems

    

    tmp_dir = os.path.join(outdir, os.path.basename(path_ngbrs).split('.fa')[0])

    if os.path.exists(tmp_dir): shutil.rmtree(tmp_dir)

    os.makedirs(tmp_dir)

    with open(path_ngbrs) as f: defense_finder.run(f, dbtype, workers, tmp_dir)

    

    defense_finder_genes = defense_finder_posttreat.best_solution.get(tmp_dir)
    systems = build_systems(defense_finder_genes)

    genes_dict = dict()
    for item in defense_finder_genes:

        genes_dict[item['hit_id']] = item

    data = {'genes': genes_dict, 'systems': systems}

    with open(path_defin, 'wb') as f: pickle.dump(data, f)

    

def mixed_defense_mapper():

    def parse(path):

        def tokenize_line(line):

            tokens = line.split()
            tokenized_line = tokens[0:22]
            tokenized_line.append(' '.join(tokens[22:]))

            return tokenized_line

        columns = ['Target Name',
                   'accession1',
                   'tlen',
                   'query name',
                   'accession2',
                   'qlen',
                   'seq E-Value',
                   'seq score',
                   'seq bias',
                   'dom #',
                   'dom of',
                   'dom c-Evalue',
                   'dom i-Evalue',
                   'dom score',
                   'dom bias',
                   'hmm from',
                   'hmm to',
                   'ali from',
                   'ali to',
                   'env from',
                   'env to',
                   'acc',
                   'description of target']

        raw_text = read_file(path)
        lines = raw_text.splitlines()
        lines = [line for line in lines if not line.startswith('#')]

        tokenized_neigbr_lines = []

        for line in lines: tokenized_neigbr_lines.append(tokenize_line(line))

        df = pd.DataFrame(tokenized_neigbr_lines, columns=columns)

        df['accession1'] = df['accession1'].apply(lambda s: s.split('.')[0])

        return df

    def extract_all_features(path_input, path_output, pr=''):

        inp = pr + "rp62a_mixed_defense_map/RP62A.gb"
        out = pr + "rp62a_mixed_defense_map/rp62a.faa"

        path_input = inp
        path_output = out

        record: SeqRecord = SeqIO.read(path_input, "gb")  

        string = ''

        for feature in record.features:

            quals = feature.qualifiers

            if 'translation' in feature.qualifiers.keys():

                loc = LocationTuple(feature.location.start, feature.location.end,
                                    feature.location.strand).encode_location_tuple()

                if 'locus_tag' in quals.keys(): ltg = quals['locus_tag'][0]
                else: ltg = 'no_locus_tag'

                if 'protein_id' in quals.keys(): pid = quals['protein_id'][0]
                else: pid = 'no_pid'

                if 'product' in quals.keys(): prd = quals['product'][0]
                else: prd = 'no_product'

                trl = quals['translation'][0]

                s = '>{} {} {} {}\n{}\n'.format(loc, ltg, pid, prd, trl)
                string += s

        write_file(path_output, string)

    def main(pr=''):

        path_dfhmm = pr + "rp62a_mixed_defense_map/def_finder_results/defense_finder_hmmer.tsv"
        path_pfams = pr + "rp62a_mixed_defense_map/[NC_002976.3]_pfams_[na][b].txt"
        path_ngbrs = pr + "rp62a_mixed_defense_map/[NC_002976.3]_ngbrs_[na][b].fa"
        path_defrl = pr + "rp62a_mixed_defense_map/active.csv"
        path_outdf = pr + "rp62a_mixed_defense_map/output.csv"

        ccr_pfams_set = {"PF00239", "PF07508", "PF13408"}

        df_defrel = pd.read_csv(path_defrl)
        defrel_set = set(df_defrel['Pfam'].to_list())
        defrel_dic = df_defrel.set_index('Pfam').T.to_dict('list')

        df = parse(path_pfams)

        df_hmm = pd.read_table(path_dfhmm)

        rows = []

        null_value = np.nan

        for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):

            desc = seq_record.description.replace(seq_record.id, '').strip()

            
            ddf = df[df['query name'] == seq_record.id]
            pfam_set = set(ddf['accession1'].to_list())
            dpfam_set = defrel_set.intersection(pfam_set)

            len_dpfam_set = len(dpfam_set)

            ccr = 1 if ccr_pfams_set.issubset(pfam_set) else null_value
            sub_dictionary = dict((k, defrel_dic[k]) for k in dpfam_set if k in defrel_dic)

            

            ddh = df_hmm[df_hmm['hit_id'] == seq_record.id]
            hmm_set = set(ddh['gene_name'].to_list())

            len_hmm_set = len(hmm_set)

            

            if len_dpfam_set > 0 or len_hmm_set > 0: d = 1
            else: d = null_value

            

            a = len_dpfam_set
            b = len_hmm_set
            c = a + b

            

            st, ed, strand = LocationTuple.decode(seq_record.id)

            

            mobility = ""

            if "integrase" in desc: mobility += "I,"
            if "recombinase" in desc: mobility += "R,"
            if "transposase" in desc: mobility += "T,"

            

            
            rows.append([i, st, ed, strand, d, c, b, a, ccr, mobility, desc, repr(sub_dictionary), repr(hmm_set)])

        cols = ['i',
                'Start',
                'End',
                'Strand',
                'Defense Related',
                '
                '
                '
                'Contains CCR domains',
                'Mobility',
                'Description',
                'dpfam_set',
                'hmm_set']

        df = pd.DataFrame(rows, columns=cols)
        df.set_index('i', inplace=True)
        df.to_csv(path_outdf, sep=',')

    main()

class CcrDomainConfirmation:

    @staticmethod
    def ccr_domain_confirmation_1(pr):

        path = pr + "snail_analysis_epidermidis/backups/neighborhood_files_ws500_blast_version/"
        

        files = os.listdir(path)

        files = [file for file in files if "_dfcsv_" in file and not file.startswith('.')]
        pnames = []

        for file in files:

            print()
            print(file, ' --------------')

            p = path + file
            df = pd.read_csv(p)

            df = df[df['contains_ccr_domains'] == 1]
            descriptions = df['description'].to_list()

            for desc in descriptions:

                m = re.search(r'[0-9]{2}\.[0-9] +(.+)$', desc)
                pname = m.group(1)
                pnames.append(pname)

                

        print('\n\n', '#' * 80, '\n')

        print(len(pnames))
        pnames = set(pnames)
        for d in pnames: print(d)

    @staticmethod
    def ccr_domain_confirmation_2(pr):

        path = pr + "snail_analysis_epidermidis/neighborhood_files_ws500/"

        files = os.listdir(path)

        files = [file for file in files if "_dfcsv_" in file and not file.startswith('.')]
        pnames = []

        for file in files:

            print()
            print(file, ' --------------')

            p = path + file
            df = pd.read_csv(p)

            df = df[df['defense_related'] == 1]
            dpfam_set = df['dpfam_set'].to_list()
            descriptions = df['description'].to_list()
            ii = df['i'].to_list()

            for i, desc in enumerate(descriptions):

                if ii[i] <= 200:

                    

                    d = ast.literal_eval(dpfam_set[i])

                    msg = f"{ii[i]} {d['subtype']} {d['hit_profile_cov']} - {d['hit_seq_cov']} {desc} "
                    print(msg)

class Meme:

    def __init__(self, p_ins, p_out, p_bck, w = 18):
  
        self.p_ins = p_ins
        self.p_out = p_out
        self.p_bck = p_bck
        self.w = w

        self.p_result = p_out + "/meme.txt"

    def step1_generate_motif(self):

        print(f"Running MEME on {self.p_ins}")

        cmd = f"meme {self.p_ins} -dna " \
            f"-oc {self.p_out} " \
            f"-mod zoops -nmotifs 1 -w {self.w} " \
            f"-markov_order 0 " \
            f"-bfile {self.p_bck}"

        print(cmd)

        subprocess.check_call(cmd, shell=True)

    def step2_meme2meme_conversion(self):

        out = self.p_out + "/meme2meme.txt"

        cmd = f"meme2meme {self.p_result} > {out}"
         
        subprocess.check_call(cmd, shell=True)
    
    def step3_run_FIMO(self, target_sequence, p_fimo_out, pseudo_count, qv_thresh=0.8):

        if os.path.isfile(p_fimo_out + "/fimo.tsv"):
            print("FIMO already ran. Skipping...")

        else:

            path_motif = self.p_out + "/meme2meme.txt"

            cmd = f"fimo " \
                  f"--oc {p_fimo_out} " \
                  f"--motif-pseudo {pseudo_count} " \
                  f"--bfile motif-file " \
                  f"--thresh {qv_thresh} " \
                  f"--qv-thresh " \
                  f"--max-stored-scores 10000000 " \
                  f"{path_motif} " \
                  f"{target_sequence}"
            
            subprocess.check_call(cmd, shell=True)

        df = pd.read_csv(p_fimo_out + "/fimo.tsv", sep='\t', comment='#')
        return df

    
        
    def parse(self, mode="log-odds"):
        
        if mode == "log-odds":
            string = "log-odds matrix:"

        elif mode == "pwm":
            string = "letter-probability matrix:"

        
        
        with open(self.p_result, 'r') as file:
            lines = file.readlines()
            
        matrix = []
        parse = False

        for line in lines:

            if line.startswith(string):
                parse = True
                continue
            if line.startswith("----"):
                if parse:
                    break
                else:
                    continue
            if parse:
                row = list(map(int, line.split()))
                matrix.append(row)
        
        df = pd.DataFrame(matrix, columns=["A", "C", "G", "T"])
        return df

class MappedCsv:

    def __init__(self, path):
        self.path = path

    def parse(self):

        df = pd.read_csv(self.path)
        df_lookup = dict()

        for index in df.index:
            accession = df.loc[index]['Genome Accession']
            description = df.loc[index]['Description']
            context = df.loc[index]['Context']
            df_lookup[accession] = [description, context]

        return df, df_lookup, len(df)

    
    def write(self):
        pass

    
    @classmethod
    def create_from_from_bhl(cls, bhl):
        pass

class AccessionObj(OrderedDict):

    def __init__(self, accession=None):

        super().__init__()

        self.accession = accession
        self.valid_keys = ['ngbrs', 'allpr', 'pfams', 'defin', 
                           'alsys', 'blast', 'cored', 'coral', 
                           'dtfrm', 'dfcsv']

        for identifier in self.valid_keys:

            global object_types

            identifier_with_underscores = '_' + identifier + '_'
            Obj = object_types[identifier_with_underscores](parent=self)
            self[identifier] = Obj

    def __setitem__(self, key, value):

        if key not in self.valid_keys:
            raise KeyError("Key <{}> is not allowed".format(repr(key)))

        super().__setitem__(key, value) 

    def __str__(self):

        txt = "AccessionObj({})\n".format(self.accession)

        for k in [key for key in self.valid_keys if key in self.keys()]:

            txt += '\t{}: {}\n'.format(k, repr(self[k]))

        return txt

    def add_xObject(self, Obj):

        if Obj.accession != self.accession:

            raise ValueError('accession numbers do not match')

        Obj.parent = self

        if isinstance(Obj, xObjects.NgbrsObj): self['ngbrs'] = Obj
        if isinstance(Obj, xObjects.AllPrObj): self['allpr'] = Obj
        if isinstance(Obj, xObjects.PfamsObj): self['pfams'] = Obj
        if isinstance(Obj, xObjects.DeFinObj): self['defin'] = Obj
        if isinstance(Obj, xObjects.AlSysObj): self['alsys'] = Obj
        if isinstance(Obj, xObjects.BlastObj): self['blast'] = Obj
        if isinstance(Obj, xObjects.CoredObj): self['cored'] = Obj
        if isinstance(Obj, xObjects.CoralObj): self['coral'] = Obj
        if isinstance(Obj, xObjects.DtfrmObj): self['dtfrm'] = Obj
        if isinstance(Obj, xObjects.DfCsvObj): self['dfcsv'] = Obj

    def debug_print(self):

        print(f'Accession: {self.accession} \n' )
        for idx, Obj in self.items():
            print('\t', idx, Obj.identifier, Obj.parent.accession, Obj.path_full)

class xObjects:

    class BaseObj:

        def __init__(self, *args, **kwargs):

            if 'parent' in kwargs:

                self.parent = kwargs['parent']

            else:

                self.parent = None

            if any(args):

                self.path_full = args[0]
                self.file_name = os.path.basename(self.path_full)

                if {'acc', 'lcs', 'cir', 'identifier'}.issubset(kwargs.keys()):
                    self.accession = kwargs['acc']
                    self.location_tuple = kwargs['lcs']
                    self.cir = kwargs['cir']
                    self.identifier = kwargs['identifier']

                else:

                    a, l, c, i = FileName.parse(self.file_name)
                    self.accession = a
                    self.location_tuple = l
                    self.cir = c
                    self.identifier = i

            else:

                self.path_full = None
                self.file_name = None
                self.accession = None
                self.location_tuple = None
                self.cir = None
                self.identifier = None

        @property
        def exists(self):

            if self.path_full is None: return False
            else: return True

        def reinitialize_from_path(self, path_full):

            self.path_full = path_full
            self.file_name = os.path.basename(self.path_full)
            a, l, c, i = FileName.parse(self.file_name)
            self.accession = a
            self.location_tuple = l
            self.cir = c
            self.identifier = i

    class NgbrsObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "NgbrsObj -> {}".format(self.path_full)

        @classmethod
        def create(cls, path_file, path_genomes, accession, location, window_size, displacement):

            

            genome = gbkGenome(path_genomes, accession)
            location_tuple = LocationTuple.decode(location)

            output_features, suffix, circle_indicator = genome.get_neighbours_of_rlmh(location_tuple,
                                                                                      window_size=window_size,
                                                                                      displacement=displacement)
            file_name, _ = genome.convert2fasta(path_file,
                                                output_features,
                                                suffix,
                                                circle_indicator)

            return cls(path_file + file_name)

        def parse(self):

            return list(SeqIO.parse(self.path_full, "fasta"))

    class AllPrObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "AllPrObj -> {}".format(self.path_full)

        @classmethod
        def create(cls, path_file, path_genomes, accession, location_tuple, displacement):

            genome = gbkGenome(path_genomes, accession)

            
            output_features, suffix, circle_indicator = genome.get_neighbours_of_rlmh(location_tuple,
                                                                                      window_size=0,
                                                                                      displacement=displacement)

            file_name, _ = genome.convert2fasta(path_file, output_features, suffix, circle_indicator, identifier='allpr')

            return cls(path_file + file_name)

        def parse(self):

            return list(SeqIO.parse(self.path_full, "fasta"))

    class PfamsObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            
            

            return "PfamsObj -> {}".format(self.path_full)

        
        @classmethod
        def create(cls, path_ngbrs, path_pfams):

            command = ["hmmscan",
                       "--cut_ga",
                       "--noali",
                       "--acc",
                       "--domtblout",
                       path_pfams,
                       path_pfamdb,
                       path_ngbrs]

            proc = subprocess.Popen(command, stdout=subprocess.DEVNULL) 
            proc.wait()

            return cls(path_pfams)

        def parse(self):

            def tokenize_line(line):

                tokens = line.split()
                tokenized_line = tokens[0:22]
                tokenized_line.append(' '.join(tokens[22:]))

                return tokenized_line

            columns = ['Target Name',
                       'accession1',
                       'tlen',
                       'query name',
                       'accession2',
                       'qlen',
                       'seq E-Value',
                       'seq score',
                       'seq bias',
                       'dom #',
                       'dom of',
                       'dom c-Evalue',
                       'dom i-Evalue',
                       'dom score',
                       'dom bias',
                       'hmm from',
                       'hmm to',
                       'ali from',
                       'ali to',
                       'env from',
                       'env to',
                       'acc',
                       'description of target']

            raw_text = read_file(self.path_full)
            lines = raw_text.splitlines()
            lines = [line for line in lines if not line.startswith('#')]

            tokenized_neigbr_lines = []
            tokenized_origin_lines = []

            if self.location_tuple == (0, 0, 0):

                for line in lines:
                    tokenized_line = tokenize_line(line)
                    tokenized_neigbr_lines.append(tokenized_line)

            else:

                for line in lines:

                    tokenized_line = tokenize_line(line)

                    if self.location_tuple.encode_location_tuple() in line:
                        tokenized_origin_lines.append(tokenized_line)
                    else:
                        tokenized_neigbr_lines.append(tokenized_line)

            df_origin = pd.DataFrame(tokenized_origin_lines, columns=columns)
            df_neigbr = pd.DataFrame(tokenized_neigbr_lines, columns=columns)

            df_origin['accession1'] = df_origin['accession1'].apply(lambda s: s.split('.')[0])
            df_neigbr['accession1'] = df_neigbr['accession1'].apply(lambda s: s.split('.')[0])

            return df_origin, df_neigbr

    class DeFinObj(BaseObj):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "DeFinObj -> {}".format(self.path_full)

        @classmethod
        def create(cls, path_ngbrs, path_defin, outdir):

            run_defense_finder(path_ngbrs, path_defin, outdir) 

            return cls(path_defin)

        def parse_0(self):

            with open(self.path_full, 'rb') as f: data = pickle.load(f)

            genes = data['genes']
            systems = data['systems']

            print('---------- GENES ----------')
            for pos, gene in genes.items():

                print('\t', pos)
                for k,v in gene.items(): print('\t\t', k, ": ", v)
                print()

            print('\n\n')

            print('---------- SYSTEMS ----------')
            for system in systems:

                for k,v in system.items(): print('\t', k, ": ", v)
                print()

        def parse_1(self):

            with open(self.path_full, 'rb') as f: data = pickle.load(f)

            genes = data['genes']
            systems = data['systems']

            seqnums = []
            for pos, gene in genes.items(): seqnums.append(gene['hit_pos'])

            seqnums.sort()

            print('\t', seqnums)

        def parse(self):

            with open(self.path_full, 'rb') as f: data = pickle.load(f)

            genes = data['genes']
            systems = data['systems']

            return genes, systems

    class AlSysObj(BaseObj):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "AlSysObj -> {}".format(self.path_full)

        @classmethod
        def create(cls, path_allpr, path_alsys, outdir):

            run_defense_finder(path_allpr, path_alsys, outdir)

            return cls(path_alsys)

        def parse(self):

            with open(self.path_full, 'rb') as f: data = pickle.load(f)

            genes = data['genes']
            systems = data['systems']

            return genes, systems

    class BlastObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "BlastObj -> {}".format(self.path_full)

        def create(self):
            pass

        def parse2df(self):

            rows = []
            for line in read_file(self.path_full).splitlines():

                f = line.split()

                row = [f[0],
                       f[1],
                       float(f[2]),
                       float(f[3]),
                       float(f[4]),
                       float(f[5]),
                       float(f[6]),
                       float(f[7]),
                       float(f[8]),
                       float(f[9])]

                rows.append(row)

            columns = ['qseqid',
                       'sseqid',
                       'evalue',
                       'qlen',
                       'slen',
                       'pident',
                       'bitscore',
                       'qframe',
                       'sframe',
                       'qcovs']

            return pd.DataFrame(rows, columns = columns)

        def parse2gen(self):

            file = open(self.path_full, "r")
            for line in file:

                f = line.split()
                qseqid, sseqid, evl, pident, qcov = f[0], f[1], float(f[2]), float(f[5]), float(f[9])

                yield qseqid, sseqid, evl, pident, qcov

    class CoredObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "CoredObj -> {}".format(self.path_full)

        
        def create_with_df(self):

            if self.parent['ngbrs'].exists and self.parent['blast'].exists:

                if self.parent['cored'].exists:

                    raise ValueError('Cored already exists')

                path_ngbrs = self.parent['ngbrs'].path_full
                path_cored = FileName.convert(path_ngbrs, conv='ngbrs->cored')

                core_gene_score_dictionary = OrderedDict()
                positions = []
                for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):
                    positions.append(seq_record.id)

                for pos in positions:
                    core_gene_score_dictionary[pos] = [set(), 0]

                

                st = time.perf_counter()
                df = self.parent['blast'].parse()
                ed = time.perf_counter()
                print('Finished in {} seconds:'.format(ed-st))

                for ind in df.index:

                    sseqid = df.loc[ind]['sseqid']
                    qseqid = df.loc[ind]['qseqid']
                    qcov   = df.loc[ind]['qcovs']
                    pident = df.loc[ind]['pident']
                    evalue = df.loc[ind]['evalue']

                    if qcov > 95 and pident > 95 and evalue < 1e-6:
                        core_gene_score_dictionary[qseqid][0].add(sseqid.split('|')[0])

                

                for pos in positions:
                    core_gene_score_dictionary[pos][1] = len(core_gene_score_dictionary[pos][0])

                with open(path_cored, 'wb') as f:
                    pickle.dump(core_gene_score_dictionary, f)

                self.reinitialize_from_path(path_cored)

            else:

                raise ValueError('Ngbrs or Blast does not exist!')

        
        def create_with_gen(self):

            if self.parent['ngbrs'].exists and self.parent['blast'].exists:

                if self.parent['cored'].exists:

                    raise ValueError('Cored already exists')

                path_ngbrs = self.parent['ngbrs'].path_full
                path_cored = FileName.convert(path_ngbrs, conv='ngbrs->cored')

                positions = []
                for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):
                    positions.append(seq_record.id)

                core_gene_score_dictionary = OrderedDict()
                for pos in positions:
                    core_gene_score_dictionary[pos] = [set(), 0]

                

                for qseqid, sseqid, evalue, pident, qcov in self.parent['blast'].parse2gen():

                    if qcov > 95 and pident > 95 and evalue < 1e-6:
                        core_gene_score_dictionary[qseqid][0].add(sseqid.split('|')[0])

                

                for pos in positions:
                    core_gene_score_dictionary[pos][1] = len(core_gene_score_dictionary[pos][0])

                with open(path_cored, 'wb') as f:
                    pickle.dump(core_gene_score_dictionary, f)

                self.reinitialize_from_path(path_cored)

            else:

                raise ValueError('Ngbrs or Blast does not exist!')

        def create(self):

            
            
            pass

        def parse(self):

            with open(self.path_full, 'rb') as f: genes = pickle.load(f)
            return genes

    class CoralObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "CoralObj -> {}".format(self.path_full)

        
        def create_with_df(self):

            if self.parent['ngbrs'].exists and self.parent['blast'].exists:

                if self.parent['cored'].exists:

                    raise ValueError('Cored already exists')

                path_ngbrs = self.parent['ngbrs'].path_full
                path_cored = FileName.convert(path_ngbrs, conv='ngbrs->cored')

                core_gene_score_dictionary = OrderedDict()
                positions = []
                for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):
                    positions.append(seq_record.id)

                for pos in positions:
                    core_gene_score_dictionary[pos] = [set(), 0]

                

                st = time.perf_counter()
                df = self.parent['blast'].parse()
                ed = time.perf_counter()
                print('Finished in {} seconds:'.format(ed-st))

                for ind in df.index:

                    sseqid = df.loc[ind]['sseqid']
                    qseqid = df.loc[ind]['qseqid']
                    qcov   = df.loc[ind]['qcovs']
                    pident = df.loc[ind]['pident']
                    evalue = df.loc[ind]['evalue']

                    if qcov > 95 and pident > 95 and evalue < 1e-6:
                        core_gene_score_dictionary[qseqid][0].add(sseqid.split('|')[0])

                

                for pos in positions:
                    core_gene_score_dictionary[pos][1] = len(core_gene_score_dictionary[pos][0])

                with open(path_cored, 'wb') as f:
                    pickle.dump(core_gene_score_dictionary, f)

                self.reinitialize_from_path(path_cored)

            else:

                raise ValueError('Ngbrs or Blast does not exist!')

        
        def create_with_gen(self):

            if self.parent['ngbrs'].exists and self.parent['blast'].exists:

                if self.parent['cored'].exists:

                    raise ValueError('Cored already exists')

                path_ngbrs = self.parent['ngbrs'].path_full
                path_cored = FileName.convert(path_ngbrs, conv='ngbrs->cored')

                positions = []
                for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):
                    positions.append(seq_record.id)

                core_gene_score_dictionary = OrderedDict()
                for pos in positions:
                    core_gene_score_dictionary[pos] = [set(), 0]

                

                for qseqid, sseqid, evalue, pident, qcov in self.parent['blast'].parse2gen():

                    if qcov > 95 and pident > 95 and evalue < 1e-6:
                        core_gene_score_dictionary[qseqid][0].add(sseqid.split('|')[0])

                

                for pos in positions:
                    core_gene_score_dictionary[pos][1] = len(core_gene_score_dictionary[pos][0])

                with open(path_cored, 'wb') as f:
                    pickle.dump(core_gene_score_dictionary, f)

                self.reinitialize_from_path(path_cored)

            else:

                raise ValueError('Ngbrs or Blast does not exist!')

        def parse(self):

            with open(self.path_full, 'rb') as f: genes = pickle.load(f)
            return genes

    class DtfrmObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "DtfrmObj -> {}".format(self.path_full)

        
        def create_v0(self, df_defrel, save=True):

            

            

            loc_rlmh             = self.parent['pfams'].location_tuple.encode_location_tuple()
            df_origin, df_neigbr = self.parent['pfams'].parse()
            genes                = self.parent['cored'].parse()
            path_ngbrs           = self.parent['ngbrs'].path_full
            path_dtfrm = FileName.convert(path_ngbrs, conv='ngbrs->dtfrm')
            path_dfcsv = FileName.convert(path_ngbrs, conv='ngbrs->dfcsv')

            ccr_pfams_set = {"PF00239", "PF07508", "PF13408"}

            tcov_min = 0.8  
            qcov_min = 0.1  
            eval_max = 1e-6
            null_value = np.nan

            df = df_neigbr.astype( {"hmm to": int,
                                    "hmm from": int,
                                    "ali to": int,
                                    "ali from": int,
                                    "tlen": int,
                                    "qlen": int,
                                    "seq E-Value": float})

            df['tcov'] = (df['hmm to'] - df['hmm from']) / df['tlen']
            df['qcov'] = (df['ali to'] - df['ali from']) / df['qlen']

            defrel_set = set(df_defrel['Pfam'].to_list())
            defrel_dic = df_defrel.set_index('Pfam').T.to_dict('list')

            rows = []
            nhis = []

            

            

            for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):

                desc = seq_record.description.replace(seq_record.id, '').strip()
                ddf = df[df['query name'] == seq_record.id]
                pfam_set = set(ddf['accession1'].to_list())
                dpfam_set = defrel_set.intersection(pfam_set)

                
                len_dpfam_set = len(dpfam_set)

                if len_dpfam_set > 0:

                    if len_dpfam_set == 1:

                        pf = list(dpfam_set)[0]

                        op = ddf[(ddf['accession1'] == pf) & (ddf['tcov'] > tcov_min) & (ddf['qcov'] > qcov_min) & (
                                    ddf['seq E-Value'] < eval_max)]

                        if len(op) > 0:
                            d = 1
                        else:
                            d = null_value

                    else:

                        d = 1
                else:
                    d = null_value

                ccr = 1 if ccr_pfams_set.issubset(pfam_set) else null_value

                if seq_record.id == loc_rlmh: rlmh = 1 
                else: rlmh = null_value

                
                if "PF09848" in pfam_set:  

                    nhis.append([self.accession, desc, seq_record.id, str(pfam_set), str(len(seq_record.seq)),
                                 str(seq_record.seq)])
                    dpfam_set = {"PF09848"}
                    dus = 1
                    d = 1

                else:

                    dus = null_value

                sub_dictionary = dict((k, defrel_dic[k]) for k in dpfam_set if k in defrel_dic)

                rows.append([i, seq_record.id, d, ccr, rlmh, dus, desc, repr(sub_dictionary)])

            

            

            

            cols = ['i', 'pos', 'defense_related', 'contains_ccr_domains', 'is_rlmh', 'is_dusA', 'description', 'dpfam_set']
            df = pd.DataFrame(rows, columns=cols)
            df.set_index('pos', inplace=True)
            df['gene_conservation'] = 0

            for ind in df.index: df.at[ind, 'gene_conservation'] = genes[ind][1]

            if save:

                with open(path_dtfrm, 'wb') as f: pickle.dump(df, f)
                df.to_csv(path_dfcsv, sep=',')

            else:

                print(df)

            

        
        def create_v1(self, save=True):

            

            

            loc_rlmh             = self.parent['pfams'].location_tuple.encode_location_tuple()
            df_origin, df_neigbr = self.parent['pfams'].parse()
            genes                = self.parent['cored'].parse()
            path_ngbrs           = self.parent['ngbrs'].path_full
            path_dtfrm = FileName.convert(path_ngbrs, conv='ngbrs->dtfrm')
            path_dfcsv = FileName.convert(path_ngbrs, conv='ngbrs->dfcsv')

            ccr_pfams_set = {"PF00239", "PF07508", "PF13408"}

            tcov_min = 0.8  
            qcov_min = 0.1  
            eval_max = 1e-3  
            null_value = np.nan

            df = df_neigbr.astype( {"hmm to": int,
                                    "hmm from": int,
                                    "ali to": int,
                                    "ali from": int,
                                    "tlen": int,
                                    "qlen": int,
                                    "seq E-Value": float,
                                    "dom i-Evalue": float})

            df['tcov'] = (df['hmm to'] - df['hmm from']) / df['tlen']
            df['qcov'] = (df['ali to'] - df['ali from']) / df['qlen']
            df['pass'] = (df['tcov'] > tcov_min) & (df['qcov'] > qcov_min) & (df['dom i-Evalue'] <= eval_max)

            rows = []

            

            

            for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):

                desc = seq_record.description.replace(seq_record.id, '').strip()
                ddf = df[df['query name'] == seq_record.id]

                all_domains = ddf['Target Name'].to_list()
                def_domains, ccr_domains = [], []

                for domain in all_domains:
                    if '|' in domain: def_domains.append(domain)
                    else: ccr_domains.append(domain)

                
                len_def_domains = len(def_domains)

                if len_def_domains > 0:

                    if len_def_domains == 1:

                        op = ddf[(ddf['Target Name'] == def_domains[0]) & (ddf['pass'] == True)]

                        if len(op) > 0:d = 1
                        else:d = null_value

                    else:

                        d = 1
                else:
                    d = null_value

                
                ccr = 1 if len(set(ccr_domains)) >= 3 else null_value

                if seq_record.id == loc_rlmh: rlmh = 1 
                else: rlmh = null_value

                
                

                row = [i, seq_record.id, d, ccr, rlmh, null_value, desc, repr(def_domains)]
                
                rows.append(row)

            

            

            

            cols = ['i', 'pos', 'defense_related', 'contains_ccr_domains', 'is_rlmh', 'is_dusA', 'description', 'dpfam_set']
            df = pd.DataFrame(rows, columns=cols)
            df.set_index('pos', inplace=True)
            df['gene_conservation'] = 0

            for ind in df.index: df.at[ind, 'gene_conservation'] = genes[ind][1]

            if save:

                with open(path_dtfrm, 'wb') as f: pickle.dump(df, f)
                df.to_csv(path_dfcsv, sep=',')

            else:

                print(df)

            

        def create(self, save=True):

            loc_rlmh    = self.parent['defin'].location_tuple.encode_location_tuple()
            genes, systems = self.parent['defin'].parse()

            cored_genes = self.parent['cored'].parse()
            path_ngbrs  = self.parent['ngbrs'].path_full
            path_dtfrm  = FileName.convert(path_ngbrs, conv='ngbrs->dtfrm')
            path_dfcsv  = FileName.convert(path_ngbrs, conv='ngbrs->dfcsv')

            null_value = np.nan
            rows = []

            for i, seq_record in enumerate(SeqIO.parse(path_ngbrs, "fasta")):

                desc = seq_record.description.replace(seq_record.id, '').strip()

                mm = re.search(r"recombinase\s+family\s+protein", desc, re.IGNORECASE)
                ccr = 1 if mm is not None else null_value

                pos = seq_record.id

                if pos == loc_rlmh: rlmh = 1
                else: rlmh = null_value

                if pos in genes.keys():
                    def_related = 1
                    def_descrpt = genes[pos]
                else:
                    def_related = null_value
                    def_descrpt = null_value

                row = [i, pos, def_related, ccr, rlmh, null_value, desc, def_descrpt, cored_genes[pos][1]]
                rows.append(row)

            cols = ['i', 'pos', 'defense_related', 'contains_ccr_domains', 'is_rlmh', 'is_dusA', 'description', 'dpfam_set', 'gene_conservation']
            df = pd.DataFrame(rows, columns=cols)
            df.set_index('pos', inplace=True)

            if save:

                with open(path_dtfrm, 'wb') as f: pickle.dump(df, f)
                df.to_csv(path_dfcsv, sep=',')

            else:

                print(df)

        def parse(self):

            with open(self.path_full, 'rb') as f: df = pickle.load(f)
            return df

    class DfCsvObj(BaseObj):

        def __init__(self, *args, **kwargs):

            super().__init__(*args, **kwargs)

        def __repr__(self):

            return "DfCsvObj -> {}".format(self.path_full)

        def create(self):
            pass

        def parse(self):
            pass

    @staticmethod
    def auto_create_from_file(path_full, **kwargs):

        global object_types

        identifier = re.search("](_[npabcd][a-z]{4}_)\[", path_full).group(1)
        return object_types[identifier](path_full, **kwargs)

class FileName:

    @staticmethod
    def convert(filename, conv=''):

        global file_types

        if conv == 'get_accession':

            return filename.split(']')[0].strip('[')

        else:

            fr, to = conv.split('->')

            fr_with_underscores = "_%s_" % fr
            to_with_underscores = "_%s_" % to
            old_ext = file_types[fr_with_underscores]
            new_ext = file_types[to_with_underscores]

            assert filename.endswith(old_ext), "Old filename doesn't end with proper extension..."

            return filename.replace(fr_with_underscores, to_with_underscores).replace(old_ext, new_ext)

    @staticmethod
    def parse(filename):

        message('Parsing Filename: {}'.format(filename))  

        global file_types

        m = re.search("^\[(.+)]_([a-z]{5})_\[(.+)]\[([bonp])](\..{2,3})$", filename)

        if m is not None:

            acc = m.group(1)
            idf = m.group(2)
            lcs = m.group(3)
            cir = m.group(4)
            ext = m.group(5)

            idf_with_underscores = "_%s_" % idf

            assert (idf_with_underscores, ext) in file_types.items(), "File identifier does not match extension"

            return acc, LocationTuple.decode(lcs), cir, idf_with_underscores

        else:

            return None, None, None, None

class nbFileManager:

    def __init__(self, path, path_genomes):

        self.path = path
        self.path_genomes = path_genomes
        self.accessions = set() 
        self.AccessionObjDict = self.generate_AccessionObjDict()
        self.GenomeObjDict    = self.generate_GenomeObjDict()

    def __str__(self):

        txt = ''
        d = self.AccessionObjDict
        for k in d:

            txt += "{}:\n".format(k)

            for m in d[k]:

                txt += "\t{0: <10}: {1}\n".format(m,d[k][m])

        return txt

    def generate_AccessionObjDict(self):

        AccessionObjDict = dict()
        files = sorted(os.listdir(self.path))

        files = [f for f in files if f not in ['.DS_Store', '._.DS_Store']]

        for filename in files:

            acc, location_tuple, cir, identifier = FileName.parse(filename)

            if acc is not None:

                self.accessions.add(acc)

                if acc not in AccessionObjDict.keys():
                    AccessionObjDict[acc] = AccessionObj(acc)

                params = {'acc': acc, 'lcs': location_tuple, 'cir': cir, 'identifier': identifier}
                obj = xObjects.auto_create_from_file(self.path + filename, **params)

                AccessionObjDict[acc].add_xObject(obj)

            else: message('\tFile name invalid...')

        return AccessionObjDict
    
    def generate_GenomeObjDict(self):

        GenomeObjDict = dict()
        files = os.listdir(self.path_genomes)

        files = [f for f in files if f not in ['.DS_Store', '._.DS_Store']]

        for filename in files:

            key = filename.replace('.gb', '')

            GenomeObjDict[key] = os.path.join(self.path_genomes, filename)

        return GenomeObjDict

    def get_processing_list(self, reference_key, target_key, path_mode = "file_name_only"):

        if path_mode == "file_name_only": pre = ''
        else: pre = self.path

        todo_list = []
        for accession, AccObj in self.AccessionObjDict.items():

            ref_Obj = AccObj[reference_key]
            tar_Obj = AccObj[target_key]

            if ref_Obj.exists and not tar_Obj.exists:

                conv = '{}->{}'.format(reference_key, target_key)
                row = [accession, pre + ref_Obj.file_name, pre + FileName.convert(ref_Obj.file_name, conv)]

                todo_list.append(row)

        return todo_list

    def get_processing_list_of_acc_objects(self, target_key):

        todo_list = []
        for accession, AccObj in self.AccessionObjDict.items():

            if not AccObj[target_key].exists: todo_list.append(AccObj)

        return todo_list

    def get_path(self, accession, identifier, path_mode = "file_name_only"):

        if path_mode == "file_name_only":

            return self.AccessionObjDict[accession][identifier].file_name

        else:

            return self.path + self.AccessionObjDict[accession][identifier].file_name

    def get_xObject(self, accession, identifier):

        obj = self.AccessionObjDict[accession][identifier]
        return obj

    def debug_print(self):

        txt = ''

        for accession, AccObj in self.AccessionObjDict.items():

            txt += "{}:\n".format(accession)

            for identifier, Obj in AccObj.items():

                txt += "\t{0: <5}: {1}\n".format(identifier, Obj)
                txt += "\t{0: <5}: {1} {2}\n".format(' ', Obj.identifier, type(Obj))
                txt += "\t{0: <5}: Parent Accession - {1}\n".format(' ', Obj.parent.accession)

                txt += '\n'

        print(txt)

class DefFinderData:

    path_profiles    = "/home/baslan/.macsyfinder/data/defense-finder-models/profiles/"
    path_definitions = "/home/baslan/.macsyfinder/data/defense-finder-models/definitions/"
    path_pfamdb      = "/home/baslan/shared_folder/bioinformatics/PfamDB/PfamDB_35.0/Pfam-A.hmm" 
    path_database    = "/home/baslan/shared_folder/bioinformatics/PfamDB/DefFinderDB/"
    path_db_comps    = path_database + "comps/"
    path_db_rftrd    = path_db_comps + "refactored_hmms/"

    
    

    extra_pfams_list = ["PF00239.24", "PF07508.16", "PF13408.9"] 

    def __init__(self):

        pass

    def step1_update_profiles(self):

        cmd = "macsydata install -U --org mdmparis defense-finder-models"
        print(cmd)
        subprocess.run(cmd, shell=True)

    def step2_refactor_profiles(self):

        files = os.listdir(self.path_profiles)

        files = [file for file in files if file.endswith('.hmm')]

        for file in files:

            text_file = open(self.path_profiles + file, "r")
            profile = text_file.read()
            text_file.close()

            new_name = file.split('.')[0]

            

            profile = re.sub("^NAME(\s*)(.+)\n", lambda match: "NAME{}{}|{}\n".format(match.group(1), new_name, match.group(2)), profile, flags=re.MULTILINE)

            print(profile[0:100])
            print()

            write_file(self.path_db_rftrd + file, profile)

    def step3_create_hmm_database(self):

        

        for pfam in self.extra_pfams_list:

            self.aux_hmmfetch(self.path_pfamdb, pfam)

        

        if not os.path.isfile(self.path_db_comps + "all_minus.hmm"):

            cmd = f"cat {self.path_db_rftrd}*.hmm > {self.path_db_comps}all_minus.hmm"
            subprocess.run(cmd, shell=True)

        else:

            print('File all_minus.hmm already exists...')

        

        if not os.path.isfile(self.path_database + "all.hmm"):

            cmd = f"cat {self.path_db_comps}*.hmm > {self.path_database}all.hmm"
            subprocess.run(cmd, shell=True)

        else:

            print('File all.hmm already exists...')

        

        if not os.path.isfile(self.path_database + "all.hmm.h3i"):

            subprocess.run(f"hmmpress {self.path_database}all.hmm", shell=True)

    

    def aux_hmmfetch(self, hmm_file, key):

        out_file = self.path_db_comps + f"{key}.hmm"

        if os.path.isfile(out_file):

            print(f"File {key}.hmm already exported...")

        else:

            subprocess.run(f"hmmfetch {hmm_file} {key} > {out_file}", shell=True)

    def aux_get_profiles_report(self):

        def parse_profile(file):

            text_file = open(self.path_profiles + file, "r")
            profile = text_file.read()
            text_file.close()

            name, leng, ga_val = None, None, None

            nm = re.search(r'^NAME +(.+)\n', profile, re.MULTILINE)
            lg = re.search(r'^LENG +([0-9]+)\n', profile, re.MULTILINE)
            ga = re.search(r'^GA +(.+)\n', profile, re.MULTILINE)

            if nm is not None: name   = nm.group(1)
            if lg is not None: leng   = lg.group(1)
            if ga is not None: ga_val = ga.group(1)

            return name, leng, ga_val

        files = os.listdir(self.path_db_rftrd)

        files = [file for file in files if file.endswith('.hmm')]

        for file in files:

            name, leng, ga = parse_profile(file)

            print(50*'-')
            print(f"FILE: {file}")
            print(f">@Name: {name}, @Leng: {leng}, @GA: {ga}")
            print('\n\n')

    def aux_parsing_all_hmm(self):

        

        def parse_fields2(profile):

            name, leng, ga_val = None, None, None

            nm = re.search(r'^NAME +(.+)\n', profile, re.MULTILINE)
            lg = re.search(r'^LENG +([0-9]+)\n', profile, re.MULTILINE)
            ga = re.search(r'^GA +(.+)\n', profile, re.MULTILINE)

            if nm is not None: name   = nm.group(1)
            if lg is not None: leng   = lg.group(1)
            if ga is not None: ga_val = ga.group(1)

            return name, leng, ga_val

        text_file = open("/home/baslan/shared_folder/bioinformatics/PfamDB/DefFinderDB/all.hmm", "r")
        my_text = text_file.read()
        text_file.close()

        profiles = my_text.split('//')
        print(len(profiles))

        for profile in profiles:

            name, leng, ga = parse_fields2(profile)

            print(50*'-')
            print(f"@Name: {name}, @Leng: {leng}, @GA: {ga}")
            print('\n\n')

    @staticmethod
    def selft_test(pr):

        
        

        path_ngbrs = pr + "experiments_general/[NC_002976.3]_ngbrs_[2584175m2584655][n].fa"
        path_pfams = pr + "experiments_general/[NC_002976.3]_pfams_[2584175m2584655][n].txt"

        j = 3

        if j == 0:

            DFD = DefFinderData()
            
            DFD.step3_create_hmm_database()

        if j == 1:

            command = ["hmmscan",
                       "--cut_ga",
                       "--noali",
                       "--acc",
                       "--domtblout",
                       path_pfams,
                       "/home/baslan/shared_folder/bioinformatics/PfamDB/DefFinderDB/all.hmm",
                       path_ngbrs]

            proc = subprocess.Popen(command, stdout=subprocess.DEVNULL)  
            proc.wait()

        if j == 2:
            obj = xObjects.PfamsObj(path_pfams)
            dfo, dfn = obj.parse()

            print(dfo)
            print(dfn)

class DefPartition:

    def __init__(self, a=[], b=0, c=0): 

        self.def_systems = a
        self.n_def_systems = b
        self.n_def_genes = c

    def __getitem__(self, key):

        return getattr(self, key)

    def __str__(self):

        return "DP: {} {} {}".format(self.def_systems, 
                                     self.n_def_systems, 
                                     self.n_def_genes)

    def __add__(self, other):

        return DefPartition(self.def_systems + other.def_systems,
                            self.n_def_systems + other.n_def_systems,
                            self.n_def_genes + other.n_def_genes)

class DefenseFeature:

    def __init__(self, AccObj, path_genome, 
                 plot_params, pos_anchor = 200):

        self.accession = AccObj.accession
        self.path_genome = path_genome
        self.df = AccObj['dtfrm'].parse()

        self.all_genes, self.all_sys = AccObj['alsys'].parse()
        self.pos_anchor = pos_anchor

        normalize = lambda x: 100 * x / max_data

        cons_thresh = plot_params['th'] 
        max_data    = plot_params['n'] 
        gene_cons   = self.df['gene_conservation'].to_numpy()

        self.pos_acrb = self.get_acrb(normalize(gene_cons), 
                                      normalize(cons_thresh), 
                                      pos_anchor + 1)

    def get_nucl_dist(self, i_1, i_2):

        df = self.df
        df = df.reset_index().copy()

        lt_1 = LocationTuple.decode(df.iloc[i_1]['pos'])

        try:

            lt_2 = LocationTuple.decode(df.iloc[i_2]['pos'])

        except:

            print(f"problem -> lt_2: {self.accession} -> {df.iloc[i_2-1]['pos']}")
            lt_2 = LocationTuple.decode(df.iloc[i_2]['pos'], 
                                        ignore_fuzzy_loc=True)

        record = SeqIO.read(self.path_genome, "gb")
        dist = calculate_nucleotide_distance(lt_1, lt_2, record.seq)
        
        return dist, len(record.seq)

    def get_vector(self, mode, i):

        

        

        def map2region(i, xl, xu):

            if 0 <= i < xl:
                return 0
            elif xl <= i <= xu:
                return 1
            elif xu < i:
                return 2
            else:
                raise ValueError("Cannot determine region...")
            
        
                 
        df = self.df
        pa = self.pos_anchor
        pb = self.pos_acrb

        
        
        
        

        b = [[pa, pb],
             [pa - 199, pa + 199],
             [pa, pa + 300]]
        
        i_1 = b[i][0]
        i_2 = b[i][1]

        

        nuc_len, gnm_len = self.get_nucl_dist(i_1, i_2)

        

        df['Region'] = df['i'].apply(map2region, 
                                     xl=i_1, 
                                     xu=i_2)       
        rows = []
        for pos, gene in self.all_genes.items():

            try: r = df.loc[pos]['Region']
            except: r = 3

            rows.append({'pos': pos, 
                         'region': r,  
                         'type': gene['type'], 
                         'subtype': gene['subtype'], 
                         'gene_name': gene['gene_name']})

        df_dmap = pd.DataFrame(rows)

        

        
        
        
        
        
            
        DP = DefPartition   

        region_partitions = [DP(), DP(), DP(), DP(), DP()]

        if len(self.all_genes) != 0:

            P_all = DP() 
            for region in [0, 1, 2, 3]:

                op = df_dmap[df_dmap['region'] == region]
                def_systems = op[mode].unique().tolist()

                P = DefPartition(def_systems, len(def_systems), len(op))

                region_partitions[region] = P
                P_all += P

            region_partitions[4] = P_all
             
        return region_partitions, nuc_len, gnm_len

    

    def get_acrb(self, cons, threshold, min_accessory_length):

        for i, cons_i in enumerate(cons):

            cond1 = all(val >= threshold for val in cons[i:i + 3])
            cond2 = [val >= threshold * 0.8 for val in cons[i:i + 6]].count(True) > 5

            if  cond1 and cond2 and i > min_accessory_length:
                return i
            
        ii = len(cons) - 1
            
        print(f"\t\tNo accessory region found in {self.accession}.")
        print(f"\t\tSetting .acrb = {ii}")

        return ii

class SnailFeature:

    def __init__(self, AccObj: AccessionObj, offset, max_data, core_depth, data_length, conservation_threshold):

        self.AccObj = AccObj
        self.accession = AccObj.accession

        with open(AccObj['dtfrm'].path_full, 'rb') as f: df = pickle.load(f)

        self.df:pd.DataFrame = df.iloc[offset:offset + data_length].copy()
        self.defrelpfams = self.df['dpfam_set'].to_list()
        self.position_list = self.df.index.to_list()

        self.cons = self.df['gene_conservation'].to_numpy().astype(float)
        self.cons = 100 * self.cons / max_data
        self.dmap = self.df['defense_related'].to_numpy()
        self.dusA = self.df['is_dusA'].to_numpy()
        self.ccrX = self.df['contains_ccr_domains'].to_numpy()
        self.rlmh = self.df['is_rlmh'].to_numpy()
        self.len_ccr = 0 
        self.rlmh_loc_index = np.where(self.rlmh == 1)[0][0]

        self.accessory_region_boundary = self.get_accessory_region_boundary(100 * conservation_threshold / max_data, self.rlmh_loc_index + 1)
        self.alen = np.array([(np.nan if i!= self.accessory_region_boundary else 1) for i in range(data_length)])

        for j in range(data_length):

            if self.dmap[j] == 1 or self.rlmh[j] == 1:
                self.cons[j] = np.nan

            if j > self.accessory_region_boundary + core_depth:
                self.cons[j] = np.nan
                self.dmap[j] = np.nan
                self.dusA[j] = np.nan
                self.ccrX[j] = np.nan
                self.rlmh[j] = np.nan

            if self.rlmh_loc_index < j <= self.accessory_region_boundary:

                if self.ccrX[j] == 1: self.len_ccr +=1

    def get_accessory_region_boundary(self, threshold, min_accessory_length):

        

        for i, cons_i in enumerate(self.cons):

            cond1 = all(val >= threshold for val in self.cons[i:i + 3])
            cond2 = [val >= threshold * 0.8 for val in self.cons[i:i + 6]].count(True) > 5

            if  cond1 and cond2 and i > min_accessory_length:
                return i

        return len(self.cons)

    def get_defense_related_pfams_in_accessory_region(self):

        cum_dict = dict()
        cum_dict_counts = dict()

        a = self.rlmh_loc_index 
        b = self.accessory_region_boundary

        masked_defrelpfams = self.defrelpfams[a:b]
        masked_dmap = self.dmap[a:b]

        for i, txt in enumerate(masked_defrelpfams):

            
            d = ast.literal_eval(txt)
            valid = True if masked_dmap[i] == 1 else False
            

            if valid:

                for pfam in d:

                    def_type = d[pfam][2]
                    if def_type not in cum_dict.keys(): cum_dict[def_type] = {pfam}
                    else: cum_dict[def_type].add(pfam)

        for k in cum_dict.keys(): cum_dict_counts[k] = len(cum_dict[k])

        return cum_dict_counts

    def get_defense_signature(self, direction:int):

        def poik2(s: str):

            jobs = []
            for l in set(s):

                n = s.count(l)
                if n > 1:
                    jobs.append([l * n, l])  

            for job in jobs:
                s = s.replace(job[0], job[1])

            return s

        cdc = self.get_defense_related_pfams_in_accessory_region()

        defense_systems = []
        for def_type in cdc.keys():

            if def_type == 'RM':
                defense_systems.append('0')
            elif def_type == 'CRISPR':
                defense_systems.append('1')
            elif def_type == 'Abi':
                defense_systems.append('2')
            elif def_type == 'TA':
                defense_systems.append('3')
            elif def_type == 'BREX':
                defense_systems.append('4')
            elif "septu" in def_type.lower():
                defense_systems.append('5')
            elif def_type == 'Nhi':
                defense_systems.append('6')
            elif def_type == 'STK2':
                defense_systems.append('7')
            else:
                defense_systems.append('8')

        if direction == 0: signature = ''.join(sorted(defense_systems, reverse=False))

        else: signature = ''.join(sorted(defense_systems, reverse=True))

        signature = signature.replace('0', 'R')
        signature = signature.replace('1', 'C')
        signature = signature.replace('2', 'A')
        signature = signature.replace('3', 'T')
        signature = signature.replace('4', 'B')
        signature = signature.replace('5', 'U')
        signature = signature.replace('6', 'N')
        signature = signature.replace('7', 'S')
        signature = signature.replace('8', 'O')

        processed_signature = poik2(signature) 

        return processed_signature, cdc

    def rlmh_analysis_plot(self, path_figures, max_y_line = 89): 

        dus = None
        x = []
        y = []
        ccr_list = []

        path = "{}{}_individual_snail_plot.png".format(path_figures, self.accession)

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(x, y, 'g-')
        ax1.axhline(y=max_y_line, color='red', alpha=0.4, linewidth=2)

        dmap_array = np.array(self.dmap)
        ax2.scatter(x, np.cumsum(dmap_array), marker='|', s=20)
        

        if dus is not None:
            ax1.axvline(x=dus, color='black', alpha=1, linewidth=2)

        if len(ccr_list) > 0:
            for j in ccr_list: ax1.axvline(x=j, color='purple', alpha=1, linewidth=1)

        ax1.set_xlabel('Position')
        ax1.set_ylabel('Number of Strains Containing The Gene', color='g')
        ax2.set_ylabel('Is Defense Related', color='b')
        plt.savefig(path)
        

    def get_boundary_positions(self):

        rlmh_ind = self.rlmh_loc_index
        acrb_ind = self.accessory_region_boundary

        rlmh_pos = self.position_list[rlmh_ind]
        acrb_pos = self.position_list[acrb_ind]

        rlmh_loc = LocationTuple.decode(rlmh_pos)
        acrb_loc = LocationTuple.decode(acrb_pos)

        return rlmh_loc, acrb_loc

class SnailAnalyzerFull:

    def __init__(self, dvar):

        
        

        self.dvar = dvar
        self.plot_params  = dvar['plot_params']
        self.window_size  = dvar['win_size']
        self.numthreads   = dvar['numthreads']
        self.path_genomes = dvar['genomes']
        self.path_root    = dvar['root']
        self.path_mapped  = dvar['root'] + "mapped.csv"
        self.path_nbrhd   = dvar['root'] + "neighborhood_files_ws{}/".format(self.window_size)

        self.path_snaildt = dvar['root'] + "snail_data"
        self.path_dfntmp  = dvar['root'] + "defense_finder_results/"
        self.path_dfallr  = dvar['root'] + "defense_finder_results_allpr/"

        self.path_combined_fasta      = self.path_root + "combined.fasta"
        self.path_combined_all__fasta = self.path_root + "combined_all.fasta"

        self.path_groups  = dvar['root'] + "groups/"
        self.path_groups_all  = dvar['root'] + "groups_all/"

        ensure_directory(self.path_nbrhd)
        ensure_directory(self.path_dfntmp)
        ensure_directory(self.path_dfallr)

        ensure_directory(self.path_groups)
        ensure_directory(self.path_groups_all)

        self.displacement = 200
        self.file_manager = nbFileManager(self.path_nbrhd, self.path_genomes)
        self.all_existing_ngbrs_file_accessions = sorted(list(self.file_manager.accessions))

        self.suffix     = self.path_root.split('/')[-2].split('_')[-1]
        
        self.path_stats = self.path_root + f"defense_stats_i_{0}_{self.suffix}.txt"

    def step1_download_all_genbank(self):

        path = self.dvar['acc_list']
        path_accessions = path + "accession_list.txt"

        txt = read_file(path_accessions)
        accessions = [acc.strip() for acc in txt.splitlines() if acc != ""]

        print(len(accessions))

        download_genome(path, accessions, 'gbwithparts', 'text')

    def step2_generate_mapped_csv_from_collection(self):

        def process_accession_list(accessions_list, q):

            for accession in accessions_list:

                genome = gbkGenome(self.path_genomes, accession)
                rlmh_features, n = genome.find_protein_by_keyword(['methy', 'RlmH'], 'product')

                if len(rlmh_features) == 0:

                    prod = 'None'
                    loc = 'None'
                    Include, Note = False, 'No rlmh'

                elif rlmh_features[0].type != 'CDS' or 'translation' not in rlmh_features[0].qualifiers.keys():

                    prod = rlmh_features[0].qualifiers['product']
                    loc = str(rlmh_features[0].location)
                    Include, Note = False, 'pseudo rlmh'

                elif '>' in str(rlmh_features[0].location) or '<' in str(rlmh_features[0].location):

                    prod = rlmh_features[0].qualifiers['product']
                    loc = str(rlmh_features[0].location)
                    Include, Note = False, 'Funny Location'

                else:

                    prod = rlmh_features[0].qualifiers['product']
                    loc = str(rlmh_features[0].location)
                    Include, Note = True, 'Regular'

                row = [accession,
                       len(genome.record.features),
                       loc,
                       prod,
                       genome.record.description,
                       'chromosome',
                       Note,
                       Include]

                q.put(row)

        accessions = [file.split('.gb')[0] for file in os.listdir(self.path_genomes) if '.gb' in file]

        rows = multiprocess_with_pbar(process_accession_list,
                                         accessions,
                                         n_processes = 64,
                                         description = "Generating mapped CSV...")

        

        columns = ['Genome Accession',
                   '
                   'Feature Location',
                   'rlmh Product',
                   'Description',
                   'Context',
                   'Note',
                   'Include']

        df = pd.DataFrame(rows, columns=columns)
        df.to_csv(self.path_mapped, sep=',')

        print(df)

    def step3a_extract_rlmh_nbrhd_fasta_files(self):

        def extract_fasta_file_from_sub_mapped(df_sub_mapped, q):

            for index in df_sub_mapped.index:

                params = {'path_file': self.path_nbrhd,
                          'path_genomes': self.path_genomes,
                          'accession': df_sub_mapped.loc[index]['Genome Accession'],
                          'location': df_sub_mapped.loc[index]['Feature Location'],
                          'window_size': self.window_size,
                          'displacement': df_sub_mapped.loc[index]['Displacement']}

                obj = xObjects.NgbrsObj.create(**params)

                q.put(obj.path_full)

        df_mapped_csv, df_lookup_table, L = MappedCsv(self.path_mapped).parse()
        df_mapped_csv = df_mapped_csv[df_mapped_csv['Include'] == True]

        boolean_series = df_mapped_csv['Genome Accession'].isin(self.all_existing_ngbrs_file_accessions)
        df_mapped_csv = df_mapped_csv[~boolean_series]
        df_mapped_csv['Displacement'] = self.displacement
        processed = multiprocess_with_pbar(extract_fasta_file_from_sub_mapped,
                                         df_mapped_csv,
                                         n_processes = 64,
                                         description = "Extracting Fasta Files...")

        print('\n----------------------------')
        for p in processed: print(p)

    def step3b_extract_rlmh_all_proteins_fasta_files(self):

        processing_list = self.file_manager.get_processing_list_of_acc_objects('allpr')
        print("Number of Files to be Done: {}".format(len(processing_list)))

        for AccObj in processing_list:

            params = {'path_file': self.path_nbrhd,
                      'path_genomes': self.path_genomes,
                      'accession': AccObj.accession,
                      'location_tuple': AccObj['ngbrs'].location_tuple,
                      'displacement': self.displacement}

            AccObj['allpr'].create(**params)

    def step4a_run_defense_finder(self):

        def run_defense_finder_in_serial(processing_list):

            for accession, path_ngbrs, path_defin in processing_list:
                print(accession)
                print(path_ngbrs)
                print(path_defin)

                st = time.perf_counter()
                _ = xObjects.DeFinObj.create(path_ngbrs, path_defin, self.path_dfntmp)
                ed = time.perf_counter()
                print('Finished in:', ed - st)
                print()

        def run_defense_finder_in_batch(sub_list, q):

            for accession, path_ngbrs, path_defin in sub_list:

                _ = xObjects.DeFinObj.create(path_ngbrs, path_defin, self.path_dfntmp)
                q.put(accession)

        processing_list = self.file_manager.get_processing_list('ngbrs', 'defin', path_mode="full_path")
        print(len(processing_list))

        st = time.perf_counter()
        multiprocess_with_pbar(run_defense_finder_in_batch, processing_list, n_processes=4, description="Running Defense Finder...")
        ed = time.perf_counter()
        print('Finished in:', ed - st)

    def step4b_run_defense_finder_on_all_proteins(self):

        def run_serial(processing_list):

            for accession, path_allpr, path_alsys in processing_list:

                print(accession)
                print(path_allpr)
                print(path_alsys)
                st = time.perf_counter()
                _ = xObjects.AlSysObj.create(path_allpr, path_alsys, self.path_dfallr)
                ed = time.perf_counter()
                print('Finished in:', ed - st)
                print()

        def run_defense_finder_in_batch(sub_list, q):

            for accession, path_allpr, path_alsys in sub_list:

                _ = xObjects.AlSysObj.create(path_allpr, path_alsys, self.path_dfallr)
                q.put(accession)

        processing_list = self.file_manager.get_processing_list('allpr', 'alsys', path_mode="full_path")
        print(len(processing_list))

        st = time.perf_counter()
        multiprocess_with_pbar(run_defense_finder_in_batch, processing_list, n_processes=4, description="Running AllSystem Defense Finder...")
        ed = time.perf_counter()
        print('Finished in:', ed - st)

    def step6a_generate_levenshtein_groups_dict(self):

        

        if not os.path.isfile(self.path_combined_fasta):

            print('Creating combined fasta...')

            txt = ''
            for accession in self.all_existing_ngbrs_file_accessions:

                path = self.file_manager.get_path(accession, 'ngbrs', path_mode="full_path")
                feature_list = list(SeqIO.parse(path, "fasta")) 

                for i, seq_record in enumerate(feature_list):

                    seq_record.id = "{}|{}|{}".format(accession, seq_record.id, i)
                    txt += '>{}\n{}\n'.format(seq_record.id, seq_record.seq)

            write_file(self.path_combined_fasta, txt)

        else:

            print('Combined fasta already exists...')
            pass

        

        st = time.perf_counter()
        print('Started @ {} '.format(st))
        print()

        dynamic_acc_pos_list = list(SeqIO.parse(self.path_combined_fasta, "fasta"))

        gnum = 0

        groups_dict = dict()

        while len(dynamic_acc_pos_list) > 0:

            n_before = len(dynamic_acc_pos_list)

            SeqRec = dynamic_acc_pos_list[0]
            group_id_list = levenshtein_searcher(SeqRec.seq, dynamic_acc_pos_list)

            n = len(group_id_list)

            
            if n > 982:

                print('\t>982: ', gnum, SeqRec.id, SeqRec.seq)
                with open(self.path_groups + "lvs_group_{}.txt".format(gnum), 'wb') as f: pickle.dump(group_id_list, f)

            elif n == 0:

                print('\t=0: ', gnum, SeqRec.id, SeqRec.seq)
                group_id_list.append(SeqRec.id)

            groups_dict[gnum] = group_id_list

            dynamic_acc_pos_list = [SR for SR in dynamic_acc_pos_list if SR.id not in group_id_list] 

            n_after = len(dynamic_acc_pos_list)

            print(f"Searching {SeqRec.id}: ", gnum, n, n_before, n_after)

            gnum += 1

        
        
        with open(self.path_root + "groups_dict.dic", 'wb') as f: pickle.dump(groups_dict, f)

        ed = time.perf_counter()
        print('finished in: {} seconds'.format(ed - st))

    def step6b_generate_levenshtein_all_groups_dict(self):

        s = 1

        if s == 0:

            if not os.path.isfile(self.path_combined_all__fasta):

                print('Creating -ALL- combined fasta...')

                txt = ''
                for accession in self.all_existing_ngbrs_file_accessions:

                    path = self.file_manager.get_path(accession, 'allpr', path_mode="full_path")
                    feature_list = list(SeqIO.parse(path, "fasta")) 

                    for i, seq_record in enumerate(feature_list):

                        seq_record.id = "{}|{}|{}".format(accession, seq_record.id, i)
                        txt += '>{}\n{}\n'.format(seq_record.id, seq_record.seq)

                write_file(self.path_combined_all__fasta, txt)

            else:

                print('Combined fasta already exists...')
                pass

        if s == 1:

            st = time.perf_counter()
            print('Started @ {} '.format(st), '\n')

            dynamic_acc_pos_list = list(SeqIO.parse(self.path_combined_all__fasta, "fasta"))

            gnum = 0

            groups_dict = dict()

            while len(dynamic_acc_pos_list) > 0:

                n_before = len(dynamic_acc_pos_list)

                SeqRec = dynamic_acc_pos_list[0]
                group_id_list = levenshtein_searcher(SeqRec.seq, dynamic_acc_pos_list)

                n = len(group_id_list)

                
                if n > 89:

                    print('\t>89: ', gnum, SeqRec.id, SeqRec.seq)
                    with open(self.path_groups_all + "lvs_group_{}.txt".format(gnum), 'wb') as f: pickle.dump(group_id_list, f)

                elif n == 0:

                    print('\t=0: ', gnum, SeqRec.id, SeqRec.seq)
                    group_id_list.append(SeqRec.id)

                groups_dict[gnum] = group_id_list

                dynamic_acc_pos_list = [SR for SR in dynamic_acc_pos_list if SR.id not in group_id_list] 

                n_after = len(dynamic_acc_pos_list)

                print(f"Searching {SeqRec.id}: ", gnum, n, n_before, n_after)

                gnum += 1

            
            
            with open(self.path_root + "groups_dict_all.dic", 'wb') as f: pickle.dump(groups_dict, f)

            ed = time.perf_counter()
            print('finished in: {} seconds'.format(ed - st))

    def step7a_generate_core_dictionaries_from_lvs_groups_dict(self):

        

        genome_dict = {acc: {} for acc in self.all_existing_ngbrs_file_accessions}

        with open(self.path_root + "groups_dict.dic", 'rb') as f:
            groups_dict = pickle.load(f)

        

        for gnum, group in groups_dict.items():

            acc_set = set()
            pair_list = []

            for hit in group:
                acc, pos, _ = hit.split('|')
                pair_list.append((acc, pos))
                acc_set.add(acc)

            for acc, pos in pair_list:
                genome_dict[acc][pos] = [acc_set, len(acc_set)]

        

        for acc, AccObj in self.file_manager.AccessionObjDict.items():

            path_ngbrs = AccObj['ngbrs'].path_full
            path_cored = FileName.convert(path_ngbrs, conv='ngbrs->cored')
            print(acc)
            print('\t', path_ngbrs)
            print('\t', path_cored)
            print()

            with open(path_cored, 'wb') as f: pickle.dump(genome_dict[acc], f)

        return genome_dict

    def step7b_generate_core_dictionaries_all(self):

        

        genome_dict = {acc: {} for acc in self.all_existing_ngbrs_file_accessions}

        with open(self.path_root + "groups_dict_all.dic", 'rb') as f:
            groups_dict = pickle.load(f)

        

        for gnum, group in groups_dict.items():

            acc_set = set()
            pair_list = []

            for hit in group:
                acc, pos, _ = hit.split('|')
                pair_list.append((acc, pos))
                acc_set.add(acc)

            for acc, pos in pair_list:
                genome_dict[acc][pos] = [acc_set, len(acc_set)]

        

        for acc, AccObj in self.file_manager.AccessionObjDict.items():

            path_ngbrs = AccObj['ngbrs'].path_full
            path_coral = FileName.convert(path_ngbrs, conv='ngbrs->coral')
            print(acc)
            print('\t', path_ngbrs)
            print('\t', path_coral)
            print()

            with open(path_coral, 'wb') as f: pickle.dump(genome_dict[acc], f)

        return genome_dict

    def step8_generate_dataframes(self):

        def generate_dataframe_files(sub_plist,q):

            for AccObj in sub_plist:

                AccObj['dtfrm'].create(save=True)
                q.put(AccObj.accession)

        processing_list = self.file_manager.get_processing_list_of_acc_objects('dtfrm')
        print("Number of Files to be Done: {}".format(len(processing_list)))

        st = time.perf_counter()
        print('start time: ', st, ' | ', time.localtime())
        multiprocess_with_pbar(generate_dataframe_files, processing_list, n_processes=64, description="Generating DataFrames-parallel...")
        ed = time.perf_counter()
        print('finish time: ', ed, ' | ', time.localtime())
        print('Finished in:', ed - st, 'seconds')

    def step9_plot_snail_diagram(self):

        
        

        params = self.plot_params

        annotate_numbers = params['annotate_nums']
        annotate_defenses = params['annotate_defs']
        outer_circle = params['circle']

        n = params['n']
        m = params['m']

        conservation_threshold = params['th']
        core_depth = params['c_depth'] 
        font_size_1 = params['font_size_1']
        font_size_2 = params['font_size_2']
        res_dpi = params['res_dpi']
        fw = params['fw']
        fh = params['fh']

        def_sys_position = m-3 
        gen_num_position = m+2 

        offset = 200 - core_depth  
        max_data = n
        data = []

        print_cdc = False

        

        if not os.path.isfile(self.path_snaildt):

            print('Generating Snail Data...')

            for accession, AccObj in self.file_manager.AccessionObjDict.items():

                data.append(SnailFeature(AccObj, offset, max_data, core_depth, m, conservation_threshold))

            data.sort(key=operator.attrgetter('accessory_region_boundary', 'accession'), reverse=False)

            with open(self.path_snaildt, 'wb') as f: pickle.dump(data, f)

        else:

            print('Reading Snail Data...')

            with open(self.path_snaildt, 'rb') as f: data = pickle.load(f)

        

        full_angle = 2 * np.pi * (n - 1) / n
        r, th = np.meshgrid(np.linspace(0, m - 1, m), np.linspace(0, full_angle, n))

        tick_labels = [sf.accession for sf in data] 
        ccrs = [sf.len_ccr for sf in data]
        alns = [(sf.accessory_region_boundary - sf.rlmh_loc_index) for sf in data]

        print('Accessory Region (genes) -> min: %s, median: %s, max: %s' % (min(alns), np.median(alns), max(alns)))
        print('Accessory Region (genes) -> 5-pctl: %s, 95-pctl: %s' % (np.percentile(alns, 5), np.percentile(alns, 95)   ))

        cons = np.asarray([sf.cons for sf in data])
        dmap = np.asarray([sf.dmap for sf in data])
        ccrX = np.asarray([sf.ccrX for sf in data])
        rlmh = np.asarray([sf.rlmh for sf in data])
        alen = np.asarray([sf.alen for sf in data])

        
        
        
        

        

        dtheta = full_angle / len(th) / 2
        shading_mode = 'auto'  

        fig = plt.figure(figsize=(fw, fh), dpi=res_dpi)
        ax = fig.add_subplot(111, polar=True)

        pc = ax.pcolormesh(th + dtheta, r, cons, cmap='Blues', shading=shading_mode)  
        ax.pcolor  (th + dtheta, r, dmap, edgecolor='none', cmap=matplotlib.colors.ListedColormap(['red'])   , shading=shading_mode)
        ax.scatter (th + dtheta, r, s = ccrX    , marker = 'o', color='black' , facecolor='black', alpha=1)

        
        ax.pcolor  (th + dtheta, r, rlmh, edgecolor='none', cmap=matplotlib.colors.ListedColormap(['yellow']), shading=shading_mode)

        cbar = fig.colorbar(pc, pad=0.2)
        cbar.ax.set_ylabel('Percent Gene Conservation', rotation=90, fontsize=14)

        a = [angles[0] for angles in th]
        ax.set_xticks(a)
        ax.set_xticklabels(tick_labels, fontsize=font_size_2)

        il = 0
        label_records = []

        
        if annotate_numbers:

            for label, angle, ccr_len in zip(ax.get_xticklabels(), a, ccrs):

                label_text = label.get_text()
                il += 1
                label_records.append('({}) {}'.format(il, label_text))

                fw = 'normal'
                c = 'black'

                if 'NC_002976' in label_text: c = 'red'

                rotation = np.rad2deg(angle)

                if np.pi / 2 <= angle < 3 * np.pi / 2:

                    alignment = "right"
                    rotation = rotation + 180
                    ss = '{} ({})'.format(il, ccr_len) 

                else:

                    alignment = "left"
                    ss = '({}) {}'.format(ccr_len, il) 

                
                lab = ax.text(x=angle + dtheta,
                              y=gen_num_position,  
                              s=ss,
                              ha=alignment,  
                              va="center",  
                              color=c,
                              fontweight=fw,
                              fontsize=font_size_1,
                              rotation=rotation,
                              rotation_mode="anchor")

        
        if annotate_defenses:

            for angle, ds in zip(a, data): 

                rotation = np.rad2deg(angle)

                if np.pi / 2 <= angle < 3 * np.pi / 2:

                    alignment = "left"  
                    rotation = rotation + 180
                    
                    ss, cdc = '', ''

                else:

                    alignment = "right"  
                    
                    ss, cdc = '', ''

                lab = ax.text(x=angle + dtheta,
                              y=def_sys_position,  
                              s=ss,
                              ha=alignment,  
                              va="center",  
                              color='red',
                              fontweight='bold',
                              fontsize=font_size_1,
                              fontfamily='DejaVu Sans',
                              rotation=rotation,
                              rotation_mode="anchor")

                if print_cdc: print(cdc)

        

        write_file(self.path_root + "snail_labels.txt", ', '.join(label_records))
        ax.set_xticklabels([])
        

        if not outer_circle:
            ax.spines['polar'].set_visible(False) 

        ax.axes.get_yaxis().set_visible(False)

        

        ax.xaxis.grid(False)
        ax.yaxis.grid(False)

        
        
        fig.savefig(f"snail_diagram_{self.suffix}.png")

        

    def step10_plot_rectangular_coral(self):

        AccObj = self.file_manager.AccessionObjDict['NZ_AP017922.1'] 

        az = AccObj['coral'].path_full
        print(az)

    
        
    def parse_snail_labels(self):

        path = self.path_root + "snail_labels.txt"
        txt = read_file(path)
        labels = txt.split(', ')

        acc_order_dict = dict()
        for label in labels:

            m = re.search(r"^\(([0-9]+)\) *(.+)$", label)
            num = int(m.group(1).strip())
            acc = m.group(2).strip()
            acc_order_dict[acc] = num

        return acc_order_dict
        
    def get_diversity_dicts(self, pos_anchor, mode, i = 0):

        su = self.suffix
        pr = self.path_root
        p_csv = pr + f"defense_systems_i_{i}_{su}.csv"
        p_dat = pr + f"diversity_dicts_i_{i}_{su}.data"

        if os.path.isfile(p_dat):

            print("\tReading Diversity Dicts...")

            with open(p_dat, 'rb') as f: 

                data = pickle.load(f)
                df, defense_freq_in, defense_freq_ou = data

        else:

            print("\tGenerating Diversity Dicts...")

            FM = self.file_manager
            GD = FM.GenomeObjDict

            rows = []
            defense_freq_in = defaultdict(int)
            defense_freq_ou = defaultdict(int)

            for accession, AccObj in FM.AccessionObjDict.items():

                DFR = DefenseFeature(AccObj, 
                                    GD[accession],
                                    self.plot_params, 
                                    pos_anchor = pos_anchor)
                
                

                rp, nuc_len, gnm_len = DFR.get_vector(mode = mode, i= i)

                
                
                defs_in = sorted((rp[1])['def_systems'])
                defs_ou = sorted((rp[0] + rp[2] + rp[3])['def_systems'])

                for defense in defs_in: 
                    defense_freq_in[defense] += 1

                for defense in set(defs_ou): 
                    defense_freq_ou[defense] += 1

                n_rg1 = rp[1]['n_def_systems']
                n_all = rp[4]['n_def_systems'] + 1e-5

                dd = {'accession'       : accession,
                      'var_region_L'    : DFR.pos_acrb - DFR.pos_anchor,
                      'region defenses' : ', '.join(defs_in),
                      'region N'        : len(defs_in),
                      'other defenses'  : ', '.join(defs_ou),
                      'other N'         : len(defs_ou),
                      'ratio'           : round(n_rg1 / n_all, 2),
                      'nuc_len'         : nuc_len,
                      'gnm_len'         : gnm_len}

                rows.append(dd)

            
                
            acc_order_dict = self.parse_snail_labels()
                
            df = pd.DataFrame(rows)

            df['i'] = df['accession'].map(acc_order_dict)

            df.set_index('i', inplace=True)

            df.sort_values(['i'], ascending=True, inplace=True)

            df.to_csv(p_csv, sep=',', index=True)

            print(df)

            

            data = (df, defense_freq_in, defense_freq_ou)
            with open(p_dat, 'wb') as f: pickle.dump(data, f)
        
        return df, defense_freq_in, defense_freq_ou
        
    def statistical_analysis_per_genome(self, i, mode, thres, alpha = 0.05, pos_anchor = 200):

        def calculate_p_value_v0(row):

            
            
            

            p = row['nuc_len']/row['gnm_len']

            n_in = row['region N']
            N = n_in + row['other N']

            return 1 - binom.cdf(n_in-1, N, p) if N >= thres else np.nan

        def calculate_p_value(row):

            
            
            

            p = row['nuc_len']/row['gnm_len']

            n_in = row['region N']
            N = n_in + row['other N']

            p_val = 0
            for n in range(n_in, N+1):
                p_val += binom.pmf(n, N, p)

            return p_val if N >= thres else np.nan           

        

        df, _, _ = self.get_diversity_dicts(pos_anchor=pos_anchor, mode=mode, i=i)
        df.set_index('accession', inplace=True)

        df['p-value'] = df.apply(calculate_p_value, axis=1)
    
        df = df[df['p-value'].notna()]
        df = df.sort_values(by=['p-value'], ascending=True).copy()

        df['rank'] = df['p-value'].rank(method='first', ascending=True)
        df['q-value'] = df['rank'] * alpha / len(df)
        df['rejected'] = df['p-value'] <= df['q-value']

        M = len(df)
        n = df.loc[df['rejected']].shape[0]

        percentage = n / M

        return df, percentage, n, M
    
    def statistical_analysis_collective(self, i, mode, alpha = 0.05, pos_anchor = 200, ratio_filter = 1):

        df, _, _ = self.get_diversity_dicts(pos_anchor=pos_anchor, mode=mode, i=i)

        df['acrb_pct'] = df['nuc_len'] / df['gnm_len']
        mean_acrb_pct = df['acrb_pct'].mean()
        std_acrb_pct  = df['acrb_pct'].std()

        df['total_defenses'] = df['region N'] + df['other N']
        
        df = df[df['ratio'] <= ratio_filter]

        defenses_inside  = df['region N'].sum()
        defenses_outside = df['other N'].sum()
        defenses_total   = df['total_defenses'].sum()

        pct = mean_acrb_pct + 3 * std_acrb_pct

        expected_inside  = round(defenses_total * pct)
        expected_outside = defenses_total - expected_inside

        observed = np.array([defenses_inside, defenses_outside])
        expected = np.array([expected_inside, expected_outside])

        

        deg_freedom = len(observed) - 1

        chi = sum((observed - expected)**2 / expected)
        chi_critical = chi2.ppf(1-alpha, deg_freedom)

        p_value = chi2.sf(chi, deg_freedom) * 2

        
        
        
        
        
        
        
        
        
        

        return df, mean_acrb_pct, std_acrb_pct, expected, observed, chi, chi_critical, p_value
        
    def stats_csv(self, i, df, dfi, dfo):

        acrb_list = df['VarRegionL'].to_numpy()

        pctl_5  = np.percentile(acrb_list, 5)
        pctl_95 = np.percentile(acrb_list, 95)
        mmin = min(acrb_list)
        mmax = max(acrb_list)
        mmed = np.median(acrb_list)

        lines = ["",
                 f"{self.suffix} Data & Stats\n\n",
                 f"Accessory Region  5-percentile: {pctl_5}\n",
                 f"Accessory Region 95-percentile: {pctl_95}\n",
                 "Accessory Region:\n",
                 f" min: {mmin}\n", 
                 f" max: {mmax}\n",
                 f" med: {mmed}\n",
                 f"Defenses Inside Region 1 (defined by i={i})\n",
                 "[<i> defines what region 1 is:\n",
                 "   i = 0 -> Rlmh     <->  ACBR (variable region)\n", 
                 "   i = 1 -> Rlmh-199 <->  Rlmh+199\n", 
                 "   i = 2 -> Rlmh     <->  Rlmh+300\n",
                 repr(dict(dfi)),
                 "\n\n",
                 "Defenses Everywhere Else)\n",
                 repr(dict(dfo))]

        txt = ''.join(lines)

        write_file(self.path_stats, txt)

def PLOT_DEFENSE_DIVERSITY(dvars, mode, pos_anchor = 200, i = 0):

    global pp

    path = pp + "Defense_Diversity/"
    ensure_directory(path)

    path_csv1 = path + f"defense_diversity_all_{i}.csv"
    path_fig1 = path + f"defense_diversity_all_{i}.png"

    path_csv2 = path + f"defense_diversity_inside_{i}.csv"
    path_fig2 = path + f"defense_diversity_inside_{i}.png"

    path_stats = path + f"diversity_stats_{i}.txt"

    norm_epi = 89
    norm_aur = 982

    SA = SnailAnalyzerFull(dvars['epidermidis'])
    _, dfi, dfo = SA.get_diversity_dicts(pos_anchor, mode = mode, i=i)
    d_epi_in, d_epi_ou = dict(dfi), dict(dfo)

    SA = SnailAnalyzerFull(dvars['aureus'])
    _, dfi, dfo = SA.get_diversity_dicts(pos_anchor, mode = mode, i=i)
    d_aur_in, d_aur_ou = dict(dfi), dict(dfo)

    

    pdf_epi_in = calculate_probability_dist(d_epi_in)
    pdf_aur_in = calculate_probability_dist(d_aur_in)
    pdf_epi_ou = calculate_probability_dist(d_epi_ou)
    pdf_aur_ou = calculate_probability_dist(d_aur_ou)

    h_epi_in = calculate_shannon_entropy(pdf_epi_in)
    h_aur_in = calculate_shannon_entropy(pdf_aur_in)
    h_epi_ou = calculate_shannon_entropy(pdf_epi_ou)
    h_aur_ou = calculate_shannon_entropy(pdf_aur_ou)

    lines = [f"h_epi_in: {h_epi_in}, h_epi_ou: {h_epi_ou}",
             f"h_aur_in: {h_aur_in}, h_aur_ou: {h_aur_ou}",
             f"Diversity Ratio (Epi): {(2**h_epi_in)/(2**h_epi_ou)}",
             f"Diversity Ratio (Aur): {(2**h_aur_in)/(2**h_aur_ou)}"]
    
    txt = '\n'.join(lines)
    print(lines)
    write_file(path_stats, txt)

    
    
    

    z0 = list(d_epi_in.keys())
    z1 = list(d_aur_in.keys())
    z2 = list(d_epi_ou.keys())
    z3 = list(d_aur_ou.keys())

    x_labels = sorted(set(z0 + z1 + z2 + z3))

    df = pd.DataFrame({'x_labels': x_labels})

    df['epi_in'] = [d_epi_in.get(key, 0) for key in x_labels]
    df['aur_in'] = [d_aur_in.get(key, 0) for key in x_labels]
    df['epi_ou'] = [d_epi_ou.get(key, 0) for key in x_labels]
    df['aur_ou'] = [d_aur_ou.get(key, 0) for key in x_labels]

    df.sort_values(['epi_in', 'aur_in'], ascending=[True, False], inplace=True)

    df.to_csv(path_csv1)

    df['epi_in'] = df['epi_in'] / norm_epi
    df['aur_in'] = df['aur_in'] / norm_aur
    df['epi_ou'] = df['epi_ou'] / norm_epi
    df['aur_ou'] = df['aur_ou'] / norm_aur

    

    _, ax = plt.subplots(1, 1, figsize=(7, 7))

    bar_positions = np.arange(len(x_labels))
    bar_width = 0.40

    ax.barh(bar_positions - bar_width / 2, -df['epi_in'], height=bar_width, label='epi_in', color = 'blue')
    ax.barh(bar_positions + bar_width / 2, -df['aur_in'], height=bar_width, label='aur_in', color = 'orange')
    ax.barh(bar_positions - bar_width / 2,  df['epi_ou'], height=bar_width, label='epi_in', color = 'blue')
    ax.barh(bar_positions + bar_width / 2,  df['aur_ou'], height=bar_width, label='aur_in', color = 'orange')
    ax.axvline(x=0, color='black', label='*')
    ax.set_yticks(bar_positions)
    ax.set_yticklabels(df['x_labels'], rotation=0, fontsize=8)
    ax.set_axisbelow(True)
    ax.grid(which='major', alpha=0.7)

    plt.tight_layout()
    plt.savefig(path_fig1)

    

    
    
    

    z0 = list(d_epi_in.keys())
    z1 = list(d_aur_in.keys())

    x_labels = sorted(set(z0 + z1))

    df = pd.DataFrame({'x_labels': x_labels})

    df['epi_in'] = [d_epi_in.get(key, 0) for key in x_labels]
    df['aur_in'] = [d_aur_in.get(key, 0) for key in x_labels]

    df.sort_values(['epi_in', 'aur_in'], ascending=[False, False], inplace=True)

    df.to_csv(path_csv2)

    df['epi_in'] = df['epi_in'] / norm_epi
    df['aur_in'] = df['aur_in'] / norm_aur

    _, ax = plt.subplots(1, 1, figsize=(8, 4), dpi=150)

    bar_positions = np.arange(len(x_labels))
    bar_width = 0.40

    ax.bar(bar_positions - bar_width / 2, df['epi_in'], width=0.3, align='center', label='epi_in', color = 'blue')
    ax.bar(bar_positions + bar_width / 2, df['aur_in'], width=0.3, align='center', label='aur_in', color = 'orange')
    
    ax.set_xticks(bar_positions)
    ax.set_xticklabels(df['x_labels'], rotation=90, fontsize=8)
    ax.set_axisbelow(True)
    ax.grid(which='major', alpha=0.7)

    plt.tight_layout()
    plt.savefig(path_fig2)

    

def PLOT_DEFENSE_CONCENTRATION(dvars, mode, i = 2, plot_mode = 1):

    global pp

    path = pp + "Defense_Concentration/"
    ensure_directory(path)

    path_plot_0 = path + f"defense_concentration_i_{i}_PDF.png"
    path_plot_1 = path + f"defense_concentration_i_{i}_CDF.png"
    path_plot_2 = path + f"defense_concentration_i_{i}_PDF_seperated.png"

    

    def cumsum_exclude(a):

        reversed_data = np.flip(a)
        cumulative_sum = np.cumsum(reversed_data)
        return np.flip(cumulative_sum)
    
    def get_median_min_max(df):

        df['region_ratio'] = df['nuc_len'] / df['gnm_len']
        mmed = np.median(df['region_ratio'])
        mmin = min(df['region_ratio'])
        mmax = max(df['region_ratio'])

        return mmed, mmin, mmax
    
    bins = np.arange(0, 1.1, 0.1)

    SA = SnailAnalyzerFull(dvars['epidermidis'])
    df_epi, _, _ = SA.get_diversity_dicts(200, mode = mode, i=i)

    SA = SnailAnalyzerFull(dvars['aureus'])
    df_aur, _, _ = SA.get_diversity_dicts(200, mode = mode, i=i)

    ratios_epi = df_epi['ratio'].to_numpy()
    ratios_aur = df_aur['ratio'].to_numpy()

    
    if   plot_mode == 0:

        hist_epi, bins = np.histogram(ratios_epi, bins=bins, density=False)
        hist_aur, bins = np.histogram(ratios_aur, bins=bins, density=False)

        CDF_epi = cumsum_exclude(hist_epi) * 0.1
        CDF_aur = cumsum_exclude(hist_aur) * 0.1

        median_epi, mmin_epi, mmax_epi = get_median_min_max(df_epi)
        median_aur, mmin_aur, mmax_aur = get_median_min_max(df_aur)

        
        fg_epi = round(100*sum(hist_epi[5:]) / sum(hist_epi), 0)
        fg_aur = round(100*sum(hist_aur[5:]) / sum(hist_aur), 0)

        txt =  "Fraction of genomes with 50% or more defenses in the region\n"
        txt += f"Epi: {fg_epi}% Aur: {fg_aur}%"

        

        ep = 0 
        wd = 0.05
        N_epi = 89
        N_aur = 982
        ylim_epi =  0.5
        ylim_aur = -0.5

        _, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=250)

        ax.bar(bins[:-1] - ep, hist_epi/N_epi,  width=wd, align='center', color='blue', label='epidermidis')
        ax.bar(bins[:-1] + ep, -hist_aur/N_aur, width=wd, align='center', color='orange', label='aureus')
        ax.set_xlabel("Fraction of defenses within 300 genes downstream of rlmh")
        ax.set_ylabel("Number of Genomes")
        ax.fill_betweenx([0, ylim_epi], mmin_epi, mmax_epi, color='blue', alpha=0.2)
        line = Line2D([median_epi, median_epi], [0, ylim_epi], color='blue', linestyle='dotted', linewidth=2)
        ax.add_line(line)

        ax.fill_betweenx([0, ylim_aur], mmin_aur, mmax_aur, color='orange', alpha=0.2)
        line = Line2D([median_aur, median_aur], [0, ylim_aur], color='orange', linestyle='dotted', linewidth=2)
        ax.add_line(line)

        ax.set_ylim(ylim_aur, ylim_epi)
        ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
        ax.grid(which='major', linestyle='solid', color='grey', alpha=0.4)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(round(100*x))}%'))
        ax.set_title(txt, fontsize=9)

        ax.set_axisbelow(True)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.tight_layout()
        plt.savefig(path_plot_0)

    elif plot_mode == 1:

        ep = 0 

        hist_epi, bins = np.histogram(ratios_epi, bins=bins, density=True)
        hist_aur, bins = np.histogram(ratios_aur, bins=bins, density=True)

        CDF_epi = cumsum_exclude(hist_epi) * 0.1
        CDF_aur = cumsum_exclude(hist_aur) * 0.1

        _, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=250)

        ax.bar(bins[:-1] - ep, CDF_epi, width=0.025, align='center', color='blue', label='epidermidis')
        ax.bar(bins[:-1] + ep, -CDF_aur, width=0.025, align='center', color='orange', label='aureus')
        ax.set_xlabel("Fraction of defenses within 300 genes downstream of rlmh")
        ax.set_ylabel("Number of Genomes")
        
        ax.grid(which='major', linestyle='--', color='grey', alpha=0.5)
        ax.grid(which='major', linestyle='solid', color='grey', alpha=0.4)

        ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
        

        

        
        ax.set_axisbelow(True)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.tight_layout()
        plt.savefig(path_plot_1)

    else:

        ylim_epi = 0.5
        ylim_aur = 0.5

        _, ax = plt.subplots(2, 1, figsize=(5, 3), dpi=250)

        def asdfg(ax, hist, N, mmin, mmax, median, ylim, txt):

            ax.bar(bins[:-1], hist/N, width=wd, align='center', color='blue', label='epidermidis')

            ax.set_xlabel("Fraction of defenses within 300 genes downstream of rlmh")
            ax.set_ylabel("Number of Genomes")

            ax.fill_betweenx([0, ylim], mmin, mmax, color='blue', alpha=0.2)
            line = Line2D([median, median], [0, ylim], color='blue', linestyle='dotted', linewidth=2)
            
            ax.add_line(line)
            
            ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
            ax.grid(which='major', linestyle='solid', color='grey', alpha=0.4)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(round(100*x))}%'))

            ax.set_title(txt, fontsize=9)
            ax.set_axisbelow(True)

        asdfg(ax[0], hist_epi, N_epi, mmin_epi, mmax_epi, median_epi, ylim_epi, txt)
        asdfg(ax[1], hist_aur, N_aur, mmin_aur, mmax_aur, median_aur, ylim_aur, txt)

        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.tight_layout()
        plt.savefig(path_plot_2)

def DO_ENRICHMENT_STATS_collective(dvars, mode, i = 0):

    def main(path, plot_name, ratio_filter = 0.5, save_csv = False):

        prefix = "enrichment_testing_collective"

        _, ax = plt.subplots(4, 1, figsize=(2, 6), dpi = 250)

        cols = ['Observed', 'Expected']
        rows = ['Inside', 'Outside']

        plot_i = 0
        rows_main = []
        for s in ['epidermidis', 'aureus']:

            for i in [0, 2]:

                SA = SnailAnalyzerFull(dvars[s])
                out = SA.statistical_analysis_collective(i = i, 
                                                mode = mode,
                                                alpha = 0.05, 
                                                pos_anchor = 200, 
                                                ratio_filter = ratio_filter)

                df, mu_pct, sigma_pct, exp, obs, chi, chi_crit, p_val = out

                print("\n\n")
                print("--------------------------------------------------")
                print(f"\n\nChi-Square Test {s} | i = {i}\n")
                print("--------------------------------------------------")
                print(df)

                if save_csv:
                    df.to_csv(path + f"{prefix}_{plot_name}_{s}_{i}.csv")

                data = [[obs[0], exp[0]], [obs[1], exp[1]]]

                ax[plot_i].axis('tight')

                ax[plot_i].axis('off')

                ax[plot_i].table(cellText=data, 
                        colLabels=cols, 
                        rowLabels=rows, 
                        loc='center', 
                        cellLoc='center')  
                
                mu_pct = round(100*mu_pct, 2)
                sigma_pct = round(100*sigma_pct, 2)
                chi = round(chi, 4)
                chi_crit = round(chi_crit, 2)
                p_val = round(p_val, 2)

                ttl  = "\n------------------------------------\n"
                ttl += f"Data: {s} | i = {i}\n"
                ttl += fr"$\mu={mu_pct}\%$ | $\sigma={sigma_pct}\%$" + "\n"
                ttl += fr"$\mu + 3 \sigma = {mu_pct + 3*sigma_pct}\%$" + "\n"
                ttl += fr"$\chi^2={chi}$ | $\chi^2_{{crit}}={chi_crit}$ | $p_{{val}}={p_val}$"
                
                ax[plot_i].set_title(ttl, fontsize=4)

                r1 = [s, i, "Observed", obs[0], obs[1], mu_pct, sigma_pct, chi, chi_crit, p_val]
                r2 = [s, i, "Expected", exp[0], exp[1], mu_pct, sigma_pct, chi, chi_crit, p_val]

                rows_main.append(r1)
                rows_main.append(r2)

                plot_i += 1                    
                
        columns = ["System", "i", "Type", "Inside", 
                   "Outside", "mu_pct", "sigma_pct", 
                   "chi", "chi_crit", "p_val"]
        
        df_main = pd.DataFrame(rows_main, columns = columns)

        print(df_main)

        df_main.to_csv(path + f"{prefix}_{plot_name}_MAIN.csv")

        plt.tight_layout()

        plt.savefig(path + f"{prefix}_{plot_name}.png")

    
        
    global pp
        
    path = pp + "Enrichment_Testing_Collective/"
    ensure_directory(path)

    main(path, "full",     ratio_filter = 1.0, save_csv = True)
    

def SHANNON_DEFENSE_DIVERSITY(dvars, mode, i = 0):

    global pp

    path = pp + "Defense_Diversity/"
    ensure_directory(path)

    

    def get_system_lookup_dict():

        p = "/home/baslan/.macsyfinder/models/defense-finder-models/DefenseFinder_rules.tsv"

        system_lookup_dict = dict()

        txt = read_file(p, mode="whole")
        txt = re.sub(',', ' ', txt)
        txt = re.sub(' +' , ' ', txt)
        txt = re.sub('\t+', ' ', txt)
        data = txt.split('\n')

        data.pop(0)

        for line in data:

            row = line.split(' ')
            print(row)

            system = row[0]
            row.pop(0)

            row = [r for r in row if r != '']

            for gene in row:

                if gene not in system_lookup_dict.keys() :

                    if  gene[0].upper() == system[0].upper():

                        system_lookup_dict[gene] = system

        
        
        
        
        
        print(system_lookup_dict["Paris_I__AAA_15"])             

        return system_lookup_dict 

    def get_uniform_defense_dict(all_defenses):

        N = len(all_defenses)

        return {k: 1/N for k in all_defenses}
    
    def get_plasmid_defense_dict():

        path = "/Volume/biodata2" \
                "/sccmec_snail_analysis/snail_analysis_plasmids" \
                "/results/main_results_table.csv"
        
        df = pd.read_csv(path)
        df['raw def gene set'] = df['raw def gene set'].apply(ast.literal_eval)

        all_def_genes = set()
        for s in df['raw def gene set']:
            s = [d.split('_')[0] for d in s]
            print(s)
            all_def_genes = all_def_genes.union(s)

        gene_counts = {
        gene: sum(df['raw def gene set'].apply(lambda x: gene in x)) 
        for gene in all_def_genes}

        return gene_counts
    
    def get_epiderm_defense_dict():

        mode = "new"

        if mode == "old":

            d_epi_in =  {'Cas': 7, 'RM': 35, 'SEFIR': 3, 
                        'Nhi': 3, 'Stk2': 2, 
                        'AbiJ': 11, 'AbiR': 2, 
                        'Gabija': 5, 'CBASS': 2, 
                        'AbiG': 2, 'AbiH': 1, 
                        'RloC': 1, 'Rst_PARIS': 1}  
                    
            d_epi_ou =  {'Abi2': 62, 'RM': 19, 
                        'Kiwa': 2, 'AbiD': 2, 'RloC': 2, 
                        'Gabija': 3, 'AbiQ': 2, 'AbiH': 1, 
                        'PD-T4-9': 1, 'AbiK': 2, 'AbiG': 1}
            
            d_epi = Counter(d_epi_in) + Counter(d_epi_ou)

            return dict(d_epi)
        
        if mode == "new":

            dvars = load_dvars(pr)
            dvar = dvars['epidermidis']
            DD = DefMapDist(dvar)
            adg = DD.generate_defense_finder_raw_genes_df()

            return adg
    
    def get_aureuss_defense_dict():

        mode = "new"

        if mode == "old":

            d_aur_in =  {'Stk2': 86, 'RloC': 126, 
                        'RM': 241, 'Gabija': 11, 
                        'PD-T7-2': 39, 'Septu': 1, 
                        'AbiJ': 54, 'Pycsar': 48, 
                        'Cas': 9, 'Nhi': 41, 
                        'Lamassu-Fam': 49, 'CBASS': 14, 
                        'Bunzi': 6, 'Gao_Iet': 3, 
                        'AbiH': 4, 'Dodola': 16, 
                        'Hachiman': 5, 'PD-T4-8': 1, 'AbiA': 1}

            d_aur_ou =  {'Retron': 190, 'Abi2': 824, 
                        'RosmerTA': 554, 'AbiJ': 58, 
                        'Gabija': 77, 'AbiD': 347, 
                        'Dodola': 260, 'PD-Lambda-1': 94, 
                        'Kiwa': 6, 'AVAST': 221, 
                        'RM': 206, 'AbiZ': 33, 'Thoeris': 61, 
                        'RloC': 2, 'AbiK': 10, 'ShosTA': 21, 
                        'Shango': 25, 'PD-T4-9': 13, 
                        'PD-Lambda-5': 9, 'SoFIC': 12, 
                        'PD-T7-2': 2, 'AbiH': 1, 'AbiP2': 8}   
                    
            d_aur = Counter(d_aur_in) + Counter(d_aur_ou)

            return dict(d_aur)
        
        if mode == "new":

            dvars = load_dvars(pr)
            dvar = dvars['aureus']
            DD = DefMapDist(dvar)
            adg = DD.generate_defense_finder_raw_genes_df()

            return adg

    def test_distribution_normality(pdf):

        return sum(pdf.values())

    j = 0

    if j == 0:

        def plot_dict(d_dict, path_plot):

            

            plt.rcParams['figure.dpi'] = 600

            
            df = pd.DataFrame(d_dict)

            
            df.plot.bar(stacked=True)

            
            plt.xticks(fontsize=6)

            plt.tight_layout()

            plt.savefig(path_plot)   

        

        d_epi = get_epiderm_defense_dict()
        d_aur = get_aureuss_defense_dict()
        d_plm = get_plasmid_defense_dict()

        all_defenses = set()
        for dictionary in [d_epi, d_aur, d_plm]:

            for key in dictionary.keys(): all_defenses.add(key)

        d_unf = get_uniform_defense_dict(all_defenses)

        for key in all_defenses: 

            d_epi.setdefault(key, 0)
            d_aur.setdefault(key, 0)
            d_plm.setdefault(key, 0)

        
            
        pdf_epi = calculate_probability_dist(d_epi)
        pdf_aur = calculate_probability_dist(d_aur)
        pdf_plm = calculate_probability_dist(d_plm)
        pdf_unf = calculate_probability_dist(d_unf)

        
        
        
        
        
        

        

        h_epi = calculate_shannon_entropy(pdf_epi)
        h_aur = calculate_shannon_entropy(pdf_aur)
        h_plm = calculate_shannon_entropy(pdf_plm)
        h_unf = calculate_shannon_entropy(pdf_unf)

        

        kl_epi = calculate_kl_divergence(pdf_epi, pdf_unf)
        kl_aur = calculate_kl_divergence(pdf_aur, pdf_unf)
        kl_plm = calculate_kl_divergence(pdf_plm, pdf_unf)
        kl_unf = calculate_kl_divergence(pdf_unf, pdf_unf)

        

        
        
        
        
        
        
        

        data = [[h_epi, kl_epi], 
                [h_aur, kl_aur], 
                [h_plm, kl_plm]]
        
        df = pd.DataFrame(data, columns=['H', 'KL'], 
                            index=['Epidermidis', 'Aureus', 'Plasmids'])
        
        
        df['H']  = df['H'].apply(lambda x: round(x, 2))
        df['KL'] = df['KL'].apply(lambda x: round(x, 2))

        print(df)

        d_dict = {'Epidermidis': pdf_epi,
                    'Aureus'     : pdf_aur,
                    'Plasmids'   : pdf_plm,
                    'Uniform'    : pdf_unf}

        d_dict = {'Epidermidis': pdf_epi,
                    'Aureus'     : pdf_aur,
                    'Plasmids'   : pdf_plm}
        
        path_plot = path + "mixed_dist.png"
        plot_dict(d_dict, path_plot)  

    if j == 1: 
        
        p_plm = get_plasmid_defense_dict()
        for k,v in p_plm.items():
            print(f"{k}:".ljust(10), v)

    if j == 3:

        print("section: ", u,j)

        
            
        

        dvars = load_dvars(pr)
        dvar = dvars['epidermidis'] 
        DD = DefMapDist(dvar)
        adg = DD.generate_defense_finder_raw_genes_df()

        
        print()
        for k,v in adg.items(): print(k, v)

def PLOT_ENRICHMENT_TESTING(dvars, mode, plot_type='bar', i=0):

    for s in ['epidermidis', 'aureus']: 

        SA = SnailAnalyzerFull(dvars[s])

        thresholds = [1, 2, 3, 4, 5]
        percentages = []
        nn = []
        MM = []
        dd = []

        for th in thresholds:

            df, pct, n, M = SA.statistical_analysis_per_genome(
                            i = i, 
                            mode  = mode,
                            thres = th,
                            alpha = 0.05, 
                            pos_anchor = 200)
            
            percentages.append(pct)
            nn.append(n)
            MM.append(M)
            dd.append(M-n)

        _, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=250)

        if plot_type == 'bar':
            
            ax.bar(thresholds, 
                   nn, 
                   width=0.5, 
                   label='enriched', 
                   color='blue')

            ax.bar(thresholds, 
                   dd, 
                   width=0.5, 
                   bottom=nn, 
                   label='not enriched', 
                   color='lightblue')
            
            ax.set_xticks(thresholds)

            ax.legend(loc='upper right', fontsize=8)
            plt.tight_layout()
            plt.savefig("enrichment_testing_bar.png")

        else:

            ax2=ax.twinx()

            ax.plot(thresholds, percentages, marker='o', 
                    color='black', label='Percentage')
            
            ax2.plot(thresholds, MM, marker='o', 
                    color='red',  label='
            
            ax.set_xlabel("Minimum Number of Defense Systems Required")
            ax.set_ylabel("Percentage of Genomes Enriched")
            ax2.set_ylabel("Number of Surviving Genomes")

            ax.set_ylim(0, 1.1)
            

            ax.set_xticks(thresholds)

            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()   
            ax.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=8)

            plt.tight_layout()
            plt.savefig(f"enrichment_testing_line_{s}.png")

def DO_ENRICHMENT_STATS_per_genome(dvars, mode):

    

    thres = 2

    txt = ""
    for s in ['epidermidis', 'aureus']:

        SA = SnailAnalyzerFull(dvars[s])

        for i in [0, 2]:

            df, percentage, n, M = SA.statistical_analysis_per_genome(i = i, 
                                    mode = mode,
                                    thres = thres,
                                    alpha = 0.05, 
                                    pos_anchor = 200)
            
            pct = round(100*percentage, 2)
            
            df.to_csv(f"enrichment_testing_results_i_{i}_{s}.csv")
            txt += f"{s} - {i}: {pct}%\n"

    write_file("enrichment_testing_summary.txt", txt)

def PLOT_ENRICHMENT_TESTING_pie(dvars, mode, plot_type='bar', i=0):

    from matplotlib import patches

    def autopct_format(values):
        def my_format(pct):
            total = sum(values)
            val = int(round(pct*total/100.0))
            return '{v:d} ({p:.2f}%)'.format(p=pct,v=val)
        return my_format

    _, ax = plt.subplots(2, 1, figsize=(5, 8), dpi=250)

    for k, s in enumerate(['epidermidis', 'aureus']): 

        SA = SnailAnalyzerFull(dvars[s])

        df, pct, n, M = SA.statistical_analysis_per_genome(
                        i = i, 
                        mode  = mode,
                        thres = 0,
                        alpha = 0.05, 
                        pos_anchor = 200)

        df['total N'] = df['region N'] + df['other N']

        df0 = df[df['total N'] >= 2]
        df1 = df0[df0['rejected'] == True]

        a = len(df)
        b = len(df0)
        c = len(df1)

        title = ''
        title += f"Enrichment Testing: {s}\n"
        title += f"Total Genomes: {a}\n"
        title += f"Genomes with >= 2 Defenses: {b}\n"
        title += f"Enriched Genomes: {c}\n"

        labels   = '<= 1', 'Not Enriched', 'Enriched'
        colors   = ['white', 'lightblue', 'blue']
        sizes = [a-b, b-c, c]     

        ax[k].pie(sizes, labels=labels, 
               autopct='', 
               startangle=90, colors=colors, 
               wedgeprops=dict(edgecolor='black')) 

        arc = patches.Arc((0, 0), 2.1, 2.1, theta1=90+360*(a-c)/a, theta2=90, edgecolor='black', lw=4)

        ax[k].set_title(title, fontsize=8)
        ax[k].add_patch(arc)
        ax[k].axis('equal')  

    plt.tight_layout()
    plt.savefig("enrichment_testing.png")

def PLOT_ENRICHMENT_TESTING_scatter_v0(dvars, mode, i = 0):

    from dataclasses import dataclass

    s =  'aureus'

    @dataclass
    class Point:

        accession: str
        enriched: bool
        ratio_nuc: float
        ratio_def: float
        total_def: int

    SA = SnailAnalyzerFull(dvars[s])

    df, _, _, _ = SA.statistical_analysis_per_genome(
                  i = i, 
                  mode = mode,
                  thres = 0,
                  alpha = 0.05, 
                  pos_anchor = 200)
    
    Data = []
    for accession, row in df.iterrows():
        
        total_def = row['region N'] + row['other N']

        enriched = row['rejected']
        ratio_nuc = row['nuc_len'] / row['gnm_len']
        if total_def == 0: 
            ratio_def = -0.25
        else:
            ratio_def = row['region N'] / total_def

        P = Point(accession, enriched, ratio_nuc, ratio_def, total_def)
        Data.append(P)

    _, ax = plt.subplots(1, 1, figsize=(8, 4), dpi=250)

    for point in Data:
        color = 'red' if point.enriched else 'blue'
        if point.total_def == 0: color = 'black'
        ax.scatter(point.ratio_nuc, point.ratio_def, color=color, marker='.')

    ax.set_xlim(0, 0.06)
    ax.set_xlabel('ratio_nuc')
    ax.set_ylabel('ratio_def')
    plt.tight_layout()
    plt.savefig(f"scatter_plot_{s}.png")

def PLOT_ENRICHMENT_TESTING_scatter(dvars, mode, i = 0):

    from dataclasses import dataclass

    s =  'epidermidis'

    @dataclass
    class Point:

        accession: str
        enriched: bool
        ratio_nuc: float
        ratio_def: float
        total_def: int

    SA = SnailAnalyzerFull(dvars[s])

    df, _, _, _ = SA.statistical_analysis_per_genome(
                  i = i, 
                  mode = mode,
                  thres = 1,
                  alpha = 0.05, 
                  pos_anchor = 200)
    
    Data = []
    for accession, row in df.iterrows():
        
        total_def = row['region N'] + row['other N']

        if total_def > 0: 

            enriched = row['rejected']
            ratio_nuc = row['nuc_len'] / row['gnm_len']
            ratio_def = row['region N'] / total_def

            if ratio_def == 0: ratio_def = random.uniform(-0.25, 0)

            P = Point(accession, enriched, ratio_nuc, ratio_def, total_def)
            Data.append(P)

    _, ax = plt.subplots(1, 1, figsize=(4, 8), dpi=250)

    for point in Data:
        color = 'red' if point.enriched else 'blue'
        ax.scatter(point.ratio_nuc, point.ratio_def, color=color, marker='.')

    
    ax.set_xlabel('ratio_nuc')
    ax.set_ylabel('ratio_def')
    plt.tight_layout()
    plt.savefig(f"scatter_plot_{s}.png")

class SeqX:

    def __init__(self, path_genome):

        self.record = SeqIO.read(path_genome, "gb")
        self.seq = self.record.seq
        self.seq_len = len(self.seq)
        self.pos = list(range(self.seq_len))

    def realign_and_reorient(self, rlmh_loc):

        a, b, strand = rlmh_loc

        if strand == 1: 

            self.shift(a)

        else:

            self.shift(b)
            self.reverse_complement()

    def get_anuclen(self, rlmh_loc, acrb_loc):

        f = calculate_nucleotide_distance
        span = f(rlmh_loc, acrb_loc, self.record.seq)
        
        return span  

    
    def position_is_intergenic(self, x0, x1):

        
        
        

        is_intergenic = True
        container = None
        for feature in self.record.features:

            proper = feature.type == "CDS"

            if (proper and x0 in feature) or (proper and x1 in feature):
                is_intergenic = False
                container = feature

        return is_intergenic, container   

    def shift(self, n):

        self.seq = self.seq[n:] + self.seq[0:n]
        self.pos = self.pos[n:] + self.pos[0:n]

    def reverse_complement(self):

        self.seq = self.seq.reverse_complement()
        self.pos = list(reversed(self.pos))

    def old2new(self, x):

        return self.pos.index(x)

    def new2old(self, x):

        return self.pos[x]

    def get_segment_with_old_coordinates(self, a, b):

        a_new = self.old2new(a)
        b_new = self.old2new(b)

        if a_new > b_new:
            a_new, b_new = b_new, a_new

        return self.seq[a_new:b_new]
    
    def save_seq_as_fasta(self, fasta_id, path_save):

        SecRec = SeqRecord(self.seq, id=fasta_id, name= "...")
        SeqIO.write(SecRec, path_save, "fasta")

    @staticmethod
    def selftest2():

        
        

        s = "ATGCTAGCTACGAATCAGGATCAATTGCACAAGCTAACTAGAGTTCCATGGCATTTGCACCACGCT"
        a, b, c, d = 17, 38, 5, 14

        t0 = time.perf_counter()

        seqx = SeqX(s)
        seqx.shift(b)
        seqx.reverse_complement()

        print('org rlmh:', s[a: b])
        print('org seq :', s)
        print('new seq :', seqx.seq)

        print('new rlmh:', seqx.get_segment_with_old_coordinates(b-5-1, b-1))

        print(seqx.seq_len)
        print(seqx.new2old(0))
        print(seqx.new2old(seqx.old2new(38)))

        t1 = time.perf_counter()
        print('secs', round(t1-t0, 3))

class CutSite:

    def __init__(self, start, stop, score, pvalue, qvalue, strand, HitSeq: Seq, is_accessory: bool, is_intergenic: bool):

        self.start    = start
        self.stop     = stop
        self.score    = score
        self.pvalue   = pvalue
        self.qvalue   = qvalue
        self.strand   = strand
        self.HitSeq   = HitSeq
        self.is_accessory = is_accessory
        self.is_intergenic = is_intergenic

        pid  = f"{hash(str(HitSeq))}"
        desc = f"pos: {start}, strand: {strand}, " \
               f"score: {score}, pvalue: {pvalue}, " \
               f"qvalue: {qvalue}, " \
               f"is_accessory: {is_accessory}"

        self.SecRec = SeqRecord(HitSeq, id=pid, name= "...", description=desc)

    def __str__(self):

        output = f"Seq: {self.HitSeq}, " \
                 f"Score: {self.score}, " \
                 f"p-value: {self.pvalue}, " \
                 f"q-value: {self.qvalue}, " \
                 f"Start: {self.start}, " \
                 f"Stop: {self.stop}, " \
                 f"Strand: {self.strand}, " \
                 f"is_accessory: {self.is_accessory}, " \
                 f"is_intergenic: {self.is_intergenic}"
        
        return output

class CutSiteList(list):

    def __init__(self, *args, eps=200000):

        super().__init__(args)
        self.eps = eps
        self.num_cassettes = None
        self.median_dist = None
        self.pairs = []
        self.distances = []

    def __repr__(self):
        
        txt = "\n\n"
        for cs in self:
            txt += f"{cs}\n"
        return txt
    
    def identify_pairs(self):

        

        df: pd.DataFrame = self.df

        

        df_fwd = df[df['strand'] == 1].copy()
        df_rev = df[df['strand'] == -1].copy()

        df_fwd.sort_values('start', inplace=True)
        df_rev.sort_values('start', inplace=True)

        df_fwd['dist'] = df_fwd['start'].diff()
        df_rev['dist'] = df_rev['start'].diff()

        df_fwd['is_pair'] = df_fwd['dist'] <= self.eps
        df_rev['is_pair'] = df_rev['dist'] <= self.eps

        

        distances_fwd = df_fwd[df_fwd['is_pair']]['dist'].to_numpy()
        distances_rev = df_rev[df_rev['is_pair']]['dist'].to_numpy()

        self.distances = np.concatenate([distances_fwd, distances_rev])
        self.median_dist = median(self.distances) if any(self.distances) else -self.eps
        self.num_cassettes = len(self.distances)

        print("\nFWD:")
        print(df_fwd)
        print("\n\n")
        print("\nREV:")
        print(df_rev)
        print("\n\n")
        print(self.distances)
        print("\n\n")
        print(self.median_dist)
        print("\n\n")
        
    def get_stats(self, mode = 'pvalue'):

        v = [getattr(cs, mode) for cs in self]

        return min(v), max(v), median(v)
    
    @property
    def get_SecRecs(self):

        return [cs.SecRec for cs in self]
    
    @property
    def get_HitStrings(self):

        return repr([str(cs.HitSeq) for cs in self])
    
    @property
    def get_len(self):

        return len(self)
    
    @property
    def df(self):

        rows = []
        cols = ['start', 
                'stop', 
                'score', 
                'pvalue', 
                'qvalue', 
                'strand', 
                'HitSeq', 
                'is_accessory', 
                'is_intergenic']
        
        for cs in self: 
            
            rows.append([cs.start, 
                         cs.stop, 
                         cs.score, 
                         cs.pvalue, 
                         cs.qvalue, 
                         cs.strand, 
                         str(cs.HitSeq), 
                         cs.is_accessory, 
                         cs.is_intergenic])
            
        return pd.DataFrame(rows, columns=cols).sort_values('qvalue')
    
class CutSiteMap:

    def __init__(self, AccObj, plot_params, path_genomes, 
                 path_sites, dist_data, threshold, 
                 putative_cassette_length):
        
        self.AccObj       = AccObj
        self.plot_params  = plot_params
        self.path_genomes = path_genomes
        self.path_sites   = path_sites
        self.dist_data    = dist_data
        self.threshold    = threshold

        self.putative_cassette_length = putative_cassette_length
        
    def analyze(self, pr, alpha=5e-6):

        def calculate_pvalue(score, score_vec, CDF):
        
            idx = np.searchsorted(score_vec, score)
            return 1 - CDF[idx]            

        alpha     = 1e-7
        self.threshold = 5

        eps        = self.putative_cassette_length
        core_depth = self.plot_params['c_depth']
        max_data   = self.plot_params['n']

        sf = SnailFeature(self.AccObj, 
                               200 - core_depth, 
                               max_data,
                               core_depth, 
                               self.plot_params['m'],
                               self.plot_params['th'])
        
        
        rlmh_loc, acrb_loc = sf.get_boundary_positions()

        myseqx        = SeqX(f"{self.path_genomes}{self.AccObj.accession}.gb")
        anuclen       = myseqx.get_anuclen(rlmh_loc, acrb_loc)
        myseqx.realign_and_reorient(rlmh_loc)

        
        

        pssm      = self.dist_data['pssm']
        CDF       = self.dist_data['CDF']
        score_vec = self.dist_data['x']

        motif_length = pssm.length
        search_res = pssm.search(myseqx.seq, threshold=self.threshold)
        
        M = len(myseqx.seq)
        csl = []
        for x, score in search_res:

            if x < 0:

                strand = -1
                x = M + x
                x0, x1 = int(x), int(x + motif_length)
                HitSeq = myseqx.seq[x0: x1].reverse_complement()

            else: 
                
                strand = 1
                x0, x1 = int(x), int(x + motif_length)
                HitSeq = myseqx.seq[x0: x1]

            csl.append(CutSite(start=x0, 
                               stop = x1,
                               score=score,  
                               pvalue=calculate_pvalue(score, score_vec, CDF),
                               qvalue=None,
                               strand=strand, 
                               HitSeq=HitSeq, 
                               is_accessory=x1 <= anuclen, 
                               is_intergenic = False))

        
            
        
        

        csl.sort(key=lambda x: x.pvalue)
        n = len(csl)
        i = np.arange(1, n+1)
        c = alpha * i / n

        for idx, obj in enumerate(csl):
            obj.i = i[idx]
            obj.c = c[idx]
            obj.rejected = obj.pvalue < obj.c

        
        
        rows_bhc = []
        for cs in csl:
            rows_bhc.append([str(cs.HitSeq), cs.score, cs.pvalue, cs.i, cs.c, cs.rejected])

        df_bhc = pd.DataFrame(rows_bhc, columns=['HitSeq', 'Score', 'p-value', 'i', 'c', 'rejected'])
        print(df_bhc)

        

        fig, ax = plt.subplots(3, 1, figsize=(5, 5))
        ax2 = ax[0].twinx()

        col1 = 'p-value'
        col2 = 'Score'

        ax[0].plot(df_bhc['i'], df_bhc[col1], label=col1, color='black')
        ax2.plot(df_bhc['i'], df_bhc[col2], label=col2, color='orange')
        ax[0].set_xlabel('rank')
        ax[0].set_ylabel(col1)
        ax2.set_ylabel(col2)
        
        ax[1].plot(self.dist_data['x'][:-1], self.dist_data['PVL'], label='p-value', color='black')
        ax[1].set_xlabel('Score')
        ax[1].set_ylabel('p-value')
        ax[1].set_xlim(5, 20)
        ax[1].set_ylim(0, 1e-5)

        ax[2].scatter(df_bhc[col2], df_bhc[col1], label=col1, color='black', s=0.1)

        ax[0].legend()
        plt.tight_layout()
        plt.savefig('asdf.png')
        
        
                    
        csl_i = CutSiteList(eps = eps)
        csl_o = CutSiteList(eps = eps)

        for cs in csl:

            if cs.is_accessory: csl_i.append(cs)
            else: csl_o.append(cs)

        csl_i.identify_pairs()
        csl_o.identify_pairs()

        res = {'accession': self.AccObj.accession,
               'rlmh_translation': myseqx.seq[0:150].translate(),
               'rlmh_loc': rlmh_loc,
               'acrb_loc': acrb_loc,
               'anuclen': anuclen, 
               'ACRB': sf.accessory_region_boundary,
               'eps': eps,
               'cls_i_len': csl_i.get_len,
               'csl_i_hits': csl_i.get_HitStrings,
               'csl_i': csl_i, 
               'cls_i_num_cassettes': csl_i.num_cassettes, 
               'csl_i_median_dist': csl_i.median_dist,
               'csl_i_median_dist/eps': csl_i.median_dist/eps,
               'cls_o_len': csl_o.get_len,
               'csl_o_hits': csl_o.get_HitStrings,
               'csl_o': csl_o, 
               'cls_o_num_cassettes': csl_o.num_cassettes,
               'csl_o_median_dist': csl_o.median_dist,
               'csl_o_median_dist/eps': csl_o.median_dist/eps,
               }
        
        for k,v in res.items(): 

            if not isinstance(v, str): v = repr(v)
            print(f"{k.ljust(22).upper()}: {v}")
            
        row = [res['accession'],
               res['eps'],
               res['cls_i_len'],
               res['cls_o_len'],
               res['csl_i_median_dist'],
               res['cls_i_num_cassettes'],
               res['csl_o_median_dist'],
               res['cls_o_num_cassettes'],
               res['mdfb'], 
               res['anuclen'],
               res['ACRB'],
               res['csl_i_hits'],
               res['csl_o_hits']]
        
        return row, csl_i.get_SecRecs, csl_o.get_SecRecs
    
class CutSiteAnalyzer:

    def __init__(self, dvar, mode="meme", iteration = 0):

        self.dvar = dvar
        self.plot_params  = dvar['plot_params']
        self.window_size  = dvar['win_size']
        self.numthreads   = dvar['numthreads']
        self.path_genomes = dvar['genomes']
        self.path_root    = dvar['root']
        self.path_nbrhd   = dvar['root'] + "neighborhood_files_ws{}/".format(self.window_size)
        self.path_snaildt = dvar['root'] + "snail_data"
        self.iteration    = iteration

        self.path_combined_fasta      = self.path_root + "combined.fasta"

        self.file_manager = nbFileManager(self.path_nbrhd, self.path_genomes)
        self.accessions = sorted(list(self.file_manager.accessions))

        self.csl_i_collection = []
        self.csl_o_collection = []
        self.putative_cassette_length = 200000
        self.cut_site_score_threshold = 10 

        self.site_dic = {'GAAGCATATCATAAATGA': 'A',
                         'GAAGCATATCATAAATAA': 'B',
                         'GAAGGGTATCGTAAGTGA': 'C',
                         'GAAGCGTATCGTAAGTGA': 'D'}

        if mode == "custom":

            self.dist_data = self.generate_bg_Dist(pseudocounts=0.1)
            self.analyze_method = "analyze"
            self.path_root = self.path_root + "cut_site_analysis/custom/"

        elif mode == "meme":

            self.dist_data = None
            self.analyze_method = "analyze_with_meme"
            self.path_root = self.path_root + f"cut_site_analysis/meme/iteration{iteration}/"
        
        else: raise ValueError("Invalid mode")

        self.path_sites   = self.path_root + "CcrAB_DRs.fa"
        self.path_dfi     = self.path_root + "cut_site_df_i.data"
        self.path_dfo     = self.path_root + "cut_site_df_o.data"
        self.path_plot    = self.path_root + 'cut_site_plot.png'
        self.path_outside = self.path_root + 'cut_site_list_outside_for_logo.fasta'
        self.path_inside  = self.path_root + 'cut_site_list_inside_for_logo.fasta'

    def analyze(self):

        c1 = os.path.isfile(self.path_dfi)
        c2 = os.path.isfile(self.path_dfo)

        if not (c1 and c2):

            AO = self.file_manager.AccessionObjDict
            rows_i = []
            rows_o = []
            for accession in AO.keys(): 

                AccObj = AO[accession]

                CSM = CutSiteMap(AccObj, 
                                self.dvar['plot_params'], 
                                self.path_genomes, 
                                self.path_sites,
                                self.dist_data,
                                self.cut_site_score_threshold,
                                self.putative_cassette_length)
                
                stats_i, stats_o = CSM.analyze(self.path_root) 

                rows_i.append(stats_i)
                rows_o.append(stats_o)

                self.csl_i_collection.extend(stats_i['seq_sites']) 
                self.csl_o_collection.extend(stats_o['seq_sites'])

            df_i = pd.DataFrame(rows_i)
            df_o = pd.DataFrame(rows_o)
    
            with open(self.path_dfi, 'wb') as f: pickle.dump(df_i, f)
            with open(self.path_dfo, 'wb') as f: pickle.dump(df_o, f)

            SeqIO.write([SeqRecord(S, id='na ') for S in self.csl_i_collection], 
                         self.path_inside, "fasta")

            SeqIO.write([SeqRecord(S, id='na ') for S in self.csl_o_collection], 
                         self.path_outside, "fasta")
            
        else:
                
            with open(self.path_dfi, 'rb') as f: df_i = pickle.load(f)
            with open(self.path_dfo, 'rb') as f: df_o = pickle.load(f)

        df_i.drop(['seq_sites'], axis=1, inplace=True)
        df_o.drop(['seq_sites'], axis=1, inplace=True)
            
        print(df_i)
        print('\n\n\n\n')
        print(df_o)
        
    def generate_bg_dict(self):

        path_bg_dict = self.path_root + "cut_site_bg_dict.data"

        if os.path.isfile(path_bg_dict):
                
            with open(path_bg_dict, 'rb') as f: bg_dict = pickle.load(f)
        
        else:

            items = self.file_manager.AccessionObjDict.items()

            bg_dict = {}
            for accession, _ in items:

                p = f"{self.path_genomes}{accession}.gb"
                grec = SeqIO.read(p, "gb")
                bg = BackgroundModel.zeroth_order(grec.seq)
                bg_dict[accession] = bg

            with open(path_bg_dict, 'wb') as f: pickle.dump(bg_dict, f)

        
            
        df_bg = pd.DataFrame(bg_dict.values(), index=bg_dict.keys())

        bg_average_dict = {'A': df_bg['A'].mean(), 
                           'C': df_bg['C'].mean(), 
                           'G': df_bg['G'].mean(), 
                           'T': df_bg['T'].mean()}

        

        return bg_average_dict
        
    def generate_bg_Dist(self, pseudocounts = 0.01):

        def generate_plot(data, path_save):

            _, ax = plt.subplots(nrows = 2, ncols = 1, figsize=(10, 5))
            ax2 = ax[0].twinx()

            ax[0].plot(data['x'][:-1], data['pdf'])
            ax2.plot(data['x'][:-1], data['CDF'])
            ax[0].set_xlabel('Score')
            ax[0].set_ylabel('pdf')
            ax2.set_ylabel('CDF')
            ax[0].set_title('Background Score Distribution')
            ax[1].plot(data['x'][:-1], data['PVL'])
            ax[1].set_xlabel('Score')
            ax[1].set_ylabel('P-Value')
            ax[1].set_title('Background Score P-Value Distribution')
            ax[1].set_xlim([5, 20])
            ax[1].set_ylim([0, 1e-5])
            plt.tight_layout()
            plt.savefig(path_save, dpi=300)

        path_fpssm      = self.path_root + "cut_site_bg_fpssm.data"
        path_all_scores = self.path_root + "cut_site_bg_all_scores.data"
        path_data       = self.path_root + "cut_site_bg_data.data"
        path_dist_plots = self.path_root + 'cut_site_bg_dist_plots.png'

        if not os.path.isfile(path_all_scores):

            print('Generating Background Scores')

            bg_avg_dict = self.generate_bg_dict()

            ins   = [S.seq for S in SeqIO.parse(self.path_sites, "fasta")]
            motif = motifs.create(ins)
            pwm   = motif.counts.normalize(pseudocounts=pseudocounts)
            fpssm = pwm.log_odds(bg_avg_dict)
            rpssm = fpssm.reverse_complement()
            items = self.file_manager.AccessionObjDict.items()

            print()
            print(f"motif consensus: {motif.consensus}")
            print(f"motif length: {motif.length}")
            print()

            all_scores = np.array([])
            for accession, _ in items:

                p = f"{self.path_genomes}{accession}.gb"
                print(p)

                grec = SeqIO.read(p, "gb")
                scores_fwd = fpssm.calculate(grec.seq)
                scores_rev = rpssm.calculate(grec.seq)

                all_scores = np.concatenate((all_scores, scores_fwd))
                all_scores = np.concatenate((all_scores, scores_rev))

            with open(path_fpssm, 'wb') as      f: pickle.dump(fpssm, f)
            with open(path_all_scores, 'wb') as f: pickle.dump(all_scores, f)

        if not os.path.isfile(path_data):

            print('Generating Background Score Distribution')

            with open(path_fpssm, 'rb') as      f: fpssm = pickle.load(f)
            with open(path_all_scores, 'rb') as f: all_scores = pickle.load(f)

            
            all_scores = all_scores[~np.isnan(all_scores)] 

            mmin       = int(all_scores.min()-1)
            mmax       = int(all_scores.max()+1)
            dx         = 0.1
            x          = np.arange(mmin, mmax, dx)
            pdf, _     = np.histogram(all_scores, bins=x,density=True)
            CDF        = np.cumsum(pdf) * dx
            PVL        = [(1 - CDF[i]) for i, x0 in enumerate(x[:-1])]

            data       = {'min' : mmin, 
                          'max' : mmax, 
                          'x'   : x,
                          'dx'  : dx,
                          'pdf' : pdf, 
                          'CDF' : CDF,
                          'PVL' : PVL,
                          'pssm': fpssm}
            
            with open(path_data, 'wb') as f: pickle.dump(data, f)
            generate_plot(data, path_dist_plots)

            return data
        
        with open(path_data, 'rb') as f: data = pickle.load(f)

        return data
    
class CutSiteAnalyzer_NEW:

    def __init__(self, dvar, folder_name):

        self.plot_params   = dvar['plot_params']
        self.window_size   = dvar['win_size']
        self.numthreads    = dvar['numthreads']
        self.path_genomes  = dvar['genomes']
        self.path_nbrhd    = dvar['root'] + f"neighborhood_files_ws{self.window_size}/"
        self.path_snaildt  = dvar['root'] + "snail_data"
        self.pr            = dvar['root'] + f"cut_site_analysis/meme/{folder_name}/"

        with open(self.pr + 'inputs_sheet.json', 'r') as f:

            self.params = json.load(f)
            print(self.params)

        self.file_manager = nbFileManager(self.path_nbrhd, self.path_genomes)

        self.path_dfi       = self.pr + "df_i.data"
        self.path_dfo       = self.pr + "df_o.data"
        self.path_dfi_csv   = self.pr + "df_i.csv"
        self.path_dfo_csv   = self.pr + "df_o.csv"
        self.path_plot      = self.pr + "plot.png"
        self.path_outside   = self.pr + "outside.fasta"
        self.path_inside    = self.pr + "inside.fasta"
        self.path_instances = self.pr + "all.fasta"
        self.labels         = self.pr + "sites_labels.txt"
        
        self.path_logo_i    = self.pr + "logo_inside.png"
        self.path_logo_o    = self.pr + "logo_outside.png"
        self.path_logo_all  = self.pr + "logo_all.png"

    def step0_create_FIMO_motif(self):

        
        

        w = 18
        pr = self.pr
        p_ins = pr + "motif_inputs/all.fasta" 
        p_bck = pr + "motif_inputs/zero-th.txt"
        p_out = pr + "motif/"

        

        cmd0 = f"meme {p_ins} -dna " \
               f"-oc {p_out} " \
               f"-mod zoops -nmotifs 1 -w {w} " \
               f"-markov_order 0 " \
               f"-bfile {p_bck}"

        subprocess.check_call(cmd0, shell=True)

        
        
        out  = p_out + "/meme2meme.txt"
        res  = p_out + "/meme.txt"
        cmd1 = f"meme2meme {res} > {out}"
         
        subprocess.check_call(cmd1, shell=True)

    def step_1_analyze(self):

        def get_minimal_df(df):

            dfc = df.copy()
            col = 'seq_sites'
            drop_list = ['distances', 'median_dist', 'ACRB', 'anuclen', 'eps']
            
            dfc[col] = dfc[col].apply(lambda x: [str(y) for y in x])
            dfc.drop(drop_list, axis=1, inplace=True)

            cols = list(dfc.columns)
            cols.remove(col)
            cols.append(col)
            dfc = dfc[cols]

            

            return dfc

        def get_stats(df, accession, eps, anuclen, ACRB):

            df = df.sort_values(['strand', 'start'])
            df['dist'] = df.groupby('strand')['start'].diff()
            df['is_pair'] = df['dist'] <= eps  

            distances = df[df['is_pair']]['dist'].to_numpy()
            median_dist = median(distances) if any(distances) else -eps
            hits = [Seq(x) for x in df['matched_sequence']]

            output = {'accession'    : accession,
                    'eps'          : eps,
                    'anuclen'      : anuclen,
                    'ACRB'         : ACRB,
                    'num_sites'    : len(df),
                    'seq_sites'    : hits,
                    'distances'    : distances,
                    'median_dist'  : median_dist,
                    'num_cassettes': len(distances),
                    'max_pvalue'   : df['p-value'].max(),
                    'min_pvalue'   : df['p-value'].min(),
                    'max_qvalue'   : df['q-value'].max(),
                    'min_qvalue'   : df['q-value'].min()}
            
            print()
            print('---------------------')
            print(df)
            print()

            return output

        pseudo = self.params['pseudo']
        thres1 = self.params['thres1'] 
        key    = self.params['thres2'] 
        val    = self.params['value2'] 
        eps    = self.params['eps'] 

        pr  = self.pr
        AO  = self.file_manager.AccessionObjDict
        ppr = self.plot_params
        pg  = self.path_genomes

        cd       = ppr['c_depth']
        max_data = ppr['n']

        rows_i = []
        rows_o = []
        csl_i  = []
        csl_o  = []

        for accession in AO.keys(): 

            ensure_directory(pr + "analysis/")
            ensure_directory(pr + f"analysis/{accession}/")
            
            p_vars     = pr + f"analysis/{accession}/{accession}.vars"
            p_fasta    = pr + f"analysis/{accession}/{accession}.fasta"
            
            path_motif = pr + "motif//meme2meme.txt"
            path_fimo  = pr + f"analysis/{accession}/{accession}.fimo_out/"
            path_tsv   = path_fimo + "/fimo.tsv"
            
            cmd = f"fimo " \
                    f"--oc {path_fimo} " \
                    f"--motif-pseudo {pseudo} " \
                    f"--bfile motif-file " \
                    f"--thresh {thres1} " \
                    f"--qv-thresh " \
                    f"--max-stored-scores 10000000 " \
                    f"{path_motif} " \
                    f"{p_fasta}"
            
            
            

            if not os.path.isfile(p_vars):

                sf = SnailFeature(AO[accession], 200-cd, max_data, cd, ppr['m'], ppr['th'])
                rlmh_loc, acrb_loc = sf.get_boundary_positions()
                ACRB               = sf.accessory_region_boundary

                myseqx  = SeqX(f"{pg}{accession}.gb")
                anuclen = myseqx.get_anuclen(rlmh_loc, acrb_loc)
                myseqx.realign_and_reorient(rlmh_loc)
                myseqx.save_seq_as_fasta(accession, p_fasta)
                
                

                subprocess.check_call(cmd, shell=True)
                df = pd.read_csv(path_tsv, sep='\t', comment='#')

                
                
                variabless = {'rlmh_loc': rlmh_loc,
                                'acrb_loc': acrb_loc,
                                'ACRB': ACRB,
                                'anuclen': anuclen,
                                'myseqx': myseqx,
                                'df': df}
                
                with open(p_vars, 'wb') as f: pickle.dump(variabless, f)

            else:

                with open(p_vars, 'rb') as f: variabless = pickle.load(f)
                rlmh_loc = variabless['rlmh_loc']
                acrb_loc = variabless['acrb_loc']
                ACRB     = variabless['ACRB']
                anuclen  = variabless['anuclen']
                myseqx   = variabless['myseqx']
                df       = variabless['df']

            mkitil = myseqx.seq[0:150].translate()
            print(mkitil)

            dlist = ['motif_id', 'motif_alt_id', 'sequence_name']
            df.drop(dlist, axis=1, inplace=True)

            df = df[df[key] <= val].copy()

            df['start'] = df['start'] - 1
            df['is_accessory'] = df['stop'] <= anuclen

            dfi = df[ df['is_accessory']].copy()
            dfo = df[~df['is_accessory']].copy()

            stats_i = get_stats(dfi, accession, eps, anuclen, ACRB)
            stats_o = get_stats(dfo, accession, eps, anuclen, ACRB)

            
            

            rows_i.append(stats_i)
            rows_o.append(stats_o)

            csl_i.extend(stats_i['seq_sites'])
            csl_o.extend(stats_o['seq_sites'])

        df_i = pd.DataFrame(rows_i)
        df_o = pd.DataFrame(rows_o)

        df_i.sort_values(by=['ACRB', 'accession'], inplace=True)
        df_i.insert(0, 'i', range(1, 1 + len(df_i)))
        df_o.sort_values(by=['ACRB', 'accession'], inplace=True)
        df_o.insert(0, 'i', range(1, 1 + len(df_o)))

        
        
        

        with open(self.path_dfi, 'wb') as f: pickle.dump(df_i, f)
        with open(self.path_dfo, 'wb') as f: pickle.dump(df_o, f)

        get_minimal_df(df_i).to_csv(self.path_dfi_csv, index=False)
        get_minimal_df(df_o).to_csv(self.path_dfo_csv, index=False)

        

        records_i = [SeqRecord(S, id=f'inside_{i}', description='') for i, S in enumerate(csl_i)]
        SeqIO.write(records_i, self.path_inside, "fasta")

        records_o = [SeqRecord(S, id=f'outside_{i}', description='') for i, S in enumerate(csl_o)]
        SeqIO.write(records_o, self.path_outside, "fasta")

        records = records_i + records_o
        SeqIO.write(records, self.path_instances, "fasta")

        print(f"
        print(f"

        
        
        

        cmds = [f"weblogo -F png --resolution 1200 -c classic < {self.path_inside} > {self.path_logo_i}",
                f"weblogo -F png --resolution 1200 -c classic < {self.path_outside} > {self.path_logo_o}",
                f"weblogo -F png --resolution 1200 -c classic < {self.path_instances} > {self.path_logo_all}"]
                
        for cmd in cmds: subprocess.check_call(cmd, shell=True)

        
        
        
        
        
        
  
    def step_2_plot(self):

        max_val = 22
        max_val_ext = max_val + 5

        df_i: pd.DataFrame
        df_o: pd.DataFrame

        with open(self.path_dfi, 'rb') as f: df_i = pickle.load(f)
        with open(self.path_dfo, 'rb') as f: df_o = pickle.load(f)

        n_rows = len(df_i.index)
        n_cols = len(df_i.columns)
        dummy_row = [0 for _ in range(n_cols)]

        df_i.loc[n_rows] = dummy_row
        df_o.loc[n_rows] = dummy_row
        
        N = len(df_i.index)

        sites_i = df_i['num_sites'].to_list()
        sites_o = df_o['num_sites'].to_list()

        print(f'Totals -> Inside Sites: {sum(sites_i)}, Outside Sites: {sum(sites_o)}')

        width = 2 * np.pi / N
        dtheta = width/2

        indexes = list(range(0, N))
        angles = np.array([element * width for element in indexes])

        plt.figure(figsize=(20, 10), dpi=400)
        ax = plt.subplot(111, polar=True)

        bar1_params = {
            'x': angles + dtheta,
            'height': sites_i,
            'width': width,
            'bottom': 0,
            'linewidth': 0.5,
            'edgecolor': "#A6CDDE",
            'color': "#A6CDDE"
        }

        bar2_params = {
            'x': angles + dtheta,
            'height': sites_o,
            'width': width,
            'bottom': sites_i,
            'linewidth': 0.5,
            'edgecolor': "#05437E",
            'color': "#05437E"
        }

        bar1 = ax.bar(**bar1_params)
        bar2 = ax.bar(**bar2_params)
        bars: BarContainer = bar1 + bar2
        
        il=0
        label_records = []

        zz = zip(angles, df_o["accession"])
        for angle, label in zz:

            il += 1
            label_records.append('({}) {}'.format(il, label))

            c = 'black'
            if isinstance(label, str):
                if 'NC_002976' in label: c = 'red'

            rotation = np.rad2deg(angle)

            if np.pi / 2 <= angle < 3 * np.pi / 2:
                alignment = "right"
                rotation = rotation + 180

            else:
                alignment = "left"

            ax.text(
                x=angle + dtheta,
                y=max_val_ext + 0.5,
                s=il, 
                ha=alignment,
                va='center',
                color = c,
                fontweight='normal',
                fontsize=13,
                rotation=rotation,
                rotation_mode="anchor")

        zz = zip(angles, sites_i, sites_o, df_o['num_cassettes'], df_o['median_dist'])
        for angle, si, so, tt, medd in zz:

            ss = ''
            rotation = np.rad2deg(angle)

            if np.pi / 2 <= angle < 3 * np.pi / 2:

                alignment = "right"
                rotation = rotation + 180
                
                if medd != 0: ss = (u"\u25CF" + '') * tt

            else:

                alignment = "left"
                
                if medd != 0: ss = (u"\u25CF" + '') * tt

            ax.text(
                x=angle + dtheta,
                y= si + so + 0.5, 
                s= ss,
                ha=alignment,
                va='center',
                color = 'red',
                fontsize=7,
                fontweight = 'bold',
                rotation=rotation,
                rotation_mode="anchor")

        write_file(self.labels, ', '.join(label_records))

        ax.set_xticks(angles)
        ax.xaxis.set_ticklabels([])

        ax.set_ylim([0, max_val_ext])

        major_yticks = np.arange(0, max_val_ext, 5)
        minor_yticks = np.arange(0, max_val_ext, 1)

        ax.set_yticks(minor_yticks, minor=True)
        ax.set_yticks(major_yticks)

        ax.set_rlabel_position(0)
        ax.tick_params(axis='y', which='major', labelsize=12, labelcolor = 'purple')
        for label in ax.get_yticklabels(): label.set_fontweight('bold')

        ax.tick_params(axis='y', which='minor', labelsize=0)

        ax.grid(which='minor', alpha=0.4)
        ax.grid(which='major', alpha=1, linewidth = 1)

        plt.savefig(self.path_plot, dpi=300)
        

    def create_next_iteration(self):

        pr = self.pr
        pr = pr[:-2] + str(int(pr[-2]) + 1) + '/'
        print(pr)

        ensure_directory(pr)
        ensure_directory(pr + 'analysis/')
        ensure_directory(pr + 'motif_inputs/')

        shutil.copy(self.pr + 'motif_inputs/zero-th.txt', pr + 'motif_inputs/zero-th.txt')
        shutil.copy(self.pr + 'all.fasta', pr + 'motif_inputs/all.fasta')
 
def parse_meme_motif(p):

    with open(p, 'r') as file:
        
        lines = file.readlines()
        
    matrix = []
    parse = False

    for line in lines:

        if line.startswith("letter-probability matrix:"):
            parse = True
            continue
        if line.startswith("\n"):
            if parse:
                break
            else:
                continue
        if parse:
            row = list(map(float, line.split()))
            matrix.append(row)
    
    df = pd.DataFrame(matrix, columns=["A", "C", "G", "T"]).T

    return df
    
class MyMotif:

    def __init__(self, df, bg, normalized = False):

        mm =self.convert_df_to_matrix(df)

        if not normalized:
            self.nr = self.normalize_matrix(mm)
        else:
            self.nr = mm

        self.bg = bg

    def relative_entropy(self):

        
        
        

        p = self.nr + 1e-32
        b = self.bg

        M = p.shape[0] 
        N = p.shape[1] 

        H = np.zeros(N) 

        for j in range(N): 

            for i in range(M): 

                H[j] += p[i,j] * np.log2(p[i,j]/b[i])

        return H
        
    def pssm(self):
            
        
        

        p = self.nr + 1e-32
        b = self.bg

        M = p.shape[0] 
        N = p.shape[1]

        pssm = np.zeros((M, N)) 

        for j in range(N): 
            for i in range(M):  
                pssm[i,j] = np.log2(p[i,j]/b[i])
        
        return pssm

    def calculate_score1(self, sequence):

        def letter_to_index(letter):
                
            if letter == 'A': return 0
            if letter == 'C': return 1
            if letter == 'G': return 2
            if letter == 'T': return 3

            return -1

        pssm = self.pssm()
        score = 0

        for j, letter in enumerate(sequence):
            score += pssm[letter_to_index(letter), j]
        
        return score
    
    def pretty_print_pssm(self):

        pssm = self.pssm()

        df = pd.DataFrame(pssm, columns=[i for i in range(pssm.shape[1])], index=['A', 'C', 'G', 'T'])
        df = df.round(3)
        print(df.T)

    @staticmethod
    def convert_df_to_matrix(df):

        matrix = np.zeros(df.shape)
        for i, row in enumerate(df.values): 
            matrix[i] = row

        return matrix

    @staticmethod
    def normalize_matrix(matrix):

        return matrix / matrix.sum(axis=0)
    
    @staticmethod
    def get_background_frequencies(mm):

        ss = mm.sum(axis=1)
        ss = ss / ss.sum()

        return ss
    
    @staticmethod
    def self_test1():

        

        counts = [[3, 7, 0, 2, 1], 
                  [0, 0, 5, 2, 6], 
                  [0, 0, 0, 3, 0], 
                  [4, 0, 2, 0, 0]]
        
        df = pd.DataFrame(counts, columns=[1, 2, 3, 4, 5], index=["A", "C", "G", "T"])

        MYM = MyMotif(df, np.array([0.25, 0.25, 0.25, 0.25]), normalized=False)
        H = MYM.relative_entropy()
        print(H)

    @staticmethod
    def self_test2():

        
        
        
        
        

        counts = [[3, 7, 0, 2, 1], 
                  [0, 0, 5, 2, 6], 
                  [0, 0, 0, 3, 0], 
                  [4, 0, 2, 0, 0]]
        
        df = pd.DataFrame(counts, columns=[1, 2, 3, 4, 5], index=["A", "C", "G", "T"])

        MYM = MyMotif(df, np.array([0.25, 0.25, 0.25, 0.25]), normalized=False)
        PSSM = MYM.pssm()
        print(PSSM)
 
def calculate_probability_dist(d):

    total_sum = sum(d.values())
    pdf = {k: v / total_sum for k, v in d.items()}
    return pdf

def calculate_shannon_entropy(pdf):

    return -sum([p * math.log2(p) for p in pdf.values() if p > 0])

def calculate_kl_divergence(pdf, pdf_uniform):

    return sum([p * math.log2(p / pdf_uniform[k]) for k, p in pdf.items() if p > 0])

def probability_sampler_v0(lst, sub_size):

    total_subsections = len(lst) - sub_size + 1

    rows = []
    for i in range(total_subsections):

        subsection = lst[i:i+sub_size]
        has_one = 1 in subsection
        rows.append([i, subsection, has_one])

    df = pd.DataFrame(rows, columns=['index', 'subsection', 'has_one'])

    probability = df['has_one'].sum() / df.shape[0]

    return probability

def probability_sampler(lst, sub_size):

    lst = np.array(lst)
    window = np.ones(sub_size)
    convolved = np.convolve(lst, window, mode='valid')
    has_one = convolved > 0
    probability = has_one.sum() / len(has_one)
    
    return probability

def count_islands(array):

    array = np.insert(array, 0, 0)
    return np.sum(np.diff(array) > 0)

def generate_def_gene_prob_distributions_with_filter(defense_maps_dict, 
                                                     sizes, 
                                                     snail_data_dict=None, 
                                                     region_filter=False):

    def accessory_region_filter(direction, a, b, c, d, dmap):

        
        if direction == 1: S, E = d, a
        else:              S, E = b, c

        if S >= E:
            return np.concatenate((dmap[S+1:], dmap[0:E]))
        
        return dmap[S+1:E]

    x = []
    y = []

    for sub_size in sizes:

        sample = []
        for acc in defense_maps_dict.keys():

            rlmh_loc, acrb_loc = snail_data_dict[acc].get_boundary_positions()

            direction = rlmh_loc[2]
            a         = rlmh_loc[0]
            b         = rlmh_loc[1]
            c         = acrb_loc[0]
            d         = acrb_loc[1]

            dmap = defense_maps_dict[acc]
            len1 = len(dmap)
            num_islands1 = count_islands(dmap)

            if region_filter:

                dmap = accessory_region_filter(direction, a, b, c, d, dmap)

            maxp = len(dmap) - sub_size + 1
            len2 = len(dmap)
            num_islands2 = count_islands(dmap)

            print(f"-> {acc}: ", len1-len2, num_islands1, num_islands2, direction, rlmh_loc, acrb_loc)
            
            for i in range(100):

                start = np.random.randint(0, maxp)
                end = start + sub_size
                random_section = dmap[start:end]
                sample.append(count_islands(random_section))

        x.append(sub_size)
        y.append(sample)

    return x, y 

class PlasmidAnalyzer:

    class GenBank:

        def __init__(self, fd):

            self.record: SeqRecord = SeqIO.read(fd['genbank'], "gb") 
            self.accession = fd['acc']
            self.source = self.record.annotations['source']
            self.description = self.record.description

            self.length = len(self.record.seq)

            features = list(self.record.features)
            assert features[0].type == "source"
            self.features = features

        def get_plasmid_info(self):

            keys = self.features[0].qualifiers.keys()

            if 'plasmid' in keys:

                is_plasmid = True
                plasmid = self.features[0].qualifiers['plasmid'][0]

            else:

                is_plasmid = False
                plasmid = None

            return is_plasmid, plasmid
        
        def get_row(self):

            row = []

            strain = self.source.replace('Staphylococcus', '_')
            strain = re.sub(r'[Ss]taphylococcus *', '', self.source)

            plasmid_info = self.get_plasmid_info()

            complete = 'complete sequence' in self.description

            row = {"accession": self.accession, 
                "description": self.description, 
                "species": strain, 
                "seq length": self.length, 
                "is plasmid": plasmid_info[0], 
                "plasmid name": plasmid_info[1],
                "complete sequence": complete}

            return row
        
        def get_coding_list(self):

            coding_list = []
            for feature in self.features:

                c1 = feature.type == "CDS"
                c2 = 'translation' in feature.qualifiers.keys()
                if c1 and c2:
                    coding_list.append(feature)

            return coding_list

        def convert2fasta(self, path):

            string = ''

            coding_list = self.get_coding_list()

            for feature in coding_list:

                quals = feature.qualifiers

                loc = LocationTuple(feature.location.start, 
                                    feature.location.end, 
                                    feature.location.strand).encode_location_tuple()

                if 'locus_tag' in quals.keys(): 
                    ltg = quals['locus_tag'][0]
                else: 
                    ltg = 'no_locus_tag'

                if 'protein_id' in quals.keys(): 
                    pid = quals['protein_id'][0]
                else: 
                    pid = 'no_pid'

                if 'product' in quals.keys(): 
                    prd = quals['product'][0]
                else: 
                    prd = 'no_product'

                trl = quals['translation'][0]

                s = '>{} {} {} {}\n{}\n'.format(loc, ltg, pid, prd, trl)
                string += s

            write_file(path, string)

    def __init__(self):

        root = "/Volume/biodata2/"

        self.pg        = root + "/genbank_collection/genbank_plasmids/"
        self.pr        = root + "/sccmec_snail_analysis/snail_analysis_plasmids/"
        
        self.p_results = self.pr + "results/"
        self.p_csv     = self.p_results + "main_results_table.csv"

    def step0_download_genomes(self):

        accessions = self.get_accessions()
        download_genome(self.pr + "genomes/", accessions, 'gbwithparts', 'text')

    def step1_export_protein_fasta(self):

        for acc in self.get_accessions():
        
            fd   = self.fdict(acc)
            GB   = self.GenBank(fd)
            GB.convert2fasta(fd['fasta'])

    def step2_run_defense_finder(self):

        for acc in self.get_accessions():
        
            fd   = self.fdict(acc)
            self.run_defense_finder_PLASMID(fd)

    def step3_generate_df(self):

        rows = []
        for acc in self.get_accessions():

            print(acc)
        
            fd   = self.fdict(acc)

            GB   = self.GenBank(fd)
            r    = GB.get_row()
            cdl  = GB.get_coding_list()

            DFE  = DefenseFinderHitsExtractor(fd['deftp'], fd['rawhmm'])
            df   = DFE.export_defense_finder_hmmer_hits()

            df['gene_name']   = df['gene_name'].str.split("__").str[0]
            defense_gene_set  = set(df['gene_name'].tolist())
            num_raw_def_genes = len(set(df['hit_id'].tolist()))

            g, s    = self.get_genes_and_systems(fd)
            systems = ", ".join([system['subtype'] for system in s])

            

            rows.append({"acc"                : acc, 
                        "total CDS count"    : len(cdl),
                        "raw def gene count" : num_raw_def_genes, 
                        "raw def gene set"   : defense_gene_set,
                        "systems"            : systems,
                        "description"        : r['description'],
                        "species"            : r['species'],
                        "seq nucl. length"   : r['seq length'],
                        "is plasmid"         : r['is plasmid'],
                        "plasmid name"       : r['plasmid name']})
            
        
        
        df = pd.DataFrame(rows)
        df.set_index('acc', inplace=True)
        df.sort_values(by=['raw def gene count'], inplace=True, ascending=False)
        df.to_csv(self.p_csv, index=True)   
        print(df)

    

    def step3_generate_length_histogram(self):

        print("Reading csv...")
        df = pd.read_csv(self.p_csv)

        length_list = df['seq length'].tolist()
        _, ax = plt.subplots(1, 1, figsize=(7, 7))
        ax.hist(length_list, bins=200)
        ax.set_xlabel("Sequence Length")
        ax.set_ylabel("Frequency")
        plt.savefig(self.p_results + "length_hist.png")

    def step4_generate_species_histogram(self):

        print("Reading csv...")

        df = pd.read_csv(self.p_csv)
        species = df['species'].tolist()

        _, ax = plt.subplots(1, 1, figsize=(7, 7))
        plt.hist(species, bins=100)
        plt.xlabel("Species")
        plt.ylabel("Frequency")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.xticks(fontsize=6)
        plt.savefig(self.p_results + "species_hist.png")

    def step5_plot_length_vs_def_gene_count(self):

        df = pd.read_csv(self.p_csv)

        _, ax = plt.subplots(1, 1, figsize=(7, 7))
        ax.scatter(df['seq nucl. length'], df['raw def gene count'], s=1)
        ax.set_xlabel("Sequence Length")
        ax.set_ylabel("Defense Gene Count")
        plt.savefig(self.p_results + "length_vs_def_gene_count.png")

    
    def step6_plot_stacked_hist(self):

        def thousands(x, pos):
            'The two args are the value and tick position'
            return '%1.0fK' % (x * 1e-3)
        
        formatter = FuncFormatter(thousands)

        df = pd.read_csv(self.p_csv)
        

        df_x = df[df['raw def gene count'] > 0]
        df_n = df[df['raw def gene count'] == 0]

        len_x = df_x['seq nucl. length']
        len_n = df_n['seq nucl. length']

        
        

        b = np.arange(0, 100000, 1000)

        hist_x, bb = np.histogram(len_x, bins=b)
        hist_n, bb = np.histogram(len_n, bins=b)

        _, ax = plt.subplots(1, 1, figsize=(7, 7))

        ax.bar(bb[:-1], 
            hist_x, 
            width=600, 
            label='defense genes', 
            color='red')

        ax.bar(bb[:-1], 
            hist_n, 
            width=600, 
            bottom=hist_x, 
            label='no defense genes', 
            color='blue')
        
        ax.set_xlim(0, 85000)
        ax.set_ylim(1e-1, 220)
        
        

        ax.xaxis.set_major_locator(plt.MultipleLocator(5000))
        ax.xaxis.set_major_formatter(formatter)

        ax.set_xlabel("Sequence Length")
        ax.set_ylabel("Number of Genomes")

        ax.grid(which='major', linestyle='--', color='grey', alpha=0.5)

        plt.xticks(fontsize=9)
        plt.xticks(rotation=90) 
        plt.tight_layout()
        plt.savefig(self.p_results + "length_hist_stacked.png")    

    def step7_plot_defense_gene_diversity(self):

        print("PLOT: Diversity of defense genes")

        df = pd.read_csv(self.p_csv)
        df['raw def gene set'] = df['raw def gene set'].apply(ast.literal_eval)

        all_def_genes = set()
        for s in df['raw def gene set']:
            all_def_genes = all_def_genes.union(s)

        gene_counts = {
            gene: sum(df['raw def gene set'].apply(lambda x: gene in x)) 
            for gene in all_def_genes}

        df_counts = pd.DataFrame.from_dict(gene_counts, orient='index', columns=['count'])
        df_counts.sort_values(by=['count'], inplace=True, ascending=False)
        print(df_counts)   

        

        counts_dict = {}

        for index, row in df_counts.iterrows():

            key = index.split('_')[0]

            if key in counts_dict:

                counts_dict[key] += row['count']

            else:

                counts_dict[key] = row['count']

        print()

        df_counts = pd.DataFrame.from_dict(counts_dict, orient='index', columns=['count'])
        df_counts.sort_values(by=['count'], inplace=True, ascending=False)
        print(df_counts) 

        

        _, ax = plt.subplots(1, 1, figsize=(7, 4))
        ax.bar(df_counts.index, df_counts['count'], width=0.5, align='center', color='blue')
        ax.set_xlabel("Defense Gene")
        ax.set_ylabel("Number of Genomes")

        ax.grid(which='major', linestyle='--', color='grey', alpha=0.5)

        ax.set_axisbelow(True)

        plt.xticks(fontsize=9)
        plt.xticks(rotation=90)            
        plt.tight_layout()
        plt.savefig(self.p_results + "defense_gene_diversity.png")  

    
        
    def get_accessions(self):

        p_list = self.pr + "accessions/accessions.txt"
        p_remv = self.pr + "accessions/remove.txt"

        with open(p_list) as f:

            accessions = f.readlines()

        accessions = [x.strip() for x in accessions]

        with open(p_remv) as f:

            remove = f.readlines()

        remove = [x.strip() for x in remove if not x.startswith('#')]

        accessions = [x for x in accessions if x not in remove]

        return accessions   

    def fdict(self, acc):

        return  {"acc"     : acc,
                 "genbank" : self.pg +              f"{acc}.gb",
                 "fasta"   : self.pr + "output/"  + f"{acc}.fa",
                 "defin"   : self.pr + "output/"  + f"{acc}.pkl",
                 "rawhmm"  : self.pr + "output/"  + f"{acc}.df",
                 "idx"     : self.pr + "output/"  + f"{acc}.fa.idx",
                 "deftp"   : self.pr + "defense_finder_results_allpr/" + f"{acc}/"}

    @staticmethod
    def run_defense_finder_PLASMID(fd, dbtype='ordered_replicon' , workers=32):

        path_fasta = fd['fasta']
        path_deftp = fd['deftp']
        
        print("\nRunning defense finder on:")
        print(f"\t> {path_fasta}")

        if not os.path.exists(path_deftp): 

            os.makedirs(path_deftp)

            with open(path_fasta) as f: 
                defense_finder.run(f, dbtype, workers, path_deftp)

        else:

            print("\t> defense finder already ran...")

    @staticmethod
    def get_genes_and_systems(fd):
            
        #----------------------------------------------------------

        def build_systems(genes):

            system_groups = [list(it) for k, it in groupby(genes, lambda val: val['sys_id'])]
            systems = []

            for system_group in system_groups:

                item = {}
                first_item = system_group[0]
                last_item = system_group[-1]
                item['sys_id'] = first_item['sys_id']
                item['sys_beg'] = first_item['hit_id']
                item['sys_end'] = last_item['hit_id']
                item['type'] = first_item['type']
                item['subtype'] = first_item['subtype']
                item['protein_in_syst'] = reduce(lambda acc, s: acc + ',' + s['hit_id'], system_group, '')[1:]
                item['genes_count'] = len(system_group)
                item['name_of_profiles_in_sys'] = reduce(lambda acc, s: acc + ',' + s['gene_name'], system_group, '')[1:]
                systems.append(item)

            return systems

        

        path_defin = fd['defin']
        path_deftp = fd['deftp']

        if not os.path.isfile(path_defin):

            defense_finder_genes = defense_finder_posttreat.best_solution.get(path_deftp)
            systems = build_systems(defense_finder_genes)

            genes_dict = dict()
            for item in defense_finder_genes:

                genes_dict[item['hit_id']] = item

            data = {'genes': genes_dict, 'systems': systems}

            with open(path_defin, 'wb') as f: pickle.dump(data, f)

        else:

            with open(path_defin, 'rb') as f: data = pickle.load(f)

        return data['genes'], data['systems']    

class DefMapDist:

    def __init__(self, dvar):

        self.dvar   = dvar
        self.path_d = dvar['root'] + "defense_finder_results_allpr/"
        self.path_o = dvar['root'] + "defense_finder_raw_genes_df/"
        self.path_m = dvar['root'] + "defense_maps_dict.pkl"
        self.path_g = dvar['genomes']

    def generate_defense_maps_dict(self):

        def set_indices_to_one(array, a, b):
            array[a:b+1] = 1
            return array
        
        defense_maps_dict = {}
        for path in os.listdir(self.path_d):

            acc = self.parse_filename(path)

            rec = SeqIO.read(self.path_g + f"{acc}.gb", "gb")
            N   = len(rec.seq)  

            DFE = DefenseFinderHitsExtractor(tmp_dir  = self.path_d + path, 
                                             path_out = self.path_o + f"{acc}.df")

            df = DFE.export_defense_finder_hmmer_hits()
            
            f = lambda x: LocationTuple.decode(x, True)
            df['decoded'] = df['hit_id'].map(f)

            defense_map = np.zeros(N, dtype=int)

            for i, row in df.iterrows():

                a, b, strand = row['decoded']
                defense_map = set_indices_to_one(defense_map, a, b)

            defense_maps_dict[acc] = defense_map
            
        with open(self.path_m, 'wb') as f: pickle.dump(defense_maps_dict, f)

    def generate_defense_finder_raw_genes_df(self):

        all_def_genes = []
        for path in os.listdir(self.path_d):

            
            acc = self.parse_filename(path)

            DFE = DefenseFinderHitsExtractor(tmp_dir  = self.path_d + path, 
                                             path_out = self.path_o + f"{acc}.df")
            
            df = DFE.export_defense_finder_hmmer_hits()
            df['gene_name']   = df['gene_name'].str.split("_").str[0]
            incrementl =  list(set(df['gene_name'].tolist()))

            
            
            all_def_genes.extend(incrementl)

        all_def_genes = Counter(all_def_genes)

        return all_def_genes

    
    
    @staticmethod
    def parse_filename(txt):

        return txt.split(']_')[0][1:]
    
    
    @staticmethod
    def find_overlapping_ranges(ranges):
        
        sorted_ranges = sorted(ranges, key=lambda x: x[0])

        overlapping_ranges = []
        for i in range(1, len(sorted_ranges)):
            
            if sorted_ranges[i][0] <= sorted_ranges[i-1][1]:
                overlapping_ranges.append((sorted_ranges[i-1], sorted_ranges[i]))

        return overlapping_ranges    
        
class DefenseFinderHitsExtractor:

    def __init__(self, tmp_dir, path_out):

        self.tmp_dir = tmp_dir
        self.path_out = path_out

        self.hmmer_keys = ["hit_id", 
                           "replicon_name", 
                           "position_hit", "hit_sequence_length", "gene_name", "i_eval", 
                           "score", "profile_coverage", "sequence_coverage", 
                           "begin", 
                           "end"]

    def get_hit_sort_attr(self, hit):
        return hit['hit_id']

    def remove_duplicates(self, hmmer_hits):
        return list({v['hit_id']:v for v in hmmer_hits}.values())

    def parse_hmmer_results_file(self, path):

        
        
        
        

        with open(path) as tsv_file:
            tsv = csv.reader(tsv_file, delimiter='\t')
            data = [row for row in tsv if not row[0].startswith('#')]
        return [dict(zip(self.hmmer_keys, l)) for l in data if l]

    def get_hmmer_paths(self):
        family_dirs = os.listdir(self.tmp_dir)
        family_dirs = [f for f in family_dirs if not f.startswith('.')]
        files = []
        for family_dir in family_dirs:
            hmmer_results_dir = os.path.join(self.tmp_dir, family_dir, 'hmmer_results')
            with os.scandir(hmmer_results_dir) as it:
                for entry in it:

                    c1 = entry.name.endswith('extract')
                    c2 = entry.is_file()
                    c3 = not entry.name.startswith('.')
                    if c1 and c2 and c3:
                        files.append(entry.path)

        return files

    def hmmer_to_list(self, hmmer_hits):
        out = [self.hmmer_keys]
        for s in hmmer_hits:
            out.append([s[key] for key in self.hmmer_keys])
        return out

    def export_defense_finder_hmmer_hits(self):

        if os.path.isfile(self.path_out):

            

            with open(self.path_out, 'rb') as f: df = pickle.load(f)
            return df

        paths = self.get_hmmer_paths()
        hmmer_hits = []
        for path in paths:
            d = self.parse_hmmer_results_file(path)
            hmmer_hits.extend(self.remove_duplicates(d))
        sorted_hmmer_hits = sorted(hmmer_hits, key=self.get_hit_sort_attr)
        hmmer_hits_list = self.hmmer_to_list(sorted_hmmer_hits)

        columns = hmmer_hits_list.pop(0)
        df = pd.DataFrame(hmmer_hits_list, columns=columns)
        with open(self.path_out, 'wb') as f: pickle.dump(df, f)
        return df
     
message_enabled = False

file_types = {"_ngbrs_": ".fa",
              "_allpr_": ".fa",
              "_pfams_": ".txt",
              "_defin_": ".dic",
              "_alsys_": ".dic",
              "_blast_": ".txt",
              "_cored_": ".dic",
              "_coral_": ".dic",
              "_dtfrm_": ".df",
              "_dfcsv_": ".csv"}

object_types = {"_ngbrs_": xObjects.NgbrsObj,
                "_allpr_": xObjects.AllPrObj,
                "_pfams_": xObjects.PfamsObj,
                "_defin_": xObjects.DeFinObj,
                "_alsys_": xObjects.AlSysObj,
                "_blast_": xObjects.BlastObj,
                "_cored_": xObjects.CoredObj,
                "_coral_": xObjects.CoralObj,
                "_dtfrm_": xObjects.DtfrmObj,
                "_dfcsv_": xObjects.DfCsvObj}

if __name__ == "__main__":

    u = 5

    if u == -4:

        path = pr + "snail_analysis_epidermidis/neighborhood_files_ws500/[NC_002976.3]_ngbrs_[2584175m2584655][n].fa"

        for i, seq_record in enumerate(SeqIO.parse(path, "fasta")):

            desc = seq_record.description.replace(seq_record.id, '').strip()

            mm = re.search(r"recombinase\s+family\s+protein", desc, re.IGNORECASE)
            if mm is not None: print(desc, mm)

    if u == -1:

        import functools
        import datetime
        import inspect

        enable = True

        class LogFile:
            def __init__(self):
                self.columns = ['Time', 'CalledBy', 'Function', 'Arguments', 'Output']
                self.df = pd.DataFrame(columns=self.columns)
                self.empty = 1

            def add(self, called_by, method_name, arguments, output):
                self.empty = 0
                tm = datetime.datetime.now().strftime("%H:%M:%S.%f")
                rows = [[tm, called_by, method_name, arguments, output], ]
                df2 = pd.DataFrame(rows, columns=self.columns)
                self.df = self.df.append(df2, ignore_index=True)

            def __str__(self):
                if self.empty: return ''
                log_text = 'LOG FILE: \n\n' + self.df.to_string()
                return log_text

        log_file = LogFile()

        def logger(f):
            @functools.wraps(f)
            def with_logger(*args, **kwargs):
                output = f(*args, **kwargs)
                global log_file
                global enable
                if enable:
                    curframe = inspect.currentframe()
                    calframe = inspect.getouterframes(curframe, 2)
                    argument = ''
                    if f.__name__ == 'self_test': argument = args[1]
                    line = (calframe[1][3], f.__name__, argument, '')
                    log_file.add(*line)
                return output

            return with_logger

        @logger
        def myfunct1(x):

            return x ** 2

        @logger
        def myfunct2(x):

            y = 2 * x
            return myfunct1(y)

        a = myfunct2(5)
        print(log_file)

    if u == 1:

        params = {'path_file': pr + "aureus_snail_analysis_TEST/",
                    'path_genomes': pr + "genbank_collection/genbank_aureus_collection/",
                    'accession': 'NZ_CP033112.1',
                    'location': "[33658:34138](+)",
                    'window_size': 20,
                    'displacement': 10}

        
        
        
        
        
        

        

    if u == 2:

        j = 2

        path_ngbrs = pr + "aureus_snail_analysis_TEST/neighborhood_files_ws500/[NC_002745.2]_ngbrs_[33691p34171][p].fa"
        path_pfams = pr + "aureus_snail_analysis_TEST/neighborhood_files_ws500/[NC_002745.2]_pfams_[33691p34171][p].txt"
        path_blast = pr + "aureus_snail_analysis_TEST/neighborhood_files_ws500/[NC_002745.2]_blast_[33691p34171][p].txt"
        path_cored = pr + "aureus_snail_analysis_TEST/neighborhood_files_ws500/[NC_002745.2]_cored_[33691p34171][p].dic"

        if j == 0:

            Obj = xObjects.BlastObj(path_blast)

            st = time.perf_counter()
            df = Obj.parse2df()
            ed = time.perf_counter()
            print('Finished in {} seconds:'.format(ed - st))

            for ind in df.index[0:20]:
                print(df.loc[ind]['sseqid'])

        if j == 1:

            Obj = xObjects.BlastObj(path_blast)

            line_gen = Obj.parse2gen()

            for n in range(15):
                print(next(line_gen))

        if j == 2:
            Obj = xObjects.PfamsObj(path_pfams)
            df0, df1 = Obj.parse()

            print(df0)
            print()
            print(df1)

    if u == 3:

        pr = pr + "snail_analysis_epidermidis/neighborhood_files_ws500"
        pg = pg + "genbank_epidermidis_collection/"

        FM = nbFileManager(pr, pg)

        j = 2

        if j == 0:

            print(FM)

            for accession, AccObj in FM.AccessionObjDict.items():

                print('\n' + 100 * '-' + '\n')
                print(accession)
                

                for idx, obj in AccObj.items():
                    print('\t', idx, )
                    print('\t', type(obj), obj.parent.accession, obj.identifier)

                
                

            
            
            
            
            
            
                    
        if j == 1:

            for k, v in FM.GenomeObjDict.items():

                print(k)
                print("\t", v, "\n")

        if j == 2:

            set_a = set(FM.GenomeObjDict.keys())
            set_b = set(FM.AccessionObjDict.keys())

            print(set_a == set_b)

    if u == 5:

        j = 9.2
        dvars = load_dvars(pr)
        def_mode = "subtype"

        
        if j == 0:

            print("debug: File Manager Debug Print")
            SA = SnailAnalyzerFull(dvars['epidermidis_test'])
            SA.file_manager.debug_print()

        
        if j == 1:

            SA = SnailAnalyzerFull(dvars['test'])

            def report_on_cored(obj):

                genes = obj.parse()

                accessions_bucket = set()
                rows = []

                for pos, item in genes.items():
                    accessions_bucket = accessions_bucket.union(item[0])
                    rows.append("\t{}:{} - {}".format(pos, len(item[0]), item[1]))

                print()
                print('\tTotal number of all unique accessions in Cored: {}'.format(len(accessions_bucket)))
                print()

                for row in rows[200:220]: print(row)

            def report_on_cored2(obj):

                genes = obj.parse()

                accessions_bucket = set()
                rows = []

                for pos, item in genes.items():

                    if item[1] == 501:
                        accessions_bucket = accessions_bucket.union(item[0])

                print()
                print('\tTotal number of unique accessions in 500 genes: {}'.format(len(accessions_bucket)))
                print()

            accessions = SA.all_existing_ngbrs_file_accessions

            for accession in accessions:

                print("\nAccession: {}".format(accession))

                obj = SA.file_manager.get_xObject(accession, 'cored')

                if obj.exists:

                    report_on_cored2(obj)

                else:

                    print("\n\t\tcored doesn't exist")

        
        if j == 2:

            path = pr + "snail_analysis_epidermidis/blast_groups/group_1226.txt"

            with open(path, 'rb') as f: data = pickle.load(f)

            aa = [ll.split('|')[0] for ll in data]
            aa = set(aa)

            print(len(data), len(aa))

        
        if j == 3:

            def compare(acc1, acc2, SA):

                genes1 = SA.file_manager.get_xObject(acc1, 'cored').parse()
                genes2 = SA.file_manager.get_xObject(acc2, 'cored').parse()

                m = len(genes1)
                n = len(genes2)

                c = np.zeros((m, n))

                x, y = np.meshgrid(np.linspace(0, m - 1, m), np.linspace(0, n - 1, n))

                for i, (pos1, item1) in enumerate(genes1.items()):

                    for j, (pos2, item2) in enumerate(genes2.items()):

                        if item1[1] == 500 and item2[1] == 500:

                            c[i, j] = len(item1[0].intersection(item2[0]))

                        else:

                            c[i, j] = np.nan

                return x, y, c

            
            def analyzer_coreds():

                SA = SnailAnalyzerFull(dvars['test'])

                accessions = SA.all_existing_ngbrs_file_accessions

                acc1, acc2 = accessions[0], accessions[1]

                x, y, c = compare(acc1, acc2, SA)

                fig = plt.figure(dpi=600)
                ax = fig.add_subplot(111)

                pc = ax.pcolormesh(x, y, c, edgecolor='none', cmap='Blues', shading='auto')
                pc.set_clim(vmin=0, vmax=500)
                cbar = fig.colorbar(pc, pad=0.2)
                cbar.ax.set_ylabel('Intersection Size', rotation=90, fontsize=10)

                plt.show()

        
        if j == 4:

            print('File Renaming')
            path = pr + "epidermidis_snail_analysis/neighborhood_files_ws500/"

            ff = []
            for filename in os.listdir(path):

                if "]_cored_[" in filename:
                    new_filename = filename.replace("]_coredict_[", "]_cored_[")
                    print(filename, ' --> ', new_filename)
                    ff.append(filename)
                    

            print(len(ff))

        
        if j == 5:

            
            
            

            genome_dict = {}  
            
            

            print(genome_dict)

            for acc, AccObj in SA.file_manager.AccessionObjDict.items():

                features = AccObj['ngbrs'].parse()

                print(f'\nAccession: {acc}')

                for Seq in features:

                    

                    if Seq.id not in genome_dict[acc].keys(): print('\t Not Found:', acc, Seq.id)

        
        if j == 6:

            dynamic_acc_pos_list = list(SeqIO.parse(SA.path_combined_fasta, "fasta"))
            query = "AFKLKPDCHCTSKYLNNLIEQDHRHIKVRKTRYQSINTAKNTLKGIECIYALYKKNRRSLQICGFSPCHEISIMLAS"  
            hits = levenshtein_searcher(query, dynamic_acc_pos_list)
            for h in hits: print(h)

        
        if j == 7:

            
            dvar = dvars['aureus']
            SA = SnailAnalyzerFull(dvar)

            

            

            

            

            

            

            

            

            

            

            

            

            

        if j == 9.0:

            

            for s in ['epidermidis', 'aureus']:

                for i in [0, 2]:

                    dvar = dvars[s]
                    SA = SnailAnalyzerFull(dvar)
                    SA.get_diversity_dicts(200, mode = def_mode, i=i)

        if j == 9.1:

            
            
            DO_ENRICHMENT_STATS_per_genome(dvars, mode=def_mode)
            PLOT_ENRICHMENT_TESTING_pie(dvars, mode=def_mode, plot_type='bar', i=0)
            PLOT_ENRICHMENT_TESTING_scatter(dvars, mode = def_mode)
            PLOT_ENRICHMENT_TESTING(dvars, mode=def_mode, plot_type='line', i=0)

        if j == 9.2:

            print("Hobolob")

            DO_ENRICHMENT_STATS_collective(dvars, mode=def_mode)

        if j == 9.3:

            
         
            PLOT_DEFENSE_DIVERSITY(dvars, mode = def_mode, pos_anchor=200, i = 0)

        if j == 9.4:

            SHANNON_DEFENSE_DIVERSITY(dvars, mode=def_mode, i = 0)

        if j == 9.5:

            PLOT_DEFENSE_CONCENTRATION(dvars, mode = def_mode, i = 2, plot_mode = 0)

    if u == 6:

        j = 1

        
        if j == 0:

            j = 2

            pinp = pr + "test/input/[NC_002976.3]_ngbrs_[2584175m2584655][n].fa"
            pout = pr + "test/output/"

            if j == 0:

                st = time.perf_counter()
                
                ed = time.perf_counter()
                print(f"finished in {ed-st} seconds... ")

            if j == 1:

                pp = pr + "test/iop_systems.pickle"

                with open(pp, 'rb') as f: systems = pickle.load(f)

                for system in systems:

                    print('----')

                    for k,v in system.items():
                        print('\t', k, ": ", v)

            if j == 2:

                pp = pr + "test/pp_iop_genes.pickle"

                with open(pp, 'rb') as f: systems = pickle.load(f)

                print(type(systems))

                for system in systems:

                    print('----')

                    for k,v in system.items():
                        print('\t', k, ": ", v)

        
        if j == 1:

            p = "snail_analysis_epidermidis/neighborhood_files_ws500"
            path_allsys = pr + f"{p}/[NC_002976.3]_alsys_[2584175m2584655][n].dic"
            path_dtfrm  = pr + f"{p}/[NC_002976.3]_dtfrm_[2584175m2584655][n].df"
            path_genome = pg + "genbank_epidermidis_collection/NC_002976.3.gb"

            obj0 = xObjects.AlSysObj(path_allsys)
            obj1 = xObjects.DtfrmObj(path_dtfrm)
            AccObj = AccessionObj('NC_002976.3')

            AccObj.add_xObject(obj0)
            AccObj.add_xObject(obj1)
            

            k = 1

            if k == 0:

                plot_params = dict()
                plot_params['th'] = 85
                plot_params['n']  = 89

                p = DefenseFeature(AccObj, path_genome, plot_params, 200)
                print(p.pos_anchor)
                print(p.pos_acrb)
                print()
                print(p.path_genome)
                print()
                rp, nuc_len, gnm_len = p.get_vector(mode = 'type', i=0)

                for pp in rp:
                    print(pp)

            if k == 1:

                genes, systems = AccObj['alsys'].parse()

                for pos, data in genes.items():

                    print(pos)
                    print()
                    for key in data:
                        print('\t', key, ': ', data[key])

    if u == 12:

        j = 3

        dvars = load_dvars(pr)
        dvar = dvars['epidermidis']

        
        if j == -1:

            k = 1

            pr = pr + "snail_analysis_epidermidis/neighborhood_files_ws500/"
            pg = pg + "genbank_epidermidis_collection/"

            FM = nbFileManager(pr, pg)
            print(FM)

            if k == 0:

                print("Fixing CSL and experiments...")

                for accession, AccObj in FM.AccessionObjDict.items():

                    print('\n' + 100 * '-' + '\n')
                    print(accession)
                    print(AccObj)
                    
                    
                    

            if k == 1:

                print("Fixing CSL and experiments...")

                pp = {'n': 89, 'm': 355, 'th': 85, 'c_depth': 100, 
                    'font_size_1': 8, 'font_size_2': 8, 'res_dpi': 200, 
                    'fw': 16, 'fh': 12, 'annotate_nums': True, 
                    'annotate_defs': False, 'circle': True}
                
                AccObj = FM.AccessionObjDict['NC_002976.3']

                print(AccObj)

                sf = SnailFeature(AccObj, 
                                    200-pp['c_depth'], 
                                    pp['n'], 
                                    pp['c_depth'], 
                                    pp['m'], 
                                    pp['th'])
        
                rlmh_loc,acrb_loc = sf.get_boundary_positions()

                print(AccObj.accession, rlmh_loc, acrb_loc)

            if k == 2:

                AccObj = FM.AccessionObjDict['NC_002976.3']
                print(AccObj)

                with open(AccObj['dtfrm'].path_full, 'rb') as f: df = pickle.load(f)

                print(df)

            
            if k == 3:

                rlmh = Seq("ATCATC")
                mmmm = str(rlmh.reverse_complement())

                fwd = [rlmh, Seq(f"AACGTG{rlmh}CCGT"), (6, 12,  1)]
                rev = [mmmm, Seq(f"AACGTG{mmmm}CCGT"), (6, 12, -1)]

                s = rev

                print(s[0])
                print(s[1])

                myseq = SeqX(s[1])
                myseq.realign_and_reorient(s[2])
                print(myseq.seq)

            
            if k == 4:

                for accession, AccObj in FM.AccessionObjDict.items():

                    rlmh_raw = AccObj['defin'].location_tuple.encode_location_tuple()
                    rlmh_loc = LocationTuple.decode(rlmh_raw)
                    path_genome = pg + f"genbank_epidermidis_collection/{AccObj.accession}.gb"

                    
                    
                    
                    
                    
                    

                    myseqx = SeqX(SeqIO.read(path_genome, "gb").seq)
                    myseqx.realign_and_reorient(rlmh_loc)
                    print()
                    print(myseqx.seq[0:210].translate(), f" | {accession}")

                    

        
        if j == 1:

            k = -5

            if k == 0:

                p = "/Volume/biodata2/sccmec_snail_analysis" \
                    "/snail_analysis_epidermidis/cut_site_analysis/meme" \
                    "/iteration_1/motif/meme2meme.txt"

                df = parse_meme_motif(p)
                b0 = [0.33994502393617976, 0.16014331538707866, 0.1605355701496629, 0.33853455549359557]
                b1 = [0.25, 0.25, 0.25, 0.25]
                b2 = [0.49, 0.01, 0.01, 0.49]
                b3 = [0.01, 0.49, 0.49, 0.01]

                M0 = MyMotif(df, np.array(b0), normalized=True)
                M1 = MyMotif(df, np.array(b1), normalized=True)

                a0 = sum(M0.relative_entropy())
                a1 = sum(M1.relative_entropy())

                p0 = M0.pssm()
                p1 = M1.pssm()

                print(M0.nr)
                print()
                
                
                
                
                
                print(p0)
                print()
                print(M0.pretty_print_pssm())
                
                
                

            if k == 1:

                pr = 0.25 * np.ones((4, 18))
                print(pr)
                cols = [i for i in range(pr.shape[1])]

                df = pd.DataFrame(pr, columns=cols, index=["A", "C", "G", "T"])

                b = [0.25, 0.25, 0.25, 0.25]
                b = [0.33994502393617976, 0.16014331538707866, 0.1605355701496629, 0.33853455549359557]

                M0 = MyMotif(df, np.array(b), normalized=True)
                print(M0.pssm())

        
        if j == 2:

            CSA = CutSiteAnalyzer(dvars['epidermidis'], mode='custom')

            k = 4

            if k == -1:

                print('Building multi-genome background model')
                dist_data = CSA.generate_bg_Dist()
                print(dist_data['pssm'])

            if k == 0:

                CSA.analyze()

            if k == 1:

                CSA.plot()

            if k == 2:

                AccObj    = CSA.file_manager.AccessionObjDict["NC_002976.3"]
                dist_data = CSA.generate_bg_Dist(pseudocounts=0.1)
                row, cslic, csloc = CSA.analyze_single(AccObj, dist_data)

                
                
                
                
                
                
                
                

            
            if k == 5:

                path = CSA.path_root + 'cut_site_df'
                with open(path, 'rb') as f: df1 = pickle.load(f)

                path = CSA.path_root + 'cut_site_df_old_ro1'
                with open(path, 'rb') as f: df2 = pickle.load(f)

                def compared_dfs_v0(df1, df2):

                    df1.set_index('Accession', inplace=True)
                    df2.set_index('Accession', inplace=True)

                    rows = []
                    for i in df1.index:

                        row1 = df1.loc[i].tolist()
                        row2 = df2.loc[i].tolist()
                        row1.append(1)
                        row2.append(2)
                        row1.append(i)
                        row2.append(i)

                        rows.append(row1)
                        rows.append(row2)

                    cols = df1.columns.tolist()
                    cols.append('DF #')
                    cols.append('Accession')

                    df = pd.DataFrame(rows, columns = cols)
                    df.set_index('Accession', inplace=True)
                    print(df)

                def compared_dfs_v1(df1, df2):

                    df1.set_index('Accession', inplace=True)
                    df2.set_index('Accession', inplace=True)

                    rows = []
                    for i in df1.index:
                        row1 = np.array(df1.loc[i].tolist())
                        row2 = np.array(df2.loc[i].tolist())
                        row = row1 - row2
                        row = np.append(row, i)
                        rows.append(row)

                    cols = df1.columns.tolist()
                    
                    cols.append('Accession')

                    df = pd.DataFrame(rows, columns=cols)
                    df.set_index('Accession', inplace=True)
                    print(df)

                compared_dfs_v0(df1, df2)

        
        if j == 3:

            CSA = CutSiteAnalyzer_NEW(dvars['epidermidis'], "no_pseudo_iteration_1")

            k = 2

            if k == 0:
            
                CSA.step0_create_FIMO_motif()

            if k == 1:
                    
                CSA.step_1_analyze()

            if k == 2:

                CSA.step_2_plot()

            if k == 3:

                CSA.create_next_iteration()

    if u == 13:

        u = 1

        if u == 1:

            path_tmp = "/Volume/biodata2/sccmec_snail_analysis" \
                       "/snail_analysis_plasmids/defin_original_output/NZ_CP030247.1/"
        

            path_out = "/Volume/biodata2/sccmec_snail_analysis" \
                        "/snail_analysis_plasmids/NZ_CP030247.1_raw_def_genes.df"

            DFE = DefenseFinderHitsExtractor(path_tmp, path_out)
            df = DFE.export_defense_finder_hmmer_hits()
            print(df)
    
    if u == 15:

        path_m      = "/Volume/biodata2/sccmec_snail_analysis" \
                      "/snail_analysis_epidermidis/defense_maps_dict.pkl"

        path_snaild = "/Volume/biodata2" \
                      "/sccmec_snail_analysis/snail_analysis_epidermidis/snail_data"

        with open(path_m, 'rb') as      f: defense_maps_dict = pickle.load(f)

        with open(path_snaild, 'rb') as f: snail_data = pickle.load(f)

        print('Data Loaded...')

        sf: SnailFeature
        snail_data_dict = {sf.accession: sf for sf in snail_data}

        j = 3

        if j == -1:

            
            
            

            
            

            print(len(defense_maps_dict))
            probabilities = []

            for acc, defense_map in defense_maps_dict.items():

                num_defenses = count_islands(defense_map)
                p = probability_sampler(defense_map, 1000)

                if p < 0.95:

                    print(f"{acc}: -> num_def: {num_defenses}, p: {p}")
                    probabilities.append(p)

                else:

                    print(f"{acc}: -> IGNORED -> P = {p}")

            
                
            probabilities = np.array(probabilities)
            prob = probabilities.mean()
            print(prob)
            print()
            print(probabilities)

        if j == 0:

            sizes = np.arange(5000, 60000, 5000)
            f    = generate_def_gene_prob_distributions_with_filter
            x, y = f(defense_maps_dict, sizes)

            _, ax = plt.subplots()

            k = 3

            if k == 0:

                medianprops = dict(linestyle='-', linewidth=15, color='black') 
                ax.boxplot(y, positions=x, widths=1000, showfliers=False, medianprops=medianprops)
                plt.show()

            if k == 1:

                means = []
                var = []

                for sample in y:
                    means.append(np.mean(sample))
                    var.append(np.var(sample))  

                ax.scatter(means, var)
                
                lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
                        np.max([ax.get_xlim(), ax.get_ylim()])]

                ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)

                plt.show()

            if k == 2:

                means = []
                for sample in y: means.append(np.mean(sample))
                ax.scatter(x, means)
                plt.show()

        if j == 3:

            sizes = np.arange(5000, 60000, 5000)
            sizes = np.array([20000, 40000])

            f    = generate_def_gene_prob_distributions_with_filter
            x, y = f(defense_maps_dict, sizes, snail_data_dict, region_filter=True)

            _, ax = plt.subplots(nrows=len(x), ncols=1, sharex=True, figsize=(8, 18))

            vals = np.arange(0, 8, 1)

            for i, sample in enumerate(y):

                hist, vals = np.histogram(sample, bins=vals, density=True)
                ax[i].bar(vals[:-1], hist, width=1, align='edge')
                ax[i].set_ylabel(f"size: {x[i]}")

            
            plt.savefig("dist.png")

        if j == 4:

            print(defense_maps_dict.keys())

            print()

            sf = snail_data_dict["NC_002976.3"]
            bp = sf.get_boundary_positions()
            print(type(bp[0]))

