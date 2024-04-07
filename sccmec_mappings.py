import sys
sys.path.append("/Volume/biocode2/")

import os.path
import time
import shutil
import re
import random
import subprocess
from math import floor, ceil
from statistics import median, stdev, mean
import pandas as pd
import numpy as np
import scipy.stats
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from main_classes.anticrispr import MapperAN1
from main_classes.sjsMapper import IGV_Mapper, dvarGen, BreSeq
from main_classes.main_module import ensure_directory

import scipy
import matplotlib
import matplotlib.pyplot as plt 

class DataMat:

    def __init__(self, d):

        name = d["Name"]

        self.datasets      = d["DataSets"]

        self.path_root     = "/Volume/biodata2/sccmec_mappings/"

        self.path_mappings = self.path_root + f"mappings_{name}/"

        self.path_queries  = self.path_root + f"QUERIES_{name}/"

        self.path_results  = self.path_root + f"results_{name}/"

        self.path_df       = self.path_results + "counts.csv"

        self.path_ar       = self.path_results + "alignment_rates.csv"

        self.path_plots    = self.path_results  
           
    def gd(self, print_dvar=False, forceDvarGen=False):

        if forceDvarGen: print("Forcing dvarGen...")

        for dataset in self.datasets:

            dvars = dvarGen(self.path_mappings + f"{dataset}/", 
                            self.path_queries, 
                            force=forceDvarGen)
                
            for name in dvars:

                if print_dvar: print(dvars[name])
                
                yield name, dvars[name], MapperAN1(dvars[name]), dataset

class Cassette:

    def __init__(self, name, train, x, y, coordinate_system = 'snapgene'):

        if coordinate_system == "snapgene":

            x = x-1

        else:

            pass

        self.name, self.train, self.x, self.y = name, train, x, y

    def get_cassette(self):

        return self.train[self.x:self.y]

    def get_junction(self, N, m=0, n=0):

        x, y = self.x + m, self.y + n
        return self.train[y - N:y] + self.train[x:x + N]

    def get_scar(self, N, m=0, n=0):

        x, y = self.x + m, self.y + n
        return self.train[x - N:x] + self.train[y:y + N]

    def write2file(self, path, N, m=0, n=0):

        for segment_type in ["scar", "junction"]:
            sequence = getattr(self, f"get_{segment_type}")(N, m, n)
            record = SeqRecord(Seq(sequence), id=f"{self.name}_{segment_type}", description="")
            SeqIO.write(record, f"{path}qry_{self.name}_{segment_type}.fasta", "fasta")

    def readFromFile(self, path):

        record_dict = {}

        r_scar = SeqIO.read(f"{path}qry_{self.name}_scar.fasta"    , "fasta")
        r_junc = SeqIO.read(f"{path}qry_{self.name}_junction.fasta", "fasta")

        record_dict['scar']     = r_scar
        record_dict['junction'] = r_junc

        return record_dict

class CassetteGenerator:

    def __init__(self, path_ref, path_write, cassette_dict, L = 3000):

        self.path_ref = path_ref
        self.path_write = path_write
        self.cassette_dict = cassette_dict
        self.L = L

    def generate_cassettes(self):

        train_record = SeqIO.read(self.path_ref, "fasta")

        for name,coords in self.cassette_dict.items():

            cc = Cassette(name, train_record.seq, 
                          coords[0], coords[1],
                          coordinate_system = 'snapgene')
            
            cc.write2file(self.path_write, int(self.L/2))

    def verify_cassettes(self):
        
        train_record = SeqIO.read(self.path_ref, "fasta")

        for name,coords in self.cassette_dict.items():

            cc = Cassette(name, train_record.seq, 
                          coords[0], coords[1],
                          coordinate_system = 'snapgene')
            
            record_dict = cc.readFromFile(self.path_write)

            print(f"\n\n{name}:")
            print("------------------------------------------------------")

            for qtype, rec in record_dict.items():

                ls, rs = self._split_string(rec.seq, n=20)

                print(f"\t{qtype.ljust(10)}: {ls} --- {rs}")

    @staticmethod
    def _split_string(s, n=0):

        length = len(s)
        midpoint = length // 2
        left = s[:midpoint]
        right = s[midpoint:]

        if n > 0:

            ls = left[-n:]
            rs = right[:n]
            return ls, rs
        
        else:

            return left, right
        
    @staticmethod
    def generate_random_cassettes(k=10, L=150000):

        
        start_interval = 1100000
        stop_interval  = 1500000

        for i in range(1, k):

            while True:

                start = random.randint(start_interval, stop_interval)
                stop  = random.randint(start_interval, stop_interval)

                if L*2 > stop - start > L: break

            print(f'"rnd{i}": ({start}, {stop}),')

class Processor:

    def __init__(self, datamat = None):
            
        self.datamat = datamat

    def step0_do_the_mappings(self, forceDvarGen):

        IGVM: MapperAN1

        print("Doing the mappings...")

        dit = self.datamat.gd(print_dvar=False, 
                              forceDvarGen=forceDvarGen)

        for name, dvar, IGVM, dataset in dit:

            print(f"\n\n>>>>>>>>>>>>> {name}|{dataset} >>>>>>>>>>>>>\n")
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")

            print(dvar)
            IGVM.execute_all()
            IGVM.calculate_coverage()

    def step0_create_the_dataframe(self, tol = 30, mapq = 0):

        def extract_core_names(name):

            m = re.search(r"^qry_(.+)_(junction|scar)$", name)

            return m.group(1), m.group(2)

        rlmh_start, rlmh_end = 2584176, 2584655
        
        write_rows = True

        do_rlmh = False; write_mows = False

        rows = list()
        mows = list()

        for name, dvar, IGVM, dataset in self.datamat.gd(print_dvar=False):

            if name.startswith("qry_"):

                core_name, qtype = extract_core_names(name)

                c, mq, nq      = self.count_reads_at_center(IGVM.path_sortd, 
                                                            tolerance = tol, 
                                                            mapq_threshold = mapq)
                
                allreads_count = 0 
                row = [dataset, core_name, 
                        qtype, c, allreads_count, 
                        mq, nq]
                
                rows.append(row)

            if name.startswith("RP62A") and do_rlmh:

                rate, TOTAL    = IGVM.read_overall_alignment_rate()
                allreads_count = 0 
                bamfile = pysam.AlignmentFile(IGVM.path_sortd, "rb") 
                out = bamfile.get_index_statistics()
                
                n = 0
                for read in bamfile.fetch(reference=out[0].contig, 
                                                start=rlmh_start, 
                                                end=rlmh_end):
                    
                    c1 = not read.is_unmapped
                    c2 = read.mapping_quality >= 0
                    if c1 and c2: n += 1
                    else: pass

                bamfile.close()

                row = [dataset, "rlmh", 
                        "N/A", n, allreads_count, 
                        "N/A", "N/A"]
                
                mows.append([dataset, rate, TOTAL])
                rows.append(row)
            
        if write_mows:
                
            columns = ['Dataset', 'Rate', 'Total']
            
            df = pd.DataFrame(mows, columns = columns)
            df.sort_values(by=['Dataset'], inplace=True)
            df.to_csv(self.datamat.path_ar, index=False)

        if write_rows:

            columns = ['Dataset', 'Name', 'Type', 'Counts', 'All Reads',
                        'Median Quality', '95th Percentile']
            
            df = pd.DataFrame(rows, columns = columns)
            df.sort_values(by=['Dataset','Name'], inplace=True)
            df.to_csv(self.datamat.path_df, index=False)
    
    def step1_RP62A_coverage_plots(self):

        for name, dvar, IGVM, dataset in self.datamat.gd(print_dvar=False):
                                    
            if name == "RP62A":

                print(dvar)

                path_save = self.datamat.path_plots + f"cov_{dataset}_{name}.png"
                
                IGVM.parse_and_plot_coverage_regular(path_save=path_save, dtype="RPM")

    def step1_assemble_unmapped_reads(self):

        print("Assembling Unmapped Reads...")
        
        for name, dvar, IGVM, dataset in self.datamat.gd(print_dvar=False):
                                    
            if name == "RP62A":

                print(dvar)
                IGVM.export_unmapped_read_pairs_as_fastq()
                IGVM.assemble_unmapped_read_pairs()

    
    def step1_plot_v1(self, epsilon = 0.05):

        

        def plot(qtype, epsilon, df, anove_dict=None):

            df = df.stack().reset_index()

            df.columns = ['Name', 'Dataset', 'Value']
            df['Dataset'] = df['Dataset'].str.extract('(alone|andhra|southeast)')

            categories = df['Dataset'].unique()
            df['Dataset'] = df['Dataset'].map({category: i for i, category in enumerate(categories)})

            names = df['Name'].unique()
            df['Name'] = df['Name'].map({name: i*epsilon for i, name in enumerate(names)})

            print("Processed df\n", df, "\n\n")

            colors = matplotlib.colormaps['tab10']

            plt.figure(figsize=(10, 6))
            for i, name in enumerate(names):
                plt.scatter(df[df['Name'] == i*epsilon]['Dataset'] + 
                            df[df['Name'] == i*epsilon]['Name'], 
                            df[df['Name'] == i*epsilon]['Value'],
                            color=colors(i), label=name)

            plt.xticks(np.arange(len(categories)), categories)
            plt.title(f"ANOVA Plot for {qtype}s")
            plt.ylim(0, 60)
            plt.legend()

            rlmh_present = False

            if 'rlmh' in df['Name'].unique(): rlmh_present = True

            plt.savefig(self.datamat.path_plots + 
                        f"anova_{qtype}s_rlmh_present_{rlmh_present}.png")

        dfs_n, dfj_n = self.get_normalized_dfs()

        for qtype, df in [('scar', dfs_n), ('junction', dfj_n)]: 
            
            plot(qtype, epsilon, df)

    def step1_plot_v2(self):

        x_dist = 2

        dfs_n, dfj_n = self.get_normalized_dfs()

        new_order = ['all', 
                     'crispr',
                     'crisprnhi',
                     'nhi',
                     'nhimec', 
                     'scc_mec',
                     'cassette1', 
                     'cassette2', 
                     'rnd1', 
                     'rnd2', 
                     'rnd3', 
                     'rnd4', 
                     'rnd5']    

        dfs_n.set_index('Name', inplace=True)
        dfj_n.set_index('Name', inplace=True)

        dfs_n = dfs_n.reindex(new_order)
        dfj_n = dfj_n.reindex(new_order)

        dfs_n.reset_index(inplace=True)
        dfj_n.reset_index(inplace=True)

        dfs_n['x'] = dfs_n.index * 10 + x_dist
        dfj_n['x'] = dfj_n.index * 10 - x_dist

        dfs_n = dfs_n.round(3)
        dfj_n = dfj_n.round(3)

        print(dfs_n, "x\n\n")
        print(dfj_n, "x\n\n")

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 8), dpi=400)

        for i in dfs_n.index:

            vs = [dfs_n.iloc[i][col] for col in dfs_n.columns if col not in ['Name', 'x']]
            xs = [dfs_n.iloc[i]['x'] for a in vs]

            vj = [dfj_n.iloc[i][col] for col in dfj_n.columns if col not in ['Name', 'x']]
            xj = [dfj_n.iloc[i]['x'] for a in vj]

            pa = False
            my_alpha = 0.3
            my_dot_size = 40
            color_scar = 'red'
            color_junc = 'blue'

            ax.scatter(xs, vs, marker='o', color=color_scar, alpha=my_alpha, s=my_dot_size)
            ax.scatter(xj, vj, marker='o', color=color_junc, alpha=my_alpha, s=my_dot_size)

            medianprops = dict(linestyle='-', linewidth=2.5, color='black') 

            box_s = ax.boxplot(vs, positions=[dfs_n.iloc[i]['x']], widths=1, showfliers=False, zorder=2, medianprops=medianprops, patch_artist=pa)
            box_j = ax.boxplot(vj, positions=[dfj_n.iloc[i]['x']], widths=1, showfliers=False, zorder=2, medianprops=medianprops, patch_artist=pa)

            plt.savefig("fig.png")

        xticks =  (dfj_n['x'] + x_dist).to_list()
        xlabels = dfs_n['Name'].to_list()
        xlabels = [x.replace("cassette", "distal_") for x in xlabels]

        ax.set_xlim(left=-10, right=130)

        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
        ax.set_xlabel("Cassettes")
        ax.set_ylabel("Normalized Counts")

        for tick in ax.get_xticklabels(): tick.set_rotation(90)

        ax.grid(which='major', linestyle=':', linewidth='0.5', color='black')

        plt.tight_layout()
        plt.savefig(self.datamat.path_plots + "cvcv.png")
 
    def get_normalized_dfs(self):

        path_scar = self.datamat.path_results + "scar_norm.csv"
        path_junc = self.datamat.path_results + "junction_norm.csv"

        if os.path.isfile(path_scar) and os.path.isfile(path_junc):

            dfs_n = pd.read_csv(path_scar)
            dfj_n = pd.read_csv(path_junc)

        else:

            dfs   = self.transform_dataframe_type2(qtype="scar", mode=1)
            dfs_n = self.deseq2_norm(dfs)
            dfs_n.to_csv(path_scar)

            dfj   = self.transform_dataframe_type2(qtype="junction", mode=1)
            dfj_n = self.deseq2_norm(dfj)
            dfj_n.to_csv(path_junc)

        return dfs_n, dfj_n
        
    def transform_dataframe_type1(self, qtype = "junction", mode = 1):
        
        if mode == 0:

            df = (pd.read_csv(self.datamat.path_df)
                    .assign(norm_counts=lambda df: 1e6 * df['Counts'] / df['All Reads'])
                    .pivot(index=['Dataset', 'Name'], columns='Type', values='norm_counts')
                    .fillna(0)
                    .reset_index()
                    .query('Dataset.str.contains("control") == False')
                    )
                    
        if mode == 1:

            df = (pd.read_csv(self.datamat.path_df)
                    .pivot(index=['Dataset', 'Name'], columns='Type', values='Counts')
                    .fillna(0)
                    .reset_index()
                    .query('Dataset.str.contains("control") == False')
                    )

        return df

    def transform_dataframe_type2(self, qtype = "junction", mode = 1):

        df = pd.read_csv(self.datamat.path_df)

        if mode == 0:

            df['norm_counts'] = 1e6 * df['Counts'] / df['All Reads']
            val = 'norm_counts'

        if mode == 1:

            val = 'Counts'

        df = (df.pivot(index=['Name', 'Type'], columns='Dataset', values=val)
            .fillna(0)
            .reset_index()
            .query(f"Type == '{qtype}' or Name == 'rlmh'")
            .drop(['Type'], axis=1)
            .set_index("Name"))
        
        if "rlmh" in df.index:
        
            df.loc['rlmh'] = df.loc['rlmh'] / 10
            
        return df

    def step1_create_pie_charts(self):

        def plot(df_scar, df_junc, ds):

            fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

            zipped = zip([df_scar, df_junc], [0, 1], ['Scar', 'Junction'])

            for (df, i, dsType) in zipped:
                
                
                ax[i].pie(df['norm_counts'], 
                          labels=df['Name'], 
                          autopct='%1.1f%%', 
                          startangle=90)

                ax[i].set_title(dsType)
                ax[i].axis('equal')

            plt.tight_layout()
            plt.savefig(f"{self.datamat.path_root}PLOTS/{ds}.png")

        

        df = pd.read_csv(self.datamat.path_df)
        df['norm_counts'] = 1e6 * df['Counts'] / df['All Reads']
        
    @staticmethod
    def count_reads_at_center(path_sortd, tolerance, mapq_threshold):

        counts = 0
        qualities = list()
        lenghts = list()

        read: pysam.libcalignedsegment.AlignedSegment 
        bamfile: pysam.libcalignmentfile.AlignmentFile 
        bamfile = pysam.AlignmentFile(path_sortd, "rb") 

        out = bamfile.get_index_statistics()
        contig = out[0].contig 

        hrl = bamfile.lengths[0] // 2
        center_left  = hrl - 1
        center_right = hrl

        for read in bamfile.fetch(reference=contig, start=center_left, end=center_right):

            if (not read.is_unmapped) and (read.mapping_quality >= mapq_threshold):

                if read.reference_start < center_left-tolerance and \
                   read.reference_end   > center_right+tolerance:

                    qualities.append(read.mapping_quality)
                    lenghts.append(read.query_length) 
                    
        bamfile.close()

        counts = len(qualities)

        if counts > 0:
            return counts, median(qualities), np.percentile(qualities, 95)
        else:
            return counts, None, None
        
    @staticmethod
    def count_all_reads(path_sortd):

        result = subprocess.run(f"samtools view --threads 16 -c {path_sortd}", 
                                capture_output=True, 
                                text=True, 
                                check=True, 
                                shell=True)

        return int(result.stdout)
    
    @staticmethod
    def deseq2_norm(df_input):

        df = df_input.copy()
            
        with np.errstate(divide="ignore"): log_counts = np.log(df)
            
        logmeans = log_counts.mean(axis=1)

        
        filtered_genes = ~np.isinf(logmeans)

        
        log_ratios = log_counts.loc[filtered_genes, :] - logmeans[filtered_genes].values[:, np.newaxis]

        
        log_medians = np.median(log_ratios, axis=0)

        
        size_factors = np.exp(log_medians)
        deseq2_counts = df / size_factors

        return deseq2_counts
    
    @staticmethod
    def get_anova_pvalues(df):

        data_dict = dict()

        for name in df.index:

            data = list()
            for dset in ['alone', 'andhra', 'southeast']:

                data.append(df.filter(like=dset).loc[name].to_list())

            data_dict[name] = data
        
        anova_dict = dict()

        for name, data in data_dict.items():
                
            F, p = scipy.stats.f_oneway(*data)
            anova_dict[name] = {'F': F, 'p': p}

        return anova_dict
        
    @staticmethod
    def deseq2_norm_copilot_version(df_input):

        

        df = df_input.copy()

        
        gm = scipy.stats.gmean(df)

        
        log_ratios = np.log(df / gm)

        
        size_factors = np.exp(np.median(log_ratios, axis=0))

        
        deseq2_counts = df / size_factors

        return deseq2_counts

if __name__ == "__main__":

    d1 = {"Name": "MAIN",
          
          "DataSets": ['andhra1', 'andhra2', 'andhra3',
                       'southeast1', 'southeast2', 'southeast3'],

          "Cassettes": {"crispr"    : (2493653, 2519738),
                        "nhi"       : (2519739, 2536572),
                        "scc_mec"   : (2536573, 2584193),
                        "all"       : (2493653, 2584193),
                        "nhimec"    : (2519739, 2584193),
                        "crisprnhi" : (2493653, 2536572),
                        "cassette1" : (1114161, 1259476),
                        "cassette2" : (1271306, 1309937),
                        "rnd1"      : (1190020, 1341940),
                        "rnd2"      : (1138053, 1312854),
                        "rnd3"      : (1105250, 1306664),
                        "rnd4"      : (1302176, 1470876),
                        "rnd5"      : (1170359, 1356097),}}
    
    j = 1

    
    if j == 0:

        c = {"rnd4"      : (1302176, 1470876),
             "rnd5"      : (1170359, 1356097)}

        print(c)

        p_ref = "/Volume/biodata2/sccmec_mappings/QUERIES_MAIN/RP62A.fasta"
        p_wrt = "/Volume/biodata2/sccmec_mappings/"

        CS = CassetteGenerator(p_ref, p_wrt, c)
        CS.generate_cassettes()
        CS.verify_cassettes()

    if j == 1:

        data_mat = DataMat(d1)

        P = Processor(data_mat)
        
        P.step1_plot_v2()

    if j == 2:

        for ss in sys.path: print(ss)

    
    if j == 3:

        import main_classes.assembly_pipeline as ap

        p = "/Volume/biodata2" \
            "/sccmec_mappings/mappings/alone3/mappings/RP62A/unmapped_assembly/"

        dvars = ap.dvarGenAssembly(p)

        for name in dvars:

            SA = ap.SpadesAssembler(dvars[name], output_index=0)
            
            
            SA.generate_summary(kind="scaffolds", L_cutoff=200, C_cutoff=2, sort_by_col='Length')
            

    
    if j == 4:

        dsets = ["alone1", "alone2", "alone3", 
                 "andhra1", "andhra2", "andhra3",
                 "southeast1", "southeast2", "southeast3"]
        
        for dset in dsets:

            p_root = "/Volume/biodata2" \
                     "/sccmec_mappings/"

            p_sour1 = p_root + f"mappings/{dset}" \
                      f"/mappings/RP62A/unmapped_assembly/assembly_results"
            
            p_sour2 = p_root + f"mappings/{dset}" \
                      f"/mappings/RP62A/unmapped_assembly/assembly" \
                      f"/spades_output_0/assembly_graph_with_scaffolds.gfa"
            
            p_dest = p_root + f"Exported_Unmapped_Results/{dset}/" \
                     f"assembly_results/"
            
            dryrun = False

            if dryrun:

                print("---------------------------------------------")
                print(p_sour2)
                print(os.path.isfile(p_sour2))
                print(p_dest)
                print()

            else:

                if not os.path.isdir(p_dest) or True:

                    shutil.copy2(p_sour2, p_dest)

                else:
                    print('\talready copied...')

    
    if j == 5:

        input_file = "/Volume/biodata2/analysis_QUO1006861/mappings_RP62A_cassettes/queries/qry_rlmh.fasta"
        output_file = "/Volume/biodata2/analysis_QUO1006861/mappings_RP62A_cassettes/queries/qry_rlmh_new.fasta"

        with open(output_file, "w") as out_handle:
            for record in SeqIO.parse(input_file, "fasta"):
                record.seq = record.seq.upper()
                SeqIO.write(record, out_handle, "fasta")

    if j == 6:

        pass

    if j == 7:

        path = "/Volume/biodata2" \
               "/sccmec_mappings/mappings_MAIN/alone1" \
               "/mappings/qry_cassette1_junction/bowtie2out_sorted.bam"
        
        for n in [0, 1, 5, 10, 20, 30]:

            c, mq, nq = Processor.count_reads_at_center(path, tolerance=n, mapq_threshold=0)
            print(f"{n} -> {c} | {mq} | {nq}")

    if j == 8:

        CassetteGenerator.generate_random_cassettes()
