
# SYSTEM REQUIREMENTS AND SOFTWARE DEPENDENCIES

Any reasonably modern computer with a python environment should be able to run these scripts (sccmec_main.py, sccmec_mappings.py). 

The main script depends on several standard and third-party packages. The primary dependencies include:

- **Standard libraries**: `os`, `re`, `ast`, `time`, `json`, `random`, `csv`, `shutil`, `pickle`, `math`, `operator`, `statistics`, `collections`, `itertools`, `functools`, `subprocess`, `multiprocessing`, and `dataclasses`.
- **Third-party packages**: `pandas`, `numpy`, `Levenshtein`, `matplotlib`, `tqdm`, `scipy`, and `biopython`.

Other software tools required (e.g. hmmer, samtools, bowtie2, etc) are listed in the main manuscript along with their versions.

Furthermore, macsyfinder and defense-finder can be found here: 

- https://github.com/gem-pasteur/macsyfinder
- https://github.com/mdmparis/defense-finder

# INSTALLATION GUIDE

No installation necessary. The script can be downloaded and run on any python environment that has the necessary python packages pre-installed. However, mdm_df.zip must be unzipped and placed in the same directory as the script.

# DEMO

The main script(s) are a collection of custom functions and classes that are meant to be used in an exploratory fashion. The concept of small dataset for demo is inapplicable. They are written in a clear, standalone, linear fashion that makes the flow of the code itself the best instructions of how to use it.

For demo, the user can download one or both the genome file collections (full genbank files) that are listed in the manuscript (also here on github: list_aureus.txt, list_epidermidis.txt) and follow the instructions below to analyze them. The expected outputs should match the results in the manuscript. 

Some steps of the pipeline can take a few hours depending on the system.

# INSTRUCTIONS FOR USE

The `SnailAnalyzerFull` class is the core of the pipeline, orchestrating a multi-step analysis that generates the polar plot and other results. Once properly initialized as an instance of the class, it is meant to be run sequentially, step by step, as the names of the methods suggest.

1. First generate a json file that looks like the example at the end of this README file
2. Create a "accession_list.txt" file (i.e. list_aureus.txt, list_epidermidis.txt. They need to be renamed to accession_list.txt) that contains the accession numbers of the genomes you want to analyze, and place it in the directory specified in the json file.

An example would be to append the following at the end of the main script, and run it by enabling the desired steps (i.e. uncommenting them):

```python

if __name__ == "__main__":

    pr = "/path/to/json_file.json"
    dvars = load_dvars(pr)
    dvar = dvars['epidermidis'] # or dvars['aureus']

    SA = SnailAnalyzerFull(dvar)
    SA.step1_download_all_genbank()
    # SA.step2_generate_mapped_csv_from_collection()
    # SA.step3a_extract_rlmh_nbrhd_fasta_files()
    # SA.step3b_extract_rlmh_all_proteins_fasta_files()
    # SA.step4a_run_defense_finder()
    # SA.step4b_run_defense_finder_on_all_proteins()
    # SA.step6a_generate_levenshtein_groups_dict()
    # SA.step6b_generate_levenshtein_all_groups_dict()
    # SA.step7a_generate_core_dictionaries_from_lvs_groups_dict()
    # SA.step7b_generate_core_dictionaries_all()
    # SA.step8_generate_dataframes()
    # SA.step9_plot_snail_diagram()
```

# Appendix: Example JSON file

```json
{
    "epidermidis": {
        "genomes": "/biodata2/genbank_collection/genbank_epidermidis_collection/",
        "root": "/biodata2/sccmec_snail_analysis/snail_analysis_epidermidis/",
        "acc_list": "/biodata2/sccmec_snail_analysis/genbank_collection/genbank_aureus_collection/",
        "win_size": 500,
        "numthreads": 64,
        "plot_params": {
            "n": 89,
            "m": 355,
            "th": 85,
            "c_depth": 100,
            "font_size_1": 8,
            "font_size_2": 8,
            "res_dpi": 500,
            "fw": 16,
            "fh": 12,
            "annotate_nums": true,
            "annotate_defs": false,
            "circle": true
        }
    },
    "aureus": {
        "genomes": "/biodata2/genbank_collection/genbank_aureus_collection/",
        "root": "/biodata2/sccmec_snail_analysis/snail_analysis_aureus/",
        "acc_list": "/biodata2/sccmec_snail_analysis/genbank_collection/genbank_aureus_collection/",
        "win_size": 500,
        "numthreads": 64,
        "plot_params": {
            "n": 982,
            "m": 355,
            "th": 930,
            "c_depth": 100,
            "font_size_1": 2,
            "font_size_2": 2,
            "res_dpi": 500,
            "fw": 16,
            "fh": 12,
            "annotate_nums": true,
            "annotate_defs": false,
            "circle": false
        }
    }
}
```
