# !/usr/bin/python

# IMPORT
from Bio import SeqIO
import os
from os import listdir
from os.path import isfile, isdir, join
import pandas as pd
import numpy as np
import multiprocessing as mp
import argparse
import sys

# OBJECT
class blast_tab(object):
    def __init__(self, filepath, ref_prefix):
        self.filepath = filepath
        self.ref_prefix = ref_prefix
        if os.path.isfile(self.filepath):
            self.parse()
    def parse(self, filepath = None, ref_prefix= None):
        if filepath is None:
            filepath = self.filepath
        if ref_prefix is None:
            ref_prefix = self.ref_prefix
        ref_dict = {}
        for i in ref_prefix:
            ref_dict[i] = [] # Store identity and alignment length
        f = open(filepath,"r")
        lines = [i.split("\t") for i in f.readlines()]
        f.close()
        for line in lines:
            (frag_id, ref_id, identity, align) = line[:4]
            identity = float(identity)/100
            align_pct = float(int(align))/1020
            if identity >= 0.3 and align_pct >= 0.7:
                ref_dict[ref_id].append(identity)
        ANI_dict = {}
        for each_ref in ref_prefix:
            ANI_dict[each_ref] = np.mean(ref_dict[each_ref],dtype=np.float64)
        df = pd.DataFrame.from_dict(ANI_dict,orient='index')
        df.columns = ['ANI']
        new_df = df.sort_index(by=['ANI'],ascending=False)
        self.ANI_table = new_df


# FUNCTION
def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="FastANI"
    )
    parser.add_argument("-i", dest="reference_folder", help="Folder where the FASTA files are")
    parser.add_argument("-o", dest="work_dir", help="Output dir")
    args = parser.parse_args()
    return args

def split_query_seq(query_seq, frag_size=1020):
    """
    input: object from concatenated query genome string
    output: a list of fragmented sequences
    If the last fragment is not exactly 1020 bp, it will be discarded
    """
    total_size = len(query_seq)
    fragments = []
    append_fragments = fragments.append
    for i in range(0, total_size, frag_size):
        append_fragments(query_seq[i:i + frag_size])
    if len(fragments[-1]) != frag_size:
        fragments = fragments[:-1]
    else:
        fragments = fragments
    return fragments


def concate_reference_files(reference_folder,work_dir):
    """
    Write reference files into one
    :param reference_folder:
    :return: A mapping file for multi-processing, [each_query_filepath, ref_prefix]
    """
    reference = []
    append_reference = reference.append
    reference_files = [join(reference_folder,file) for file in listdir(reference_folder) if isfile(join(reference_folder, file))]
    concat_ref_file = open(join(workdir,'concat_ref.fasta'), 'w')
    for file in reference_files:
        ref_prefix = "".join(file.split("/")[-1].split(".")[:-1])
        append_reference(ref_prefix)
        concat_ref_file.write(">{0}\n".format(ref_prefix))
        each_ref = open(file, "r")
        records = list(SeqIO.parse(each_ref, "fasta"))
        each_ref.close()
        seq = ""
        for record in records:
            seq_modified = str(record.seq).replace("N", "")
            seq = seq + seq_modified
        concat_ref_file.write(seq + "\n")
    concat_ref_file.close()
    mapping = [[i, reference] for i in reference_files]
    print mapping
    return mapping


def makeblastdb(workdir):
    cmd = "makeblastdb -dbtype nucl -in concat_ref.fasta -title ref_genome " \
          "-out {0}ref_genome_blastdb".format(join(workdir,""))
    os.system(cmd)


def parse_blast_tab_get_best(ANI_table):
    top_hit = ANI_table.index[0]
    if top_hit != prefix_query:
        best_match = top_hit
        best_ANI = ANI_table.loc[best_match,'ANI']
    else:
        best_match = ANI_table.index[1]
        best_ANI = ANI_table.loc[best_match,'ANI']
    best_ANI = "{0:.5f}%".format(best_ANI*100)
    return best_match, best_ANI

def single_blast_run(each_mapping):
    """
    Process single file against reference database
    :param file: file path
    :param work_dir:
    :return:
    """
    (query_filepath,ref_prefix) = each_mapping
    prefix_query = "".join(query_filepath.split("/")[-1].split(".")[:-1])
    handler_NewGenome = open(query_filepath,"r")
    records_NewGenome = list(SeqIO.parse(handler_NewGenome,"fasta"))
    handler_NewGenome.close()
    concatenated_seq_NewGenome = ""
    for record in records_NewGenome:
        concatenated_seq_NewGenome = concatenated_seq_NewGenome + str(record.seq).replace("N","")
    fragments_NewGenome = split_query_seq(query_seq=concatenated_seq_NewGenome,frag_size=1020)
    handler_fragments_NewGenome = open((prefix_query + "_query.fna"), "w")
    for i in range(len(fragments_NewGenome)):
        handler_fragments_NewGenome.write(">fragment_{0}\n".format(i))
        handler_fragments_NewGenome.write(fragments_NewGenome[i]+"\n")
    handler_fragments_NewGenome.close()
    blastall_cmd = "blastall -p blastn -o {0}_result.tab -i {1} -d {2} " \
                   "-X 150 -q -1 -F F -e 1e-15 " \
                   "-b 1 -v 1 -m 8 -a 4" \
        .format(prefix_query, prefix_query + "_query.fna", "ref_genome_blastdb")
    os.system(blastall_cmd)
    FilePath_blast_tab = prefix_query + "_result.tab"
    ANI_table = blast_tab(FilePath_blast_tab, ref_prefix).ANI_table
    best_match, best_ANI = parse_blast_tab_get_best(ANI_table=ANI_table)
    ANI_table.columns = [prefix_query]
    with open("best_match.tab","a") as best_table:
        best_table.write("{0}\t{1}\t{2}\n".format(prefix_query, best_match, best_ANI))
    ANI_table.to_csv("{0}_ANI.csv".format(prefix_query))


def FastANI(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()
    reference_folder = args.reference_folder
    work_dir = args.work_dir
    if isdir(work_dir):
        os.system("rm -rf {0}".format(work_dir))
        os.mkdir(work_dir)
    else:
        os.mkdir(work_dir)
    # Prepare blast db
    mapping = concate_reference_files(reference_folder=reference_folder, work_dir=work_dir)
    makeblastdb(work_dir=work_dir)
    os.chdir(work_dir)
    # Multi-processing
    pool = mp.Pool(processes=4)
    pool.map(single_blast_run, mapping)
    ref_prefix = mapping[0][1]
    df = pd.DataFrame(index=ref_prefix)
    pairwise_ANI_files = [file for file in listdir("./") if file.endswith(".csv")]
    for file in pairwise_ANI_files:
        each_df = pd.DataFrame.from_csv(file, header=0, index_col=0)
        colname = each_df.columns[0]
        df[colname] = each_df
    df.to_csv("pairwise_ANI.csv")

if __name__ == "__main__":
    FastANI()

