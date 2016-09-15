#!/usr/bin/python3

# !/usr/bin/python

# IMPORT
from Bio import SeqIO
import os
from os import listdir
from os.path import isfile, isdir, join
import pandas as pd
import numpy as np

# OBJECT
class blast_tab(object):
    def __init__(self, filepath, ref_prefix):
        self.filepath = filepath
        self.ref_prefix = ref_prefix
        if os.path.isfile(self.filepath):
            self.parse()
    def parse(self, filepath = self.filepath, ref_prefix=self.prefix):
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


def concate_reference_files(reference_folder, work_dir):
    reference = []
    append_reference = reference.append
    reference_files = [file for file in listdir(reference_folder) if isfile(join(reference_folder, file))]
    concat_ref_file = open(join(work_dir,'concat_ref.fasta'), 'w')
    for file in reference_files:
        ref_prefix = "".join(file.split(".")[:-1])
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
    return ref_prefix


def makeblastdb(work_dir):
    cmd = "makeblastdb -dbtype nucl -in {0}/concat_ref.fasta -title ref_genome " \
          "-out {0}/ref_genome_blastdb".format(work_dir)
    os.system(cmd)


def parse_blast_tab(prefix_query, ref_prefix, work_dir):
    FilePath_blast_tab = join(work_dir, "result.tab")
    ANI_table = blast_tab(FilePath_blast_tab,ref_prefix).ANI_table
    top_hit = ANI_table.index[0]
    if top_hit != prefix_query:
        best_match = top_hit
        best_ANI = ANI_table.loc[best_match,'ANI']
    else:
        best_match = ANI_table.index[1]
        best_ANI = ANI_table.loc[best_match,'ANI']
    return best_match, best_ANI




def FastANI(query_file_path, reference_folder, work_dir):
    FilePath_NewGenome = query_file_path
    prefix_query = "".join(query_file_path.split("/")[-1].split(".")[-1])
    handler_NewGenome = open(FilePath_NewGenome, "r")
    records_NewGenome = list(SeqIO.parse(handler_NewGenome, "fasta"))
    handler_NewGenome.close()
    concatenated_seq_NewGenome = ""
    for record in records_NewGenome:
        concatenated_seq_NewGenome = concatenated_seq_NewGenome + str(record.seq).replace("N", "")
    fragments_NewGenome = split_query_seq(query_seq=concatenated_seq_NewGenome, frag_size=1020)
    handler_fragments_NewGenome = open(join(work_dir, prefix_query, "query.fna"), "w")
    for i in range(len(fragments_NewGenome)):
        handler_fragments_NewGenome.write(">fragment_{0}\n".format(i))
        handler_fragments_NewGenome.write(fragments_NewGenome[i] + "\n")
    handler_fragments_NewGenome.close()
    ref_prefix = concate_reference_files(reference_folder=reference_folder, work_dir=work_dir)
    makeblastdb(work_dir=work_dir)
    blastall_cmd = "blastall -p blastn -o {0}/result.tab -i {1} -d {2} " \
                   "-X 150 -q -1 -F F -e 1e-15 " \
                   "-b 1 -v 1 -m 8 -a 4"\
        .format(work_dir, join(work_dir,prefix_query,"_query.fna"), join(work_dir,"ref_genome_blastdb"))
    os.system(blastall_cmd)
    best_match, best_ANI = parse_blast_tab(prefix_query=prefix_query, ref_prefix=ref_prefix, work_dir=work_dir)


