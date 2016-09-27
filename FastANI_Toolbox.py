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
from datetime import datetime

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
            ref_dict[i] = {} # Store identity and alignment length
        f = open(filepath,"r")
        lines = [i.strip().split("\t") for i in f.readlines()]
        f.close()
        for line in lines:
            best_hit = True
            (frag_id, ref_id, identity, align, mismatchCount, gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, evalue, bitScore) = line
            identity = float(identity)/100
            align = int(align) - int(gapOpenCount)
            align_pct = float(align)/1020
            if identity >= 0.3 and align_pct >= 0.7:
                if frag_id not in ref_dict[ref_id]:
                    ref_dict[ref_id][frag_id] = [identity, align]
                else:
                    continue
        ANI_dict = {}
        align_dict = {}
        for each_ref in ref_prefix:
            fragments = ref_dict[each_ref].keys()
            identities = [ref_dict[each_ref][fragment][0] for fragment in fragments]
            align_pcts = [ref_dict[each_ref][fragment][1] for fragment in fragments]
            align_dict[each_ref] = len(align_pcts)
            if len(identities) > 0:
                ANI_dict[each_ref] = np.mean(identities,dtype=np.float64)
        df_ANI = pd.DataFrame.from_dict(ANI_dict,orient='index')
        df_align = pd.DataFrame.from_dict(align_dict,orient='index')
        df_align.columns = ['ALIGNED']
        df_ANI.columns = ['ANI']

        self.ANI_table = df_ANI
        self.cov_table = df_align


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
    fragments_count = 0
    for i in range(0, total_size, frag_size):
        fragments_count += 1
        append_fragments(query_seq[i:i + frag_size])
    if len(fragments[-1]) != frag_size:
        fragments_count -= 1
        fragments = fragments[:-1]
    return fragments, fragments_count


def concate_reference_files(reference_folder,work_dir):
    """
    Write reference files into one
    :param reference_folder:
    :return: A mapping file for multi-processing, [each_query_filepath, ref_prefix]
    """
    reference = []
    append_reference = reference.append
    reference_files = [join(reference_folder,file) for file in listdir(reference_folder) if isfile(join(reference_folder, file))]
    concat_ref_file = open(join(work_dir,'concat_ref.fasta'), 'w')
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
    return mapping


def makeblastdb(work_dir):
    cmd = "makeblastdb -dbtype nucl -in {0}concat_ref.fasta -title ref_genome " \
          "-out {0}ref_genome_blastdb".format(join(work_dir,""))
    os.system(cmd)



def parse_blast_tab_get_best(ANI_table,prefix_query):
    sorted_ANI_table = ANI_table.sort_index(by=['ANI'], ascending=False)
    top_hit = sorted_ANI_table.index[0]
    best_match = sorted_ANI_table.index[1]
    best_ANI = sorted_ANI_table.loc[best_match,'ANI']
    best_ANI = "{0:.5f}%".format(best_ANI*100)
    return best_match, best_ANI

def single_blast_run(each_mapping):
    """
    Process single file against reference database
    :param file: file path
    :param work_dir:
    :return:
    """
    try:
        (query_filepath,ref_prefix) = each_mapping
        prefix_query = "".join(query_filepath.split("/")[-1].split(".")[:-1])
        handler_NewGenome = open(query_filepath,"r")
        print "Reading {0} genome".format(prefix_query)
        records_NewGenome = list(SeqIO.parse(handler_NewGenome,"fasta"))
        handler_NewGenome.close()
        concatenated_seq_NewGenome = ""
        print "Concatenating {0} sequences".format(prefix_query)
        for record in records_NewGenome:
            concatenated_seq_NewGenome = concatenated_seq_NewGenome + str(record.seq).replace("N","")
        print "Cutting {0} genome into consecutive 1020 bp fragments".format(prefix_query)
        fragments_NewGenome, num_fragments = split_query_seq(query_seq=concatenated_seq_NewGenome,frag_size=1020)
        num_fragments_recorder = open("num_fragments.tab", "a")
        num_fragments_recorder.write("{0}\t{1}\n".format(prefix_query,num_fragments))
        num_fragments_recorder.close()
        handler_fragments_NewGenome = open((prefix_query + "_query.fna"), "w")
        print "Writing {0} fragments into a new FASTA file".format(prefix_query)
        for i in range(len(fragments_NewGenome)):
            handler_fragments_NewGenome.write(">fragment_{0}\n".format(i))
            handler_fragments_NewGenome.write(fragments_NewGenome[i]+"\n")
        handler_fragments_NewGenome.close()
        print "Blasting {0} against the database".format(prefix_query)
        blastall_cmd = "blastall -p blastn -o {0}_result.tab -i {1} -d {2} " \
                       "-X 150 -q -1 -F F -e 1e-15 " \
                       "-m 8" \
            .format(prefix_query, prefix_query + "_query.fna", "ref_genome_blastdb")
        os.system(blastall_cmd)
        FilePath_blast_tab = prefix_query + "_result.tab"
        print "Parsing blast result of {0}".format(prefix_query)
        blast_out_obj = blast_tab(FilePath_blast_tab, ref_prefix)
        ANI_table = blast_out_obj.ANI_table
        cov_table = blast_out_obj.cov_table.fillna(value=0)
        print "Retrieving the best match of {0}".format(prefix_query)
        best_match, best_ANI = parse_blast_tab_get_best(ANI_table=ANI_table,prefix_query=prefix_query)
        ANI_table.columns = [prefix_query]
        cov_table.columns = [prefix_query]
        print "Writing best match result of {0}".format(prefix_query)
        with open("best_match.tab","a") as best_table:
            best_table.write("{0}\t{1}\t{2}\n".format(prefix_query, best_match, best_ANI))
        ANI_table.to_csv("{0}_ANI.csv".format(prefix_query))
        cov_table.to_csv("{0}_aligned.csv".format(prefix_query))
    except:
        pass


def FastANI(argv=None):
    print datetime.now()
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
    pool = mp.Pool(processes=8)
    pool.map(single_blast_run, mapping)
    ref_prefix = mapping[0][1]
    df_ANI = pd.DataFrame(index=ref_prefix)
    pairwise_ANI_files = [file for file in listdir("./") if file.endswith("_ANI.csv")]
    print "Concatenating pairwise ANI result"
    for file in pairwise_ANI_files:
        each_df = pd.DataFrame.from_csv(file, header=0, index_col=0)
        colname = each_df.columns[0]
        df_ANI[colname] = each_df
    df_ANI.to_csv("pairwise_ANI.csv")
    print "Concatenating coverage result"
    num_fragments = {}
    num_fragments_recorder = open("num_fragments.tab","r")
    lines = [i.strip().split("\t") for i in num_fragments_recorder.readlines()]
    num_fragments_recorder.close()
    for each_line in lines:
        num_fragments[each_line[0]] = float(each_line[1])
    cov_files = [file for file in listdir("./") if file.endswith("_aligned.csv")]
    df_cov = pd.DataFrame(index=ref_prefix)
    for file in cov_files:
        each_df = pd.DataFrame.from_csv(file,header=0, index_col=0)
        colname = each_df.columns[0]
        df_cov[colname] = each_df
        df_cov[colname] = df_cov[colname]/num_fragments[colname]
    df_cov.to_csv("total_aligned.csv")
    calibrated_ANI = pd.DataFrame(index=ref_prefix)
    for each_col in ref_prefix:
        for row in ref_prefix:
            calibrated_ANI.loc[row, each_col] = df_ANI.loc[row, each_col] * df_cov.loc[row, each_col]
    calibrated_ANI.to_csv("calibrated_ANI.csv")
    print datetime.now()


if __name__ == "__main__":
    FastANI()

