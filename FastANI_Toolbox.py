#!/usr/bin/python3

# !/usr/bin/python

# IMPORT
from Bio import SeqIO
import os
from os import listdir
from os.path import isfile, isdir, join


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


def concate_reference_files(top10, cursor):
    concat_ref_file = open('concat_ref.fasta', 'w')
    for each_genome_ID in top10:
        concat_ref_file.write(">{0}\n".format(each_genome_ID))
        cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(each_genome_ID))
        filepath = c.fetchone()[0]
        each_ref = open(filepath, "r")
        records = list(SeqIO.parse(each_ref, "fasta"))
        each_ref.close()
        seq = ""
        for record in records:
            seq_modified = str(record.seq).replace("N", "")
            seq = seq + seq_modified
        concat_ref_file.write(seq + "\n")
    concat_ref_file.close()


def makeblastdb():
    cmd = "makeblastdb -dbtype nucl -in concat_ref.fasta -title ref_genome -out ref_genome_blastdb"
    os.system(cmd)


def run_blastn(fragment):
    return


def FastANI(new_Geonme_ID, top10_list, cursor):
    cursor.execute("SELECT FilePath FROM Genome WHERE Genome_ID={0}".format(
        new_Geonme_ID))  # Get the sequence of the new genome by Genome_ID
    FilePath_NewGenome = cursor.fetchone()[0]
    handler_NewGenome = open(FilePath_NewGenome, "r")
    records_NewGenome = list(SeqIO.parse(handler_NewGenome, "fasta"))
    handler_NewGenome.close()
    concatenated_seq_NewGenome = ""
    for record in records_NewGenome:
        concatenated_seq_NewGenome = concatenated_seq_NewGenome + str(record.seq).replace("N", "")
    fragments_NewGenome = split_query_seq(query_seq=concatenated_seq_NewGenome, frag_size=1020)
    handler_fragments_NewGenome = open("query.fna", "w")
    for i in range(len(fragments_NewGenome)):
        handler_fragments_NewGenome.write(">fragment_{0}\n".format(i))
        handler_fragments_NewGenome.write(fragments_NewGenome[i] + "\n")
    handler_fragments_NewGenome.close()
    concate_reference_files(top10=top10_list, cursor=cursor)
    makeblastdb()
    blastall_cmd = "blastall -p blastn -o {0}.blast_tab -i {1} -d {2} " \
                   "-X 150 -q -1 -F F -e 1e-15 " \
                   "-b 1 -v 1 -m 8 -a 4".format("test", "query.fna", "ref_genome_blastdb")
    os.system(blastall_cmd)


