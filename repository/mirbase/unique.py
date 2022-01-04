import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='Combine redundant sequences with a cumulative header')
parser.add_argument('-i', '--input', nargs=1, required=True, help='input file path')
parser.add_argument('-o', '--output', nargs=1, required=True, help='output file path')
args = parser.parse_args()

sequences = {}

# @todo: maybe it would be better to use accession IDs and not miR names

# go through all sequences, store them as key of a dict and store all IDs in a list at this
# position. So we get a dict of sequences with each a list of their IDs
for record in SeqIO.parse(args.input[0], "fasta"):
    if str(record.seq) in sequences:
        sequences[str(record.seq)].append(record.id)
    else:
        sequences[str(record.seq)] = [record.id]


# create a seq object of each sequence and join the IDs to one string (| separated)
uniqueSequences = []
for sequence in sequences:
    seq = SeqRecord(Seq(sequence),
        id='|'.join(sequences[sequence]),
        description="")
    uniqueSequences.append(seq)

# write all sequences to the given output file
SeqIO.write(uniqueSequences, args.output[0], "fasta")
