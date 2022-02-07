from Bio import AlignIO
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Make summary on substitutions in a reference sequence from an alignmnet.')
parser.add_argument('-a', '--aln', action='store', type=str, help='Alignment file in FASTA format.')
parser.add_argument('-o', '--out', action='store', type=str, help='Output file.')
parser.add_argument('-e', '--exclude_file', action='store', type=str, default='-', help='File with IDs of genes that should be excluded from the analysis. One ID per line.')
args = parser.parse_args()

# Read IDs that need to be excluded
exclude_id = []
if args.exclude_file != '-':
	with open(args.exclude_file) as exclude_file:
	    exclude_id = [line.strip() for line in exclude_file]

# Open alignment file
alignments = AlignIO.parse(args.aln, 'fasta')
alignments_parsed = []

# Exclude specified IDs
for aln in alignments:
    for record in aln:
        if record.id not in exclude_id:
            alignments_parsed.append(record)

alignments_parsed = AlignIO.MultipleSeqAlignment(alignments_parsed)
alignments_length = len(alignments_parsed._records[0].seq)

aa_occurences = pd.DataFrame(
    0,
    index=range(alignments_length),
    columns=[
        'G',
        'A',
        'L',
        'M',
        'F',
        'W',
        'K',
        'Q',
        'E',
        'S',
        'P',
        'V',
        'I',
        'C',
        'Y',
        'H',
        'R',
        'N',
        'D',
        'T',
        '-'
    ]
)

# Count substitutions
print('Count occurences')
for i in range(alignments_length):
    column = alignments_parsed[:, i]
    aa_set = set(column)
    aa_counts = [column.count(a_acid) for a_acid in aa_set]
    aa_dict = dict(zip(aa_set, aa_counts))
    for aa in aa_dict:
        aa_occurences.loc[i, aa] = aa_dict.get(aa)

print('Write tsv')
aa_occurences.transpose().to_csv(args.out, sep='\t')

