def main():
    from Bio import AlignIO
    import pandas as pd
    import argparse
    from os import path
    import re

    parser = argparse.ArgumentParser(description='Make summary on substitutions and deletions in alignmnet.')
    parser.add_argument('-a', '--aln', action='store', type=str, help='Alignment file in FASTA format.')
    parser.add_argument('-o', '--out_dir', action='store', default='.', type=str, help='Output dir.')
    parser.add_argument('-p', '--prefix', action='store', default='', type=str, help='Prefix for output files. Default: no prefix.')
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
    alignments_length = len(alignments_parsed[0].seq)

    accepted_aa = [
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

    aa_occurences = pd.DataFrame(
        0,
        index=range(alignments_length),
        columns=accepted_aa
    )

    # Count substitutions
    print('Count substitutions')
    for position in range(alignments_length):
        column = alignments_parsed[:, position]
        aa_set = set(column)
        aa_counts = [column.count(a_acid) for a_acid in aa_set]
        aa_dict = dict(zip(aa_set, aa_counts))
        for aa in aa_dict:
            aa_occurences.loc[position, aa] = aa_dict.get(aa)

    # Count deletions
    deletions = {}

    print('Count deletions')
    # Go through each record of alignment
    for record in alignments_parsed:
        # Count gaps (deletion events) with regexp
        deletion_match = [_ for _ in re.finditer('-+', str(record.seq))]
        if len(deletion_match):
            # Start poditions of deletions
            starts = [_.start() for _ in deletion_match]
            # Lengths of deletions
            del_len = [len(_.group()) for _ in deletion_match]
            for del_event in tuple(zip(starts, del_len)):
                # setdefault will set 1 count if event was not observed else count will be increased by 1
                if deletions.setdefault(del_event, 1) > 1:
                    deletions[del_event] += 1

    print('Write otput')
    prefix = ''
    if args.prefix != '':
        prefix = args.prefix + '_'

    # Output substitutions
    # Change from 0-based to 1-based positions
    aa_occurences.index += 1
    # In case there will be other letters except accepted letters additional processing of data frame is needed
    aa_occurences = aa_occurences.fillna(0)
    aa_occurences = aa_occurences[aa_occurences.columns].astype('int64')
    aa_occurences.transpose().to_csv(path.join(args.out_dir, prefix + 'stat_variant.tsv'), sep='\t')

    # Output deletions
    len_del = len(deletions)

    # Data Frame for output
    deletion_df = pd.DataFrame({
        'Start_pos': [0] * len_del,
        'Length': [0] * len_del,
        'Count': [0] * len_del
    })

    index = 0
    for event in deletions.items():
        # +1 change 0-based positions ot 1-based positions
        deletion_df.iloc[index, ] = [event[0][0] + 1, event[0][1], event[1]]
        index += 1

    deletion_df.to_csv(path.join(args.out_dir, prefix + 'stat_deletions.tsv'), sep='\t', index=False)


# Execute
if __name__ == '__main__':
    main()
