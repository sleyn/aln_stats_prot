def main():
    from Bio import AlignIO
    import pandas as pd
    import argparse
    import os
    import re

    parser = argparse.ArgumentParser(description='Make summary on substitutions and deletions in alignmnet.')
    parser.add_argument('-a', '--aln', action='store', type=str, help='Alignment file in FASTA format.')
    parser.add_argument('-o', '--out_dir', action='store', default='.', type=str, help='Output dir.')
    parser.add_argument('-p', '--prefix', action='store', default='', type=str,
                        help='Prefix for output files. Default: no prefix.')
    parser.add_argument('-r', '--reference', action='store', default='', type=str,
                        help='ID of the reference protein. If set all variants will be attributed to the positions '
                             'in the reference protein.')
    parser.add_argument('-e', '--exclude_file', action='store', type=str, default='-',
                        help='File with IDs of genes that should be excluded from the analysis. One ID per line. '
                             'Works only if "--include_file" is not set.')
    parser.add_argument('-i', '--include_file', action='store', type=str, default='-',
                        help='File with IDs of genes that should be kept from the analysis. One ID per line. '
                             'Overrides "--exclude_file" option.')
    args = parser.parse_args()

    # Read IDs that need to be excluded and included
    exclude_id = []
    include_id = []

    # Flag if include file is set
    set_include = False

    if args.include_file != '-':
        set_include = True
        with open(args.include_file) as include_file:
            include_id = [line.strip() for line in include_file]
    elif args.exclude_file != '-':
        with open(args.exclude_file) as exclude_file:
            exclude_id = [line.strip() for line in exclude_file]


    # Open alignment file
    alignments = AlignIO.parse(args.aln, 'fasta')
    alignments_parsed = []

    # Information about reference
    ref_length = 0
    ref_aln_seq = ''

    # Exclude or include specified IDs
    # Check if reference sequence was found
    # False - not found
    # True - was found
    ref_check = False

    for aln in alignments:
        for record in aln:
            if set_include:
                if record.id in include_id:
                    alignments_parsed.append(record)
                    if record.id == args.reference:
                        ref_length = len(record.seq.ungap('-'))
                        ref_aln_seq = record.seq
                        ref_check = True
            elif record.id not in exclude_id:
                alignments_parsed.append(record)

    # Check if reference protein is found in the alignmnet
    if ref_check != bool(args.reference):
        print('Reference is set but not found.')
        exit(1)

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

    if args.reference == '':
        index_list = range(alignments_length)
    else:
        index_list = range(ref_length)

    aa_occurences = pd.DataFrame(
        0,
        index=index_list,
        columns=accepted_aa
    )

    # Count substitutions
    print('Count substitutions')
    ref_position = 0
    for position in range(alignments_length):
        # Skip potion if it is a gap in the reference
        if ref_check:
            if ref_aln_seq[position] == '-':
                continue
        column = alignments_parsed[:, position]
        aa_set = set(column)
        aa_counts = [column.count(a_acid) for a_acid in aa_set]
        aa_dict = dict(zip(aa_set, aa_counts))
        for aa in aa_dict:
            aa_occurences.loc[ref_position, aa] = aa_dict.get(aa)
        ref_position += 1

    # Count deletions
    deletions = {}

    print('Count deletions')
    # Go through each record of alignment
    for record in alignments_parsed:
        # If reference was set remove reference positions with gaps in alignment
        if ref_check:
            record_seq = ''
            record_seq_original = str(record.seq)
            for pos in range(len(ref_aln_seq)):
                if ref_aln_seq[pos] != '-':
                    record_seq += record_seq_original[pos]
        else:
            record_seq = str(record.seq)

        # Count gaps (deletion events) with regexp
        deletion_match = [_ for _ in re.finditer('-+', record_seq)]
        if len(deletion_match):
            # Start positions of deletions
            starts = [_.start() for _ in deletion_match]
            # Lengths of deletions
            del_len = [len(_.group()) for _ in deletion_match]
            for del_event in tuple(zip(starts, del_len)):
                # setdefault will set 1 count if event was not observed else count will be increased by 1
                if del_event in deletions:
                    deletions[del_event] += 1
                else:
                    deletions[del_event] = 1

    print('Write output')
    # Create output directory if it not exists
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    prefix = ''
    if args.prefix != '':
        prefix = args.prefix + '_'

    # Output substitutions
    # Change from 0-based to 1-based positions
    aa_occurences.index += 1
    # In case there will be other letters except accepted letters additional processing of data frame is needed
    aa_occurences = aa_occurences.fillna(0)
    aa_occurences = aa_occurences[aa_occurences.columns].astype('int64')
    aa_occurences.transpose().to_csv(os.path.join(args.out_dir, prefix + 'stat_variant.tsv'), sep='\t')

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

    deletion_df.to_csv(os.path.join(args.out_dir, prefix + 'stat_deletions.tsv'), sep='\t', index=False)


# Execute
if __name__ == '__main__':
    main()
