# Protein alignment statiscitcs
Script producing protein alignment variants statistics.

## Usage

```
usage: analyze_aln.py [-h] [-a ALN] [-o OUT_DIR] [-p PREFIX] [-r REFERENCE] [-e EXCLUDE_FILE] [-i INCLUDE_FILE]

Make summary on substitutions and deletions in alignment.

optional arguments:
  -h, --help            show this help message and exit
  -a ALN, --aln ALN     Alignment file in FASTA format.
  -o OUT_DIR, --out_dir OUT_DIR
                        Output dir.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files. Default: no prefix.
  -r REFERENCE, --reference REFERENCE
                        ID of the reference protein. If set all variants will be attributed 
                        to the positions in the reference protein.
  -e EXCLUDE_FILE, --exclude_file EXCLUDE_FILE
                        File with IDs of genes that should be excluded from the analysis. 
                        One ID per line. Works only if "--include_file" is not set.
  -i INCLUDE_FILE, --include_file INCLUDE_FILE
                        File with IDs of genes that should be kept from the analysis. 
                        One ID per line. Overrides "--exclude_file" option.
```

Alignment should be in the Fasta format.

To remove some sequences from the statistics computing a file with ID that should be excluded could be specified. The file should have one ID per line.

Example:
```
b1234
b1235
b1346
```

## Output

Script outputs two files:
* `[OUT_DIR]/[PREFIX]_stat_variant.tsv` - table of amino acid counts for each alignment position. Columns are 1-based alignment positions and rows are one-letter amino acid codes or gaps. Each cell contains count of the corresponding letter in a corresponding position. If unusual letters will be observed in the alignment (letters out of amino acid alphabet like 'X' or 'Z') they will be preseved for alignment debugging purpose. 
* `[OUT_DIR]/[PREFIX]_stat_deletions.tsv` - table of deletion counts. Table has three columns:
  * Start_pos - position of alignment where deletion starts
  * Length - length of the deletion
  * Count - number of alignments where the deletion event was observed.
