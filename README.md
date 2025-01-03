# Sequence Logo Generator for Splice Sites

This script generates sequence logos to visualize the information content and base distribution at splice sites (5' and 3') from a given set of sequences. The logos use custom images for nucleotide bases (A, C, G, T) to represent their relative contributions.

## Features
- Reads input sequences from a FASTA file.
- Calculates relative frequencies, Shannon entropy, and information content for splice site bases.
- Generates sequence logos for 5' and 3' splice sites using base images.
- Outputs a high-resolution PNG figure.

## Usage
### Input Files
- **Sequences File** (`--sequenceFile`): Path to a FASTA file containing splice site sequences.
- **A, C, G, T Files** (`--Afile`, `--Cfile`, `--Gfile`, `--Tfile`): Paths to PNG images for nucleotides A, C, G, and T.

### Output
- **Output File** (`--outFile`): The filename for the generated sequence logo (default: `logo.png`).

### Command-Line Arguments
```bash
python script_name.py -i <sequence_file> -A <A_file> -C <C_file> -G <G_file> -T <T_file> -o <output_file>
