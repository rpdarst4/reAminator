reAminator is a collection of command-line scripts to process long-read (200+ bp), high-throughput bisulfite genomic sequencing. reAminator comprises 6 Python scripts:

barcode.py  -- label reads by barcode
bsAlign.py  -- align reads to a library of references
bsDraw.py   -- extract aligned data from bsAlign.py output,
               filtering for length and percent deamination
meMapper.py -- Python wrapper to call MethylMapper from bsDraw.py
meTools.py  -- search for methylation patterns of interest
omelet.py   -- draw PNG images of methylation patterns

REQUIREMENTS

reAminator requires Python 2.7 with the BioPython module, a local installation of BLAST+, and FASTools. These may be obtained at

http://biopython.org/wiki/Main_Page
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
http://genome.ufl.edu/rivalab/fastools/

Optional dependencies are MethylMapper and the Python Imaging LIbrary:

http://www.pythonware.com/products/pil/

Input files for barcode.py, bsAlign.py, and bsDraw.py must be in FASTA format.

INSTALLATION

Open the file "reAminator.cfg" in a text editor and update the paths to local BLAST+ and FASTools. Optional: update path to MethylMapper, and set whether Python Imaging Library is present (for omelet.py).

USAGE

usage: barcode.py [-h] [-d DIRECTORY] [-n] [-p] reads codes

positional arguments:
  reads                 FASTA-format file of 454 reads
  codes                 FASTA-format file of barcodes

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        output file path
  -n, --names           track barcodes by ID (default by sequence)
  -p, --perfect         no mismatches



usage: bsAlign.py [-h] [-dest DEST] [-mask] [-pair PAIR] [-1] [-2] seqs refs

positional arguments:
  seqs            FASTA-format file of read sequences
  refs            FASTA-format file of reference sequences

optional arguments:
  -h, --help      show this help message and exit
  -dest DEST      output file name
  -mask           enable lowercase masking in references
  -pair PAIR      2nd FASTA-format file for paired-end sequencing
  -1, --a1_or_b1  convert BS-read G to A, i.e.seq. primer = a1 or b1
  -2, --a2_or_b2  convert BS-read C to T, i.e.seq. primer = a2 or b2



usage: bsDraw.py [-h] [-dest DEST] [-codes [CODES [CODES ...]]]
                 [-refs [REFS [REFS ...]]] [-strand STRAND]
                 [-bisulfite BISULFITE] [-length LENGTH] [-uniques]
                 [-Sites SITES] [-TSS] [-Weights WEIGHTS] [-Omelet]
                 [-window WINDOW]
                 [alignments [alignments ...]]

positional arguments:
  alignments            bsBlast output table file/s

optional arguments:
  -h, --help            show this help message and exit
  -dest DEST            output directory
  -codes [CODES [CODES ...]]
                        barcode seqs (leave blank for all)
  -refs [REFS [REFS ...]]
                        ref. IDs (leave blank for all)
  -strand STRAND        strands ("a" = C to T, "b" = G to A)
  -bisulfite BISULFITE  minim. percent deaminated (0-100)
  -length LENGTH        minim. bp length (default) or percent
  -uniques              use only unique non-deamination patterns
  -Sites SITES          enter alternate methylation dictionary (e.g.
                        CG=1,CC=1) NOTE: not compatible with MethylMapper
  -TSS                  read TSS distance after "+" symbol in ref. id (e.g.
                        >YFG+400)
  -Weights WEIGHTS      ratio of CG to GC weight for MethylMapper (e.g. 50,50)
  -Omelet               make omelet .png files
  -window WINDOW        basepair width of omelet .png files (or zero for
                        automatic)

EXAMPLE OF USE

To label read sequences in EXAMPLE/reads.fa with barcodes from EXAMPLE/barcodes.fa, type at the command line (omit the prompt "$ "):

	$ python barcode.py EXAMPLE/reads.fa EXAMPLE/barcodes.fa

To align the barcoded reads to sequences from EXAMPLE/references.fa:

	$ python bsAlign.py EXAMPLE/reads.barcodes.fa EXAMPLE/references.fa

To extract all FASTA alignment files from the database EXAMPLE/reads.references.db:

	$ python bsDraw.py EXAMPLE/reads.references.db

To extract all FASTA alignment files, and produce omelet.py PNG images:

	$ python bsDraw.py EXAMPLE/reads.references.db -O

(c) 2013, Russell Darst, University of FLorida 