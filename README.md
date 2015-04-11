# MaP2D

Analysis of 2D signal in mutational profiling sequencing data

Ensure that you have the following installed:

* Novobarcode, part of the Novoalign software package, which is freely available for educational and not-for-profit use. Download the latest version of Novoalign at http://www.novocraft.com/support/download/
* ShapeMapper, software developed by the Weeks lab at UNC Chapel Hill for 1D analysis of mutational profiling data. Available at http://www.chem.unc.edu/rna/software.html
* nwalign, software for Needleman-Wunsch sequence alignment developed by brentp. Available at https://pypi.python.org/pypi/nwalign/
* numpy
* The RDATkit for handling RDAT data files. Available at https://github.com/hitrace/rdatkit

Add the MaP2D folder to your PATH.

To test on example sequencing data, go to the Tutorial folder and run:

    MaP2D.py sequence.fa RTBbarcodes.fa Sample1_S1_L001_R1_001_sub.fastq Sample1_S1_L001_R2_001_sub.fastq --offset 72 [--config example.cfg]
* 'sequence.fa' [required] is a fasta-formatted file with the name and sequence of the RNA.
* 'RTBbarcodes.fa' [required] is a Novobarcode-formatted file with the names and sequences of the barcodes in the RTB primers.
* The read 1 and read 2 FASTQs are required inputs.
* 'example.cfg' [optional] is a config file in the format required by ShapeMapper. If a config file is provided, ShapeMapper analysis will be performed to get 1D reactivity profiles. The correct names of the conditions in the pilot experiment used for this tutorial are already input in the [alignments] and [profiles] sections of the config file. These sections must be edited for different experiments. NOTE: In the tutorial, ShapeMapper encountered an error with failed quality filtering; for now, do not input a config file.

The ExampleResults archive contains the expected output of running the analysis on the example sequencing data (a subset of Ann Kladwang's pilot MaP experiment) in the Tutorial folder.