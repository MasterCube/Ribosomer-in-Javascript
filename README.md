PhD thesis in computational biology.
************************************
This repository is one of the two legacies of my PhD thesis. The PhD thesis has two deliverables that are open to 
the scientific community studying translation by ribosomes. The second deliverable attached to this GitHub repository 
is an instance of an agent-based model (ABM) of protein synthesis by a pool of ribosomes acting on a set of transcripts.
This second instance was developed in Javascript and used 'p5.js processing' to increase html interactivity and 
enforce dynamical displays on html pages.

To have an overall introduction to the functions and motivations of this simulation tool, the reader is 
invited to read the PhD thesis attached in this repository and entitled 'Protein synthesis by ribosomes, an 
agent-based model of mRNA translation rate incorporating tRNA modifications effects'.

This repository contains all the directories and files that are necessary to conduct protein synthesis simulation by a pool 
of ribosomes on any transcriptome set that the user can provide as data input.

This instance of the agent-based model is meant to be used as a didactical tool for training purposes. The idea is to get familiar
with the protein elongation dynamics and see what are the effects of a change in factors affecting protein elongation by ribosomes.

MyData directory:
----------------
The directory 'MyData' contains files that are used as input of the program simulating protein synthesis by ribosomes.

  # transcriptome data:
  --------------------
  The transcriptome data is a fasta file. The file CDSfasta01.txt contains a list of arbitrarily chosen transcripts (mRNA) 
  sequences with their metadata that can be provided by the user. This file is used as input data in the 'Ribosomer' agent-based model to provide the 
  codons sequences of each transcript that will be used for the protein synthesis simulations.

  # dataJSONyeast.json:
  ---------------------
  This .json file provides the kinetic parameters that are used throughout a simulation for the elongation process 
  of protein synthesis by ribosomes at codon resolution. An elongation cycle of a ribosome on a given codon has 
  3 substeps: 
  1) accommodation of a loaded-tRNA and proofreading (of anticodon to codon) at the A-site.
  2) peptide bond formation of the amino-acid at the A-site with the carboxy-terminal amino-acid at the P-site.
  3) transocation of the ribosome to the next codon.

  Each of this sub-step can be approximated by a first-order kinetics law. The resulting queueing time of a ribosome 
  on a codon is equivalent to a global queueing time that is the convolution product of three exponentially distributed
  queueing times -one for each of the three sub-steps occuring sequentially. The three sub-steps complete a full 
  elongation cycle of the ribosome on a given codon [Joiret, M., F. Kerff, F. Rapino, P. Close, and L. Geris (2023a).
  “A simple geometrical model of the electrostatic environment around the catalytic center of the ribosome and 
  its significance for the elongation cycle kinetics”. In: Computational and Structural Biotechnology Journal 21, 
  pp. 3768–3795. doi: 10.1016/j.csbj.2023.07.016].

  Each sub-step has a rate that is inherently dependent on the nature of the codon. The 'dataJSONyeast.json' file is
  a dictionnary with 61 keys - the 61 sense codons, and 3 rates per key - one rate for each of the elongation sub-step.
  The calibration of this dictionary is species dependent and this one is for yeast - Saccharomyces cerevisiae. 
  The rate values for each key were derived from published metadata analysis [Dana, A. (Oct. 2014). 
  “Properties and determinants of codon decoding time distributions”. In: BMC Genomics 15, S13].
  [Dana, A. and T. Tuller (Nov. 2012). “Determinants of Translation Elongation Speed and Ribosomal Profiling Biases 
  in Mouse Embryonic Stem Cells”. In: PLoS computational biology 8, e1002755].
  [Dana, A. and T. Tuller (2014). “The effect of tRNA levels on decoding times of mRNA codons”. 
  In: Nucleic Acids Research 42, pp. 9171 –9181].

  # tunnelElectrostaticsAF.json:
  -----------------------------
  This .json file provides the profile of the axial forces applied in the charged amino acids of the peptide 
  nascent chain embedded in the ribosome exit tunnel. The profile was calibrated from the electrostatic 
  potential profile modeled in [Joiret, M., F. Kerff, F. Rapino, P. Close, and L. Geris (2022). 
  “Ribosome exit tunnel electrostatics”. In: Physical Review E 105.1, p. 014409. doi: 10.1103/PhysRevE. 105.014409].

MyInput directory:
-----------------
The directory 'MyInput' contains instances of input files that are used to provide information on 
RNA-seq relative abundance and on the relative fold-change in initiation rates of each transcript 
in the transcriptomic set of interest. These input pieces of information are used throughout the simulations runs. 
Different input files can be used to compare the effects of a change in transcripts relative abundance or 
the effects of a change in initiation rates while keeping the transcripts abundance constant.
