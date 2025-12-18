# Package index

## All functions

- [`BuildMDP()`](https://bowmanr.github.io/scDNA/reference/BuildMDP.md)
  : The main function builds an adjacency list of the theoretical order
  of mutations. The same MDP mask can be reused different mutations
- [`GetMembership()`](https://bowmanr.github.io/scDNA/reference/GetMembership.md)
  : Get membership quality of cells vs empty
- [`GetStainIndex()`](https://bowmanr.github.io/scDNA/reference/GetStainIndex.md)
  : Get Stain Index for cells
- [`add_cell_annotation()`](https://bowmanr.github.io/scDNA/reference/add_cell_annotation.md)
  : Title
- [`annotate_variants()`](https://bowmanr.github.io/scDNA/reference/annotate_variants.md)
  : Annotate variants of interest This function takes in h5 files to
  extract DNA information through given TXDB files (primarily hg19) The
  gene names are mapped for their specific nucleotide positions. The
  starting and ending positions for the chromosome are listed. The
  consequence of the mutation such as synonymous, nonsynonymous, as well
  as if this is a coding region, on the exon boundry or intronic. Short
  amino acid changes are also labeled. The variance matrix also extracts
  the amplicon that the variant is found on.
- [`attach_weights()`](https://bowmanr.github.io/scDNA/reference/attach_weights.md)
  : Assigns observed counts as weights to the Markov Decision Process
- [`cell_confidence_labeling()`](https://bowmanr.github.io/scDNA/reference/cell_confidence_labeling.md)
  : cell_confidence_labeling This function provides additional labeling
  to he quality of the cells based on dna and protein data.
- [`clonograph()`](https://bowmanr.github.io/scDNA/reference/clonograph.md)
  : Plotting clonographs
- [`compare_VAFs()`](https://bowmanr.github.io/scDNA/reference/compare_VAFs.md)
  : Title
- [`compute_clone_statistics()`](https://bowmanr.github.io/scDNA/reference/compute_clone_statistics.md)
  : Title
- [`demultiplex_samples()`](https://bowmanr.github.io/scDNA/reference/demultiplex_samples.md)
  : Demultiplex wrapper function to handle splitting samples
- [`enumerate_clones()`](https://bowmanr.github.io/scDNA/reference/enumerate_clones.md)
  : Enumerate clones
- [`extract_droplet_size()`](https://bowmanr.github.io/scDNA/reference/extract_droplet_size.md)
  : Extract protein library size, dna library size, and amplicon size
  for all droplets.
- [`fcs_export()`](https://bowmanr.github.io/scDNA/reference/fcs_export.md)
  : Title
- [`gemerate_txdb()`](https://bowmanr.github.io/scDNA/reference/gemerate_txdb.md)
  : Title
- [`get_own_path()`](https://bowmanr.github.io/scDNA/reference/get_own_path.md)
  : A way to navigate the MDP for any given starting root node to leaf
  node.
- [`impute_cluster()`](https://bowmanr.github.io/scDNA/reference/impute_cluster.md)
  : Title
- [`loom_to_sce()`](https://bowmanr.github.io/scDNA/reference/loom_to_sce.md)
  : Title
- [`match_clonal_graph()`](https://bowmanr.github.io/scDNA/reference/match_clonal_graph.md)
  : Title
- [`mdp_Q_learning_with_linklist()`](https://bowmanr.github.io/scDNA/reference/mdp_Q_learning_with_linklist.md)
  : This file run Reinforcement Learning (model-free Q-learning) to
  evaluate the most likely mutation paths from the data
- [`normalize_protein_data()`](https://bowmanr.github.io/scDNA/reference/normalize_protein_data.md)
  : Normalize Protein Data
- [`optimize_matrix()`](https://bowmanr.github.io/scDNA/reference/optimize_matrix.md)
  : Title
- [`quality_output()`](https://bowmanr.github.io/scDNA/reference/quality_output.md)
  : Produce long form quality metrics for each cell-variant pair for
  read depth, allele frequency, and genotype quality
- [`readDNA_CN_H5()`](https://bowmanr.github.io/scDNA/reference/readDNA_CN_H5.md)
  : This function generates the Copy Number by determining the ploidy of
  each mutation
- [`select_clones()`](https://bowmanr.github.io/scDNA/reference/select_clones.md)
  : Select clones of interest on the basis QC metrics
- [`tabulate_mutations()`](https://bowmanr.github.io/scDNA/reference/tabulate_mutations.md)
  : Title
- [`tapestri_h5_to_sce()`](https://bowmanr.github.io/scDNA/reference/tapestri_h5_to_sce.md)
  : Import Tapestri H5 data and extract genotype matrix
- [`trajectory_analysis()`](https://bowmanr.github.io/scDNA/reference/trajectory_analysis.md)
  : Run Trajectory Analysis after extraction from SingleCellExperiment
  object
- [`trajectory_of_interest_BSCITE_format()`](https://bowmanr.github.io/scDNA/reference/trajectory_of_interest_BSCITE_format.md)
  : BSCITE Style single trajectory Plot
- [`trajectory_of_interest_figure()`](https://bowmanr.github.io/scDNA/reference/trajectory_of_interest_figure.md)
  : Single Trajectory visualization of interest.
- [`variant_ID()`](https://bowmanr.github.io/scDNA/reference/variant_ID.md)
  : Variant identification and frequency tallies The DNA variants for
  each cell are pulled from the specified H5 files. Each variant for
  each cell is genotyped to be WildType, Heterozygous, Homozygous or
  Missing. The genotyping rate is determined by taking WT+Het+Hom over
  total cell calls (including missing). The VAF is determined by number
  of allele copies we see in a weighted sum. A filter is applied to both
  of these calculatins to include or exclude variants of interest. These
  are then annotated to include the variant information such as gene
  name, nucleotide location, and short amino acid changes.
- [`visualize_WT_dominant_clone()`](https://bowmanr.github.io/scDNA/reference/visualize_WT_dominant_clone.md)
  : Title
- [`visualize_all_WT_dominant_clone()`](https://bowmanr.github.io/scDNA/reference/visualize_all_WT_dominant_clone.md)
  : Title
- [`visualize_any_optimal_path()`](https://bowmanr.github.io/scDNA/reference/visualize_any_optimal_path.md)
  : Title
- [`visualize_full_network()`](https://bowmanr.github.io/scDNA/reference/visualize_full_network.md)
  : Visualize Full Network
