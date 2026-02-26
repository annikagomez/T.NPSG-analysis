# T.NPSG-analysis
Code and data used to perform analysis and generate figures in "A globally distributed phage-plasmid infects marine N2-fixing cyanobacterium Trichodesmium". 

## Code files
- `fig_1_genome_vizualization.ipynb`: Code to process phage-plasmid genome annotations, create files needed for NCBI GenBank deposit, and generate Figure 1a and Supplementary Figures 2 and 3.
- `fig_1_orthogroup_vizualization.ipynb`: Code to process OrthoFinder output and generate Figure 1b.
- `fig_2_metaT.ipynb`: Code to process metatranscriptome mapping, assembly gene annotations, and DESeq results and generate Figures 2a and 2b and Supplementary Figures 4 and 5.
- `fig_3_methylation.ipynb`: Code to process MicrobeMod output and generate Figures 3a and 3b and Supplementary Figure 6.
- `fig_4_metagenome_mapping.ipynb`: Code to process mapping results from global metagenomes and generate Figure 4a and Supplementary Figures 7 and 8.
- `fig_4_model_output.ipynb`: Code to generate Figure 4b. 
- `deseq.R`: R script to perform DESeq analysis.

## Source data
- annotation_files/
  - `unified_annotations_*.csv`: Manually curated annotations of T.NPSG.1 (col) and T.NPSG.2 (slick)
  - `*_genes.bed`: BED files for both genomes
  - `*_unified.faa`: Fasta file with predicted protein sequences in both genomes
  - `feature_table_*.txt`: Feature annotation table geneated for NCBI GenBank deposit
  - `*_phold.gbk`: Genbank file produced by Pharokka/Phold pipeine
  - `tricho_phage_NPSG_*.gbk`: Genbank file updated with manual edits to annotations; source data for Figure 1a and Supplementary Figure 2.
  - `col_v_slick_phage_blastn.out`: T.NPSG.1 vs T.NPSG.2 blastn output; source data for Supplementary Figure 3
- mapping_files/
  - `FYP_Tricho_Metadata.csv`: Information on all publicly available picked colony metagenomes
  - `tricho_abundance_tara.csv`: Abundance of Trichodesmium MAGs in all Tara Oceans metagenomes (From Delmont, et al. 2021)
  - `tara_metadata_*.csv`: Metadata for stations, samples, and ENA from Tara Oceans samples
- metaT_files/
  - `*.tab`: Outputs from BWA mapping metatranscriptome reads to hybrid slick assembly + T.NPSG.2 genome
  - `all_sig_tery.csv`: All Trichodesmium erythraeum genes which were differentially expressed between the slick and colonies based on DESeq2 analysis with eggnog-mapper annotations; Source data for Supplementary Figure 5
  - `bacterial_phylum_relative_expression.csv`: Sum of TPMs for transcripts assigned to major bacterial phyla by GTDB
  - `cyanobacterial_orders_relative_expression.csv`: Sum of TPMs for transcripts assigned to major cyanobacterial orders by GTDB; source data for Supplementary Figure 4
  - `gtdb_taxonomy_result.tsv`: Output from MMSeqs2 annotation of hybrid slick assembly using GTDB; source data for Supplementary Figure 4
  - `hybrid_assembly_prots.bed`: Location of all protein coding region within the hybrid slick assembly
  - `phage_tpm.csv`: TPM by sample of each T.NPSG.2 CDS
  - `slick_metadata.csv`: Maps sample number to slick or colony for metaT samples
  - `t_e_deseq_resdf.csv`: DESeq2 output
  - `tery_raw_gene_counts.csv`: Raw number of reads mapped per coding region for Trichodesmium erythraeum contigs; input for DESeq analysis
  - `tery_tpm.csv`: TPM by sample of each Trichodesmium erythraeum CDS
- methylation_files/
  - `*_ont.cov`: Coverage of the T.NPSG.2 and Trichodesmium erythraeum genomes by nanopore read mapping
  - `slick_*_rm.rm.genes.tsv`: MicrobeMod annotation of RM-related genes in T.NPSG.2 and Trichodesmium erythraeum genome
  - `tery_IMS101_ref.gbk`: T. erythraeum IMS101 genome (NC_008312.1)
  - `*_methylated_sites.tsv`: Methylation calls - MicrobeMod output
  - `*_motifs.tsv`: Methylated motifs - MicrobeMod output
- model_output/
  - netcdf files produced by global ecosystem model  
- orthogroup_files/
  - `Orthogroups.tsv`: OrthoFinder output
  -  `genome_type.csv`: Key mapping genome name to category for visualization
- output_figures/
  - Vector files used to create main and supplementary figures
- output_tables/
  - Supplementary tables resulting from analysis/processing by the code in this repository
