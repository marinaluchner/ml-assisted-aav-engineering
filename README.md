# Machine learning-assisted directed evolution yields dramatic improvement on novel AAV engineering task
---
<p align="center">
  <img src="screening_workflow.png" width="600">
</p>

This repository provides source code to reproduce, verify and extend the study "Machine learning-assisted directed evolution yields dramatic improvement on novel AAV engineering task". Source data files necessary to reproduce the manuscript figures are provided on NCBI SRA (https://www.ncbi.nlm.nih.gov/bioproject/1473000).

## Abstract

Directed evolution enables the discovery of protein mutants with improved fitness through iterative rounds of selection and has been widely applied to adeno-associated virus (AAV) capsid engineering. Machine learning (ML) can augment this process by prioritising mutants for experimental validation, but whether its benefits outweigh the added cost of ML-designed library construction remains unclear. Here, we address this question by considering improvements in manufacturing efficiency of AAV capsids in the context of exosomal encapsulation, which is a technique for improving AAV immunogenicity. We generated a directed evolution dataset comprising 53,974 mutants across three rounds of selection and an independent assessment dataset of 472 ML-designed mutants with detailed profiling. Directed evolution alone yielded a 3-fold improvement, while ML-assisted directed evolution yielded a 41-fold improvement in manufacturing efficiency, demonstrating that ML-assisted directed evolution outperforms directed evolution alone. We show that data quality – specifically the number of selection rounds and the statistical power of mutant counts – has a greater impact on model performance than model architecture. Finally, we provide practical guidance on experimental design and implementation of ML tuning strategies that can augment model performance in protein engineering applications.  

---

## Requirements
This software was tested under Windows Subsystem for Linux 2 (WSL2; Linux kernel 6.6.114.1) on a workstation equipped with a 13th Generation Intel Core i7-13700H processor and 16GB RAM. No GPU acceleration was used.

The requirements to run the data processing and associated analysis are stored in data_processing.yaml. Meanwhile, the requirements for model training, analysis and visualization can be found in the modelling.yaml. Conda environment can be created from these files using the commands below.
```bash
conda env create -f data_processing.yaml
conda env create -f modelling.yaml
```
Installing the data processing environment from the yaml file can take several minutes, while the modelling environment can take up to an hour.

## Example data processing

<p align="center">
  <img src="data_processing.png" width="700">
</p>

Example input to run NGS data processing and calculate the log enrichment score S are included in example/raw_NGS_files. Small test files (15,145 reads, example/raw_NGS_files/sample_01/example_L01_R1.fq) can be run from the root directory via the following commands:
```bash
conda activate data_prep_env
bash 1_combine_FASTQ_files.bash example/raw_NGS_files example/outputs
python 2_align_and_merge_paired_end_reads.py -i example/outputs -o example/outputs/sample_01
bash 3_filter_and_trim_reads.bash example/outputs example/outputs/sample_01
python 4_calculate_mutant_counts.py -i example/outputs -o example/outputs
python 5_calculate_log_enrichment_score.py -pre_encapsulation_i example/outputs/sample_01/sample_01_example_variant_count.xlsx -post_encapsulation_i example/outputs/sample_01/sample_01_example_variant_count.xlsx -o example/outputs/sample_01
python 6_combine_datasets_across_replicates.py -i example/outputs/sample_01/enrichment_score_sample_01_example_variant_count_sample_01_example_variant_count.xlsx example/outputs/sample_01/enrichment_score_sample_01_example_variant_count_sample_01_example_variant_count.xlsx -o example/outputs/sample_01
python 7_translate_nucleotide_sequence.py -i example/outputs/sample_01 -o example/outputs/sample_01
```
Alignment and merging are the most time-consuming steps and take ~5 min to run with 4 threads. The expected outputs are included in example/outputs/enrichment_scores.csv.

## Example model training

Example configuration and input files to train and compare sequence-to-function models are included in example/sequence-to-function/. Model training and evaluation can be tested by running the following command:

```bash
conda activate modelling_env
python make_figure_2d.py -i example/sequence-to-function/enrichment_score_with_amino_acid_sequence_threshold.xlsx -o example/sequence-to-function
```
Training and evaluating the models takes ~15 min without GPU usage. The expected output will be stored in example/sequence-to-function/figure_2d.png.

## Full reproduction instructions

For full reproduction of results, download the source data (https://www.ncbi.nlm.nih.gov/bioproject/1473000) and repeat the steps outlined under example data processing. Next, run the visualization Python scripts and Jupyter notebooks using the resulting merged_enrichment_score_with_amino_acid_sequences.xlsx files as input. 

## Reproduction with pre-processed files

Since raw NGS files encompass 41.9 GB of compressed data, downloading and processing can be time-intensive. Hence, we provide input data in a more pre-processed manner to reproduce, verify and extend the findings of our study. Follow the commands below and find the output graphs under plots/.

```bash
conda activate modelling_env
jupyter execute make_figure_2a_and_supplementary_figure_S1a.ipynb
jupyter execute make_figure_2b_and_supplementary_figure_S1b.ipynb
python make_figure_2c.py
python make_figure_2d.py
conda activate data_prep_env
jupyter execute make_figure_3a_and_supplementary_figure_S6.ipynb
jupyter execute make_figure_3b.ipynb
conda activate modelling_env
jupyter execute make_figure_4a.ipynb
jupyter execute make_figure_4b.ipynb
jupyter execute make_figure_4c.ipynb
jupyter execute make_figure_4d.ipynb
jupyter execute make_supplementary_figure_S2a_and_S3a.ipynb
jupyter execute make_supplementary_figure_S2b_and_S3b.ipynb
jupyter execute make_supplementary_figure_S4.ipynb
conda activate data_pred_env
python make_supplementary_figure_S7.py
```

## Contact

Marina Luchner, Department of Engineering Science, University of Oxford

For questions regarding the code, data or analyses presented in this repository, please contact:
[marina.luchner@eng.ox.ac.uk](mailto:marina.luchner@eng.ox.ac.uk)

Bug reports and feature requests are welcome through the GitHub Issues page.
