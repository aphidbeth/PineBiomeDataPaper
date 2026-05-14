# PineBiomeDataPaper

Scripts and files associated with the PineBiome project data paper "Metabarcode and transcriptome datasets of Pinus sylvestris to assess fungal phyllosphere and disease dynamics". [![DOI](https://zenodo.org/badge/1237974155.svg)](https://doi.org/10.5281/zenodo.20179422)

Information on the background of the study can be found at https://www.hutton.ac.uk/project/pinebiome/.

# Folder layout

There are two main folders: 
1. Metabarcoding => scripts, plots and data related to the ITS metabarcoding of the tree needles.
2. RNA_data => scripts and data related to the transcriptomic libraries of the tree needles. 

Minimally processed OTU tables can be found in /Metabarcoding/Data/OTU_tables.
These are available as tsv files:
* [T1_ITS_min_clean_OTU_table.tsv](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/OTU_tables/T1_ITS_min_clean_OTU_table.tsv)
* [T2_ITS_min_clean_OTU_table.tsv](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/OTU_tables/T2_ITS_min_clean_OTU_table.tsv)
* [T3_ITS_min_clean_OTU_table.tsv](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/OTU_tables/T3_ITS_min_clean_OTU_table.tsv)

or can be imported directly into R from the rds objects: 

* [T1_ITS_min_clean_OTU.rds](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/OTU_tables/T1_ITS_min_clean_OTU.rds)
* [T2_ITS_min_clean_OTU.rds](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/OTU_tables/T2_ITS_min_clean_OTU.rds)
* [T3_ITS_min_clean_OTU.rds](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/OTU_tables/T3_ITS_min_clean_OTU.rds)

OTU are named by the feature hashes from the qiime pipeline. Relevant taxonomic information for these is included in the rds object or can be cross referenced with the "seqs-and-taxonomy" tables:
* [seqs-and-taxonomy_T1_ITS.tsv](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/Sequences/seqs-and-taxonomy_T1_ITS.tsv)
* [seqs-and-taxonomy_T2_ITS.tsv](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/Sequences/seqs-and-taxonomy_T2_ITS.tsv)
* [seqs-and-taxonomy_T3_ITS.tsv](https://github.com/aphidbeth/PineBiomeDataPaper/blob/main/Metabarcoding/Data/Sequences/seqs-and-taxonomy_T3_ITS.tsv)


## License

### Data
All proccessed datasets (including OTU tables) are released under CC BY 4.0.

### Code
All scripts are released under the MIT License.
