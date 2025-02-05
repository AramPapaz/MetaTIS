# MetaTIS: A tool to predict eukaryotic translation initiation sites.

![Image](https://github.com/user-attachments/assets/d96ecb2c-ce25-4d22-a6fa-4e0cc4b8dece)

When scanning mRNA sequences, ribosomes start translation typically at an AUG codon flanked by a so-called Kozak region. However, the first AUG of an mRNA is not always effective in initiating translation. In rare cases, also near cognate codon sequences may be recognized as start sites. The ribosome profiling technique, where translation elongation is stalled by chemical agents, is able to identify actively used translation start sites. In the past, various bioinformatics classifiers have been trained on such data to predict putative translation initiation sites from mRNA sequence features. Here, we formulated a stacking approach that can differentiate between false and true translation initiation sites. This was trained on experimental data for translation initiation in HEK293 cells produced by the so-called TISCA protocol. Our classifier gave a good overall performance on its own test set (accuracy 0.93) as well as multiple external validation sets. Moreover, it was able to predict almost quantitatively whether overlapping open-reading frames suppress translation from the main ORF for 11 genes in HeLa cells as validated by experimental luciferase assays. The MetaTIS tool is publicly available as a [webserver](https://service.bioinformatik.uni-saarland.de/metatis/)

## Data & Models
Some of the training data which include the positive, downstream negative, and upstream negative samples are stored on [zenodo](https://doi.org/10.5281/zenodo.14809153). Additionally, you can find the KmersERF and FlanksERF models there. The models are composed of 40 random forest classifiers which are stored as a dictionary. scikit-learn version 1.5.1 was used to create these models.

Here we provide a usage script which illustrates how the MetaTIS webserver functions. The ATGEff and NearCognateEff text files are the Noderer et al. translation efficiency values. The RBps pickle file represents the 9 position weight matrices of the RNA binding proteins being used. The Flanks_cols, Kmers_cols, and meta_cols text files contain the features utilized by FlanksERF, KmersERF, and MetaTIS respectively. 






