# MetaTIS: A tool to predict eukaryotic translation initiation sites.

![Image](https://github.com/user-attachments/assets/d96ecb2c-ce25-4d22-a6fa-4e0cc4b8dece)

As the ribosome makes its foray unto a sequence of mRNA, ribosomes typically commence translation at a methionine-encoding AUG codon flanked by a so-called Kozak region, a short nucleic acid motif serving as an initiation site in many eukaryotes. Though, the characteristic AUG start codon of an mRNA is not always effective in initiating translation. In more seldom cases, near-cognate codon sequences may also be recognized as start sites. Ribosome profiling techniques, characterized by the stymieing of mRNA-ribosome complex translation function via chemical treatment, are able to elucidate active translation start sites. Historically, several bioinformatics classifiers have been trained on translation start site data, gleaned from ribosomal profiling, to predict putative translation initiation sites from mRNA sequence features. A stacking approach was formulated for the MetaTIS tool that can differentiate spurious and true translation initiation sites. The tool was trained on experimental data for translation initiation in HEK293 cells produced by the TISCA protocol, a method allowing for accurate translation initiation site identification. Our classifier delivers a notable ROC-AUC of 0.93 while performing on its own test set, as well as multiple external validation sets. Moreover, it was able to almost quantitively predict whether overlapping open-reading frames suppress translation from the main ORF for 11 genes in HeLa cells, as validated by experimental luciferase assays. The MetaTIS tool is publicly available as a [webserver](https://service.bioinformatik.uni-saarland.de/metatis/)

## Data & Models
Some of the training data which include the positive, downstream negative, and upstream negative samples are stored on [zenodo](https://doi.org/10.5281/zenodo.14809153). Additionally, you can find the KmersERF and FlanksERF models there. The models are composed of 40 random forest classifiers which are stored as a dictionary. scikit-learn version 1.5.1 was used to create these models.

Here we provide a usage script which illustrates how the MetaTIS webserver functions. The ATGEff and NearCognateEff text files are the Noderer et al. translation efficiency values. The RBps pickle file represents the 9 position weight matrices of the RNA binding proteins being used. The Flanks_cols, Kmers_cols, and meta_cols text files contain the features utilized by FlanksERF, KmersERF, and MetaTIS respectively. 






