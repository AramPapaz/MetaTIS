# MetaTIS
Ribosomes typically commence translation at a methionine-encoding AUG codon flanked by a so-called Kozak region, a short nucleic acid motif that serves as an initiation site in humans. Though, the characteristic AUG start codon of an mRNA is not always effective in initiating translation. Seldomly, near-cognate codon sequences may also be recognized as start sites. Several types of ribosomal profiling techniques have been developed that elucidate active translation initiation sites (TIS). Based on data gathered from these techniques, machine learning models have been developed to predict translation start sites using mRNA sequence features. Here, a meta-model termed MetaTIS was implemented by combining outputs of genomic and protein language models fine-tuned on five different TIS datasets plus the Ensembl annotations. The model proficiently differentiates between spurious and true TIS in four distinct test sets, for both canonical and noncanonical instances. While analysing one of the base models with integrated gradients, it was found that the model considered the importance of the Kozak sequence context and the presence of an upstream open reading frame (uORF) for classifying a position as an actual TIS. We further demonstrated how MetaTIS predictions can be used to detect important uORFs in three cancer types. The tool is also available as a [webserver](https://service2.bioinformatik.uni-saarland.de/metatis/).

## Data & Models
The *Data* folder contains all the TIS reported from five distinct studies. In the *Seq* column, *S* represents the position of the start codons. The train, validation, and test datasets for NT and ESM2 are available on [Hugging Face](https://huggingface.co/datasets/Arampapaz/MetaTIS_Train_Val_Test).  
The *Models* folder contains the calibrated meta-model plus the NT and ESM2 optimal checkpoints having the lowest loss on the validation set. *usage.ipynb* shows how these models can be used to obtain TIS predictions.
## Python and package versions
*python* ==  3.12.4  
*scikit-learn* == 1.7.2  
*pandas* ==2.2.2  
*numpy* == 2.0.0  
*torch* == 2.8.0  
*transformers* == 4.55.3  
*peft* == 0.17.0  
