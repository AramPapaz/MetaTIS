# MetaTIS
Ribosomes typically commence translation at a methionine-encoding AUG codon flanked by a so-called Kozak region, a short nucleic acid motif that serves as an initiation site in humans. Though, the characteristic AUG start codon of an mRNA is not always effective in initiating translation. Seldomly, near-cognate codon sequences may also be recognized as start sites. Several types of ribosomal profiling techniques have been developed that elucidate active translation initiation sites (TIS). Based on data gathered from these techniques, machine learning models have been developed to predict translation start sites using mRNA sequence features. Here, a meta-model termed MetaTIS was implemented by combining outputs of genomic and protein language models fine-tuned on five different TIS datasets plus the Ensembl annotations. The model proficiently differentiates between spurious and true TIS in four distinct test sets, for both canonical and noncanonical instances. While analysing one of the base models with integrated gradients, it was found that the model considered the importance of the Kozak sequence context and the presence of an upstream open reading frame (uORF) for classifying a position as an actual TIS. We further demonstrated how MetaTIS predictions can be used to detect important uORFs in three cancer types. 


## Python and package versions
*python* ==  3.12.4  
*scikit-learn* == 1.7.2  
*pandas* ==2.2.2  
*numpy* == 2.0.0  
*torch* == 2.8.0  
*transformers* == 4.55.3  
*peft* == 0.17.0  
