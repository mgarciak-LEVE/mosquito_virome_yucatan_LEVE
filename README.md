# Characterization of *Aedes taeniorhynchus* and *Aedes serratus* viromes from Yucatán.

## Introduction.

Of all arboviruses that cause human diseases, approximately 90% are transmitted by mosquitoes (Öhlund et al., 2019). However, most of the mosquito RNA virus repertoire can only replicate within the mosquito and does not maintain viral cycles between a vector and a vertebrate animal; these are therefore known as insect-specific viruses (ISVs). Consequently, although not all viruses naturally harbored by different mosquito species are necessarily pathogens and do not represent a public health problem per se, they can influence pathogen maintenance, modulate transmission, vector competence, and affect viral emergence (Lewis et al., 2023; Liu et al., 2023; Olmo et al., 2023).

The concept of a core microbiome aligns with the identification of common factors shared among microbiomes of various individuals, implying a broad range of shared characteristics among individuals in a population (Neu et al., 2021). The viruses that compose this microbiome have a commensal role, characterized by the absence of clinical symptoms in the host (Haynes & Rohwer, 2010). Thus, the term core virome was specifically coined to describe a set of viruses shared among the majority of individuals in a mosquito population (Shi et al., 2019). Furthermore, the virus-mosquito relationship is determined by extrinsic and intrinsic factors; these conditions interact and contribute to virome dynamics. For example, metagenomic studies reveal that they exhibit: (1) seasonal variations with greater viral diversity and abundance in warm months (Liu et al., 2023), (2) differences between urban and wild habitats, associated with changes in host availability and ecological conditions (Wang et al., 2024), (3) different abundance and composition at the viral family level depending on vector distributions (Pettersson et al., 2019), and (4) geographical variability related to the vector's taxon, where viromes are generally specific at the genus level and those occupying a larger geographical range tend to harbor more diverse viromes (Shi et al., 2017).

In recent years, the development of experimental techniques and changes in methodological approaches, such as the popularization of metagenomic studies to analyze mosquitoes, has led to an increase in the discovery of viral species (Rosenberg et al., 2013). For example, metatranscriptomic studies have identified more than 40 new viruses in wild mosquitoes of the genera Aedes, Culex, and Culiseta (Batson et al., 2021). However, the mere detection of viral sequences at the total RNA level does not allow discrimination between active infections, residual viral material, or Endogenous Viral Elements (EVEs). In mosquitoes, both ISV and arbovirus infections activate the RNA interference (RNAi) pathway, which are classified according to their size, biogenesis, and associated proteins into: microRNAs (~20-22 nt), small interfering RNAs (siRNAs, ~20-22 nt), and PIWI-interacting RNAs (piRNAs, ~25-30 nt) (Tikhe & Dimopoulos, 2021; Trammell & Goodman, 2021). In this context, non-retroviral EVEs acquire particular relevance, as in mosquitoes of the genus Aedes, they are enriched in piRNA clusters and are transcribed generating antisense piRNAs against cognate viruses. In the case of siRNAs, these are consistently found covering viral genomes, allowing the reconstruction of viral sequences that infect mosquitoes (Palatini et al., 2022).

Collectively, this project aims to characterize viromes based on total RNA and small RNA data in mosquitoes belonging to the species Ae. serratus and Ae. taeniorhynchus in order to identify possible patterns and biogeographic, ecological, and/or taxonomic factors that shape their composition and dynamics.
Objectives

## General Objective.
Characterize viral diversity in Ae. serratus and Ae. taeniorhynchus mosquitoes from the Mexican Republic collected in the Yucatán Peninsula and contrasting environments (conserved, diversified, and urban regions) using total RNA and small RNAs.

### Specific Objectives.

- Characterize the virome of mosquito species collected in the Yucatán Peninsula from 4 habitat types and perform genetic and virological characterization of detected viruses.
- Identify and analyze possible novel viral sequences.
- Determine if there are possible patterns and biogeographic, ecological, and/or taxonomic factors associated with the composition of their viromes.

## Workflow.

The proposed workflow for this project consists of several stages: sequence cleaning to ensure adequate quality...




## Prerequisites.
- Quality control: FastQC (Andrews, 2010) and MultiQC (Ewels et al., 2016).
- Sequence trimming: Trimmomatic (Bolger et al., 2014).
- Generation of super-reference from different mosquito genomes and alignment for host read depletion: STAR (Dobin et al., 2013).
- Assembly of reads that did not align to the host: MEGAHIT (Li et al., 2015) and SPAdes (Prijbelski et al., 2020).
- Annotation of sequences and viral proteins: BLAST (Altschul et al., 1990) and Diamond (Buchfink et al., 2015).

## Directory Structure
The code follows the assumption that the directory follows the structure outlined in the project. If you modify the directory structure, make sure to update the paths in the code accordingly.

```
mosquito_virome_yucatan_LEVE/
│
├── data/
│   ├── raw/
│   |   ├── total_RNA/
│   |   │   └── [samples]
│   |   ├── small_RNA/
│   |   │   └── [samples]/
│   |   └── metadata/
│
├── references/
│   ├── mosquito_genomes/
│   │   └── aedes_super_index/
│   └── databases/
│       ├── BLAST/
│       └── DIAMOND/
│
├── results/
│   ├── untrimmed_qc/
│   │   └── fastqc/
│   │   └── multiqc/
│   ├── trimmed/
│   ├── trimmed_qc/
│   │   └── fastqc/
│   │   └── multiqc/
│   ├── aligned/
│   └── assembly/
│       ├── statistics/
│       ├── rnaSPAdes/
│       │   └── [samples]/
│       ├── metaSPAdes/
│       │   └── [samples]/
│       ├── MEGAhit/
│       │   └── [samples]/
│
├── logs/
│   ├── trimming/
│   ├── mapping/
│   ├── assembly/
│   └── blast/
│
├── docs/
│   └── aedes_genomes_specs/
│
└── scripts/
    ├── aedes_reference_genomes/
    ├── databases/
    ├── pipeline_whole/
    └── individual_analyses/
```


## References
1. Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Retrieved from http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
2. Batson, J., Dudas, G., Haas-Stapleton, E., Kistler, A. L., Li, L. M., Logan, P., Ratnasiri, K., & Retallack, H. (2021). Single mosquito metatranscriptomics identifies vectors, emerging pathogens and reservoirs in one assay. ELife, 10, e68353. https://doi.org/10.7554/eLife.68353
3.  Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
4.  Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635
5.  Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
6.  Genome [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; 2004 – [cited 2025/08/25]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/
7.  Haynes, M., & Rohwer, F. (2010). The Human Virome. In Metagenomics of the Human Body (pp. 63–77). Springer. https://doi.org/10.1007/978-1-4419-7089-3_4
8.  Lewis, J., Gallichotte, E. N., Randall, J., Glass, A., Foy, B. D., Ebel, G. D., & Kading, R. C. (2023). Intrinsic factors driving mosquito vector competence and viral evolution: a review. Frontiers in Cellular and Infection Microbiology, 13. https://doi.org/10.3389/fcimb.2023.1330600
9.  Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033
10. Liu, Q., Cui, F., Liu, X., Fu, Y., Fang, W., Kang, X., Lu, H., Li, S., Liu, B., Guo, W., Xia, Q., Kang, L., & Jiang, F. (2023). Association of virome dynamics with mosquito species and environmental factors. Microbiome, 11(1). https://doi.org/10.1186/s40168-023-01556-4
11. Neu, A. T., Allen, E. E., & Roy, K. (2021). Defining and quantifying the core microbiome: Challenges and prospects. Proceedings of the National Academy of Sciences, 118(51). https://doi.org/10.1073/pnas.2104429118
12. Öhlund, P., Lundén, H., & Blomström, A.-L. (2019). Insect-specific virus evolution and potential effects on vector competence. Virus Genes, 55(2), 127–137. https://doi.org/10.1007/s11262-018-01629-9
13. Olmo, R. P., Todjro, Y. M. H., Aguiar, E. R. G. R., de Almeida, J. P. P., Ferreira, F. V., Armache, J. N., de Faria, I. J. S., Ferreira, A. G. A., Amadou, S. C. G., Silva, A. T. S., de Souza, K. P. R., Vilela, A. P. P., Babarit, A., Tan, C. H., Diallo, M., Gaye, A., Paupy, C., Obame-Nkoghe, J., Visser, T. M., & Koenraadt, C. J. M. (2023). Mosquito vector competence for dengue is modulated by insect-specific viruses. Nature Microbiology, 8(1), 135–149. https://doi.org/10.1038/s41564-022-01289-4
14. Palatini, U., Alfano, N., Carballar-Lejarazu, R., Chen, X.-G., Delatte, H., & Bonizzoni, M. (2022). Virome and nrEVEome diversity of Aedes albopictus mosquitoes from La Reunion Island and China. Virology Journal, 19(1). https://doi.org/10.1186/s12985-022-01918-8
15. Pettersson, J. H.-O. ., Shi, M., Eden, J.-S., Holmes, E. C., & Hesson, J. C. (2019). Meta-Transcriptomic Comparison of the RNA Viromes of the Mosquito Vectors Culex pipiens and Culex torrentium in Northern Europe. Viruses, 11(11), 1033. https://doi.org/10.3390/v11111033
16. Prijbelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes de novo assembler. Current Protocols in Bioinformatics, 70(1), e102. https://doi.org/10.1002/cpbi.102
17. Shi, C., Beller, L., Deboutte, W., Yinda, K. C., Delang, L., Vega-Rúa, A., Failloux, A.-B., & Matthijnssens, J. (2019). Stable distinct core eukaryotic viromes in different mosquito species from Guadeloupe, using single mosquito viral metagenomics. Microbiome, 7(1). https://doi.org/10.1186/s40168-019-0734-2
18. Shi, M., Neville, P., Nicholson, J., Eden, J.-S., Imrie, A., & Holmes, E. C. (2017). High-Resolution Metatranscriptomics Reveals the Ecological Dynamics of Mosquito-Associated RNA Viruses in Western Australia. Journal of Virology, 91(17). https://doi.org/10.1128/jvi.00680-17
19. Tikhe, C. V., & Dimopoulos, G. (2021). Mosquito antiviral immune pathways. Developmental & Comparative Immunology, 116, 103964. https://doi.org/10.1016/j.dci.2020.103964
20. Trammell, C. E., & Goodman, A. G. (2021). Host Factors That Control Mosquito-Borne Viral Infections in Humans and Their Vector. Viruses, 13(5), 748. https://doi.org/10.3390/v13050748
21. World Health Organization. (2017). Global vector control response 2017–2030. https://www.who.int/publications/i/item/9789241512978