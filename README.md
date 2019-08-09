# A continually expanding collection of microbiome analysis tools

These notes are not intended to be comprehensive. They include notes about methods, packages and tools I am learing and/or would like to explore. For a comprehensive overview of the subject, consider [other bioinformatics resources](https://github.com/mdozmorov/Bioinformatics_notes) and [collections of links to various resources](https://github.com/mdozmorov/MDmisc_notes). Issues with suggestions and pull requests are welcome!

# Table of content

* [Pipelines](#pipelines)
* [Downstream analysis](#downstream-analysis)
  * [Integrative analysis](#integrative-analysis)
* [Taxonomy](#taxonomy)
* [Phylogenetics](#phylogenetics)
* [Differential analysis](#differential-analysis)
* [Data](#data)
* [Misc](#misc)

## Pipelines

- ATLAS - Three commands to start analysing your metagenome data. Documentation, https://metagenome-atlas.readthedocs.io/en/latest/, GitHub, https://github.com/metagenome-atlas/atlas

- `F1000_workflow` - Microbiome workflow. RSV instead of OTU. Data preprocessing from raw reads. DADA2 pipeline, ASV summary tables using RDP (Greengenes and SILVA are available), phylogenetic tree reconstruction (pangorn). phyloseq downstream analysis, from filtering to agglomeration, transformation, various ordination visualizations (from PCoA, DPCoA, rank PCA, to CCA), supervised learning, graph-based visualization and testing, multi-omics analyses. https://github.com/spholmes/F1000_workflow
    - Callahan, Ben J., Kris Sankaran, Julia A. Fukuyama, Paul J. McMurdie, and Susan P. Holmes. “Bioconductor Workflow for Microbiome Data Analysis: From Raw Reads to Community Analyses.” F1000Research 5 (2016): 1492. https://doi.org/10.12688/f1000research.8986.2.
    
- `bioBakery` - an environment for metagenomics analysis. VM running on Vagrantr/VirtualBox, Docker image, Google Cloud and Amazon Machine Image. Homebrew/Linuxbrew installation. AnADAMA2 controls the workflows. Wiki, https://bitbucket.org/biobakery/biobakery/wiki/Home, workflows and tutorials, http://huttenhower.sph.harvard.edu/biobakery_workflows
    - McIver, Lauren J, Galeb Abu-Ali, Eric A Franzosa, Randall Schwager, Xochitl C Morgan, Levi Waldron, Nicola Segata, and Curtis Huttenhower. “BioBakery: A Meta’omic Analysis Environment.” Edited by John Hancock. Bioinformatics 34, no. 7 (April 1, 2018): 1235–37. https://doi.org/10.1093/bioinformatics/btx754.

- `DADA2` - resolves sequencing errors and reconstructs sequences for finer-resolution clustering. Complete pipeline to process PI FASTQ into merged, denoised, chimera-free, error-corrected sample sequences. The error model quantifies the rate $\lambda_{ij}$ at which an amplicon read with sequence $i$ is produced from sample sequence $j$ as a function of sequence composition and quality, Poisson distribution. The NCBI RefSeq 16S rrna database (RefSeq) and the Genome Taxonomy Database (GTDB) are both now available to use with dada2's assignTaxonomy function! https://zenodo.org/record/2541239#.XEyoLc9Kjfa. DADA2 page: https://github.com/benjjneb/dada2. A DADA2 workflow for Big Data, https://benjjneb.github.io/dada2/bigdata.html
    - Callahan, Benjamin J., Paul J. McMurdie, Michael J. Rosen, Andrew W. Han, Amy Jo A. Johnson, and Susan P. Holmes. “DADA2: High-Resolution Sample Inference from Illumina Amplicon Data.” Nature Methods 13, no. 7 (2016): 581–83. https://doi.org/10.1038/nmeth.3869.

- `Deblur` - resolves Illumina sequencing errors and creates sub-operational taxonomic unit (sOTU) clusters. Operates on individual samples. Plugin for QIIME2 exists. Competing methods - DADA2, UNOISE2. Methods in the supplementary text S1. https://github.com/biocore/deblur
    - Amir, Amnon, Daniel McDonald, Jose A. Navas-Molina, Evguenia Kopylova, James T. Morton, Zhenjiang Zech Xu, Eric P. Kightley, et al. “Deblur Rapidly Resolves Single-Nucleotide Community Sequence Patterns.” Edited by Jack A. Gilbert. MSystems 2, no. 2 (April 25, 2017). https://doi.org/10.1128/mSystems.00191-16.

- `HiMAP` - high-resolution microbial analysis pipeline for 16S data analysis. Wraps many DADA2 functions. Comparison with DADA2, QIIME, detects more species. https://github.com/taolonglab/himap
    - Segota, Igor, and Tao Long. “A High-Resolution Pipeline for 16S-Sequencing Identifies Bacterial Strains in Human Microbiome.” BioRxiv, March 4, 2019. https://doi.org/10.1101/565572.

- `HUMAnN2`: The HMP Unified Metabolic Analysis Network 2 - functional profiling and pathway reconstruction of metagenomes. Tiered approach: 1) Screening for known species with MetaPhlAn2; 2) mapping against pangenomes; 3) mapping against protein sequences. These mappings can help to assign metabolic and functional annotations. http://huttenhower.sph.harvard.edu/humann2
    - Franzosa, Eric A., Lauren J. McIver, Gholamali Rahnavard, Luke R. Thompson, Melanie Schirmer, George Weingart, Karen Schwarzberg Lipson, et al. “Species-Level Functional Profiling of Metagenomes and Metatranscriptomes.” Nature Methods 15, no. 11 (November 2018): 962–68. https://doi.org/10.1038/s41592-018-0176-y.

- `microbial-rnaseq` - microbial composition from host's RNA-seq data, https://github.com/FredHutch/microbial-rnaseq

- `Microbiome Helper` - wrapper scripts and tutorials for metagenomics analysis. https://github.com/LangilleLab/microbiome_helper/wiki.
    - Comeau, André M., Gavin M. Douglas, and Morgan G. I. Langille. “Microbiome Helper: A Custom and Streamlined Workflow for Microbiome Research.” Edited by Jonathan Eisen. MSystems 2, no. 1 (February 28, 2017). https://doi.org/10.1128/mSystems.00127-16.

- `SqueezeMeta` - a pipeline for metagenomics/metatranscriptomics for co-assembly (SPAdes, Canu), gene and rRNA prediction (prodigal, RDP classifier), binning, gene abundance estimation, taxonomic annotation (fast LCA). Support for MinION nanopore sequencing data (long, error-prone reads). Table 1 - comparison with other pipelines. https://github.com/jtamames/SqueezeMeta

## Downstream analysis

- `Calour` - Heatmap-based visual exploration of microbiome data. Input - sOTU table (Deblur-processed) and phenodata. Normalization,sorting, filtering, interface with annotation databases, machine learning methods from scikit-learn. Python, Jupyter notebooks. http://biocore.github.io/calour/
    - Xu, Z.Z., Amir, A., Sanders, J., Zhu, Q., Morton, J.T., Bletz, M.C., Tripathi, A., Huang, S., McDonald, D., Jiang, L., et al. (2019). Calour: an Interactive, Microbe-Centric Analysis Tool. MSystems 4, e00269-18.

- `microbiome` R package with rich set of functions for microbiome analysis, visualization, statistical analysis. Supports phyloseq objects. Leo Lahti, Sudarshan Shetty et al. (2017). Tools for microbiome analysis in R. Version 1.5.23. URL: http://microbiome.github.com/microbiome

- `microbiomeSeq` - An R package for microbial community analysis in an environmental context. GitHub, https://github.com/umerijaz/microbiomeSeq, and tutorial, http://userweb.eng.gla.ac.uk/umer.ijaz/projects/microbiomeSeq_Tutorial.html

- `Metaviz` - visual exploratory data analysis of annotated microbiome data. Java/D3 implementation. Imports metagenomeSeq object, works with phyloseq objects. Web interface with 33 demo datasets, http://metaviz.cbcb.umd.edu/, `metavizr` R package, https://www.bioconductor.org/packages/release/bioc/html/metavizr.html. Docker https://epiviz.github.io/tutorials/metaviz/usingDocker/. GitHub, https://github.com/epiviz/metavizr. Documentation, https://epiviz.github.io/tutorials/metaviz/
    - Wagner, Justin, Florin Chelaru, Jayaram Kancherla, Joseph N Paulson, Alexander Zhang, Victor Felix, Anup Mahurkar, Niklas Elmqvist, and Héctor Corrada Bravo. “Metaviz: Interactive Statistical and Visual Analysis of Metagenomic Data.” Nucleic Acids Research 46, no. 6 (April 6, 2018): 2777–87. https://doi.org/10.1093/nar/gky136.
 
- `phyloseq` R package for import of the most OTU clustering data formats, preprocessing (normalization, standartization, subsampling, filtering), visualization (various definitions of distance, dimensionality reduction methods), and analysis (comparative) of microbiome data. phyloseq-class with four components (otu_table, sample_data, tax_table, phy_tree). Plotting functions using ggplot2 graphics. http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html, http://joey711.github.io/phyloseq/, https://github.com/joey711/phyloseq
    - McMurdie, Paul J., and Susan Holmes. “Phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data.” Edited by Michael Watson. PLoS ONE 8, no. 4 (April 22, 2013): e61217. https://doi.org/10.1371/journal.pone.0061217.

- Shiny-phyloseq app. http://joey711.github.io/shiny-phyloseq/
    - McMurdie, P. J., and S. Holmes. “Shiny-Phyloseq: Web Application for Interactive Microbiome Analysis with Provenance Tracking.” Bioinformatics 31, no. 2 (January 15, 2015): 282–83. https://doi.org/10.1093/bioinformatics/btu616.

- `WHAM!` - data exploration, clustering, visualization, differential expression analysis (ANOVA-like). Input format uses bioBakery pipeline output. Shiny app, https://ruggleslab.shinyapps.io/wham_v1/, GitHub, https://github.com/ruggleslab/jukebox/tree/master/wham_v1
    - Devlin, Joseph C., Thomas Battaglia, Martin J. Blaser, and Kelly V. Ruggles. “WHAM!: A Web-Based Visualization Suite for User-Defined Analysis of Metagenomic Shotgun Sequencing Data.” BMC Genomics 19, no. 1 (June 25, 2018): 493. https://doi.org/10.1186/s12864-018-4870-z.

### Integrative analysis

- `MMCA` - microbiome and metabolome correlation analysis pipeline. Five correlation methods (Spearman default), unsupervised (CCA, O2PLS) and supervised (PCA, PLS-DA, OPLS-DA, RF) analyses, network (WGCNA) analysis, KEGG enrichment in modules (Tax4Fun2). Report generation. http://mmca.met-bioinformatics.cn/
    - Ni, Yan, Gang Yu, Yongqiong Deng, Xiaojiao Zheng, Tianlu Chen, Junfeng Fu, and Wei Jia. “MMCA: A Web-Based Server for the Microbiome and Metabolome Correlation Analysis.” BioRxiv, January 1, 2019, 678813. https://doi.org/10.1101/678813.


## Taxonomy

- `Centrifuge` - microbial classification using BWT and FM index. Compression of genome sequences before indexing, then progressive exact matching of k-mers. Memory-efficient and faster than Kraken. https://github.com/infphilo/centrifuge
    - Kim, Daehwan, Li Song, Florian P. Breitwieser, and Steven L. Salzberg. “Centrifuge: Rapid and Sensitive Classification of Metagenomic Sequences.” Genome Research 26, no. 12 (2016): 1721–29. https://doi.org/10.1101/gr.210641.116.

- `Kraken` - assigning taxonomic labels to metagenomic DNA sequences. Exact matching of k-mers (31bp) against a database (different versions for memory considerations). Their own optimized algorithm for k-mer match search. https://ccb.jhu.edu/software/kraken2/
    - Wood, Derrick E., and Steven L. Salzberg. “Kraken: Ultrafast Metagenomic Sequence Classification Using Exact Alignments.” Genome Biology 15, no. 3 (March 3, 2014): R46. https://doi.org/10.1186/gb-2014-15-3-r46.

- `KrakenUniq` - Extension of the original k-mer-based classification with a HyperLogLog algorithm for assessing the coverage of unique k-mers (cardinality). Better handling of false positives. https://github.com/fbreitwieser/krakenuniq
    - Breitwieser, F. P., D. N. Baker, and S. L. Salzberg. “KrakenUniq: Confident and Fast Metagenomics Classification Using Unique k-Mer Counts.” Genome Biology 19, no. 1 (December 2018). https://doi.org/10.1186/s13059-018-1568-0.

- `metagenomeFeatures` - R package for annotating OTUs with Greengene IDs (v.13.8), RDP and SILVA (in future?). https://bioconductor.org/packages/release/bioc/html/metagenomeFeatures.html
    - Olson, Nathan D, Nidhi Shah, Jayaram Kancherla, Justin Wagner, Joseph N Paulson, and Hector Corrada Bravo. “MetagenomeFeatures: An R Package for Working with 16S RRNA Reference Databases and Marker-Gene Survey Feature Data.” Edited by Janet Kelso. Bioinformatics, March 1, 2019. https://doi.org/10.1093/bioinformatics/btz136.


## Phylogenetics

- `iTOL` - display and annotation of phylogenetic trees, https://itol.embl.de/

- `ggtree` R package for phylogenetic tree visualization, coloring, and annotation. Support for multiple file formats. https://github.com/GuangchuangYu/ggtree
    - Yu, Guangchuang, David K. Smith, Huachen Zhu, Yi Guan, and Tommy Tsan-Yuk Lam. “Ggtree: An R Package for Visualization and Annotation of Phylogenetic Trees with Their Covariates and Other Associated Data.” Edited by Greg McInerny. Methods in Ecology and Evolution 8, no. 1 (January 2017): 28–36. https://doi.org/10.1111/2041-210X.12628.

- `phyloT` - generates phylogenetic trees from based on the NCBI taxonomy. Input: NCBI scientific names and more, output: tree in Newick and other formats. Results can be visualized in iTOL, interactive Tree Of Life. https://phylot.biobyte.de/

## Differential analysis

- `ALDEx2` R package - a compositional data analysis tool that uses Bayesian methods to infer technical and statistical errors. Works with RNA-seq, microbiome, and other compositional data. Distinction between absolute counts and compositional data. Counts are converted to probabilities by Monte Carlo sampling (128 by default) from the Dirichlet distribution with a uniform prior. Centered log-ratio transformation, clr - divide by the geometric mean. https://bioconductor.org/packages/release/bioc/html/ALDEx2.html
    - Fernandes, Andrew D., Jennifer Ns Reid, Jean M. Macklaim, Thomas A. McMurrough, David R. Edgell, and Gregory B. Gloor. “Unifying the Analysis of High-Throughput Sequencing Datasets: Characterizing RNA-Seq, 16S RRNA Gene Sequencing and Selective Growth Experiments by Compositional Data Analysis.” Microbiome 2 (2014): 15. https://doi.org/10.1186/2049-2618-2-15.

- `metagenomeSeq` R package. Differential microbial abundance analysis. New normalization - Cumulative-sum scaling (CSS) - raw counts are divided by the cumulative sum of counts up to a percentile determined using a data-driven approach, e.g., the 75th percentile of each sample’s nonzero count distribution. Zero-inflated Gaussian (ZIG) distribution mixture model that accounts for biases in differential abundance testing resulting from undersampling of the microbial community, https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html
    - Paulson, Joseph N, O Colin Stine, Héctor Corrada Bravo, and Mihai Pop. “Differential Abundance Analysis for Microbial Marker-Gene Surveys.” Nature Methods 10, no. 12 (September 29, 2013): 1200–1202. https://doi.org/10.1038/nmeth.2658.


## Data

- American Gut project paper (http://americangut.org/). In collaboration with Earth Microbiome Project (http://www.earthmicrobiome.org/), Human Food Project (http://humanfoodproject.com/). Many results and statistical analyses in Methods. Data and results page, http://americangut.org/resources/, Jupyter analysis notebooks, https://github.com/knightlab-analyses/american-gut-analyses
    - McDonald, Daniel, Embriette Hyde, Justine W. Debelius, James T. Morton, Antonio Gonzalez, Gail Ackermann, Alexander A. Aksenov, et al. “American Gut: An Open Platform for Citizen Science Microbiome Research.” Edited by Casey S. Greene. MSystems 3, no. 3 (May 15, 2018). https://doi.org/10.1128/mSystems.00031-18.

- Earth Microbiome Project (EMP), ~100 studies, meta-analysis. 16S V4 sequencing. Deblur to summarize the data at the ASV level. Greengenes and SILVA rRNA gene databases. Analysis, alpha-, beta-, Shannon-, chao2, faith_pd diversity, PCoA. Data on Qiita. Code and data on https://github.com/biocore/emp. Unfiltered and filtered OTU tables at ftp://ftp.microbio.me/emp/release1/otu_tables
    - Thompson, Luke R., Jon G. Sanders, Daniel McDonald, Amnon Amir, Joshua Ladau, Kenneth J. Locey, Robert J. Prill, et al. “A Communal Catalogue Reveals Earth’s Multiscale Microbial Diversity.” Nature 551, no. 7681 (23 2017): 457–63. https://doi.org/10.1038/nature24621.


- MicrobiomeHD database - 28 case-control studies of gut microbiome, ten diseases (IBD, Crohn's, etc.), 16S sequencing. 16S processing pipeline: https://github.com/thomasgurry/amplicon_sequencing_pipeline, Python code to reproduce all analyses: https://github.com/cduvallet/microbiomeHD, raw sequencing data, metadata, OTU tables https://zenodo.org/record/1146764#.XDQHec9KjfY
    - Duvallet, Claire, Sean M. Gibbons, Thomas Gurry, Rafael A. Irizarry, and Eric J. Alm. “Meta-Analysis of Gut Microbiome Studies Identifies Disease-Specific and Shared Responses.” Nature Communications 8, no. 1 (December 2017). https://doi.org/10.1038/s41467-017-01973-8.

- Database for Preterm Birth Research. Other databases - dbPTB, GEneSTATION. Different studies. Predominantly microbiome, but also CyTOF, RNA-Seq, cell-free DNA and RNA sequencing, and genotyping. Open access. http://www.immport.org/resources/mod
    - Sirota, Marina, Cristel G. Thomas, Rebecca Liu, Maya Zuhl, Payal Banerjee, Ronald J. Wong, Cecele C. Quaintance, et al. “Enabling Precision Medicine in Neonatology, an Integrated Repository for Preterm Birth Research.” Scientific Data 5 (November 6, 2018): 180219. https://doi.org/10.1038/sdata.2018.219. 

- LogMIPE study (Landscape of gut Microbiome - Pan-India Exploration) - FASTQ files from 1004 subjects from 14 geographical locations. Data, https://www.ebi.ac.uk/ena/data/view/PRJEB25642, and scripts for processing them using QIIME or Mothur pipelines, https://github.com/anirbanbhaduri/LogMPIE
    - Dubey, Ashok Kumar, Niyati Uppadhyaya, Pravin Nilawe, Neeraj Chauhan, Santosh Kumar, Urmila Anurag Gupta, and Anirban Bhaduri. “LogMPIE, Pan-India Profiling of the Human Gut Microbiome Using 16S RRNA Sequencing.” Scientific Data 5 (October 30, 2018): 180232. https://doi.org/10.1038/sdata.2018.232.

- Qiita - database and analysis platform for meta-analysis of microbiome studies. redbiom plugin allows for searching. Web site, https://qiita.ucsd.edu/, GitHub, https://github.com/biocore/qiita
    - Gonzalez, Antonio, Jose A. Navas-Molina, Tomasz Kosciolek, Daniel McDonald, Yoshiki Vázquez-Baeza, Gail Ackermann, Jeff DeReus, et al. “Qiita: Rapid, Web-Enabled Microbiome Meta-Analysis.” Nature Methods 15, no. 10 (October 2018): 796–98. https://doi.org/10.1038/s41592-018-0141-9.


## Misc

- Links to microbiome-related resources by Guangchuang Yu, https://github.com/GuangchuangYu/bookmarks

- Metagenomics blog https://microbe.land/ and a list of metagenomics-related methods, tools, and many more at the Google Doc [Metagenomics - Tools, Methods and Madness](https://docs.google.com/document/d/e/2PACX-1vQbLMrFcpFh8asvZsUv95wQWwTzQYBgtadDiVKffSA33Oi_vZNdi0czrEPUL1seOZLd1HaqWs29H6hp/pub)

- `VAMB` - Binning metagenomic sequences using variational autoencoders for representing the data for clustering, visualization. Uses tetranucleotide (ATCG) frequency of 4-mers, and their abundance (sequencing depth). Application to human gut data, benchmarking against other methods. https://github.com/jakobnissen/vamb
     - Nissen, Jakob Nybo, Casper Kaae Sonderby, Jose Juan Almagro Armenteros, Christopher Heje Groenbech, Henrik Bjorn Nielsen, Thomas Nordahl Petersen, Ole Winther, and Simon Rasmussen. “Binning Microbial Genomes Using Deep Learning.” BioRxiv, January 1, 2018, 490078. https://doi.org/10.1101/490078.

- MetaMap - microbial reads in human RNA-seq data. Unmapped reads ran through CLARK-S to get OTU tables. Data http://gigadb.org/dataset/100456.
    - Simon, L M, S Karg, A J Westermann, M Engel, A H A Elbehery, B Hense, M Heinig, L Deng, and F J Theis. “MetaMap: An Atlas of Metatranscriptomic Reads in Human Disease-Related RNA-Seq Data.” GigaScience 7, no. 6 (June 1, 2018). https://doi.org/10.1093/gigascience/giy070.

- MetaMap web app. Seven visualization methods, interactive links show them. GitHub https://github.com/theislab/MetaMap, and web-tool http://146.107.176.18:3838/MetaMap/R/
    - Simon, Lukas, George Tsitsiridis, Philipp Angerer, and Fabian Theis. “MetaMap, an Interactive Webtool for the Exploration of Metatranscriptomic Reads in Human Disease-Related RNA-Seq Data,” September 30, 2018. https://doi.org/10.1101/425439.

- Microbial species behavior in soil. Good visualization, how-to collected in https://github.com/visualization-bioinformatics/cool.visualization.techniques
    - George, Paul B. L., Delphine Lallias, Simon Creer, Fiona M. Seaton, John G. Kenny, Richard M. Eccles, Robert I. Griffiths, et al. “Divergent National-Scale Trends of Microbial and Animal Biodiversity Revealed across Diverse Temperate Soil Ecosystems.” Nature Communications 10, no. 1 (March 7, 2019): 1107. https://doi.org/10.1038/s41467-019-09031-1.

- `PathSeq` - GATK tool to detect microbial DNA in human sequencing. Reads that do not align to human are progressively aligned to microbial genomes, and a table of detected microbial organisms is generated. https://gatkforums.broadinstitute.org/gatk/discussion/23205/cross-species-contamination-identification-with-pathseq


# Local files and folders

- `LogMPIE` - Microbiome Processing Pipeline, Python scripts for running Mothur, QIIME. https://github.com/anirbanbhaduri/LogMPIE

