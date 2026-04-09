# Longitudinal clinical monitoring reveals widespread effects of medications on gut microbiome dynamics

Ashwin Chetty<sup>1</sup>, Matthew A. Odenwald<sup>2</sup>, Ramanujam Ramaswamy<sup>3</sup>, Nicholas Dylla<sup>3</sup>, Huaiying Lin<sup>3</sup>, Victoria Burgo<sup>3</sup>, Ashley M. Sidebottom<sup>3</sup>, Eric G. Pamer<sup>3,4</sup>, Ran Blekhman<sup>5</sup>

<sup>1</sup>Committee on Genetics, Genomics and Systems Biology, University of Chicago, Chicago, IL, USA

<sup>2</sup>Department of Medicine, Section of Gastroenterology, Hepatology, and Nutrition, University of Chicago, Chicago, IL, USA

<sup>3</sup>Duchossois Family Institute, University of Chicago, Chicago, IL, USA 

<sup>4</sup>Department of Medicine, Section of Infectious Diseases and Global Health, University of Chicago, Chicago, IL, USA

<sup>5</sup>Section of Genetic Medicine, Department of Medicine, The University of Chicago, Chicago, IL, USA

### Abstract
Medications exert strong effects on the gut microbiome, impacting host physiology and affecting health outcomes. Studies have characterized associations between medications and the gut microbiome, primarily through in vitro analyses and, when conducted in humans, either using cross-sectional designs or focusing on specific medication classes. However, how medications directly reshape the microbiome within individual patients over time, especially in clinical settings, remains poorly understood.  Here, we leverage longitudinal stool metagenomics, metabolomics, and electronic health records from 3,469 samples representing 1,122 patients at the University of Chicago Medical Center between 2020 and 2024. Using a longitudinal design where each patient serves as their own control, we compare patient microbiome composition before and after starting a medication while controlling for the effects of medical procedures, diagnoses, demographic variables, and other drugs taken simultaneously. In total, we identify 36,637 associations relating 138 of the most commonly given medications, such as oxycodone and pantoprazole, with 106 microbial genera, 204 microbial species, 627 microbial pathways, and 35 gut metabolites. We find, for example, that oral prednisone is associated with decreases in Parabacteroides and Enterobacter abundance and increases in microbial ergosterol and long-chain fatty acid biosynthesis pathways. Furthermore, the selective serotonin reuptake inhibitor sertraline is associated with increases in Alistipes abundance, increases in microbial pathways related to dopamine degradation, and increases in fecal concentration of the related metabolite kynurenic acid. Critically, we discover that medication effects are time-dependent, with microbiome recovery to baseline occurring over days to weeks. Our findings represent the first systematic characterization of medication-microbiome dynamics within individual patients in a clinical setting, establishing temporal patterns and mechanistic pathways that inform microbiome-guided therapeutic strategies. 

### Description

This page contains code to run the models and generate the figures in this paper. Code, output files, and selected de-identified data files are available on Zenodo. 
To generate the figures, run:
`./build_figures.py`



### Dependencies

| Package  | Version |
| ------------- |:-------------:|
| numpy | 1.23.5 |
| pandas | 1.5.3 |
| matplotlib | 3.6.3 |
| sklearn | 1.2.1 |
| skbio | 0.5.8 |
| seaborn | 0.12.2 |
| scipy | 1.10.0 |
| statsmodels | 0.13.5 |
| matplotlib | 3.6.3 |
| networkx | 3.0 |

