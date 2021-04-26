# Pipeline to analyze proteomic dataset for the protocol paper "Mapping proximity networks in yeast mitochondria using BioID coupled to proteomics (DOI: 10.1016/j.xpro.2020.100219)"

### How to run two-group comparison(using LIMMA) code?
* Implementation of LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group comparision in a proteomic experiment [1].
* The pipeline is implemented in R programming language and needed libraries are installed automatically.
* Make sure you have mitochondrial proteins (identified with uniprot id) in your proteingroups otherwise code will terminate with a message. See Mitochondrial_Proteins folder for a list of experimentally determined mitochondrial proteins.
* Run the limma_main.R code by clicking **Source** in RStudio and select a MaxQuant outputted proteingroups.txt file with either iBAQ/LFQ values (see  Example/MQ_output/proteinGroups_example.txt). limma_helper_functions.R, the file containing the helper functions must be in the same directory as the limma_main.R.
* There are two modes: either use full data or remove exclusive proteins before analysis. Proteins that have zero intensities in all the replicates of one group are defined as exclusive proteins. For example  a set of three 'exclusively enriched' proteins in bait 'Cbp1' is defined as follows

 |Uniprot | Symbol  |Cbp1_1  |  Cbp1_2  |  Cbp1_3  |  contrl_1 |   contrl_2  |  contrl_3|
 |--------|---------|--------|----------|----------|-----------|-------------|----------|
 |P25554  |SGF29    |0       |     0    |    0     |   2810900 |    4903800  |     0    |
 
* Proteins that do not participate in two group comparison are saved in the results folder as tsv files with their iBAQ/LFQ intensities.
* The code is interactive and will help users to provide the correct names as input for bait and control as it appears in the proteingroups file columns. Part of the name (case insensitive) is also accepted. It will also ask users if they want to median normalize their data prior to two-group comparison and if chosen to do so, it will boxplot (in RStudio 'Plots' panel) the data before and after normalization. It will also plot the data distribution (histogram) before and after missing value imputation.           
* After successful run, it will create a volcano plot in html format and a tsv file containing final data inside a folder called "Results_timestamp" with the current system "timestamp" in the same directory where the limma_main.R file is present. You can find two such sample result folders in Example/Limma_output folder. 
* You can view the plot in any browser and save it as png by clicking camera icon in plot
* In addition, the code will save 'exclusively enriched' proteins (if any) in control and bait replicates with corresponding LFQ/iBAQ values in the "Results_timestamp" folder.
* Tested on Windows10 64 bit, Rstudio v1.3.959, R v4.0.2 and using the following packages

| Package  | Version |
| ------------- | ------------- |
| dplyr         | >= 0.8.3  |
| stringr       | >= 0.4.0  |
| MASS       | >= 7.3-5.4 |
| plotly       | >= 4.9.0  |
| htmlwidgets       | >= 0.3 |
| limma       | >= 3.42.0  |
| qvalue       | >= 2.8.0  |

* We recommend using the same versions of the packages (atleast for plotly)

### How to generate edge table (for Cytoscape) from LIMMA output?
*	If you run multiple two-group comparisons, you will end up creating multiple "Results_timestamp" folders. Inside each "Results_timestamp" folder there will be one "timestamp_final_data.tsv" file.
*	Create a folder in your preferred location and copy those "timestamp_final_data.tsv" into that folder.
*	Then run the "automate_nw_tab_gen.R" code in RStudio and select the folder where you copied those .tsv files (see Example/Limma_output/Final_tsvs folder).
*	The code will first ask users to specify a log fold change cutoff. Then for every .tsv file in that folder it will print a list of iBAQ/LFQ column names and ask you to enter a bait name from that list. Here assumption is that iBAQ/LFQ columns will contain bait names.
*	Using those information, it will create a virtual list of bait-prey interactions that are >= cutoff and < pval(0.05).
*	After that it will print a message if that file is processed successfully.
*	Finally, when all .tsv files are processed successfully, code will concatenate those virtual lists into a single one, print top 10 rows that list and save it as "Links.xlsx" in the same directory where you copied .tsv files earlier (See Example/Limma_output/Final_tsvs/Links.xlsx).
*	Then you can load "Links.xlsx" to visualize a bait-prey interaction network.
*	Open Cytoscape and click File &#8594; Import &#8594; Network from File. and load "Links.xlsx".
*	Open Cytoscape and click File &#8594; Import &#8594; Network from File. and load "Links.xlsx".
*	For 'Bait' select Meaning: 'Source Node' and Data type: 'String'. For 'Prey' select Meaning: 'Target Node' and Data type: 'String'. For 'Fold change' select Meaning: 'Edge attribute' and Data type: 'Floating point'. Click 'OK'.
*	Click Edit &#8594; Remove Self-Loops. to remove self-interactions involving baits.
*	The network displayed represents Bait-prey undirected network which can be customized depending on the user preference.
*	Tested on Windows 10 (64 bit) using Cytoscape 3.7.2.

### Cite this repo
Salvatori, R., Aftab, W., Forne, I., Imhof, A., Ott, M. and Singh, A.P., 2020. Mapping protein networks in yeast mitochondria using proximity-dependent biotin identification coupled to proteomics. [STAR protocols](https://star-protocols.cell.com/protocols/327), p.100219.

#### Reference:
[1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.
