# This is a pipeline to analyze proteiomic data for the protocol paper "Mapping proximity networks in yeast mitochondria using BioID coupled to proteomics"

### How to run two-group comparison(using LIMMA) code?
* Implementation of LIMMA (Linear Models for Microarray Data), an empirical Bayes method for two group comparision in a proteomic experiment [1].
* The pipeline is implemented in R programming language and needed libraries are installed automatically.
* Run the limma_main.R code and select a MaxQuant outputted proteingroups.txt file with either iBAQ/LFQ values (see Example/MQ_output/proteinGroups Mrpl4, cbp3, cbs1.txt file).
* Make sure you have mitochondrial proteins (identified with uniprot id) in your proteingroups otherwise code will terminate with a message.
* The helper functions must be in the same directory as the main.
* After successful run it will create volcano plots both in html format and a tsv file containing final data inside a folder called  "Results_timestamp" with the current "timestamp" (in the same directory where the limma_main.R file is present).           
* You can also save the plot as png from the browser.
* In addition, the code will create files containng exclusive proteins present (if any) in control and bait/treatment replicates with LFQ/iBAQ values, for example an exclusive protein in control is defined as follows, i.e. it is not qunatified in the bait but in one or more control replictes.

 |Uniprot | Symbol  |treat_1 |  treat_2 |  treat_3 |  contrl_1 |   contrl_2  |  contrl_3|
 |--------|---------|--------|----------|----------|-----------|-------------|----------|
 |P25554  |SGF29    |0       |     0    |    0     |   2810900 |    4903800  |     0    |
* The code will help you to provide correct names as input for bait/treatment and control as it appears in the proteingroups file. Part of the name (case insensitive) is also acceptrd.
* Tested on Windows10 64 bit, Rstudio v1.3.959, R v4.0.2 and using the following packages

| Package  | Version |
| ------------- | ------------- |
| dplyr         | 0.8.3  |
| stringr       | 0.4.0  |
| MASS       | 7.3-5.4 |
| plotly       | 4.9.0  |
| htmlwidgets       | 0.3 |
| limma       | 3.42.0  |
| qvalue       | 2.8.0  |

* We recommend using the same versions of the packages (atleast for plotly)

### How to generate edge table (for Cytoscape) from LIMMA output?
* If you run multiple two-group comparision, you will end up creating multiple "Results_timestamp" folders.
* Inside each "Results_timestamp" folder there will be one "timestamp_final_data.tsv" file.
* Create a folder in your preffered location and copy those "timestamp_final_data.tsv" into that folder.
* Then run the "automate_nw_tab_gen.R" code in RStudio and select the folder where you copied those .tsv files(see Example/Limma_output folder).
* For every .tsv file in that folder and each bait involing a two-group comparison, code will first ask users to specify a log fold change cutoff, then it will print a list of column names and ask you to enter a bait name
* Using those information, it will create a vitual list of bait-prey interactions that are >= cutoff and < pval(0.05)
* After that it will print a message if that file is processed successfully
* Finally, when all .tsv files are processed successfully, code will concatenate those vitual lists into a single one, print top 10 rows that list and save it as "Links.xlsx" in the same directory where you copied .tsv files earlier.
* Then you can load "Links.xlsx" to visualize a bait-prey interaction network.



#### Reference:
[1] Kammers, Kai, et al. "Detecting significant changes in protein abundance." EuPA open proteomics 7 (2015): 11-19.
