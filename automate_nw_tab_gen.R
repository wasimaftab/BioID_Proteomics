## clear screen
cat('\014')

## clear environment
rm(list = ls())

## load lib(s)
library(rstudioapi)
library(stringr)
## select a directory
select_folder <- function(my_caption) {
    
    ##select dir using rstudioapi package
    selectDirectory(caption = my_caption)
}

## Change to source dir
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

## automate network file gen.
tsv_data_dir <- select_folder("Select folder containing Limma final data files(tsv)")
files <- list.files(path = tsv_data_dir, pattern = '\\.tsv$')
nw_tab <- data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE)
fc <- as.double(readline('Enter the fold change threshold = '))
pval <- 0.05

for (i in 1:length(files)) {
    d <- read.table(file = paste0(tsv_data_dir, '/', files[i]), sep = '\t', header = TRUE)
    idx_pval <- which(d$p.mod < pval)
    idx_high_fc <- which(d$logFC >= fc)
    idx_desired <- intersect(idx_high_fc, idx_pval)
    idx <- grep('iBAQ', colnames(d), ignore.case = FALSE)
    print(colnames(d)[idx])
    bait_name <- readline('Check the list above and enter a bait name = ')
    nw_tab <- rbind(nw_tab, cbind(rep(bait_name, length(idx_desired)),
                                  d$Symbol[idx_desired],
                                  d$logFC[idx_desired]))
    print(paste('Successfully processed', files[i]))
    cat('\n')
}
colnames(nw_tab) <- c('source', 'target', 'fc')

## convert fc column to numeric
nw_tab$fc <- as.numeric(nw_tab$fc)

## use stringr::str_to_sentence to make first letter capital
nw_tab$source <- str_to_sentence(nw_tab$source)
nw_tab$target <- str_to_sentence(nw_tab$target)

print('Network table created whose top few rows are as follows:')
cat('\n')
print(head(nw_tab))

# ## write network table as excel file
writexl::write_xlsx(nw_tab, path = paste0(tsv_data_dir, '/', 'Links.xlsx'), col_names = TRUE)
