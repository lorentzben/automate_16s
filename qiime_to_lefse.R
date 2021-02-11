#! /usr/bin/env Rscript --vanilla
require(dplyr)
require(tibble)
require(qiime2R)
require(phyloseq)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "lefse_formatted.txt"
}
ioi <- args[1]
phyloseq_to_lefs <- function(physeq){
  # aggregate to genus level
  ps_lefse <- physeq %>% tax_glom(taxrank = 'Genus', NArm = F)
  
  # extract taxonomic data from phyloseq object and then stored in a vector called lefse_tax
  lefse_tax <- ps_lefse %>% tax_table %>% data.frame(stringsAsFactors=FALSE)
  lefse_tax <- replace(lefse_tax, is.na(lefse_tax), 'Unclassified')
  lefse_tax <- lefse_tax %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% mutate(id = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "|")) %>% ungroup %>% pull(id)
  
  # extract otu abundance matrix from phyloeq object and annotated with tax information
  lefse_matrix <- otu_table(ps_lefse) %>% data.frame(stringsAsFactors = F) %>% t %>% data.frame
  colnames(lefse_matrix) <- lefse_tax 
  # I added this line to fix an error of '-' in sample names being turned into '.'
  rownames(lefse_matrix) <- ps_lefse %>% sample_data %>% data.frame(stringsAsFactors = F) %>% rownames
  
  # extract sample matadata and order sample same in lefse_matrix
  lefse_sample <- sample_data(ps_lefse)
  
  # convert factor in lefse_sample into character in order to combine sample and abundance data
  lefse_sample_isfactor <- sapply(lefse_sample, is.factor)
  lefse_sample[,lefse_sample_isfactor] <- lefse_sample[,lefse_sample_isfactor] %>% lapply(as.character)
  lefse_sample <- lefse_sample %>% data.frame
  
  lefse_table <- full_join(rownames_to_column(lefse_sample), rownames_to_column(lefse_matrix), by = ("rowname" = "rowname")) %>% t
  
  return(lefse_table)
}

cycle_1 <- qza_to_phyloseq("table-dada2.qza","rooted-tree.qza","taxonomy.qza","metadata.tsv")

# modifications to select item of interest and remove the rest of the metadata
new_samp_2 <- data.frame(sample_data(cycle_1))
new_samp_2 <- new_samp_2 %>% 
  select(as.name(ioi))

# generates phyloseq object with sample data of only item of interest so that phyloseq to lefs can run
cycle_2 <- phyloseq(cycle_1@otu_table, cycle_1@tax_table, cycle_1@phy_tree, sample_data(new_samp_2))

# transforms phyloseq object into lefse input object 
mm <- phyloseq_to_lefs(cycle_2)

write.table(mm, args[2] , sep="\t", row.names = T, col.names = F, quote = F)
