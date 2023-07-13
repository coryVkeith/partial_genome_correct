#!/bin/env Rscript

# Loads necessary libraries
library(ggplot2)
library(gridExtra)
library(plyr)
library(Rcpp)

# Gathers arguments and passes only the last arguments. Trailingonly allows the positional arguments passed from the bash script to be set as 1 and 2.
args <- commandArgs(trailingOnly = TRUE)

# Sets passed arguments to new variables in R.
HELPER_DIR <- args[1]
SMPLE <- args[2]

# Calls files in designated folder to iterate through in following code. Matches the regular pattern ".bam.coverage". 
file_path <- gsub(" ","", paste(HELPER_DIR,"/",SMPLE,"/bowtie2/"))
files <- list.files(path=file_path, pattern=".coverage", full.names=FALSE, recursive=FALSE)
setwd(file_path)

# Sets a list of plots that is the length of files in the folder. Plots for each file generated in the following code will be saved to this list.
plist <- vector("list", length(files))

# Sets the names of the plist to the names of the files.
names(plist) = files

#Iterates through files
for( i in files ){
    # Strips filenames of only informative parts and creates a variable for the plot name.
    #samp_suff <- gsub("\\.*REF_*","",i)
    #samp_suff <- gsub("(.*)_REF_.*","\\1",i) 
    #print(samp_suff)
    # Opens file and creates a dataframe
    df=read.delim(i, header=F)
    names(df) <- c("name", "Position", "Depth")
    df$name <- "names(i)"
 
    # Creates a plot of the open coverage file.
    p <- ggplot(df, aes(x=Position, y=Depth, fill = name)) + geom_area()
    p <- p + ggtitle(names(i)) + theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill="white", colour="gray53"))
    p <- p + scale_fill_manual(values="black") + guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill="white", colour="gray53"), text = element_text(size=20))
    
    # Adds plot to plist outside of for loop.
    plist[[i]] = p
}

# Sets a dataframe as the plist, list of plots, generated above.
df = ldply(plist, function(x) x$data)
names(df)[1] <- "files"

# Plots plist.
p = ggplot(df, aes(x=Position, y=Depth, fill = name)) + geom_area()
p = p + ggtitle(samp_suff) + theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill="white", colour="gray53"))
p = p + scale_fill_manual(values="black") + guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill="white", colour="gray53"), text = element_text(size=20))
p = p + facet_wrap(~files, scales="free")

#Writes the file
svg_file <- gsub(" ","", paste(samp_suff,"_combined.svg"))
svg(filename = svg_file, width = 24, height = 8)
p1 <- grid.arrange(p, ncol =1)
print(p1)
dev.off()

