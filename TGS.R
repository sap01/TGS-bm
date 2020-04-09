#!/usr/bin/env Rscript
## Goal of this program: Implement the TGS Algorithm.
##
## Author: Saptarshi Pyne (saptarshipyne01@gmail.com)
## Last modified on: Mar 6, 2018
##
## The coding practices followed are documented at:
## https://github.com/sap01/coding_practices/blob/master/R_coding_practices.md
##
## How to execute this script:
## For Unix-alike OSes:
## Let us assume that this script is inside directory '/home/saptarshi/R/R-3.3.2/projects/repoTGS' and
## the Rscript file is inside directory '/home/saptarshi/R/R-3.3.2/bin'.
## Then, execute this script using the following commands:
## $ cd /home/saptarshi/R/R-3.3.2/projects/repoTGS/  
## $ nohup time /home/saptarshi/R/R-3.3.2/bin/Rscript /home/saptarshi/R/R-3.3.2/projects/repoTGS/TGS.R input.json &
## where 'asset/input.json' contains the user-defined parameters. A file 
## named 'nohup.out' will be generated inside 
## '/home/saptarshi/R/R-3.3.2/projects/repoTGS/'.
##
## For Windows OSes:
## Let us assume that this script is inside directory 'D:\R\R-3.3.2\projects\repotgsr' and
## the 'Rscript.exe' file is inside directory 'C:\Program Files\R\R-3.3.1\bin'.
## Then, execute this script using the following commands (the '>' symbol
## represents the DOS command prompt):
## >cd "D:\R\R-3.3.2\projects\repotgsr"
## >"C:\Program Files\R\R-3.3.1\bin\Rscript.exe" TGS.R input.json
## where '/repotgsr/asset/input.json' contains the user-defined parameters.
##
## Input: A time series gene expression dataset with multiple time series.
##
## Output: Time-varying Gene Regulatory Networks and a corresponding rolled up network.

## Remove all objects in the current workspace
rm(list = ls())

##------------------------------------------------------------
## Begin: Load the Required Libraries
##------------------------------------------------------------
## For reading from and writing to '.json' files
library(rjson)

##------------------------------------------------------------
## End: Load the Required Libraries
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Read User-defined input Params
##------------------------------------------------------------
input.args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 1)
{
  stop("Exactly one input file must be supplied.", call.=FALSE)
}

input.params <- rjson::fromJSON(file = paste(getwd(), 'asset', input.args, sep = '/'))
rm(input.args)

## Input file for time-series gene expression data
input.data.filename <- input.params$input.data.filename
input.data.filename <- paste(getwd(), 'asset', input.data.filename, sep = '/')

## Number of time points (T)
num.timepts <- input.params$num.timepts

## True rolled network file.
## If 'true.net.filename' is an empty string, then the true rolled network
## is not known a prior. Otherwise, it is known a prior and would be used
## to evaluate performance metrics of the predicted rolled network.
true.net.filename <- input.params$true.net.filename
if (true.net.filename != '')
{
  true.net.filename <- paste(getwd(), 'asset', true.net.filename, sep = '/')
}

## Input file for Wild Type (WT) values of the genes.
## If 'input.wt.data.filename' is an empty string, then
## the WT values are not required for further computation. Otherwise, 
## WT values are required for further computation.
input.wt.data.filename <- input.params$input.wt.data.filename
if (input.wt.data.filename != '')
{
  input.wt.data.filename <- paste(getwd(), 'asset', input.wt.data.filename, sep = '/')
}

## is.discrete is true or false, implying whether the input data file is already
## discretized or not.
is.discrete <- input.params$is.discrete

## Number of discrete levels in the input data or
## number of discrete levels in which the data needs to be discretized.
num.discr.levels <- input.params$num.discr.levels

## Name of the discretization algorithm to be applied
discr.algo <- input.params$discr.algo

mi.estimator <- input.params$mi.estimator

apply.aracne <- input.params$apply.aracne

## Name of the CLR algorithm to be applied
clr.algo <- input.params$clr.algo

## The maximum number of regulators a gene can have
max.fanin <- input.params$max.fanin

## allow.self.loop takes values true or false, dependending on whether 
## to allow self loops in the predicted rolled network or not.
allow.self.loop <- input.params$allow.self.loop

rm(input.params)
##------------------------------------------------------------
## End: Read User-defined input Params
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Create the output directory
##------------------------------------------------------------
init.path <- getwd()

## Output directory name
output.dirname <- paste('output', format(Sys.time(), "%Y%m%d%H%M%S"), sep = '')

if(.Platform$OS.type == 'windows') {
  
  if(! output.dirname %in% shell("ls asset" , intern = TRUE)) {
    
    ## Output directory name for Windows OSes
    output.dirname <- paste('asset', output.dirname, sep = '/')
    output.dirname <- paste(init.path, output.dirname, sep = '/')
    
    ## Convert directory path to canonical form for the Windows OS.
    ## It raises the warning if the directory does not exist, which
    ## is expected. Therefore, please ignore the warning.
    output.dirname <- normalizePath(output.dirname, winslash = '\\', mustWork = NA)
    
    shell(paste('mkdir ', output.dirname, sep = ''), intern = TRUE, mustWork =NA)
  }
} else {
  ## .Platform$OS.type != 'windows'
  
  if(.Platform$OS.type == 'unix') {
    if(! output.dirname %in% system("ls asset" , intern = TRUE))
    {
      output.dirname <- paste('asset', output.dirname, sep = '/')
      output.dirname <- paste(init.path, output.dirname, sep = '/')
      
      system(paste('mkdir ', output.dirname, sep = ''))
    }
  } 
} 
##------------------------------------------------------------
## End: Create the output directory
##------------------------------------------------------------

# ##------------------------------------------------------------
# ## Begin: Create the output directory
# ##------------------------------------------------------------
# init.path <- getwd()
# 
# ## Output directory name
# output.dirname <- paste('asset/output', format(Sys.time(), "%Y%m%d%H%M%S"), sep = '')
# output.dirname <- paste(init.path, output.dirname, sep = '/')
# 
# if(.Platform$OS.type == "unix") {
#   if(! output.dirname %in% system("ls" ,intern=TRUE))
#   {
#     system(paste('mkdir ', output.dirname, sep = ''))
#   }
# } else{# if(.Platform$OS.type == "unix"){
#   shell(paste('mkdir ', output.dirname, sep = ''), intern = TRUE, mustWork =NA)
# }
# ##------------------------------------------------------------
# ## End: Create the output directory
# ##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Load the Required External Functions
##------------------------------------------------------------
source(paste(init.path, 'discretizeData.R', sep = '/'))
source(paste(init.path, 'compute_cmi.R', sep = '/'))
source(paste(init.path, 'learn_mi_net_struct.R', sep = '/'))

# Parallel programming with degree of parallelism 1 and Serial Programming
source(paste(init.path, 'learnDbnStruct3dParDeg1.R', sep = '/'))

source(paste(init.path, 'calcPerfDiNet.R', sep = '/'))
# source(paste(init.path, 'learnCmiNetStruct.R', sep = '/'))
source(paste(init.path, 'CompareNet.R', sep = '/'))
source(paste(init.path, 'rollDbn.R', sep = '/'))
source(paste(init.path, 'RDataToCytoscape.R', sep = '/'))
##------------------------------------------------------------
## End: Load the Required External Functions
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Main program
##------------------------------------------------------------

## Print the output dir name in 'nohup.out'
print('The output directory name is:')
print(output.dirname)
print('') ## to append a blank line

## Save console output in a file named 'output.txt' inside the output directory.
output.filename <- paste(output.dirname, 'output.txt', sep = '/')
output.file.conn <- file(output.filename, open = "wt")
sink(output.file.conn)

##------------------------------------------------------------
## Begin: Read input data file
##------------------------------------------------------------

## Begin: Find file extension of the input data file. Only '.tsv' and '.RData'
## are allowed.
## Split the string at every '.' and consider the last substring as the 
## file extension.
input.data.filename.ext <- unlist(strsplit(input.data.filename, '[.]'))
## End: Find file extension of the input data file. Only '.tsv' and '.RData'
## are allowed.

## Initialize input data
input.data <- NULL
if (input.data.filename.ext[length(input.data.filename.ext)] == 'tsv') {
  input.data <- read.table(input.data.filename, header = TRUE, sep="\t")
  
  timepts.names <- input.data[1:num.timepts, 1]
  
  ## Remove first col i.e. the time point names
  input.data <- input.data[, -1]
  
} else if (input.data.filename.ext[length(input.data.filename.ext)] == 'RData') {
  ## Loads an object named input.data
  load(input.data.filename)
  
  timepts.names <- 1:num.timepts
}

## Begin: Replace original node names with {v1, v2, ..., vV}
## V = total number of nodes in the input data.
## The replacement is crucial as unexpected characters in the node names
## may generate an error in further computation.

## Save the original node names in case they are required in future
orig.node.names <- colnames(input.data)

node.names <- c()
for (col.idx in 1:ncol(input.data))
{
  new.node.name <- paste('v', as.character(col.idx), sep = '')
  node.names <- c(node.names, new.node.name)
}
rm(col.idx)
colnames(input.data) <- node.names
## End: Replace original node names with {v1, v2, ..., vV}
## V = total number of node in the input data.

num.nodes <- ncol(input.data)

## Max fanin must be restricted to 14. Because, it is empirically observed that bnstruct::learn.network()
## function can learn a BN with upto 15 nodes without segmentation fault, given a 32 GB main memory. A max 
## fanin restriction of 14 limits the number of nodes in the to-be-learnt BN to 15 (1 target node and a 
## maximum of 14 candidate parents).
max.fanin <- min(num.nodes, 14)

num.samples.per.timept <- (nrow(input.data) / num.timepts)

## If input data is already discretized
if (is.discrete)
{
  input.data.discr <- input.data
  
} else
{
  if (discr.algo == '')
  {
    stop('Please specify the value of discr.algo.')
    
  } else if (discr.algo == 'discretizeData.2L.wt.l')
  {
    input.data.discr <- discretizeData.2L.wt.l(input.data, input.wt.data.filename)
    
  } else if (discr.algo == 'discretizeData.2L.Tesla')
  {
    input.data.discr <- discretizeData.2L.Tesla(input.data)
  }
  
  save(input.data.discr, file = paste(output.dirname, 'input.data.discr.RData', sep = '/'))
}
##------------------------------------------------------------
## End: Read input data
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Create a 3D array using discretized input data.
## Here,
## dim1 = time points,
## dim2 = nodes,
## dim3 = time series.
## It is useful for some CLR algos and the BN struct learning
## algos.
##------------------------------------------------------------
input.data.discr.matrix <- data.matrix(input.data.discr)

input.data.discr.3D <- array(NA, c(num.timepts, num.nodes, num.samples.per.timept),
                             dimnames = c(list(timepts.names), list(node.names),
                                          list(1:num.samples.per.timept)))

for (sample.idx in 1:num.samples.per.timept) {
  start.row.idx <- (1 + (num.timepts * (sample.idx - 1)))
  end.row.idx <- (num.timepts * sample.idx)
  input.data.discr.3D[ , , sample.idx] <- input.data.discr.matrix[start.row.idx:end.row.idx, ]
}
rm(sample.idx)

rm(input.data.discr.matrix)
##------------------------------------------------------------
## End: Create a 3D array using discretized input data.
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Learn MI Net Structure
##------------------------------------------------------------

start.time <- proc.time() # start the timer

##------------------------------------------------------------
## Begin: Initialize mutual information network
##------------------------------------------------------------

## Initialize mutual info matrix for static CLR net 
mi.net.adj.matrix <- NULL

## Initialize mutual info matrix for time-varying CLR nets
mi.net.adj.matrix.list <- NULL

mut.info.matrix <- NULL
if (clr.algo == 'CLR') {
  
  # Initialize mutual information matrix with zeroes
  mut.info.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                            dimnames = c(list(node.names), list(node.names)))
  
  if (mi.estimator == 'mi.pca.cmi') {
    ## Build mutual information matrix
    for (col.idx in 1:(num.nodes - 1)) {
      for (col.idx.2 in (col.idx + 1):num.nodes) {
        
        ## 'compute_cmi.R' 
        mut.info <- ComputeCmiPcaCmi(input.data.discr[, col.idx], input.data.discr[, col.idx.2])
        
        mut.info.matrix[col.idx, col.idx.2] <- mut.info
        mut.info.matrix[col.idx.2, col.idx] <- mut.info
      }
      rm(col.idx.2)
    }
    rm(col.idx)
    
  } else if (mi.estimator == 'mi.empirical') {
    mut.info.matrix <- minet::build.mim(input.data.discr,
                                        estimator = 'mi.empirical',
                                        disc = 'none')
    
  } else if (mi.estimator == 'mi.mm') {
    mut.info.matrix <- minet::build.mim(input.data.discr, 
                                        estimator = 'mi.mm', 
                                        disc = 'none')
  }
  
  if (apply.aracne == TRUE) {
    
    mut.info.matrix.pre.aracne <- mut.info.matrix
    
    mut.info.matrix <- minet::aracne(mut.info.matrix)
    
    mut.info.matrix.post.aracne <- mut.info.matrix
    
    elapsed.time <- (proc.time() - start.time)
    writeLines('elapsed.time just after the ARACNE step= \n')
    print(elapsed.time)
    rm(elapsed.time)
    
    save(mut.info.matrix.pre.aracne, 
         file = paste(output.dirname, 'mut.info.matrix.pre.aracne.RData', sep = '/'))
    
    save(mut.info.matrix.post.aracne, 
         file = paste(output.dirname, 'mut.info.matrix.post.aracne.RData', sep = '/'))
    
    rm(mut.info.matrix.pre.aracne, mut.info.matrix.post.aracne)
    
  } else {
    ## apply.aracne == FALSE
    
    # writeLines('mut.info.matrix = \n')
    # print(mut.info.matrix)
    save(mut.info.matrix, file = paste(output.dirname, 'mut.info.matrix.RData', sep = '/'))
  }
} else {
  ## clr.algo != 'CLR'
  
  rm(mut.info.matrix)
} 

if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1') | (clr.algo == 'spearman')) {
  
  ## CLR net is not time-varying
  rm(mi.net.adj.matrix.list)
  
  mi.net.adj.matrix <- matrix(0, nrow = num.nodes, ncol = num.nodes, 
                              dimnames = c(list(node.names), list(node.names)))
  
} else if (clr.algo == 'CLR3') {
  
  ## CLR net is not static
  rm(mi.net.adj.matrix)
  
  ## Pre-allocate an empty list of length = number of time intervals
  num.time.ivals <- (num.timepts - 1)
  mi.net.adj.matrix.list <- vector(mode = 'list', length = num.time.ivals)
  rm(num.time.ivals)
}
##------------------------------------------------------------
## End: Initialize mutual information network
##------------------------------------------------------------

# entropy.matrix <- computEntropy(input.data.discr) #----Verify the name

## Initialize filename where 'mi.net.adj.matrix.list' is to be saved
## in case 'clr.algo == CLR3'
mi.net.adj.matrix.list.filename <- NULL
if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
  rm(mi.net.adj.matrix.list.filename)
}

## source('learn_mi_net_struct.R')
if (clr.algo == 'CLR') {
  # mi.net.adj.matrix <- LearnMiNetStructZstat(mut.info.matrix, mi.net.adj.matrix, entropy.matrix, alpha)
  # mi.net.adj.matrix <- LearnMiNetStructClr(mut.info.matrix, mi.net.adj.matrix, num.nodes)
  mi.net.adj.matrix <- LearnClrNetMfi(mut.info.matrix, mi.net.adj.matrix, num.nodes, max.fanin, output.dirname)
  
} else if (clr.algo == 'CLR2') {
  mi.net.adj.matrix <- LearnClr2NetMfi(input.data.discr, num.nodes, node.names, num.timepts, 
                                       max.fanin, output.dirname, mi.net.adj.matrix)
  
} else if (clr.algo == 'CLR2.1') {
  mi.net.adj.matrix <- LearnClrNetMfiVer2.1(input.data.discr, num.nodes, node.names, num.timepts, 
                                            max.fanin, output.dirname, mi.net.adj.matrix)
  
} else if (clr.algo == 'CLR3') {
  mi.net.adj.matrix.list <- LearnClr3NetMfi(input.data.discr.3D, num.nodes, node.names, num.timepts, 
                                            max.fanin, mi.net.adj.matrix.list)
  
  ## Since 'mi.net.adj.matrix.list' is very large, save it in a specific file
  ## and remove it. Then load it when necessary. No need to retain it in the
  ## workspace even when not required.
  mi.net.adj.matrix.list.filename <- paste(output.dirname, 'mi.net.adj.matrix.list.RData', sep = '/')
  save(mi.net.adj.matrix.list, file = mi.net.adj.matrix.list.filename)
  rm(mi.net.adj.matrix.list)
}

elapsed.time <- (proc.time() - start.time) # Check time taken by CLR
writeLines('elapsed.time just after the CLR step= \n')
print(elapsed.time)
rm(elapsed.time)

if (clr.algo == 'CLR') {
  rm(mut.info.matrix)
}

if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
  # writeLines('\n mi.net.adj.matrix = \n')
  # print(mi.net.adj.matrix)
  save(mi.net.adj.matrix, file = paste(output.dirname, 'mi.net.adj.matrix.RData', sep = '/'))
  ## Check max number of nbrs for a node in the MI net
  # writeLines('\n Max number of nbrs = \n')
  # print(max(colSums(mi.net.adj.matrix)))
  
  ## Identify which nodes have the max number of nbrs
  # writeLines('\n Nodes having max number of nbrs = \n')
  # print(colnames(mi.net.adj.matrix[, colSums(mi.net.adj.matrix) == max(colSums(mi.net.adj.matrix))]))
}

# stop('MI net struct learning completed.')

# Rgraphviz::plot(as(mi.net.adj.matrix, 'graphNEL'))
##------------------------------------------------------------
## End: Learn MI Net Structure
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Learn Network Structures
##------------------------------------------------------------
unrolled.DBN.adj.matrix.list <- NULL

if ((clr.algo == 'CLR') | (clr.algo == 'CLR2') | (clr.algo == 'CLR2.1')) {
  unrolled.DBN.adj.matrix.list <- learnDbnStructMo1Layer3dParDeg1_v2(input.data.discr.3D, mi.net.adj.matrix, 
                                                                     num.discr.levels, num.nodes, num.timepts, 
                                                                     max.fanin, node.names, clr.algo)
  rm(mi.net.adj.matrix)
  
} else if (clr.algo == 'CLR3') {
  
  num.time.ivals <- (num.timepts - 1)
  unrolled.DBN.adj.matrix.list <- vector(mode = 'list', length = num.time.ivals)
  # rm(num.time.ivals)
  
  ## Adjacency matrix for a time-interval-specific DBN
  time.ival.spec.dbn.adj.matrix <- matrix(0, nrow = num.nodes, 
                                          ncol = num.nodes, 
                                          dimnames = c(list(node.names), 
                                                       list(node.names)))  
  for (time.ival.idx in 1:num.time.ivals) {
    unrolled.DBN.adj.matrix.list[[time.ival.idx]] <- time.ival.spec.dbn.adj.matrix
  }
  rm(time.ival.idx)
  
  rm(num.time.ivals, time.ival.spec.dbn.adj.matrix)
  
  unrolled.DBN.adj.matrix.list <- LearnDbnStructMo1Clr3Ser(input.data.discr.3D, mi.net.adj.matrix.list.filename, 
                                                           num.discr.levels, num.nodes, num.timepts, max.fanin, 
                                                           node.names, unrolled.DBN.adj.matrix.list)
  rm(mi.net.adj.matrix.list.filename)
}

save(unrolled.DBN.adj.matrix.list, file = paste(output.dirname, 'unrolled.DBN.adj.matrix.list.RData', sep = '/'))
rm(input.data.discr.3D)

## Learn the rolled DBN adj matrix
## source(paste(init.path, 'rollDbn.R', sep = '/'))
# rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, roll.method, allow.self.loop)
# rolled.DBN.adj.matrix <- rollDbn(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix, 'any', FALSE)
rolled.DBN.adj.matrix <- rollDbn_v2(num.nodes, node.names, num.timepts, unrolled.DBN.adj.matrix.list, 
                                    'any', allow.self.loop)
di.net.adj.matrix <- rolled.DBN.adj.matrix
rm(rolled.DBN.adj.matrix)

# writeLines('\n di.net.adj.matrix = \n')
# print(di.net.adj.matrix)
## Change the node names back to the original node names
rownames(di.net.adj.matrix) <- orig.node.names
colnames(di.net.adj.matrix) <- orig.node.names
save(di.net.adj.matrix, file = paste(output.dirname, 'di.net.adj.matrix.RData', sep = '/'))

## Create an '.sif' file equivalent to the directed net adjacency matrix
## that is readable in Cytoscape.
adjmxToSif(di.net.adj.matrix, output.dirname)
# rm(unrolled.DBN.adj.matrix)
rm(unrolled.DBN.adj.matrix.list)
##------------------------------------------------------------
## End: Learn Network Structures
##------------------------------------------------------------

## If the true rolled network is known a prior, then evaluate the performance
## metrics of the predicted rolled network.
if (true.net.filename != '')
{
  predicted.net.adj.matrix <- di.net.adj.matrix
  
  ## Loads R obj 'true.net.adj.matrix'
  load(true.net.filename)
  
  ## Begin: Create the format for result
  Result <- matrix(0, nrow = 1, ncol = 11)
  colnames(Result) <- list('TP', 'TN', 'FP', 'FN', 'TPR', 'FPR', 'FDR', 'PPV', 'ACC', 'MCC',  'F')
  # ## End: Create the format for result
  
  ResultVsTrue <- calcPerfDiNet(predicted.net.adj.matrix, true.net.adj.matrix, Result, num.nodes)
  rm(Result)
  writeLines('Result TGS vs True = \n')
  print(ResultVsTrue)
  rm(ResultVsTrue)
}

rm(di.net.adj.matrix)

elapsed.time <- (proc.time() - start.time) # Stop the timer
writeLines('elapsed.time = \n')
print(elapsed.time)
rm(elapsed.time)

sink()

##------------------------------------------------------------
## Begin: Save R session info in a File
##------------------------------------------------------------
sink(paste(output.dirname, 'sessionInfo.txt', sep = '/'))
sessionInfo()
sink()
##------------------------------------------------------------
## End: Save R session info in a File
##------------------------------------------------------------

##------------------------------------------------------------
## End: Main program
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: References
##------------------------------------------------------------
## 
##------------------------------------------------------------
## End: References
##------------------------------------------------------------