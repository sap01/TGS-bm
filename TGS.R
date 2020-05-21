#!/usr/bin/env Rscript
#' Benchmark the TGS algorithm
#'
#' @section Coding style:
#' \url{https://github.com/sap01/coding_practices/blob/master/R_coding_practices.md}
#'
#' @section How to execute this script:
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

library(TGS)
##------------------------------------------------------------
## End: Load the Required Libraries
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Read Command-line Arguments
##------------------------------------------------------------
input_args <- commandArgs(trailingOnly = TRUE)

if (length(input_args) != 1) {
  stop("Exactly one input file must be supplied.",
       call. = FALSE)
}

init_path <- getwd()

input_dirname <- paste(init_path, 
                       'asset',
                       sep = '/')

json_file <- paste(input_dirname,
                   input_args,
                   sep = '/')

rm(input_args)
##------------------------------------------------------------
## End: Read Command-line Arguments
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Create the output directory
##------------------------------------------------------------

## Create output dir name
output_dirname <-
  paste('output', format(Sys.time(), "%Y%m%d%H%M%S"), sep = '')

## Create output dir path
output_dirname <- paste(init_path, 
                        'asset', 
                        output_dirname, 
                        sep = '/')

if (.Platform$OS.type == 'windows') {
    
    ## Convert directory path to canonical form for Windows.
    ## It raises the warning if the directory does not exist, which
    ## is expected. Therefore, please ignore the warning.
    output_dirname <-
      normalizePath(output_dirname, 
                    winslash = '\\', 
                    mustWork = NA)
    
    shell(paste('mkdir ', output_dirname, sep = ''),
          intern = TRUE,
          mustWork = TRUE)
} else if (.Platform$OS.type == 'unix') {
      
      system(paste('mkdir ', output_dirname, sep = ''))
}

##------------------------------------------------------------
## End: Create the output directory
##------------------------------------------------------------

##------------------------------------------------------------
## Begin: Read User-defined input Params
##------------------------------------------------------------
print(json_file)

if (.Platform$OS.type == 'windows') {

  ## Convert directory path to canonical form for Windows.
  ## It raises the warning if the directory does not exist, which
  ## is expected. Therefore, please ignore the warning.
  json_file <-
    normalizePath(json_file,
                  winslash = '\\',
                  mustWork = NA)
}
print(json_file)

input_params <- rjson::fromJSON(file = json_file)
rm(json_file)

input.data.filename <- input_params$input.data.filename
num.timepts <- input_params$num.timepts
true.net.filename <- input_params$true.net.filename
input.wt.data.filename <- input_params$input.wt.data.filename
is.discrete <- input_params$is.discrete
num.discr.levels <- input_params$num.discr.levels
discr.algo <- input_params$discr.algo
mi.estimator <- input_params$mi.estimator
apply.aracne <- input_params$apply.aracne
clr.algo <- input_params$clr.algo
max.fanin <- input_params$max.fanin
allow.self.loop <- input_params$allow.self.loop
scoring.func <- input_params$scoring.func

rm(input_params)
##------------------------------------------------------------
## End: Read User-defined input Params
##------------------------------------------------------------


## Run algorithm TGS
TGS::LearnTgs(
  isfile = 1,
  json.file = json_file,
  input.dirname = input_dirname,
  input.data.filename = input.data.filename,
  num.timepts = num.timepts,
  true.net.filename = true.net.filename,
  input.wt.data.filename = input.wt.data.filename,
  is.discrete = is.discrete,
  num.discr.levels = num.discr.levels,
  discr.algo = discr.algo,
  mi.estimator = mi.estimator,
  apply.aracne = apply.aracne,
  clr.algo = clr.algo,
  max.fanin = max.fanin,
  allow.self.loop = allow.self.loop,
  scoring.func = scoring.func,
  output.dirname = output_dirname
)

## Print the output dir name in 'nohup.out'
print('The output directory name is:')
print(output_dirname)
print('') ## to append a blank line

## Save console output in a file named 'output.txt' inside the output directory.
# output.filename <- paste(output_dirname, 'output.txt', sep = '/')
# output.file.conn <- file(output.filename, open = "wt")
# sink(output.file.conn)
