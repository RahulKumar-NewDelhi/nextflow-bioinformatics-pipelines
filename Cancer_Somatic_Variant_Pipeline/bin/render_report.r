#!/usr/bin/env Rscript

library(stringr)

args = commandArgs(trailingOnly=TRUE)

render_report = function(sample_name, rmd_file, rmd_params) {
  rmarkdown::render( rmd_file, params = rmd_params, output_file = paste0(sample_name, ".html"))
}

get_params_list <- function(sample_name, args) {
  params = c("sample_name", str_replace(args[seq(1,length(args)) %% 2 == 1], "--", ""))
  value = c(sample_name, args[seq(1, length(args)) %% 2 == 0])
  names(value) <- params
  as.list(value)
}

sample_name <- args[1]
rmd_file <- args[2]
tail_args <- tail(args, -2)

names(sample_name) <- c("sample_name")
rmd_params <- if (length(tail_args) != 0 ) get_params_list(sample_name, tail_args) else as.list(sample_name)

render_report(sample_name, rmd_file, rmd_params)
