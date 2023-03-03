#!/usr/bin/env bash
#this script runs render.R for you and, if necessary, will specify the right path to pandoc. 
#Update the following line to set the PATH properly if you get an error about an older pandoc version.
#usage: this_script.sh output_format
#example: render.sh pdf


output_format=$1
export PATH=/Applications/RStudio.app/Contents/MacOS/pandoc:$PATH
Rscript --vanilla render.R $output_format
