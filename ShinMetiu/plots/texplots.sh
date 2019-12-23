#!/bin/bash

# Figure 1
name="SM_setup"
pdflatex $name.tex
rm -f "$name"-inc.pdf "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log

# Figure 2
name="mask_example"
pdflatex $name.tex
rm -f "$name"-inc.pdf "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log

# Figure 3
name="mask_uen_grid2"
pdflatex $name.tex
rm -f "$name"-inc.pdf "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log

# Figure 4
name="mask_dt"
pdflatex $name.tex
rm -f "$name"-inc.pdf "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log

# Figure 5
name="mask_gcoc_grid"
pdflatex $name.tex
rm -f "$name"-inc.pdf "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log

# Figure 6
name="linear_gcoc_grid"
pdflatex $name.tex
rm -f "$name"-inc.pdf "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log

# Figure 7
name="mask_gcoc_smalldt"
pdflatex $name.tex
rm -f "$name"-inc.pdf	 "$name".aux "$name".dvi "$name".log "$name"-inc.eps "$name".tex texput.log


