#!/bin/bash

# Exclude dat/ods/pdf resource files, tex documentation files, and
# txt files that are not CMake related, svg and png images.
wc -l `find ./ ! \( -name "*.dat" -o -name "*.ods" -o -name "*.pdf" -o -name "*.tex" -o -name "*.bib" -o \( -name "*.txt" -a ! -name "CMakeLists.txt" \) -o -name "*.png" -o -name "*.svg" \)` | sort -n

