#!/bin/csh 

set infile = "inputfile"

set record = $argv[1]

sed -i "/$record/d" $infile
