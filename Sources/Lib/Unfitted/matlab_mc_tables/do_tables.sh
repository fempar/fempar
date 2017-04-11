#!/bin/bash
MATLAB="/usr/local/MATLAB/R2015b/bin/matlab"
##MATLAB="/opt/MATLAB/R2015b/bin/matlab"
$MATLAB -nodesktop -nodisplay -nosplash -nojvm < do_tables.m
