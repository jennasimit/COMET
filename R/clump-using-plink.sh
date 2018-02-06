#!/usr/local/bin/bash

output1File=$1
RPfile=$2
chr=$3
output2File=$1"-abf"


#plink --noweb \
/software/gapi/pkg/plink/1.0.7/bin/plink --noweb \
--bfile $RPfile \
--clump $output1File \
--clump-field neg.abf \
--clump-p1 0 \
--clump-p2 0 \
--clump-r2 0.10 \
--clump-kb 500 \
--out $output2File

fout=$output2File".clumped"

R CMD BATCH '--args $fout $chr' format-plink-output.R format-plink-output.Rout


gzip $fout






