#!/bin/bash

pdffile="${1}"
shift 1
imagefiles=(${*})
infiles="${*}"
convert $infiles tmp.mng
convert tmp.mng $pdffile
rm tmp.mng

