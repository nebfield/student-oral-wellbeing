#!/usr/bin/env awk -f

BEGIN{ OFS="\t" }

function basename(file) {
    sub(".*/", "", file)
    return file
}

{
split($8,fastq,";"); print $4, "$PWD/" basename(fastq[1]), "$PWD/" basename(fastq[2])
}
