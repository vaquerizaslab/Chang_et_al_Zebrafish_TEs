#!/usr/bin/bash

datadir="/local/workdir/jnw72/Projects/drerio-tes/data"

RepeatMasker -pa 60 \
    -a \
    -s \
    -nolow \
    -gccalc \
    -gff \
    -cutoff 200 \
    -no_is \
    -species danio \
    -dir $datadir/repeatmaskerout \
    $datadir/genome/danRer11.nonalt.fa
