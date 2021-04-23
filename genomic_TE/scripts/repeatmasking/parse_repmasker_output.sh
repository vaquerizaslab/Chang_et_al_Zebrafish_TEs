#!/usr/bin/env bash

# Must be run from directory containing .out or .align files.
# The original .fa file must also be present in the directory.

parsermdir="/workdir/jnw72/Software/Parsing-RepeatMasker-Outputs"

nohup perl ${parsermdir}/parseRM.pl \
    --in danRer11.nonalt.fa.align \
    --parse \
    --fa \
    --land 60,0.25 &

