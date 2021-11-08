#!/bin/bash

binpath=$1
outdir=$2
SRR=$3

cmd=fastq-dump

cmd=$binpath$cmd
echo $cmd

$cmd -O $outdir --split-files $line