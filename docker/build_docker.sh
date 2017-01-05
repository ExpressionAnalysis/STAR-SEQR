#!/bin/bash

VERSION=$(grep version ../starseqr_utils/__init__.py | cut -d "=" -f2 | tr -d " \t\n\r\"\'"s)

cmd1='docker build -t eagenomics/star-seqr:${VERSION} --rm --file Dockerfile .'
echo $cmd1
eval $cmd1

cmd2='docker build -t eagenomics/star-seqr:latest --rm --file Dockerfile .'
echo $cmd2
eval $cmd2

