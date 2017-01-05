#!/bin/bash

## This file is intended to be run using sudo

VERSION=$(grep version ../starseqr_utils/__init__.py | cut -d "=" -f2 | tr -d " \t\n\r\"\'"s)

cmd0='sudo service docker start'
echo $cmd0
eval $cmd0

cmd1='sudo docker build --rm --no-cache -t eagenomics/starseqr:${VERSION} -t eagenomics/starseqr:latest --file Dockerfile .'
echo $cmd1
eval $cmd1


