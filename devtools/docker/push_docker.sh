#!/bin/bash

docker login --username=jasper1918

VERSION=$(grep version ../starseqr_utils/__init__.py | cut -d "=" -f2 | tr -d " \t\n\r\"\'"s)

cmd1="docker push eagenomics/starseqr:${VERSION}"
echo $cmd1
eval $cmd1

cmd2="docker push eagenomics/starseqr:latest"
echo $cmd2
eval $cmd2
