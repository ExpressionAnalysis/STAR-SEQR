#!/bin/bash
# this file is intended to be run using sudo

cmd0='service docker start'
echo $cmd0
eval $cmd0

cmd_clean='docker rm -f $(docker ps -a -q)'
echo $cmd_clean
eval $cmd_clean

cmd_clean2='docker rmi -f $(docker images -q)'
echo $cmd_clean2
eval $cmd_clean2



