#!/bin/zsh

mkdir logs
find -name '*log' | xargs cp -t logs/