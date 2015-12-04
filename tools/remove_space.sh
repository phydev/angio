#!/bin/bash
IFS="\n"
for file in *.xyz;
do
    mv "$file" "${file//[[:space:]]}"
done
