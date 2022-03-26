#!/bin/bash

in_dir=$1
out_dir=$2


if [[ -d $in_dir ]]; then
  if [[ ! -d $out_dir ]]; then
    echo "Making $out_dir"
    mkdir $out_dir
  fi
  python3 processData.py $in_dir $out_dir
else
  echo "$in_dir does not exist!"
fi
