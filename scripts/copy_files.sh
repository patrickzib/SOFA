#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    exit 1
fi

rm -rf ~/MESSI_SFA_logs/$1/$2
mkdir -p ~/MESSI_SFA_logs/$1/$2
mv ~/MESSI_logs/* ~/MESSI_SFA_logs/$1/$2
