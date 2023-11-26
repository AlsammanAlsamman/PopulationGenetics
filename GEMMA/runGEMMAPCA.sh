#!/bin/bash

scriptUsage() {
    echo "This script is used to run ****"
    echo "Usage: ./Command.sh *prameters*"
}

# if no input, print usage
if [ $# -eq 0 ]; then
    scriptUsage
    exit 1
fi

# run GEMMA PCA

