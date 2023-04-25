#!/bin/bash

cat << EOF > .config.json
{
    "DATAPATH": "$PWD/data",
    "RESULTSPATH": "$PWD/results"
}
EOF

mkdir -p data results

echo setup complete!
