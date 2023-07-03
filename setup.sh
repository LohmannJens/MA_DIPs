#!/bin/bash

cat << EOF > .config.json
{
    "DATAPATH": "$PWD/data",
    "RESULTSPATH": "$PWD/results"
}
EOF

mkdir -p data results
cd results
mkdir -p control_analysis deletion_length_and_position direct_repeats free_energy_estimations general_validation ML motif_discovery NP_density PR8 regression_length_vs_occurrence relative_occurrence_nucleotides

echo setup complete!
