#!/bin/bash
# Extract the embedded definition file from each .sif into ./defs/

# Use apptainer if available, otherwise singularity
if command -v apptainer &> /dev/null; then CT=apptainer; else CT=singularity; fi

mkdir -p defs

for sif in *.sif; do
    name="${sif%.sif}"
    if $CT inspect --deffile "$sif" > "defs/${name}.def" 2> /dev/null && [ -s "defs/${name}.def" ]; then
        echo "OK      ${name}.def"
    else
        echo "EMPTY   ${name}  (no embedded def - check manually)"
        rm -f "defs/${name}.def"
    fi
done

# Flag things to review before publishing
echo -e "\n--- Bootstrap sources ---"
grep -H -i "^Bootstrap:\|^From:" defs/*.def

echo -e "\n--- Defs with a %files section (need extra files on GitHub) ---"
grep -l "%files" defs/*.def
