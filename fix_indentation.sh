#!/usr/bin/env bash
# Usage: ./fix_indentation.sh Snakefile

# Check input
if [ -z "$1" ]; then
  echo "Usage: $0 <file>"
  exit 1
fi

file="$1"

# Replace all tab characters with 4 spaces
# The $'\t' is a Bash literal for a tab
sed -i 's/'$'\t''/    /g' "$file"
