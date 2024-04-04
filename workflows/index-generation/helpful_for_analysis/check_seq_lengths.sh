#!/bin/bash

# Check if a file name was provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 filename"
    exit 1
fi

# Initialize variables
sequence=""
sequence_name=""

# Function to print sequence length
print_length() {
    if [ -n "$sequence_name" ]; then
        echo "$sequence_name: ${#sequence}"
    fi
}

# Read the file line by line
while IFS= read -r line || [ -n "$line" ]; do
    if [[ $line == '>'* ]]; then
        # Print the length of the previous sequence
        print_length
        # Start a new sequence
        sequence_name=$line
        sequence=""
    else
        # Accumulate the sequence
        sequence+=$line
    fi
done < "$1"

# Print the length of the last sequence
print_length
