#!/usr/bin/env bash

# This script is based in : https://github.com/rrwick/Bacsort

# This script downloads assemblies from NCBI and their metadata. 
# It also uses Mash to do pairwise distances between all assemblies in each genus.


input=$1
mash_sketch_size=10000
threads=4

printf "\n"
echo "Checking if Mash is installed."
echo "......................................................................."
if command -v mash; then
    printf "\n"
    echo "Mash is installed!"
else
    echo "Error: Mash is not installed"
    exit 1
fi
printf "\n"


while IFS= read -r line; 
do
    printf "\n"
    echo "Downloading genomes from "$line
    echo "......................................................................."
    mkdir -p Genomes/$line
    cd Genomes/$line
    ncbi-genome-download --verbose --genera $line --metadata-table data.tsv --format fasta -p $threads --retries 100 bacteria

    if test -n "$(find refseq/bacteria -name '*.fna.gz' -print -quit)"
        then
        for f in refseq/bacteria/*/*.fna.gz; do
            new_name=$(echo $f | grep -oP "\w{3}_\d{9}\.\d" | head -n 1)".fna.gz"
            printf "\n"
            echo "mv "$f" "$new_name
            printf "\n"
            mv $f $new_name
        done
        rm -r refseq
        printf "\n\n"

        echo "Finding pairwise distances for "$line
        echo "......................................................................."
        mash sketch -p $threads -o mash -s $mash_sketch_size *.fna.gz
        mash dist -p $threads mash.msh mash.msh > mash_distances
        cd ../..
    else
        echo "###### WARNING ######"
        echo "No assemblies downloaded for "$line
        cd ..
        rm -r $line
        cd ..
    fi
    printf "\n"

done < "$input"
