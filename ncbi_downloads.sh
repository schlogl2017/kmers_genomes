echo "Downloading genomes from NCBI"

input="bact_taxa.txt"

while IFS= read -r line
do
  mkdir $line
  cd $line
  echo "Downloading $line genomes frinputom NCBI"
  ncbi-genome-download --genera $line bacteria -l "complete,chromosome" -F fasta --flat-output -p 4 -r 10
  cd ..
done < "$input"


