mkdir Sarani
mkdir biocomputing ; cd biocomputing
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna > wildtype.fna ; wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk > wildtype.gbk ; wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk > wildtype.gbk
mv wildtype.fna /home/sarani_admin/Sarani
rm wildtype.gbk
grep -i tatatata wildtype.fna
grep -i tatatata wildtype.fna.1 > mutant_lines.fna
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NM_001323776.1&db=nuccore&report=fasta&retmode=text" -O sequence.fasta
echo $(( $(wc -l < sequence.fasta) - 2 ))
sed 1d sequence.fasta | grep -o -i 'A' | wc -l
sed 1d sequence.fasta | grep -o -i 'G' | wc -l
sed 1d sequence.fasta | grep -o -i 'C' | wc -l
sed 1d sequence.fasta | grep -o -i 'T' | wc -l
sudo apt install bc

gc_content=$(sed 1d sequence.fasta | grep -o -i '[GC]' | wc -l)

total_bases=$(sed 1d sequence.fasta | grep -o -i '[AGTC]' | wc -l)

gc_percentage=$(echo "scale=2; $gc_content / $total_bases * 100" | bc)

echo "GC content: $gc_percentage%"

cp sequence.fasta Sarani.fasta
num_A=$(sed 1d sequence.fasta | grep -o -i 'A' | wc -l)
num_G=$(sed 1d sequence.fasta | grep -o -i 'G' | wc -l)
num_T=$(sed 1d sequence.fasta | grep -o -i 'T' | wc -l)
num_C=$(sed 1d sequence.fasta | grep -o -i 'C' | wc -l)
echo "Number of A: $num_A" >> Sarani.fasta
echo "Number of G: $num_G" >> Sarani.fasta
echo "Number of T: $num_T" >> Sarani.fasta
echo "Number of C: $num_C" >> Sarani.fasta 

git clone https://github.com/Sarani-NS/HackBio
cd HackBio/
cd Task3
mkdir -p output
cd ../
mv /c/Users/sante/Sarani/Sarani.fasta ./output/
git add output/Sarani.fasta
git commit -m "Ajout du fichier Sarani.fasta"
git push origin main

nano Sarani.sh 
