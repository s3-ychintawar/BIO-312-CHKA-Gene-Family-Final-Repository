# BIO-312-CHKA-Gene-Family-Final-Repository

# Lab: Bioinformatics Workflow for CHKA Gene Analysis
This README provides step-by-step instructions for analyzing CHKA gene homologs. The labs include homolog identification, sequence alignment, phylogenetic tree construction, and protein domain mapping.

---

## Lab 3: Homolog Identification Using BLAST
### Commands
```bash
# Navigate to the CHKA working directory
cd ~/lab03-$MYGIT/CKA

# Download the reference CHKA protein sequence
ncbi-acc-download -F fasta -m protein "NP_001268.2"

# Perform a standard BLASTP search
blastp -db ../allprotein.fas -query NP_001268.2.fa -outfmt 0 -max_hsps 1 -out CKA.blastp.typical.out

# Perform a detailed BLASTP search with custom output fields
blastp -db ../allprotein.fas -query NP_001268.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out CKA.blastp.detail.out

# Inspect the detailed BLASTP results
less -S CKA.blastp.detail.out

# Count the number of matches to Homo sapiens
grep -c H.sapiens CKA.blastp.detail.out

# Extract and count unique species
grep -o -E "^[A-Z]\.[a-z]+" CKA.blastp.detail.filtered.out | sort | uniq -c

# Filter BLASTP results by e-value (< 1e-30)
awk '{if ($6 < 1e-30) print $1}' CKA.blastp.detail.out > CKA.blastp.detail.filtered.out

# Count the number of filtered homologs
wc -l CKA.blastp.detail.filtered.out

# Summarize unique species in filtered results
grep -o -E "^[A-Z]\.[a-z]+" CKA.blastp.detail.filtered.out | sort | uniq -c
---

## Lab 4: Sequence Alignment
### Commands
```bash
# Extract homologous sequences from the proteome database
seqkit grep --pattern-file ~/lab03-$MYGIT/CKA/CKA.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/CKA/CKA.homologs.fas

# Perform multiple sequence alignment (MSA) with MUSCLE
muscle -align ~/lab04-$MYGIT/CKA/CKA.homologs.fas -output ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas

# Inspect the alignment file
less -RS ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas

# Generate alignment statistics
alignbuddy -pi ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | awk ' (NR>2) { for (i=2;i<=NF;i++) {sum+=$i;num++} } END { print(100*sum/num) } '

# Trim poorly aligned regions
alignbuddy -trm all ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | alignbuddy -al

# Generate a color-coded alignment visualization
alv -kil -w 100 ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | aha > ~/lab04-$MYGIT/CKA/CKA.homologs.al.html
---

## Lab 5: Phylogenetic Tree Construction
### Commands
```bash
# Create and navigate to the working directory for Lab 5
mkdir ~/lab05-$MYGIT/CKA
cd ~/lab05-$MYGIT/CKA

# Prepare aligned sequences for tree construction
sed 's/ /_/g' ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas

# Construct a maximum-likelihood phylogenetic tree with bootstrap analysis
iqtree -s ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas -bb 1000 -nt 2

# Reroot the tree at its midpoint
gotree reroot midpoint -i ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile

# Generate visualization files for the phylogenetic tree
nw_order -c n ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.svg -
convert ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.svg ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.pdf
---

## Lab 6: Gene-Species Tree Reconciliation
### Commands
```bash
# Create a working directory for Lab 6
mkdir ~/lab06-$MYGIT/CKA

# Copy the midpoint-rooted tree file
cp ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile

# Perform reconciliation of the gene tree with the species tree using Notung
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/CKA

# Extract reconciled species tree and visualize it
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.png | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -

# Convert reconciled tree to RecPhyloXML format
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.ntg --include.species

# Annotate the reconciled tree
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.svg

# Convert annotated tree to a high-quality PDF
convert -density 150 ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.pdf
---

## Lab 8: Protein Domain Prediction and Mapping
### Commands
```bash
# Clean homolog sequences for compatibility with RPS-BLAST
sed 's/*//' ~/lab04-$MYGIT/CKA/CKA.homologs.fas > ~/lab08-$MYGIT/CKA/CKA.homologs.fas

# Run RPS-BLAST to predict conserved domains
rpsblast -query ~/lab08-$MYGIT/CKA/CKA.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/CKA/CKA.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001

# Copy the outgroup-rooted tree file
cp ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/CKA

# Map conserved domains onto the phylogenetic tree
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/CKA/CKA.rps-blast.out ~/lab08-$MYGIT/CKA/CKA.tree.rps.pdf

# Inspect RPS-BLAST output
mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | tail -n +2 | less -S

# Count hits per sequence
cut -f 1 ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | sort | uniq -c

# Count unique domain types
cut -f 6 ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | sort | uniq -c

# Calculate domain lengths
awk '{a=$4-$3; print $1, "\t", a}' ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | sort -k2nr

# Extract specific fields from RPS-BLAST output
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/CKA/CKA.rps-blast.out
---
