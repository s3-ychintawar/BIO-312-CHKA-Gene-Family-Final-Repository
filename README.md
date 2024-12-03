# CHKA Gene Family Repository
# Welcome to the CHKA Gene Family Repository. This repository is dedicated to the analysis of the CHKA gene family, which plays a critical role in phospholipid metabolism and has implications in neurodevelopmental disorders and cancer.

# Lab Overview and Aims:
# - **Lab 3: Homolog Identification Using BLAST**
# Identify homologous sequences for the CHKA gene across proteomes using BLAST. This lab focuses on retrieving and filtering homologous sequences based on high-confidence metrics.
# - **Lab 4: Sequence Alignment**
# Align homologous sequences of the CHKA gene family to identify conserved regions and variations. The alignment provides the foundation for downstream phylogenetic analysis.
# - **Lab 5: Phylogenetic Tree Construction**
# Build and visualize a phylogenetic tree to study evolutionary relationships among CHKA homologs. This lab incorporates rooting and visual representation techniques.
# - **Lab 6: Gene-Species Tree Reconciliation**
# Reconcile gene and species trees to infer duplication and loss events, providing insights into the evolutionary history of the CHKA gene family.
# - **Lab 8: Protein Domain Prediction and Mapping**
# Predict conserved protein domains in CHKA homologs using RPS-BLAST and map them onto phylogenetic trees to analyze functional and evolutionary patterns.

---

# Lab 3: Homolog Identification Using BLAST
# Navigate to the CHKA working directory
cd ~/lab03-$MYGIT/CKA
# Sets the working directory where all BLAST-related files will be generated and stored. This step is repeated for all repositories needed for organization. 

# Download the reference CHKA protein sequence
ncbi-acc-download -F fasta -m protein "NP_001268.2"
# Retrieves the CHKA protein sequence in FASTA format from the NCBI database to use as the query for BLAST searches.

# Perform a standard BLASTP search
blastp -db ../allprotein.fas -query NP_001268.2.fa -outfmt 0 -max_hsps 1 -out CKA.blastp.typical.out
# Runs a BLASTP search to find homologous sequences in the proteome database.
# Output: CKA.blastp.typical.out

# Perform a detailed BLASTP search with custom output fields
blastp -db ../allprotein.fas -query NP_001268.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out CKA.blastp.detail.out
# Performs a more detailed BLASTP search that includes specific fields such as percent identity, alignment length, and bit score for finer analysis.
# Output: CKA.blastp.detail.out

# Inspect the detailed BLASTP results
less -S CKA.blastp.detail.out
# Opens the detailed BLASTP output in a scrollable view for verification of alignment quality and homolog detection.

# Count the number of matches to Homo sapiens
grep -c H.sapiens CKA.blastp.detail.out
# Counts the number of homologous sequences that belong to Homo sapiens within the detailed results.

# Extract and count unique species
grep -o -E "^[A-Z]\.[a-z]+" CKA.blastp.detail.filtered.out | sort | uniq -c
# Identifies and counts the unique species represented in the filtered BLAST results, summarizing homolog distribution.

# Filter BLASTP results by e-value (< 1e-30)
awk '{if ($6 < 1e-30) print $1}' CKA.blastp.detail.out > CKA.blastp.detail.filtered.out
# Selects high-confidence homologs based on their e-value, filtering out less reliable matches.
# Output: CKA.blastp.detail.filtered.out

# Count the number of filtered homologs
wc -l CKA.blastp.detail.filtered.out
# Counts the total number of homologs that passed the filtering criteria for high confidence.

# Summarize unique species in filtered results
grep -o -E "^[A-Z]\.[a-z]+" CKA.blastp.detail.filtered.out | sort | uniq -c
# Gives a breakdown of the filtered homologs by species, showing their taxonomic distribution.

# Lab 3 Outputs
# - NP_001268.2.fa: Reference CHKA protein sequence.
# - CKA.blastp.typical.out: Standard BLASTP results.
# - CKA.blastp.detail.out: Detailed BLASTP results including custom output fields.
# - CKA.blastp.detail.filtered.out: Filtered homolog sequences (e-value < 1e-30).


# Lab 4: Sequence Alignment
# Extract homologous sequences based on BLAST results
seqkit grep --pattern-file ~/lab03-$MYGIT/CKA/CKA.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/CKA/CKA.homologs.fas
#  Extracts sequences from the proteome database based on BLAST results, excluding unwanted species (e.g., "carpio").
# Output: CKA.homologs.fas

# Perform multiple sequence alignment (MSA) with MUSCLE
muscle -align ~/lab04-$MYGIT/CKA/CKA.homologs.fas -output ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas
#  Aligns the extracted homolog sequences to identify conserved regions across species using the MUSCLE tool.
# Output: CKA.homologs.al.fas

# Inspect the alignment file
less -RS ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas
# Opens the alignment file in a scrollable format, allowing manual inspection of the alignment quality.

# Generate alignment statistics
alignbuddy -pi ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | awk ' (NR>2) { for (i=2;i<=NF;i++) {sum+=$i;num++} } END { print(100*sum/num) } '
# Calculates alignment statistics such as the percentage identity between sequences to evaluate alignment quality.

# Trim poorly aligned regions
alignbuddy -trm all ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | alignbuddy -al
# Removes poorly aligned regions from the sequences to improve overall alignment quality.

# Generate a color-coded alignment visualization
alv -kil -w 100 ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | aha > ~/lab04-$MYGIT/CKA/CKA.homologs.al.html
# Produces a browser-compatible visualization of the alignment, highlighting conserved regions in color for easy review.
# Output: CKA.homologs.al.html

# Lab 4 Outputs
# - CKA.homologs.fas: Extracted homolog sequences from the proteome database.
# - CKA.homologs.al.fas: Multiple sequence alignment (FASTA format).
# - CKA.homologs.al.html: Color-coded alignment visualization in HTML format.

# Lab 5: Phylogenetic Tree Construction
# Objective: Build a phylogenetic tree of CHKA homologs to study evolutionary relationships.


# Prepare the aligned sequences for tree construction
sed 's/ /_/g' ~/lab04-$MYGIT/CKA/CKA.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas
# Cleans sequence headers by replacing spaces with underscores and removes duplicate labels to ensure compatibility with tree-building software.
# Output: CKA.homologsf.al.fas

# Construct a phylogenetic tree using IQ-TREE
iqtree -s ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas -bb 1000 -nt 2
# Generates a maximum-likelihood phylogenetic tree with 1,000 bootstrap replicates to evaluate the reliability of tree branches.
# Output: CKA.homologsf.al.fas.treefile 

# Plot the unrooted tree using an R script
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas.treefile ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas.treefile.pdf 0.4 15
# Visualizes the unrooted tree and saves it as a PDF for publication or further analysis.
# Output: CKA.homologsf.al.fas.treefile.pdf 

# Reroot the tree at its midpoint
gotree reroot midpoint -i ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile
# Reroots the phylogenetic tree at its midpoint, creating a more balanced visual representation of evolutionary relationships.
# Output: CKA.homologsf.al.mid.treefile 

# Generate a high-quality visualization of the midpoint-rooted tree
nw_order -c n ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.svg -
convert ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.svg ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.pdf

# Creates both SVG and PDF visualizations of the midpoint-rooted tree for high-quality graphical representation.
# Outputs:
# - CKA.homologsf.al.mid.treefile.svg
# - CKA.homologsf.al.mid.treefile.pdf 

# Generate a topology-only visualization of the midpoint-rooted tree
nw_order -c n ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/CKA/CKA.homologsf.al.midCl.treefile.svg -
convert ~/lab05-$MYGIT/CKA/CKA.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/CKA/CKA.homologsf.al.midCl.treefile.pdf
# Produces a simplified topology-only version of the tree for specific analytical purposes.
# Outputs:
# - CKA.homologsf.al.midCl.treefile.svg
# - CKA.homologsf.al.midCl.treefile.pdf

# Reroot the tree using specific outgroup sequences
nw_reroot ~/lab05-$MYGIT/CKA/CKA.homologsf.al.fas.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta > ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile
# Reroots the tree using specified outgroup sequences to provide biologically meaningful evolutionary insights.
# Output: CKA.homologsf.outgroupbeta.treefile 

# Visualize the outgroup-rooted tree
nw_order -c n ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile.svg -
convert ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile.svg ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile.pdf
# Creates high-quality SVG and PDF visualizations of the outgroup-rooted phylogenetic tree for analysis and publication.
# Outputs:
# - CKA.homologsf.outgroupbeta.treefile.svg 
# - CKA.homologsf.outgroupbeta.treefile.pdf

ab 5 Outputs
# - CKA.homologsf.al.fas: Cleaned and aligned sequences for tree construction.
# - CKA.homologsf.al.fas.treefile: Maximum-likelihood phylogenetic tree file in Newick format.
# - CKA.homologsf.al.mid.treefile: Midpoint-rooted tree file.
# - CKA.homologsf.al.mid.treefile.svg: Midpoint-rooted tree visualization in SVG format.
# - CKA.homologsf.al.mid.treefile.pdf: Midpoint-rooted tree visualization in PDF format.
# - CKA.homologsf.al.midCl.treefile.svg: Simplified topology-only visualization in SVG format.
# - CKA.homologsf.al.midCl.treefile.pdf: Simplified topology-only visualization in PDF format.
# - CKA.homologsf.outgroupbeta.treefile: Outgroup-rooted tree file.
# - CKA.homologsf.outgroupbeta.treefile.svg: Outgroup-rooted tree visualization in SVG format.
# - CKA.homologsf.outgroupbeta.treefile.pdf: Outgroup-rooted tree visualization in PDF format.


# Lab 6: Gene-Species Tree Reconciliation
# Objective: Reconcile the gene tree with the species tree to infer duplication and loss events.

# Create a new directory for Lab 6
mkdir ~/lab06-$MYGIT/CKA
# Initializes a dedicated directory for files and outputs related to reconciliation analysis. Ensures organization of results.

# Copy the midpoint-rooted tree file
cp ~/lab05-$MYGIT/CKA/CKA.homologsf.al.mid.treefile ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile
# Uses the midpoint-rooted tree generated in Lab 5 as input for reconciliation analysis. This provides continuity between labs.
# Output: CKA.homologsf.al.mid.treefile.

# Run Notung for reconciliation
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/CKA
# Compares the gene tree against the species tree to identify evolutionary events like duplications and losses. Outputs include reconciled trees and event annotations.
# Output: Reconciled tree files and event annotations in ~/lab06-$MYGIT/CKA

# Extract reconciled species tree information
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.png | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
# Extracts and visualizes the reconciled species tree. The `nw_display` command provides an immediate graphical representation

# Convert the reconciled tree to RecPhyloXML format
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.ntg --include.species
# Converts the reconciled tree to a standardized XML format suitable for additional visualization and downstream analysis.
# Output: CKA.homologsf.al.mid.treefile.rec.ntg.xml

# Annotate the reconciled tree
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.svg
# Annotates duplication and loss events directly onto the reconciled tree. Enhances interpretability by marking key evolutionary events.
# Output: CKA.homologsf.al.mid.treefile.rec.svg

# Convert the annotated tree to a high-quality PDF
convert -density 150 ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/CKA/CKA.homologsf.al.mid.treefile.rec.pdf
# Produces a publication-ready PDF version of the annotated tree for inclusion in reports or scientific manuscripts.
# Output: CKA.homologsf.al.mid.treefile.rec.pdf

# Lab 6 Outputs
# - Reconciled tree files: Generated by Notung during reconciliation.
# - CKA.homologsf.al.mid.treefile.rec.ntg: Raw reconciled tree file in Newick format.
# - CKA.homologsf.al.mid.treefile.rec.svg: Annotated reconciled tree visualization in SVG format.
# - CKA.homologsf.al.mid.treefile.rec.pdf: Annotated reconciled tree visualization in PDF format.


# Lab 8: Protein Domain Prediction and Mapping
# Objective: Predict conserved protein domains in CHKA homologs using RPS-BLAST and map them onto the phylogenetic tree for evolutionary insights.

# Clean homolog sequences for compatibility with RPS-BLAST
sed 's/*//' ~/lab04-$MYGIT/CKA/CKA.homologs.fas > ~/lab08-$MYGIT/CKA/CKA.homologs.fas
# Removes asterisks (*) from the sequence file to ensure compatibility with RPS-BLAST. Prepares the homolog sequences for accurate domain prediction.
# Output: CKA.homologs.fas

# Run RPS-BLAST to predict conserved domains
rpsblast -query ~/lab08-$MYGIT/CKA/CKA.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/CKA/CKA.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
# Uses the Pfam database to identify conserved protein domains in CHKA homologs. Outputs detailed information about domain matches
# Output: CKA.rps-blast.out.

# Copy the outgroup-rooted tree file
cp ~/lab05-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/CKA
# Reuses the outgroup-rooted tree from Lab 5 to provide a base for domain mapping.
# Output: CKA.homologsf.outgroupbeta.treefile

# Map conserved domains onto the phylogenetic tree
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/CKA/CKA.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/CKA/CKA.rps-blast.out ~/lab08-$MYGIT/CKA/CKA.tree.rps.pdf
# Integrates domain prediction data with the phylogenetic tree to provide a comprehensive view of evolutionary patterns.
# Output: CKA.tree.rps.pdf

# Inspect the RPS-BLAST output
mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | tail -n +2 | less -S
# Displays the RPS-BLAST results in a tabular format for manual inspection and verification.

# Count the number of hits per sequence
cut -f 1 ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | sort | uniq -c
#  Quantifies the number of domain matches for each sequence, summarizing conservation across sequences.

# Summarize unique domain types
cut -f 6 ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | sort | uniq -c
# Counts the occurrences of each unique domain type, identifying patterns of functional conservation.

# Calculate domain lengths and rank them
awk '{a=$4-$3; print $1, "\t", a}' ~/lab08-$MYGIT/CKA/CKA.rps-blast.out | sort -k2nr
#  Calculates the length of each identified domain and sorts them by size, prioritizing larger domains for further analysis.

# Extract specific fields from the RPS-BLAST output
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/CKA/CKA.rps-blast.out
# Extracts sequence IDs and their associated domain descriptions for focused downstream analysis.

# Lab 8 Outputs
# - CKA.homologs.fas: Cleaned homolog sequences.
# - CKA.rps-blast.out: Detailed domain prediction results from RPS-BLAST.
# - CKA.homologsf.outgroupbeta.treefile: Outgroup-rooted tree reused for domain mapping.
# - CKA.tree.rps.pdf: Phylogenetic tree annotated with conserved domains.

