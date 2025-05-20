# Profile Hidden Markov Model for the Kunitz-type protease inhibitor domain (Pfam: PF00014)

The principal aim of this project is to built an Hidden Markov Model for the Kunitz-type protease inhibitor domain based on structural informations.  
An Hidden Markov Model (HMM) is a statistical model used to describe systems that are probabilistic and time-dependent. The system being modeled is assumed to be a Markov process with hidden (unobservable) states, that is the probability of the next state only depends on the actual state, not on the sequence of previous states.   

## Resources required for developing this model
### To be installed on your machine:
* CD-HIT with the commands
     ```
     sudo apt upgrade                     
     sudo apt install cd-hit
* HMMER with the commands
     ```
     sudo apt upgrade                     
     sudo apt install hmmer
* [Chimera](https://www.cgl.ucsf.edu/chimera/download.html) (optional)
### Online tools
* [Protein Data Bank (PDB)](https://www.rcsb.org/) to extract protein sequences containing the Kunitz domain for training the model 
* [Uniprot](https://www.uniprot.org/) to extract a dataset of protein containig the Kunitz domain for testing the model 
* [MUSCLE](https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein) used for generating a Multiple Sequence Alignment
* [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/) used to generate a Multiple Structural Alignment
* [WebLogo](https://weblogo.berkeley.edu/logo.cgi) used to generate a logo from the Multiple Sequence Alignment
* [Skylign](https://skylign.org/) used to generate a logo from the HMM model 

## Pipeline of the project
### 1. Data collection
* Download from Uniprot the sequences of human and not-human proteins containing the Kunitz domain
  
  **Human proteins**
  
  Query:
  ```
  (ft_domain:Kunitz) AND (reviewed:true) AND (xref:pfam-PF00014) AND (organism_id:9606) 
  ```
  Output in file: human_kunitz.fasta (18 sequences)
  
  **Not human proteins**
  
  Query:
  ```
  (ft_domain:Kunitz) AND (reviewed:true) AND (xref:pfam-PF00014) NOT (organism_id:9606) 
  ```
  Output in file: not_human_kunitz.fasta (377 sequences)
  
  All kunitz sequences
  
  Query:
   ```
   (ft_domain:Kunitz) AND (reviewed:true) AND (xref:pfam-PF00014) 
   ```
  Output in file: all_kunitz.fasta (395 sequences)
* Download the entire SwissProt database of proteins (uniprot_sprot.fasta).
* Create a custom report on PDB with the query:
  ```
  Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
  ```
  containing the following information:
   - auth asym ID
   - ID
   - Sequence
   - Resolution
   - Entity ID
  
  Stored in the file rcsb_pdb_custom_report_20250506075310.csv
  
  From this report it is possible to extract the protein sequences of interest and filtering out    the ones containing multiple times the domain with the command
  ```
  cat rcsb_pdb_custom_report_20250506075310.csv |tr -d '"'| awk -F ',' '{if (length($2)>0) {name = $2}; print name,$4,$5,$6}' |grep PF00014 |awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz.fasta
  ```
  The sequences obtained from PDB are now stored in the file pdb_kunitz.fasta
## 2. Data Cleaning
Use CD-HIT to reduce redundancy with the command 
`cd-hit -i pdb_kunitz.fasta -o pdb_kunitz.clst` that will store the output in the file pdb_kunitz.clst.
Convert this output to a new fasta file (pdb_kunitz_nr.fasta) which will contain the non redundant sequences with the command ```mv pdb_kunitz.clst pdb_kunitz_nr.fasta```.
Remove a sequence that was found to be too long from the Multiple Sequence Alignment in MUSCLE with the command
```
awk '/^>/{f=($0 ~ /^>2ODY_E/)?1:0} !f' pdb_kunitz_nr.fasta > pdb_kunitz_senza_2ODY_E.fasta
```
and store it in a new file pdb_kunitz_senza_2ODY_E.fasta.  
## 3. Multiple Structural Alignment (MSA)
Perform a MSA with PDBeFold and remove also the protein 5JBT that resulted too short (pdb_kunitz_nr_23.fasta). 
Perform again the MSA with the correct ID list of 23 proteins (pdb_list.txt).
Create an empty file to store the Multiple Sequence Alignment, copy there the information and transform everything in uppercase and store it in the file pdb_uppercase.ali
```
touch pdb_correct_msa.ali
mv correct_msa.fasta pdb_correct_msa.ali
awk '/^>/ { print; next } { print toupper($0) }' pdb_correct_msa.ali > pdb_uppercase.ali
```
## 4. Hidden Markov Model build 
Build an HMM (pdb_kunitz_nr.hmm) with the command
```
hmmbuild pdb_kunitz_nr.hmm pdb_uppercase.ali
```
Create a logo in Skylign to better visualize the output ![Skylign logo](URL dell'immagine)
## 5. Database Creation 
To create a benchmark set for testing the model:  

### Creation of the positive set

First eliminate the sequences with high similarity to the PDB sequences used to train the model  
Create a BLAST database with all the sequences extracted from Uniprot 
```
makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta
```
Query the PDB sequences against the database just created and select the ones with identity grater than or equal to 95% in at least 50 residues and store it in the file redundant_ids.txt. 
```
blastp -query pdb_kunitz_nr_23.fasta -db all_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
grep -v "^#" pdb_kunitz_nr_23.blast |awk '$3 >= 95 && $4 >= 50 {print $2}' | sort -u > redundant_ids.txt
```
Get the IDs of the sequences to keep in the database 
```
grep "^>" all_kunitz.fasta | cut -d' ' -f1 | cut -c2- > all_kunitz_fasta.ids
mv redundant_ids.txt redundant.ids
comm -23 <(sort all_kunitz_fasta.ids) <(sort redundant.ids) > to_keep.ids
```
Remove the high identity sequences form the database using the python script get_seq.py and store the right sequences in a new file ok_kunitz.fasta
```
python3 get_ids.py to_keep.ids all_kunitz.fasta > ok_kunitz.fasta
```
### Creation of the negative set  
To create this set is first necessary to eliminate from the entire SwissProt database the sequences containing the kkunitz domain.  
Create a list of all the SwissProt IDs and delete from it the sequences extracted from Uniprot
```
grep ">" uniprot_sprot.fasta | cut -d "|" -f2 > sp.ids
comm -23 <(sort sp.ids) <(sort all_kunitz_fasta.ids) >negs.ids
```
Retrieve the sequences relative to the IDs with the same python script used before 
```
python3 get_seq.py negs.ids uniprot_sprot.fasta > negs.fasta
```
## 6. Model testing
Randomly divide the positive and negative sets into two subsets  
```
sort -R negs.ids > random_sp_negs.ids
sort -R to_keep.ids > random_ok_kunitz.ids
head -n 183 random_ok_kunitz.ids > pos_1.ids
tail -n 183 random_ok_kunitz.ids > pos_2.ids
head -n 286417 random_sp_negs.ids > neg_1.ids
tail -n 286417 random_sp_negs.ids > neg_2.ids
python3 get_seq.py pos_1.ids uniprot_sprot.fasta > pos_1.fasta
python3 get_seq.py pos_2.ids uniprot_sprot.fasta > pos_2.fasta
python3 get_seq.py neg_1.ids uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py neg_2.ids uniprot_sprot.fasta > neg_2.fasta
```
Perform the evaluation of the model in each set and store in an output file only the information needed
```
hmmsearch -Z 1000 --max --tblout pos_1.out pdb_kunitz_nr.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out pdb_kunitz_nr.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out pdb_kunitz_nr.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out pdb_kunitz_nr.hmm neg_2.fasta
grep -v "^#" pos_1.out | awk '{split($1,a,"\|"); print a[2],1,$5,$8}' |tr " " "\t" > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"\|"); print a[2],1,$5,$8}' |tr " " "\t" > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"\|"); print a[2],0,$5,$8}' |tr " " "\t" > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"\|"); print a[2],0,$5,$8}' |tr " " "\t" > neg_2.class
```
Readd manually some solution to te negative sets
```
comm -23 <(sort neg_1.ids) <(cut -f 1 neg_1.class |sort) |awk '{print $1"\t0\t10.0\t10.0"}' >> neg_1.class
comm -23 <(sort neg_2.ids) <(cut -f 1 neg_2.class |sort) |awk '{print $1"\t0\t10.0\t10.0"}' >> neg_2.class
```
## 7. Performance evaluation
Conctenate positive and negative in two sets
```
cat pos_1.class neg_1.class > set_1.class
cat pos_2.class neg_2.class > set_2.class
```
Evaluate the model performance in a range of E-value thresholds with the python script performance.py
```
for i in $(seq 1 12); do  python3 performance.py set_1.class 1e-$i; done
for i in $(seq 1 12); do  python3 performance.py set_2.class 1e-$i; done
```
