# Profile Hidden Markov Model for the Kunitz-type protease inhibitor domain

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
Perform a MSA with PDBeFold and remove also the protein 5JBT that resulted too short.
Perform again the MSA with the correct ID list of 23 proteins (pdb_list.txt).
Create an empty file to store the Multiple Sequence Alignment, copy there the information and transform everything in uppercase and store it in the file pdb_uppercase.ali
```
touch pdb_correct_msa.ali
mv correct_msa.fasta pdb_correct_msa.ali
awk '/^>/ { print; next } { print toupper($0) }' pdb_correct_msa.ali > pdb_uppercase.ali
```
## 4. Hidden Markov Model build 
Build an HMM with the command
```
hmmbuild pdb_kunitz_nr.hmm pdb_uppercase.ali
```
Create a logo in Skylign to better visualize the output ![Skylign logo](URL dell'immagine)
## 5. Database Creation 
To create 






