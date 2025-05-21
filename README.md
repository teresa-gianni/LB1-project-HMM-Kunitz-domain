project_pipeline.md
# Profile Hidden Markov Model for the Kunitz-type Protease Inhibitor Domain (Pfam: PF00014)

The primary goal of this project is to build a **Hidden Markov Model (HMM)** for the Kunitz-type protease inhibitor domain based on structural information. An HMM is a statistical model used to describe systems that are probabilistic and time-dependent. The system being modeled is assumed to be a Markov process with hidden (unobservable) states, meaning the probability of the next state only depends on the actual state, not on the sequence of previous states.

---

## Required Resources

### Tools to Install Locally

* **CD-HIT:**
    ```bash
    sudo apt update
    sudo apt install cd-hit
    ```
* **HMMER:**
    ```bash
    sudo apt update
    sudo apt install hmmer
    ```
* **Chimera:** (Optional) Download from [https://www.cgl.ucsf.edu/chimera/download.html](https://www.cgl.ucsf.edu/chimera/download.html)

### Online Tools

* **Protein Data Bank (PDB):** [https://www.rcsb.org/](https://www.rcsb.org/) - To extract protein sequences containing the Kunitz domain for model training.
* **UniProt:** [https://www.uniprot.org/](https://www.uniprot.org/) - To extract a dataset of proteins containing the Kunitz domain for model testing.
* **MUSCLE:** [https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein](https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein) - Used for generating a Multiple Sequence Alignment (MSA).
* **PDBeFold:** [https://www.ebi.ac.uk/msd-srv/ssm/](https://www.ebi.ac.uk/msd-srv/ssm/) - Used to generate a Multiple Structural Alignment.
* **WebLogo:** [https://weblogo.berkeley.edu/logo.cgi](https://weblogo.berkeley.edu/logo.cgi) - Used to generate a logo from the Multiple Sequence Alignment.
* **Skylign:** [https://skylign.org/](https://skylign.org/) - Used to generate a logo from the HMM model.

---

## Project Pipeline

### 1. Data Collection

This phase involves downloading protein sequences from UniProt and PDB.

* **Download from UniProt:**
    * **Human proteins:**
        ```
        (ft_domain:Kunitz) AND (reviewed:true) AND (xref:pfam-PF00014) AND (organism_id:9606)
        ```
        Output: `human_kunitz.fasta` (18 sequences)
    * **Non-human proteins:**
        ```
        (ft_domain:Kunitz) AND (reviewed:true) AND (xref:pfam-PF00014) NOT (organism_id:9606)
        ```
        Output: `nonhuman_kunitz.fasta` (377 sequences)
    * **All Kunitz sequences:**
        ```
        (ft_domain:Kunitz) AND (reviewed:true) AND (xref:pfam-PF00014)
        ```
        Output: `all_kunitz.fasta` (395 sequences)
* **Download the entire SwissProt database:**
    * `uniprot_sprot.fasta`
* **Create a Custom Report on PDB:**
    * Query:
        ```
        Data Collection Resolution <= 3.5 AND ( Identifier = "PF00014" AND Annotation Type = "Pfam" ) AND Polymer Entity Sequence Length <= 80 AND Polymer Entity Sequence Length >= 45
        ```
    * Containing the following information: `auth asym ID`, `ID`, `Sequence`, `Resolution`, `Entity ID`.
    * Generated file: `rcsb_pdb_custom_report_20250506075310.csv`
    * Extract and filter protein sequences using the command:
        ```bash
        cat rcsb_pdb_custom_report_20250506075310.csv | tr -d '"' | awk -F ',' '{if (length($2)>0) {name = $2}; print name,$4,$5,$6}' | grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz.fasta
        ```
        The sequences obtained from PDB are now in `pdb_kunitz.fasta`.

---

### 2. Data Cleaning

This phase aims to reduce redundancy in protein sequences.

* **Reduce redundancy with CD-HIT:**
    ```bash
    cd-hit -i pdb_kunitz.fasta -o pdb_kunitz.clst
    ```
    Output is saved in `pdb_kunitz.clst`.
* **Convert output to a new FASTA file:**
    ```bash
    mv pdb_kunitz.clst pdb_kunitz_nr.fasta
    ```
    The `pdb_kunitz_nr.fasta` file contains the non-redundant sequences.
* **Remove an excessively long sequence (2ODY_E):**
    ```bash
    awk '/^>/{f=($0 ~ /^>2ODY_E/)?1:0} !f' pdb_kunitz_nr.fasta > pdb_kunitz_senza_2ODY_E.fasta
    ```
    The sequence is now saved in `pdb_kunitz_senza_2ODY_E.fasta`.

---

### 3. Multiple Structural Alignment (MSA)

An MSA is performed with PDBeFold, and the file is prepared for HMMER.

* **Execute MSA with PDBeFold:**
    * Remove protein `5JBT` found to be too short from a Multiple Sequence Alignment with MUSCLE, resulting in `pdb_kunitz_nr_23.fasta`.
    * Perform again the MSA with the correct ID list of 23 proteins (`pdb_list.txt`).
    * Download in FASTA format the sequences of the 23 proteins relatives to a Multiple Sequence Alignment derived from the MSA
* **Prepare the file for HMMER:**
    * Create an empty file for the MSA and store there the sequences obtained from the previous step :
        ```bash
        touch pdb_correct_msa.ali
        mv correct_msa.fasta pdb_correct_msa.ali
        ```
    * Convert sequences to uppercase:
        ```bash
        awk '/^>/ { print; next } { print toupper($0) }' pdb_correct_msa.ali > pdb_uppercase.ali
        ```
    The `pdb_uppercase.ali` file is now ready for building the model.

---

### 4. Hidden Markov Model Build

The HMM is built using HMMER.

* **Build the HMM:**
    ```bash
    hmmbuild pdb_kunitz_nr.hmm pdb_uppercase.ali
    ```
    The model is saved in `pdb_kunitz_nr.hmm`.
* **Create a logo in Skylign:** For better output visualization.

---

### 5. Database Creation

This section describes the creation of positive and negative sets for model testing.

### Creating the Positive Set

* **Eliminate sequences with high similarity:**
    * Create a BLAST database with all UniProt sequences:
        ```bash
        makeblastdb -in all_kunitz.fasta -input_type fasta -dbtype prot -out all_kunitz.fasta
        ```
    * Query PDB sequences against the newly created database to identify those with >= 95% identity over at least 50 residues:
        ```bash
        blastp -query pdb_kunitz_nr_23.fasta -db all_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
        grep -v "^#" pdb_kunitz_nr_23.blast | awk '$3 >= 95 && $4 >= 50 {print $2}' | sort -u > redundant_ids.txt
        ```
    * Identify IDs of sequences to keep:
        ```bash
        grep "^>" all_kunitz.fasta | cut -d' ' -f1 | cut -c2- > all_kunitz_fasta.ids
        mv redundant_ids.txt redundant.ids
        comm -23 <(sort all_kunitz_fasta.ids) <(sort redundant.ids) > to_keep.ids
        ```
    * Remove high-identity sequences using the Python script `get_seq.py`:
        ```bash
        python3 get_seq.py to_keep.ids all_kunitz.fasta > ok_kunitz.fasta
        ```
        The remaining sequences are in `ok_kunitz.fasta`.

### Creating the Negative Set

* **Eliminate Kunitz domain-containing sequences from the entire SwissProt database:**
    * Create a list of all SwissProt IDs:
        ```bash
        grep ">" uniprot_sprot.fasta | cut -d "|" -f2 > sp.ids
        ```
    * Remove Kunitz sequence IDs:
        ```bash
        comm -23 <(sort sp.ids) <(sort all_kunitz_fasta.ids) > negs.ids
        ```
    * Retrieve sequences corresponding to the remaining IDs with `get_seq.py`:
        ```bash
        python3 get_seq.py negs.ids uniprot_sprot.fasta > negs.fasta
        ```

---

### 6. Model Testing

This phase involves randomly splitting the sets and running `hmmsearch`.

* **Randomly divide the sets:**
    ```bash
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
* **Evaluate the model with `hmmsearch`:**
    ```bash
    hmmsearch -Z 1000 --max --tblout pos_1.out pdb_kunitz_nr.hmm pos_1.fasta
    hmmsearch -Z 1000 --max --tblout pos_2.out pdb_kunitz_nr.hmm pos_2.fasta
    hmmsearch -Z 1000 --max --tblout neg_1.out pdb_kunitz_nr.hmm neg_1.fasta
    hmmsearch -Z 1000 --max --tblout neg_2.out pdb_kunitz_nr.hmm neg_2.fasta
    ```
* **Extract and format required information:**
    ```bash
    grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_1.class
    grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_2.class
    grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > neg_1.class
    grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > neg_2.class
    ```
* **Manually re-add some solutions to negative sets (if necessary):**
    ```bash
    comm -23 <(sort neg_1.ids) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_1.class
    comm -23 <(sort neg_2.ids) <(cut -f 1 neg_2.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_2.class
    ```

---

### 7. Performance Evaluation

This final phase calculates the model's performance.

* **Concatenate positive and negative sets:**
    ```bash
    cat pos_1.class neg_1.class > set_1.class
    cat pos_2.class neg_2.class > set_2.class
    ```
* **Evaluate model performance across a range of E-value thresholds using the script `performance.py`:**
    ```bash
    for i in $(seq 1 12); do python3 performance.py set_1.class 1e-$i; done
    for i in $(seq 1 12); do python3 performance.py set_2.class 1e-$i; done
    ```
    
