**Profile Hidden Markov Model for the Kunitz-type protease inhibitor domain**  

The principal aim of this project is to built an Hidden Markov Model for the Kunitz-type protease inhibitor domain based on structural informations.  
An Hidden Markov Model (HMM) is a statistical model used to describe systems that are probabilistic and time-dependent. The system being modeled is assumed to be a Markov process with hidden (unobservable) states, that is the probability of the next state only depends on the actual state, not on the sequence of previous states. 

**Tools required for developing this model**  
1.  Install cd-hit with the commands sudo apt upgrade
    sudo apt install cd-hit
2.                     
**Online sources required for download the files**
- [Protein Data Bank (PDB)](https://www.rcsb.org/)
- [Uniprot](https://www.uniprot.org/)

**Pipeline of the project**  
1. Download the required files from PDB and uniprot for extracting the structural information of the proteins 
2. Create a custom report on PDB in order to have an overview of data and eventually make adjustments on them
3. Use cd-hit for cleaning redundance data 
