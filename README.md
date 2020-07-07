# Conformational_Search_Workflow
Automated conformational search workflow

This workflow is designed to take any input structure and produce the lowest energy conformer. 
In short, the program leverages CREST conformational search to produce conformers, re-ranks their energies based on DFT, then fully optimizes the 10 lowest energy conformers. 

IMAGE HERE


BASIC USAGE:

    python3 autoconf.py <INPUT.xyz>    - sets up workflow structure
    
    bash INPUT/start.sh                - submits workflow
    
    
OUTPUT:

    Creates a directory named after the INPUT containing 
       - CREST, ORCA, and lowest10 directories where calculations are performed
       - start.sh used to begin calculations
       
       after running
       
       - lowest energy optimized conformer
       - INPUT-optimized-energies.txt containing energies of 10 optimized conformers
    
    
ADDITIONAL USAGE OPTIONS:



DEPENDENCIES

    CREST - 
    Gaussian - 
    ORCA - 
    

    
    
    

       
