# Conformational-Search-Workflow
Automated conformational search workflow

This workflow is designed to take any input structure and produce the lowest energy conformer. 
In short, the program leverages CREST conformational search to produce conformers, re-ranks their energies based on DFT, then fully optimizes the 10 lowest energy conformers. 

![Image of Workflow](workflow.png)

BASIC USAGE:

    python3 autoconf.py <INPUT.xyz>    - sets up workflow structure *
    
    bash INPUT/start.sh                - submits workflow
    
        *if the input file is not an xyz file, the code will attempt to use openbabel to convert to xyz
    
    
OUTPUT:

    Creates a directory named after the INPUT containing 
       - CREST, ORCA, and lowest10 directories where calculations are performed
       - start.sh used to begin calculations
       
       after running
       
       - status.txt containing calculation progress
       - lowest energy optimized conformer
       - INPUT-optimized-energies.txt containing energies of 10 optimized conformers
    
    
ADDITIONAL USAGE OPTIONS:

    >> python3 autoconf.py -h
    
    <<Usage:

       AUTOCONF: A conformational search script based on CREST
       python3 autoconf.py [inputfile] arg1 arg2...


       This code will setup the file structure and input necessary to run a CREST conformational search, followed by SP quick 5 step relaxation in ORCA, and optimization of the 10 lowest energy conformers by guassian, orca, or none. Method options can be set with variables at beginning of this script

    An input xyz should be provided, however, obabel will be used to attempt to convert other molecule file types.

    If the -c --constrain option is chosen, the file format is:

        <cat test.c

        >distance: 5, 4, auto
        >angle: 4, 3, 10, auto
        >dihedral: 5, 6, 7, 8, auto


    The type of constraint : atom#s for constraint, followed by the word "auto". This "auto" automatically sets the force constant for the constraint.



    Options:
      -h, --help            show this help message and exit
      
      -c, --constrain       If the -c --constrain option is chosen, a file with the same basename as your
                            input.xyz with the .c extension should be in the directory.        
                                Ex if input is test.xyz, constraint file should be test.c
                            
      -q CHARGE, --charge=CHARGE
                            Default: 0
                        
      -m MULTIPLICITY, --multiplicity=MULTIPLICITY
                            Default: 1
                        
      -o OPTPROGRAM, --opt=OPTPROGRAM
                            Availible programs are gaussian and orca, or none
                        
      -n, --nci             Specialized NCI run for CREST for non-covalent complexes
                        
      --sp                  Do only ORCA SP for energy rankings
      
      --nog16fail           Automatically resubmit failed gaussian jobs



DEPENDENCIES

    CREST - runs off CREST code located at /work/lopez/xtb/
    Gaussian - runs off g16 version in g16root=/work/lopez/
    ORCA - runs off ORCA_EXE=/work/lopez/orca_4_2_1_linux_x86-64_shared_openmpi216/
