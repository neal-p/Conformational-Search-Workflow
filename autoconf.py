import os
import sys
import shutil
from optparse import OptionParser

#USER PARAMETERS

ORCAmethod='B3LYP/G D3BJ 6-31G(d) def2/J NOAUTOSTART cpcm(water)'
#NOTE: optimization will be added internally, please do not request here
ORCAcores=8
ORCAmem=63
ORCApartition='lopez,short'
ORCAtime='1-00:00:00'


CRESTcores=16
CRESTmem=63
CRESTtime='1-00:00:00'
CRESTpartition='short,lopez'
#NOTE: constraints will be added separatly to this CRESTmethod line if turned on

GAUSroute='opt freq=(noraman,savenormalmodes) scrf=(iefpcm,solvent=water) empiricaldispersion=GD3BJ'
OPTcores='16'
OPTmem='63'
GAUSmethod='B3LYP'
GAUSbasis='6-31G(d)'
OPTpartition='short'
ORCAoptmethod='pbe0 opt cc-pvdz def2/J rijcosx'

#Setup

constrainhelp="""If the -c --constrain option is chosen, a file with the same basename as your input.xyz with the .c extension should be in the directory. 
    

Ex if input is test.xyz, constraint file should be test.c

"""
usage = """

       AUTOCONF: A conformational search script based on CREST
       python3 %prog [inputfile] arg1 arg2...


This code will setup the file structure and input necessary to run a CREST conformational search, followed by SP quick 5 step relaxation in ORCA, and optimization of the 10 lowest energy conformers by guassian, orca, or none. Method options can be set with variables at beginning of this script

An input xyz should be provided, however, obabel will be used to attempt to convert other molecule file types.

If the -c --constrain option is chosen, the file format is:

<cat test.c

>distance: 5, 4, auto
>angle: 4, 3, 10, auto
>dihedral: 5, 6, 7, 8, auto


The type of constraint : atom#s for constraint, followed by the word "auto". This "auto" automatically sets the force constant for the constraint.

"""
parser=OptionParser(usage=usage)
parser.add_option("-c", "--constrain", action="store_true",dest="constraint", help=constrainhelp)
parser.add_option("-q", "--charge", action="store",type="string",dest="charge",default='0',help="Default: 0")
parser.add_option("-m", "--multiplicity", action="store",type="string",dest="multiplicity", default='1', help="Default: 1")
parser.add_option("-o", "--opt", action="store", type="string", dest="optprogram", default="gaussian", help="Availible programs are gaussian and orca, or none")
parser.add_option("-n", "--nci", action="store_true", dest="nci", help="Specialized NCI run for CREST for non-covalent complexes")
parser.add_option("--sp", action="store_true", dest="SP", help="Do only ORCA SP for energy rankings")
parser.add_option( "--nog16fail",dest="g16fail",action="store_false",help="Automatically resubmit failed gaussian jobs",default="True")

(options,args)=parser.parse_args()
#infile=options.infile
if not len(sys.argv) > 1:
    print("""
    AUTOCONF SEARCH WITH CREST:

Please provide a structure file, or use -h for help

""")
    exit()
elif len(args) == 0:
    print("""
    AUTOCONF SEARCH WITH CREST:

Please provide a structure file, or use -h for help

""")
    exit()
else:
    infile=args[0]
constraint=options.constraint
g16fail=options.g16fail
charge=options.charge
mult=options.multiplicity
optprogram=options.optprogram
nci=options.nci
SP=options.SP
title=infile.split('.')[0]
ext=infile.split('.')[1]

if len(title) > 20:
    print("input file has a title longer than 20 characters. I've found that CREST can behave wierd with long strings. A truncated name will be used in the CREST directory for all input files")

if ext == 'xyz':
    if not os.path.exists('{0}.xyz'.format(title)):
        print('{0}.xyz could not be found, or is empty'.format(title))
        exit()
    print("""xyz found: {0}

Setting up calculations...

Submit autoconf workflow by executing start.sh in the {1} dir
""".format(infile,title))
else:
    os.system('obabel {0} -o xyz -O {1}.xyz'.format(infile,title))
    if os.stat('{0}.xyz'.format(title)).st_size < 1:
        print('File could not be converted by obabel, please supply xyz')
        exit()
    else:
        print("""{0} converted to {1}.xyz
Setting up calculations...

Submit autoconf workflow by executing start.sh in the {1} dir
""".format(infile,title))

if optprogram == 'gaussian':
    pass
elif optprogram == 'orca':
    pass
elif optprogram == 'none':
    pass
else:
   print("""Unrecognized optimization program. Availible options are: 
gaussian
orca
""")

if constraint == True:
    constraintfile=False
    if os.path.exists('{0}.c'.format(title)):
        print('{0}.c file found'.format(title))
        constraintfile=True
        constrainttype='c'
    elif os.path.exists('{0}.constraints'.format(title)):
        print('{0}.constraints file found'.format(title))
        constraintfile=True
        constrainttype='constraints'
    if constraintfile == False:
        print('No constraint file found, looked for {0}.c and {0}.constraints'.format(title))
        exit()


calculationinformation="""
###################################################
CREST workflow setup with the following parameters:
###################################################

constraint={0}
g16fail={1}
charge={2}
mult={3}
optprogram={4}
nci={5}
SP={6}

ORCAmethod={7}
#NOTE: optimization will be added internally, don't worry if it isn't here
ORCAcores={8}
ORCAmem={9}
ORCApartition={10}
ORCAtime={11}

CRESTcores={12}
CRESTmem={13}
CRESTtime={14}
CRESTpartition={15}

GAUSroute={16}
OPTcores={17}
OPTmem={18}
GAUSmethod={19}
GAUSbasis={20}
OPTpartition={21}
ORCAoptmethod={22}

####################################################

""".format(constraint,g16fail,charge,mult,optprogram,nci,SP,ORCAmethod,ORCAcores,ORCAmem,ORCApartition,ORCAtime,CRESTcores,CRESTmem,CRESTtime,CRESTpartition,GAUSroute,OPTcores,OPTmem,GAUSmethod,GAUSbasis,OPTpartition,ORCAoptmethod)

work=os.getcwd()

try:
    os.mkdir('{0}/{1}'.format(work,title))
except FileExistsError:
    print('re-run')
workdir='{0}/{1}'.format(work,title)

with open('{0}/status.txt'.format(workdir),'w') as status:
    status.write(calculationinformation)

try:
    os.mkdir('{0}/{1}/CREST'.format(work,title))
except FileExistsError:
    print('re-run')
CRESTdir='{0}/{1}/CREST'.format(work,title)
try:
    os.mkdir('{0}/{1}/ORCA'.format(work,title))
except FileExistsError:
    print('re-run')

shutil.copyfile('{0}/{1}.xyz'.format(work,title),'{0}/{1}.xyz'.format(CRESTdir,title))

ORCAdir='{0}/{1}/ORCA'.format(work,title)
try:
    os.mkdir('{0}/{1}/lowest10'.format(work,title))
except FileExistsError:
    print('re-run')
lowest10dir='{0}/{1}/lowest10'.format(work,title)

def CREST(title,partition,cores,mem,time,workdir,constraint,nci):
    if len(title) > 20:
        tmptitle=title[0:18]
    else:
        tmptitle=title

    CRESTmethod='/work/lopez/xtb/crest {0}.xyz -xnam /work/lopez/xtb/xtb_6.2.3/bin/xtb -g H2O -ewin 500 -T 14 -opt crude -subrmsd --verbose'.format(tmptitle)

    if constraint == True:
        CRESTmethod=CRESTmethod + ' -cinp {0}.c '.format(tmptitle)
    if nci == True:
        CRESTmethod=CRESTmethod + '-nci '
    
    sbatch="""#!/bin/bash
#SBATCH --job-name={7}-CREST
#SBATCH --output=out.o
#SBATCH --error=out.e
#SBATCH --partition={1}
#SBATCH --nodes=1
#SBATCH --ntasks={2}
#SBATCH --mail-type=END
#SBATCH --mail-user=neal.pa@husky.neu.edu
#SBATCH --mem={3}G
#SBATCH --time={4}
hostname

ulimit -s unlimited
export OMP_STACKSIZE={3}G
export OMP_NUM_THREADS={2},1

export XTBPATH="/work/lopez/xtb/"
export XTBHOME=$XTBPATH
export OMP_MAX_ACTIVE_LEVELS=1
export LD_LIBRARY_PATH=${{LD_LIBRARY_PATH}}:${{XTBHOME}}/lib
export PYTHONPATH=${{PYTHONPATH}}:${{XTBHOME}}/python
export LD_LIBRARY_PATH="/work/lopez/orca_4_2_1_linux_x86-64_shared_openmpi216/":"/work/lopez/OpenBLAS/":$LD_LIBRARY_PATH

work={5}
cd $work
time=$(date)
echo "CREST started at $time" >> ../status.txt


cp {0}.xyz {7}.xyz
cp {0}.c {7}.c
cp {7}.xyz {7}.ref

{6} > {0}.out
normaltermination=$(grep "CREST terminated normally." -c {0}.out)
if [[ $normaltermination -gt 0 ]]
    then
    cp crest_conformers.xyz ../ORCA/{0}-all.xyz
    obabel ../ORCA/{0}-all.xyz -O ../ORCA/{0}-all-sorted-conf.xyz -m
    nstruct=$(ls -la ../ORCA/{0}-all-sorted-conf*.xyz |wc -l)
    if [[ $nstruct -gt 500 ]]
        then
        nstruct=$((500))
    fi
    sed -i "s/END/$nstruct/g" ../ORCA/*sbatch
    sbatch ../ORCA/{0}-ORCA.sbatch
else
    echo "CREST FAILED" >> ../status.txt
fi
""".format(title,partition,cores,mem,time,workdir,CRESTmethod,tmptitle)
    return sbatch

def ORCA(cores,time,title,partition,mem,workdir,method,charge,mult,SP):
    if len(title) > 20:
        tmptitle=title[0:5]
    else:
        tmptitle=title

    sbatch="""#!/bin/sh
## script for ORCA-4.2.1
#SBATCH --nodes=1
#SBATCH --ntasks={0}
#SBATCH --time={1}
#SBATCH --job-name={8}-ORCA
#SBATCH --partition={3}
#SBATCH --mem={4}Gb
#SBATCH --output=%j.o
#SBATCH --error=%j.e
#SBATCH --array=1-END%100

export WORKDIR={5}
export ORCA_EXE=/work/lopez/orca_4_2_1_linux_x86-64_shared_openmpi216/
export OPENMPI=/work/lopez/openmpi-2.1.6/
export LD_LIBRARY_PATH=$OPENMPI/lib:$ORCA_EXE:$LD_LIBRARY_PATH
export PATH=$OPENMPI/bin:$PATH

cd $WORKDIR
if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
    then
    time=$(date)
    echo "ORCA started at $time" >> ../status.txt
    nstruct=$(ls -la {2}-all-sorted-conf*.xyz |wc -l)
    charge={6}
    mult={7}
    inp=$(head -n 8 {2}.inp)
    if test -f {2}-conf$nstruct.inp
        then
        rm {2}-conf*.inp
    fi 
    for ((i=1;i<=nstruct;i++))
        do
        echo -e "${{inp/CHARGE MULTIPLICITY FILE/$charge $mult {2}-all-sorted-conf$i.xyz}}" > {2}-conf$i.inp
        echo {2}-conf$i.inp >> {2}-coms.txt
    done
    sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID ../lowest10/{2}-OPT.sbatch
    sbatch --dependency=afternotok:$SLURM_ARRAY_JOB_ID ORCAfailed.sbatch
else
    sleep 120s
fi

input={2}-conf${{SLURM_ARRAY_TASK_ID}}.inp
#input=$(sed "${{SLURM_ARRAY_TASK_ID}}q;d" {2}-coms.txt)
export INPUT=${{input%.*}}

$ORCA_EXE/orca $INPUT.inp > $INPUT.out
date >> $INPUT.out

converged=$(grep "SCF NOT CONVERGED AFTER" -c $INPUT.out)
unreliablestep=$(grep "HUGE, UNRELIABLE STEP WAS ABOUT TO BE TAKEN" -c $INPUT.out)
if [[ $unreliablestep -gt 0 ]] || [[ $converged -gt 0 ]]
    then
        sed -i 's#{9}#{9} Slowconv NOSOSCF DIIS#g' $INPUT.inp
        $ORCA_EXE/orca $INPUT.inp > $INPUT.out
        date >> $INPUT.out
        converged=$(grep "SCF NOT CONVERGED AFTER" -c $INPUT.out)
        if [[ $converged -gt 0 ]]
            then
            energies=$(grep "FINAL SINGLE POINT ENERGY" -c $INPUT.out)
            if [[ $energies -gt 0 ]]
                then
                echo "$INPUT did not converge last scf, but previous energy was taken" >> /scratch/autots-errors/{2}
        exit 0
            else
                echo "FINAL SINGLE POINT ENERGY     500" >> $INPUT.out
                exit 0
fi
fi
else
    exit 0
fi

""".format(cores,time,title,partition,mem,workdir,charge,mult,tmptitle,method)

    ORCAfailedscript="""!/bin/bash
#SBATCH --job-name={0}-ORCAfailed
#SBATCH --output=resubmit.o
#SBATCH --error=resubmit.e
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00:20:00

cd {1}
echo "ORCA FAILED" >> ../status.txt
echo "Check the following:" >> ../status.tst
for i in *out
    do 
    terminated=$(grep "ORCA finished by error termination" -c $i)
    if [[ $terminated -gt 0 ]]
        then echo "    $i" >> ../status.txt
    fi
done

""".format(tmptitle,workdir)
    
    with open('{0}/ORCAfailed.sbatch'.format(workdir,title), 'w') as ORCAfailed:
        actualmem=int(mem/cores)
        ORCAfailed.write(ORCAfailedscript)
    if SP == True:
        inputfile="""!{0}
%pal nprocs {1} end
%Maxcore {2}000
*xyzfile CHARGE MULTIPLICITY FILE

""".format(method,cores,actualmem,title)
    else: 
        inputfile="""!{0} Opt
%pal nprocs {1} end
%Maxcore {2}000
%geom
MaxIter 5
end
*xyzfile CHARGE MULTIPLICITY FILE

""".format(method,cores,actualmem,title)
    
    with open('{0}/{1}.inp'.format(workdir,title), 'w') as ORCAinput:
        ORCAinput.write(inputfile)

    return sbatch

def GAUS(title,GAUScores,GAUSmem,GAUSmethod,GAUSbasis,GAUSroute,partition,lowest10dir,charge,mult,g16fail):
    if len(title) > 20:
        tmptitle=title[0:5]
    else:
        tmptitle=title
    sbatch="""#!/bin/bash
#SBATCH --job-name={5}-GAUS
#SBATCH --output=out.o
#SBATCH --error=out.e
#SBATCH --partition={1}
#SBATCH --nodes=1
#SBATCH --ntasks={2}
#SBATCH --mem={3}G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-10
hostname

work={4}
cd $work

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
    then
    if test -f complete
        then
        rm complete
        rm energies-sorted.txt
        rm energies.txt
        sed -i '1,/{6} {7}/!d' {0}*.com 
    fi
    for i in ../ORCA/*out
        do
        energy=$(grep "FINAL SINGLE POINT ENERGY" -c $i)
        if [[ $energy -lt 1 ]]
            then
            echo "$i did not complete successfully, aborting run" >> ../status.txt
            scancel $SLURM_ARRAY_JOB_ID
            exit 1234
        fi
    done   

    sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID ../{0}-EX.sbatch
    sbatch --dependency=afternotok:$SLURM_ARRAY_JOB_ID {0}-OPTfail.sbatch

    nstruct=$(ls -la ../ORCA/{0}-all-sorted-conf*.xyz |wc -l)
    for ((x=1;x<=nstruct;x++)); do  
        energy=$(tac ../ORCA/{0}-conf$x.out | grep "FINAL SINGLE POINT ENERGY" -m1)
        echo $energy >> complete
        awk "/FINAL SINGLE POINT ENERGY/{{i++}}i==$x{{print ; exit}}" complete | awk '{{ print $5}}' >> energies.txt; done
    awk '{{print NR "   "  $s}}' energies.txt > output.txt
    sed -i 's/......$//' energies.txt
    sort -k2n -g -u output.txt > energies-sorted.txt
    
    for ((x=1;x<=nstruct;x++))
        do
        next=$(( x+1 ))
        energy1=$(sed "${{x}}q;d" energies-sorted.txt | awk '{{ print $2 }}')
        energy2=$(sed "${{next}}q;d" energies-sorted.txt | awk '{{ print $2 }}')
        diff=$(bc <<< "($energy1 - $energy2)*627.5")
        cutoff='0.01'
        if [ 1 -eq "$(echo "${{diff}} < ${{cutoff}}" | bc)" ]
            then
            #sed -i "${{next}}d" energies-sorted.txt
            echo "Possible redundant conformer $next in energoes-sorted.txt" >> ../status.txt
        fi
    done
            
    lowest10=$(awk  '{{print $1}} NR==10{{exit}}' energies-sorted.txt )
    count=0
    for b in $lowest10
        do
        count=$((count+1))
        obabel ../ORCA/{0}-conf$b.out -o xyz | tail -n +3  >> {0}-conf$count.com
        obabel ../ORCA/{0}-conf$b.out -o xyz -O {0}-conf$count.xyz
        echo " " >> {0}-conf$count.com
    done
else
    sleep 120s
fi

nstruct=$(cat energies-sorted.txt |wc -l)

if [[ ${{SLURM_ARRAY_TASK_ID}} -le $nstruct ]]
    then
    input=$(sed "${{SLURM_ARRAY_TASK_ID}}q;d" {0}-coms.txt)
    export INPUT=$input
    export WORKDIR=$work
    export GAUSS_SCRDIR=$work
    export g16root=/work/lopez/
    . $g16root/g16/bsd/g16.profile

    cd $WORKDIR
    $g16root/g16/g16 $INPUT

    station=$(grep "Station" -c ${{INPUT%.*}}.log)
    if [[ $station -lt 2 ]]
        then 
        exit 1234
    fi

else
    echo "no $SLURM_ARRAY_TASK_ID generated"
    rm {0}-conf$SLURM_ARRAY_TASK_ID.com
fi
 
""".format(title,partition,GAUScores,GAUSmem,lowest10dir,tmptitle,charge,mult)

    for i in range(1,11):
        inputfile="""%chk={0}-conf{6}.chk
%nprocs={1}
%mem={2}GB
# {3}/{4} {5}

autots script

{7} {8}
""".format(title,GAUScores,GAUSmem,GAUSmethod,GAUSbasis,GAUSroute,i,charge,mult)
        
        with open('{0}/{1}-conf{2}.com'.format(lowest10dir,title,i), 'w') as com:
            com.write(inputfile)
        if i == 1:
            with open('{0}/{1}-coms.txt'.format(lowest10dir,title), 'w') as coms:
                coms.write('{0}-conf{1}.com\n'.format(title,i))
        else:
            with open('{0}/{1}-coms.txt'.format(lowest10dir,title), 'a') as coms:
                coms.write('{0}-conf{1}.com\n'.format(title,i))


    GAUSfail=r"""#!/bin/bash
#SBATCH --job-name={4}-GAUSfailed
#SBATCH --output=resubmit.o
#SBATCH --error=resubmit.e
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00:20:00
hostname

work={1}
cd $work
touch failed-script-ran
if test -f {0}-resubmit.txt
    then
    nresub=$(sed "1q;d" {0}-resubmit.txt)
else
    nresub=0
fi
nresub=$((nresub+=1))
if [[ $nresub -lt 3 ]]
    then
    rm {0}-resubmit.txt
    echo $nresub >> {0}-resubmit.txt
    for i in {0}*log
        do
        finished=$(grep 'Station' -c $i)
        echo "$i $finished" >>{0}-resublog.txt
        if [[ $finished -lt 2 ]]
            then
            convert=$(obabel $i -o xyz)
            if [[ $convert == *"0 molecules converted"* ]]
                then
                echo "${{i%.*}} could not be converted, reverting to original" >> {0}-resublog.txt
                old=$(grep "%oldchk=" -c ${{i%.*}}.com)
                sed -i '1,/{2} {3}/!d' ${{i%.*}}.com
                echo " " ${{i%.*}}.com
                if [[ $old -gt 0 ]]
                    then
                    sed -i 1d ${{i%.*}}.com
                fi
                sed -i 's/geom=check guess=read//g' ${{i%.*}}.com
                sed -i 's/opt=readfc/opt=calcfc/g' ${{i%.*}}.com
                sed -i 's/opt /opt=calcfc /g' ${{i%.*}}.com
                end=${{i#-conf*}}
                index=${{end%.*}} 
                tail -n +3 {0}-all-sorted-conf$index.xyz >> ${{i%.*}}.com
                echo " " >>  ${{i%.*}}.com
            else
            sed -i '1,/{2} {3}/!d' ${{i%.*}}.com
            termination=$(grep "Normal termination" -c $i)
            cycles=$(grep "SCF Done" -c $i)
            check=$(grep "geom=check guess=read" -c ${{i%.*}}.com)
            old=$(grep "%oldchk=" -c ${{i%.*}}.com)
            
            #Non stationary point found
                 #read in fc
            if [[ $finished -eq 1 ]] && [[ $termination -ge 2 ]]
                then
                echo "${{i%.*}} Non-stationary point found, reading previous fc" >> {0}-resublog.txt
                echo " " >> ${{i%.*}}.com
                if [[ $old -gt 0 ]]
                    then
                    sed -i 1d ${{i%.*}}.com
                fi
                sed -i 's/geom=check guess=read//g' ${{i%.*}}.com
                mv ${{i%.*}}.chk ${{i%.*}}-readingfc.chk
                sed -i 's/opt /opt=readfc /g' ${{i%.*}}.com
                sed -i 's/opt=calcfc/opt=readfc/g' ${{i%.*}}.com
                sed -i 's/opt=readfc/opt=readfc geom=check guess=read/g' ${{i%.*}}.com
                sed -i "1s/^/%oldchk=${{i%.*}}-readingfc.chk\n/" ${{i%.*}}.com
                echo " " >> ${{i%.*}}.com

             #stationary found, but didnt finish frequencies or far from starting geometry
                 #re-calculate force constants
             elif ([[ $finished -eq 1 ]] && [[ $termination -lt 2 ]]) || ([[ $cycles -gt 15 ]])
                then
                echo "${{i%.*}} failed frequencies or is far from input geometry" >> {0}-resublog.txt
                if [[ $old -gt 0 ]]
                    then
                    sed -i 1d ${{i%.*}}.com
                fi
                obabel $i -o xyz | tail -n +3 >> ${{i%.*}}.com
                echo " " >> ${{i%.*}}.com
                sed -i 's/geom=check guess=read//g' ${{i%.*}}.com
                sed -i 's/opt /opt=calcfc /g' ${{i%.*}}.com
                sed -i 's/opt=readfc/opt=calcfc/g' ${{i%.*}}.com
                echo " " >> ${{i%.*}}.com

             #no stationary point found yet
             else
                fc=$(grep "Converged?" -c $i)
  
                #if force constants were finished computing, read them
                if [[ $fc -gt 0 ]]
                    then
                    echo "${{i%.*}} did not find stationary point, reading fc" >> {0}-resublog.txt
                    echo " " >> ${{i%.*}}.com
                    if [[ $old -gt 0 ]]
                        then
                        sed -i 1d ${{i%.*}}.com
                    fi
                    sed -i 's/geom=check guess=read//g' ${{i%.*}}.com
                    mv ${{i%.*}}.chk ${{i%.*}}-readingfc.chk
                    sed -i 's/opt=calcfc/opt=readfc/g' ${{i%.*}}.com
                    sed -i 's/opt /opt=readfc /g' ${{i%.*}}.com
                    sed -i 's/opt=readfc/opt=readfc geom=check guess=read/g' ${{i%.*}}.com
                    sed -i "1s/^/%oldchk=${{i%.*}}-readingfc.chk\n/" ${{i%.*}}.com
                    echo " " >> ${{i%.*}}.com
                 #if not, take geometry and restart
                 else
                    echo "${{i%.*}} did not find stationary point and has no fc" >> {0}-resublog.txt
                    obabel $i -o xyz |tail -n +3 >> ${{i%.*}}.com
                    echo " " >> ${{i%.*}}.com
                    if [[ $old -gt 0 ]]
                        then
                        sed -i 1d ${{i%.*}}.com
                    fi
                    sed -i 's/geom=check guess=read//g' ${{i%.*}}.com
                    sed -i 's/opt=readfc/opt=calcfc/g' ${{i%.*}}.com 
                    sed -i 's/opt /opt=calcfc /g' ${{i%.*}}.com
                    echo " " >> ${{i%.*}}.com
                fi
            fi
            fi
        echo ${{i%.*}}.com >> {0}-resubmit.txt
        echo " " >> ${{i%.*}}.com
    fi
    done
    toresub=$(cat {0}-resubmit.txt |wc -l)
    currentarray=$(sed "10q;d" {0}-OPT.sbatch)
    sed -i "s/$currentarray/#SBATCH --array=2-$toresub/g" {0}-OPT.sbatch
    sed -i "s/{0}-coms.txt/{0}-resubmit.txt/g" {0}-OPT.sbatch
    ID=$(sbatch --parsable {0}-OPT.sbatch)
    sbatch --dependency=afterok:$ID ../{0}-EX.sbatch
    sbatch --dependency=afternotok:$ID {0}-OPTfail.sbatch

else
    for i in {0}*log
        do
        finished=$(grep 'Station' -c $i)
        if [[ $finished -lt 2 ]]
            then 
            echo "Check $i optimization" >> ../status.txt
        fi
    done
fi


""".format(title,lowest10dir,charge,mult,tmptitle)

    if g16fail:
        with open('{0}/{1}-OPTfail.sbatch'.format(lowest10dir,title),'w') as failed:
            failed.write(GAUSfail)
    

    return sbatch 


def ORCAopt(title,OPTcores,OPTmem,OPTmethod,OPTpartition,lowest10dir,charge,mult):
    if len(title) > 20:
        tmptitle=title[0:5]
    else:
        tmptitle=title

    sbatch="""#!/bin/bash
#SBATCH --job-name={5}-OPT
#SBATCH --output=out.o
#SBATCH --error=out.e
#SBATCH --partition={1}
#SBATCH --nodes=1
#SBATCH --ntasks={2}
#SBATCH --mem={3}G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-10
hostname

work={4}
cd $work

if [[ $SLURM_ARRAY_TASK_ID == 1 ]]
    then
    nstruct=$(ls -la ../ORCA/{0}-all-sorted-conf*.xyz |wc -l)
    
    for ((x=1;x<=nstruct;x++)); 
        do 
        energy=$(tac ../ORCA/{0}-conf$x.out | grep "FINAL SINGLE POINT ENERGY" -m1)
        echo $energy >> complete
        awk "/FINAL SINGLE POINT ENERGY/{{i++}}i==$x{{print ; exit}}" complete | awk '{{ print $5}}' >> energies.txt; done
    awk '{{print NR "   "  $s}}' energies.txt > output.txt
    sort -k2n output.txt > energies-sorted.txt
    lowest10=$(awk  '{{print $1}} NR==10{{exit}}' energies-sorted.txt )
    count=0
    for b in $lowest10
        do
        count=$((count+1))
        cp ../ORCA/{0}-all-sorted-conf$b.xyz ./{0}-conf$count.xyz
        echo ../ORCA/{0}-all-sorted-conf$b.xyz
        echo {0}-conf$count.com
    done
    sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID ../{0}-EX.sbatch 
else
    sleep 120s
fi

nstruct=$(cat energies-sorted.txt |wc -l)

if [[ ${{SLURM_ARRAY_TASK_ID}} -le $nstruct ]]
    then
    input=$(sed "${{SLURM_ARRAY_TASK_ID}}q;d" {0}-coms.txt)
    export INPUT=${{input%.*}}
    export WORKDIR={4}
    export ORCA_EXE=/work/lopez/orca_4_2_1_linux_x86-64_shared_openmpi216/
    export OPENMPI=/work/lopez/openmpi-2.1.6/
    export LD_LIBRARY_PATH=$OPENMPI/lib:$ORCA_EXE:$LD_LIBRARY_PATH
    export PATH=$OPENMPI/bin:$PATH
    $ORCA_EXE/orca $INPUT.inp > $INPUT.out
    date >> $INPUT.out
    
    #CHECK IF THE THING FINISHED
    if [[ $finished -lt 1 ]]
        then
        echo "Check $INPUT optimization" >> ../status.txt
        exit 1234
     fi 
else
    echo "no $SLURM_ARRAY_TASK_ID generated"
    rm {0}-conf$SLURM_ARRAY_TASK_ID.inp
fi

""".format(title,OPTpartition,OPTcores,OPTmem,lowest10dir,tmptitle)

    for i in range(1,11):
        inputfile="""!{6}
%pal nprocs {1} end
%geom
maxiter 500
end
%Maxcore {2}000
*xyzfile {3} {4} {0}-conf{5}.xyz

""".format(title,OPTcores, OPTmem, charge, mult, i,OPTmethod)

        with open('{0}/{1}-conf{2}.inp'.format(lowest10dir,title,i), 'w') as com:
            com.write(inputfile)
        if i == 1:
            with open('{0}/{1}-coms.txt'.format(lowest10dir,title), 'w') as coms:
                coms.write('{0}-conf{1}.inp\n'.format(title,i))
        else:
            with open('{0}/{1}-coms.txt'.format(lowest10dir,title), 'a') as coms:
                coms.write('{0}-conf{1}.inp\n'.format(title,i))


    return sbatch

def NONEopt(title,lowest10dir):
    sbatch="""#!/bin/bash
#SBATCH --job-name={0}-OPT
#SBATCH --output=out.o
#SBATCH --error=out.e
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=20:00
hostname

work={1}
cd $work

    cat ../ORCA/*out >> complete
    nstruct=$(ls -la ../ORCA/{0}-all-sorted-conf*.xyz |wc -l)
    for ((x=1;x<=nstruct;x++)); do  awk "/FINAL SINGLE POINT ENERGY/{{i++}}i==$x{{print ; exit}}" complete | awk '{{ print $5}}' >> energies.txt; done
    awk '{{print NR "   "  $s}}' energies.txt > output.txt
    sort -k2n output.txt > energies-sorted.txt
    lowest10=$(awk  '{{print $1}} NR==10{{exit}}' energies-sorted.txt )

""".format(title,OPTpartition,OPTcores,OPTmem,lowest10dir)
    return sbatch

def Start(CRESTdir,title,ORCAdir,lowest10dir):
    start="""#!/bin/bash
 
CRESTID=$(sbatch --parsable {0}/{1}-CREST.sbatch)
echo Submitted
""".format(CRESTdir,title,ORCAdir,lowest10dir)
    return start

def Extract(title,workdir,lowest10dir):
    if len(title) > 20:
        tmptitle=title[0:18]
    else:
        tmptitle=title

    sbatch="""#!/bin/bash
#SBATCH --job-name={2}-EX
#SBATCH --output=out.o
#SBATCH --error=out.e
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=20:00
hostname

work={3}
cd $work

    if test -f {0}-optimized-energies.txt
        then 
        rm {0}-optimized-energies.txt
    fi
    for i in lowest10/{0}*log;
        do
            energy=$(grep "Sum of electronic and thermal Free Energies" $i | awk '{{ print $8 }}')
            echo $i $energy >> {0}-optimized-energies.txt
     done
    if test -f "{0}-optimized-energies.txt"
        then
        lowestenergy=$(sort -u -g -k2,2 {0}-optimized-energies.txt | head -1 |awk '{{ print $1 }}')
             else
         echo "No valid structures found" >> {0}-optimized-energies.txt
     fi
     cp $lowestenergy ./{0}.log

rm lowest10/Gau*
rm ORCA/*gbw

""".format(title,lowest10dir,tmptitle,workdir)
    return sbatch

#Handle Constraints
if constraint == True:
    exclude=[]
    include=[]
    matches=["angle","distance","dihedral"]
    if constrainttype == 'c':
        shutil.copyfile('{0}/{1}.c'.format(work,title),'{0}/{1}.c'.format(CRESTdir,title))
    elif constrainttype == 'constraints':
        shutil.copyfile('{0}/{1}.constraints'.format(work,title),'{0}/{1}.c'.format(CRESTdir,title))
    shutil.copyfile('{0}/{1}.xyz'.format(work,title),'{0}/{1}.ref'.format(CRESTdir,title))
    with open('{0}/{1}.xyz'.format(work,title),'r') as xyz:
        atoms=int(xyz.readline())
    with open('{0}/{1}.c'.format(CRESTdir,title),'r') as cinpfile: 
        cinp=cinpfile.readlines()
        cinpfile.seek(0)
        cinpcontent=cinpfile.read()
    for line in cinp:
        if any(x in line for x in matches):
           for word in line.split():
               word=word.replace(",","")
               if word.isdigit():
                   exclude.append(int(word))
    missing=[]
    exclude.sort()
    numbers=exclude
    numbers.insert(0, 0) # add the minimum value on begining of the list
    numbers.append(atoms+1)  # add the maximum value at the end of the list
    for rank in range(0, len(numbers)-1):
        if numbers[rank+1] - numbers[rank] > 2:
            missing.append("%s-%s"%(numbers[rank] +1 , numbers[rank+1] - 1))
        elif numbers[rank+1] - numbers[rank] == 2:
            missing.append(str(numbers[rank]+1))
    missing=str(missing)[1:-1]
    include=missing.replace("'","")
    with open('{0}/{1}.c'.format(CRESTdir,title),'w') as cinpfile:
        cinpfile.write('$constrain\n')
        cinpfile.write(cinpcontent)
        cinpfile.write('force constant=1.0\nreference={0}.ref\n$metadyn\natoms: {1}\n$end'.format(title,include))
          
#Generate Input
with open('{0}/{1}-CREST.sbatch'.format(CRESTdir,title), 'w') as CRESTsbatch:
    CRESTsbatch.write(CREST(title,CRESTpartition,CRESTcores,CRESTmem,CRESTtime,CRESTdir,constraint,nci))
with open('{0}/{1}-ORCA.sbatch'.format(ORCAdir,title), 'w') as ORCAsbatch:
    ORCAsbatch.write(ORCA(ORCAcores,ORCAtime,title,ORCApartition,ORCAmem,ORCAdir,ORCAmethod,charge,mult,SP))
with open ('{0}/{1}-OPT.sbatch'.format(lowest10dir,title), 'w') as OPTsbatch:
    if optprogram == 'gaussian':
        OPTsbatch.write(GAUS(title,OPTcores,OPTmem,GAUSmethod,GAUSbasis,GAUSroute,OPTpartition,lowest10dir,charge,mult,g16fail))
    elif optprogram == 'orca':
        OPTsbatch.write(ORCAopt(title,OPTcores,OPTmem,ORCAoptmethod,OPTpartition,lowest10dir,charge,mult))
    elif optprogram == 'none':
        OPTsbatch.write(NONEopt(title,lowest10dir))
with open('{0}/start.sh'.format(workdir),'w') as start:
    start.write(Start(CRESTdir,title,ORCAdir,lowest10dir))
with open('{0}/{1}-EX.sbatch'.format(workdir,title),'w') as extract:
    extract.write(Extract(title,workdir,lowest10dir))

