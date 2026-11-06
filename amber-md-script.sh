#!/bin/bash

echo " please contact the PI at vaibhavadixit@gmail.com or vaibhav@niperguwahati.in
============================================================================================================================================================================================================================================= "

#source /home/$USER/.bashrc
#source $AMBERHOME/amber.sh
#module load cuda/10.1
module load amber/22
module load openbabel/3.1.1

cwd=$(pwd)
export $(xargs <$1)
echo Reading $1 file for job parameters
echo "ligand and protein used for MD are $ligand $protein " > $protein-$ligand-MD-logfile.txt
echo "$dirpath used in the calculation" >> $protein-$ligand-MD-logfile.txt
echo "$dir used in the calculation" >> $protein-$ligand-MD-logfile.txt

mkdir $dirpath/$dir

echo "created directories for oxd and red state " >> $protein-$ligand-MD-logfile.txt
echo "parameters for oxd and red states read from $oxdparm $redparm " >> $protein-$ligand-MD-logfile.txt

echo "Copying ligand and protein Changing directory to Ferrous or the reduced state i.e. $dirpath/$reddir"
#sed -i 's/UNL/LIG/g' $ligand.pdbqt

echo "parameterizing ligand and generating the protein complex"

#python oeb-pdbqt.py $ligand.oeb $ligand.pdb

#obabel $cofac.pdb -O $cofac1.mol2 -h -p 7.4 
#mv $cofac1.mol2 $cofac.pdb
#for i in $(ls $ligand*.pdb | grep -v $posenumber ); do echo removing file $i; rm $i ; done

#echo "pose number $posenumber from the vina docking output $ligand.pdbqt file was used for MD simulations " >> ../$protein-$ligand-MD-logfile.txt
#echo "pose number $posenumber from the vina docking output $ligand.pdbqt file was used for MD simulations " 
#rm $ligand"2".pdb $ligand"3".pdb $ligand"4".pdb $ligand"5".pdb $ligand"6".pdb $ligand"7".pdb $ligand"8".pdb $ligand"9".pdb

#obabel $ligand"$posenumber".pdb -O $ligand"$posenumber"-h.pdb -h

#grep ATOM $ligand"$posenumber"-h.pdb > $ligand.pdb
#echo "Give the expeced net-charge on the molecule as an integer e.g. -1, 0 or 1"
#babel test-$ligand.pdb -O test.mol2 --partialcharge gasteiger
#molchrg=$(egrep -A1000 ATOM test.mol2 | grep -v ROOT | grep LIG | sed -e "1d" | awk '{print $9}' | st | awk '{print $4}' | tail -1 | xargs printf "%1.0f")

#cp ../$ligand.mol2 .
 
$AMBERHOME/bin/antechamber -i $ligand.pdb -fi pdb -o $ligand.mol2 -fo mol2 -c bcc -s 2 -nc $molchrg
$AMBERHOME/bin/parmchk2 -i $ligand.mol2 -f mol2 -o $ligand.frcmod
$AMBERHOME/bin/antechamber -i $ligand.mol2 -fi mol2 -o $ligand.ac -fo ac 
$AMBERHOME/bin/antechamber -i $ligand.ac -fi ac -o $ligand-test.pdb -fo pdb 

$AMBERHOME/bin/antechamber -i $cofac.pdb -fi pdb -o $cofac.mol2 -fo mol2 -c bcc -s 2 -nc $molchrg
$AMBERHOME/bin/parmchk2 -i $cofac.mol2 -f mol2 -o $cofac.frcmod
$AMBERHOME/bin/antechamber -i $cofac.mol2 -fi mol2 -o $cofac.ac -fo ac 
#$AMBERHOME/bin/antechamber -i $cofac.ac -fi ac -o $cofac-test.pdb -fo pdb 

echo "Antechamber used to prepare ligand parameters and atomtypes for $ligand " >> ../$protein-$ligand-MD-logfile.txt
echo "Antechamber used to prepare ligand parameters and atomtypes for $ligand "
#cp test-$ligand.pdb test.pdb
#sed 's/[[:space:]]*[^[:space:]]*$//' $ligand-test.pdb | sed 's/[[:space:]]*[^[:space:]]*$//'  > $ligand.pdb

#sed 's/[[:space:]]*[^[:space:]]*$//' $cofac-test.pdb | sed 's/[[:space:]]*[^[:space:]]*$//'  > $cofac.pdb
#rm test.pdb
echo " "
echo "Combining $protein $ligand $cofac into a complex" | tee -a $protein-$ligand-MD-logfile.txt

cat $protein.pdb $ligand.pdb $cofac.pdb > $protein-$ligand.pdb
echo "unloading abmer22" | tee -a $protein-$ligand-MD-logfile.txt
module unload amber/22
echo "loading amber/22" | tee -a $protein-$ligand-MD-logfile.txt
module load amber/22
echo "Preparing the complex for Amber usage" | tee -a $protein-$ligand-MD-logfile.txt
$AMBERHOME/bin/pdb4amber -i $protein-$ligand.pdb -o $protein-$ligand-comp.pdb
echo "copying files $ligand $protein into $dir" | tee -a $protein-$ligand-MD-logfile.txt
cp $ligand.mol2 $dirpath/$dir
cp $ligand.frcmod $dirpath/$dir

cp $cofac.mol2 $dirpath/$dir
cp $cofac.frcmod $dirpath/$dir

cp $protein-$ligand-comp.pdb $dirpath/$dir
echo "unloading abmer18" | tee -a $protein-$ligand-MD-logfile.txt
module unload amber/22
echo "loading amber/22" | tee -a $protein-$ligand-MD-logfile.txt
module load amber/22
echo "$protein-$ligand complex formed and saved for MD simulations " | tee -a ../$protein-$ligand-MD-logfile.txt
echo "returning to the partent directory" | tee -a $protein-$ligand-MD-logfile.txt

echo "Writing tleap input file for the ferric or oxidized state" | tee -a $protein-$ligand-MD-logfile.txt
echo "
source leaprc.protein.ff14SB
source leaprc.gaff2
 " > $dirpath/$dir/$protein-$ligand-tleap.in

echo "
#LIG = loadmol2 $dirpath/$dir/$ligand.mol2
UDP = loadmol2 $dirpath/$dir/$cofac.mol2
loadamberparams $dirpath/$dir/$cofac.frcmod
#loadamberparams $dirpath/$dir/$ligand.frcmod
loadamberparams frcmod.ions1lm_126_tip3p
source leaprc.water.tip3p
mol = loadpdb $dirpath/$dir/$protein-$ligand-comp.pdb
savepdb mol $dirpath/$dir/$protein-$ligand-dry.pdb
saveamberparm mol $dirpath/$dir/$protein-$ligand-dry.prmtop $dirpath/$dir/$protein-$ligand-dry.inpcrd
charge mol
solvatebox mol TIP3PBOX 10.0
addions mol K+ 0
addions mol Cl- 0
charge mol
savepdb mol $dirpath/$dir/$protein-$ligand-solv.pdb
saveamberparm mol $dirpath/$dir/$protein-$ligand-solv.prmtop $dirpath/$dir/$protein-$ligand-solv.inpcrd
quit" >> $dirpath/$dir/$protein-$ligand-tleap.in

echo "tleap input file has been created for the ferric or oxidized state. Please check if the system charged is 0 if not make appropriate changes in the type and number of ions to be added in the addions file and rerun this main script"

echo "Writing cpptraj input files"

writemdfiles () {
	echo "Writing MD input file for a 100 ns MD production run"
	echo "Production Stage in Explicit Solvent 20 ns
&cntrl
  ntt=3,           ! Temperature scaling (=3, Langevin dynamics)
  gamma_ln=2.0,    ! Collision frequency of the Langevin dynamics in ps-1
  ntc=2,           ! SHAKE constraints (=2, hydrogen bond lengths constrained)
  ntf=2,           ! Force evaluation (=2, hydrogen bond interactions omitted)
  ntb=1,           ! Boundaries (=1, constant volume)
  cut=10.0,        ! Cutoff
  dt=0.002,        ! The time step in picoseconds
  nstlim=10000000, ! Number of MD steps to be performed
  ig=-1,           ! Random seed (=-1, get a number from current date and time)
  ntwr=5000,      ! Restart file written every ntwr steps
  ntwx=5000,      ! Trajectory file written every ntwx steps
  ntpr=5000,      ! The mdout and mdinfo files written every ntpr steps
  ioutfm=1,        ! Trajectory file format (=1, Binary NetCDF)
  iwrap=1,         ! Translate water molecules into the original simulation box
  igb=0,           ! GB model (=0, explicit solvent)
  irest=1,         ! Flag to restart the simulation
  ntx=5,           ! Initial condition (=5, coord. and veloc. read from the inpcrd file)
/
" > $dirpath/$1/prod20ns.in

#Ferrousdir=$(echo Ferrous-$protein-$ligand | sed 's/_out//g')


}

writemdfiles $dir 
#writemdfiles $reddir $dir fe2crdfe2prm fe2crdfe3prm


#echo "give leap input and output file names one after the other and press enter before entering the next file name"
leapinput=$protein-$ligand-tleap.in
leapoutput=$protein-$ligand-tleap.out
#echo "give base file names to be used for MD prmtop, inpcrd, out, rst7, mdinfo and trajectory (nc) files"
basefile=$protein-$ligand-solv
echo $leapinput $leapoutput $basefile | tee -a $protein-$ligand-MD-logfile.txt
echo "creating MD prmtop inpcrd files for the Ferric or the oxidized state in $dirpath/$dir " | tee -a $protein-$ligand-MD-logfile.txt
echo "creating gpu-job head" | tee -a $protein-$ligand-MD-logfile.txt
visdivice="0,1,2"
echo "#!/bin/bash
#SBATCH -A meity     ##account name
#SBATCH --nodes=1  ## number of nodes
#SBATCH --nodelist=node1
#SBATCH --ntasks=12
#SBATCH --error=job.%J.err    ## any error during job submission or execution
#SBATCH --output=job.%J.out   ##any output after job execution 
#SBATCH --partition=GPU_NODES    ##partition name
#SBATCH --job-name=$ligand-ugt-solv   ##job name
export CUDA_VISIBLE_DEVICES=$visdivice

hostname
module load openmpi/4.1.4
module load amber/22
" > gpu-job-head.txt


echo "creating cpu-job head " | tee -a $protein-$ligand-MD-logfile.txt
echo "#!/bin/bash
#SBATCH -A meity     ##account name
#SBATCH --nodes=1  ## number of nodes
#SBATCH --nodelist=node1
#SBATCH --ntasks=12
#SBATCH --error=job.%J.err    ## any error during job submission or execution
#SBATCH --output=job.%J.out   ##any output after job execution 
#SBATCH --partition=GPU_NODES    ##partition name
#SBATCH --job-name=$ligand-ugt-solv   ##job name


hostname
module load openmpi/4.1.4
module load amber/22

" > cpu-job-head.txt
echo "creating tleap input and sbatch file for $dir state" | tee -a $protein-$ligand-MD-logfile.txt
cat cpu-job-head.txt > $dirpath/$dir/$protein-$ligand-tleap.batch
echo "module load amber/22" >> $dirpath/$dir/$protein-$ligand-tleap.batch
echo "source $AMBERHOME/amber.sh" >> $dirpath/$dir/$protein-$ligand-tleap.batch
echo "$AMBERHOME/bin/tleap -s -f $dirpath/$dir/$leapinput > $dirpath/$dir/$leapoutput " >> $dirpath/$dir/$protein-$ligand-tleap.batch


oxdtleapjobid=$(sbatch $dirpath/$dir/$protein-$ligand-tleap.batch | awk '{print $4}')
echo "sleeping for 2m"
sleep 2m
oxdtleapjobstatus=$(sacct | grep $oxdtleapjobid | awk '{print $6}' | head -1)

#echo "creating tleap input and sbatch file for $reddir state" | tee -a $protein-$ligand-MD-logfile.txt
#cat cpu-job-head.txt > $dirpath/$reddir/$protein-$ligand-tleap.batch
#echo "creating MD prmtop inpcrd files for the Ferrous or the reduced state in $dirpath/$reddir "
#echo "module load amber/22" >> $dirpath/$reddir/$protein-$ligand-tleap.batch
#echo "source $AMBERHOME/amber.sh" >> $dirpath/$dir/$protein-$ligand-tleap.batch
#echo "$AMBERHOME/bin/tleap -s -f $dirpath/$reddir/$leapinput > $dirpath/$reddir/$leapoutput " >> $dirpath/$reddir/$protein-$ligand-tleap.batch

#redtleapjobid=$(sbatch $dirpath/$reddir/$protein-$ligand-tleap.batch | awk '{print $4}')
#echo "sleeping for 2m"
#sleep 2m
#redtleapjobstatus=$(sacct | grep $redtleapjobid | awk '{print $6}' | head -1)

echo "creating input and submitting equil jobs" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt

while [[ $redtleapjobstatus == RUNNING ]];
do
	sleep 1m
done
echo "finished equiliration of $dir structure" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt

while [[ $oxdtleapjobstatus == RUNNING ]];
do
	sleep 1m
done
echo "finished equiliration of $reddir structure" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt

echo "changing directory to $dir" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cd $dirpath/$dir
cp ../gpu-job-head.txt .
echo "creating MDequil sbatch files and submitting the job $dir state" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cat gpu-job-head.txt > $dirpath/$dir/$protein-$ligand-MDequil.batch
echo "module load amber/22" >> $dirpath/$dir/$protein-$ligand-MDequil.batch
echo "source $AMBERHOME/amber.sh" >> $dirpath/$dir/$protein-$ligand-tleap.batch

echo "/home/vaibhav/Documents/AmberMdPrep.sh -O -p $dirpath/$dir/$basefile.prmtop -c $dirpath/$dir/$basefile.inpcrd --temp 300 --ares UDP --ares LIG " >> $dirpath/$dir/$protein-$ligand-MDequil.batch
oxdequiljobid=$(sbatch $dirpath/$dir/$protein-$ligand-MDequil.batch | awk '{print $4}')
oxdequiljobstatus=$(sacct | grep $oxdequiljobid | awk '{print $6}' | head -1)

echo "sleeping for the equiljobs to finish"

while [[ $(sacct | grep $oxdequiljobid | awk '{print $6}' | head -1 ) == RUNNING ]];
do
	sleep 1m
done
echo "finished equiliration of $dir structure" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt


echo "changing directory to $dir and submitting 20ns MD jobs in triplicate (prod, prod-rp and prod-rp1) after equilibriation is success" | tee -a $dirpath/$protein-$ligand-MD-logfile.txt
cd $dirpath/$dir
sed -i 's/export CUDA_VISIBLE_DEVICES=1/export CUDA_VISIBLE_DEVICES=1/g' gpu-job-head.txt 
cat gpu-job-head.txt > $dirpath/$dir/oxdMD-rp-rp1.batch
echo "module load amber/22" >> $dirpath/$dir/oxdMD-rp-rp1.batch
echo "source $AMBERHOME/amber.sh" >> $dirpath/$dir/oxdMD-rp-rp1.batch
echo "pmemd.cuda -O -i $dirpath/$dir/prod20ns.in -p $dirpath/$dir/$basefile.prmtop -c $dirpath/$dir/final.1.ncrst -ref $dirpath/$dir/final.1.ncrst -r $dirpath/$dir/$basefile-prod.rst7 -inf $dirpath/$dir/$basefile-prod.mdinfo -o $dirpath/$dir/$basefile-prod.mdout -x $dirpath/$dir/$basefile-prod.nc 

" >> $dirpath/$dir/oxdMD-rp-rp1.batch 
oxdMDjobid=$(sbatch --dependency=afterok:$oxdequiljobid $dirpath/$dir/oxdMD-rp-rp1.batch | awk '{print $4}')
oxdMDjobstatus=$(sacct | grep $oxdMDjobid | awk '{print $6}'  | head -1)




while [[ $(sacct | grep $oxdMDjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 5m
done







