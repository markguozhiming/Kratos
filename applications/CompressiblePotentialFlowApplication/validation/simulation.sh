#!/bin/bash
#
#run this file using the command:
# bash simulation.sh
source /home/inigo/Documents/paths/salomeConverter.sh
source /home/inigo/Documents/paths/kratosMaster.sh
source /home/inigo/intel/mkl/bin/mklvars.sh intel64 lp64

#cd /home/inigo/software/kratosMaster/Kratos
GITBRANCH=$(git symbolic-ref HEAD | sed -e 's,.*/\(.*\),\1,')
#cd /home/inigo/simulations/naca0012/07_salome/03_MeshRefinement
Input_Dir=/home/inigo/simulations/naca0012/07_salome/05_MeshRefinement

#Parameters
Number_Of_Refinements=4
Number_Of_AOAS=1

Initial_AOA=5.0
AOA_Increment=1.0

Initial_Airfoil_MeshSize=1.024e-2
Airfoil_Refinement_Factor=2.0

Initial_FarField_MeshSize=2.0
FarField_Refinement_Factor=1.0

DATE=`date '+%Y%m%d_%H%M%S'`
FILE=terminal_outputs/output_terminal.txt
NAME=${FILE%.*}
EXT=${FILE#*.}
NEWFILE=${NAME}_${DATE}.${EXT}

#Set the parameters
source generate_mdpas/set_parameters.sh
cd generate_mdpas/

#Run salome: generate geometry and mesh
#python3 runSalome.py

#Convert salomes mesh into mdpa
#rm $Work_Dir/mdpas/*
#python3 use_converter.py

#Run Kratos
cd ..
rm -rf $Work_Dir/plots/cl/data/cl_*
rm -rf $Work_Dir/plots/cd/data/cd_*

rm -rf $Work_Dir/plots/cp/data/AOA*
rm $Work_Dir/plots/cp/cp_*
rm $Work_Dir/plots/cp/plots/*

rm $Work_Dir/plots/cl/figures_cl.tex
rm $Work_Dir/plots/cd/figures_cd.tex
#rm plots/cp/figures.tex

rm $Work_Dir/plots/results/*
rm -rf /media/inigo/10740FB2740F9A1C/Outputs/03_MeshRefinement/*
python3 MeshRefinement.py #> $NEWFILE

source generate_mdpas/unset_parameters.sh

#rm $Input_Dir/generate_mdpas/output_salome/*
#rm mdpas/*
rm cp*

cd $Work_Dir/plots/cl
pdflatex -interaction=batchmode main_cl.tex
cd $Work_Dir/plots/cd
pdflatex -interaction=batchmode main_cd.tex
cd $Work_Dir/aoa/
pdflatex -interaction=batchmode cl_aoa.tex
#cd ../cp/
#pdflatex -interaction=batchmode cp.tex
DIRECTORY=/media/inigo/10740FB2740F9A1C/Old_Outputs/03_MeshRefinement
mkdir -p ${DIRECTORY}_${DATE}_${GITBRANCH}
mkdir -p ${DIRECTORY}_${DATE}_${GITBRANCH}/output_gid
cd ../..
cp -r $Work_Dir/plots/ ${DIRECTORY}_${DATE}_${GITBRANCH}
cp -r /media/inigo/10740FB2740F9A1C/Outputs/03_MeshRefinement/A* ${DIRECTORY}_${DATE}_${GITBRANCH}/output_gid

