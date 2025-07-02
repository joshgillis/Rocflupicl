module purge
module load gcc/14.2.0 openmpi/5.0.7
make clean
make RFLU=1 PICL=1 FOLDER=1 SPEC=1 -j16
ls --color=auto
