\#!/bin/bash

# Template:  OpenMPI, Default (OpenIB Infiniband) Variant
# Revision:  $Id: openmpi-ib.qs 577 2015-12-21 14:39:43Z frey $
#
# Usage:
# 1. Modify "NPROC" in the -pe line to reflect the number
#    of processors desired.
# 2. Modify the value of "MY_EXE" to be your MPI program and
#    "MY_EXE_ARGS" to be the array of arguments to be passed to
#    it.
# 3. Uncomment the WANT_CPU_AFFINITY line if you want Open MPI to
#    bind workers to processor cores.  Can increase your program's
#    efficiency.
# 4. Uncomment the SHOW_MPI_DEBUGGING line if you want very verbose
#    output written to the Grid Engine output file by OpenMPI.
# 5. If you use exclusive=1, please be aware that NPROC will be
#    rounded up to a multiple of 20.  In this case, set the
#    WANT_NPROC variable to the actual core count you want.  The
#    script will "load balance" that core count across the N nodes
#    the job receives.
# 6. Jobs default to using 1 GB of system memory per slot.  If you
#    need more than that, set the m_mem_free complex.
# 




#----------------------------------
#                Farber parameters
# tells cluster to assign (20) processors
#$ -pe mpi 1
# tells cluster to allocate 1GB memory PER processor
#$ -l m_mem_free=3G
# tells cluster to give you exclusive access to node
# $ -l exclusive=1
# send script to the standby queue
#$ -l standby=1,h_rt=4:00:00
# setup messaging about (b)egin, (e)nd, (a)bort, and (s) of your program
# $ -o mpi_submit_script.sh.out
# $ -m beas
# send messages to this email address
# -M 9735252392@vtext.com
#----------------------------------

#      Load any packages you need to run it
vpkg_require openmpi/intel64
vpkg_require gnuplot/4.6
OPENMPI_FLAGS='-np 1'

#
# The MPI program to execute:
#
MY_EXE="gendomain.out"
TMP_DIR='/lustre/scratch/rohitk/FDPML_tmp'
DOMAIN="${TMP_DIR}/Domain.dat"
MASS="${TMP_DIR}/mass.dat"


cat > test_input << EOF
  &system
	PD(1) = 11, 11, 20
	crystal_coordinates = .false., 
	domain_prefix = '${DOMAIN}',
	mass_prefix = '${MASS}', 
	mass_input = .true., 
	flfrc1 = '/home/1628/QuantumEspresso/Si/results/Si_q2.fc', 
	flfrc2 = '/home/1628/QuantumEspresso/Si/results/Si_q2.fc'
  /
  &domain
	nanoparticle = .false.
	random = .false.
	c = 0.1
	Radius = 0.5
  /
EOF

mpirun -quiet ${OPENMPI_FLAGS} $MY_EXE < test_input > output_test

# check to see if it completed without errors
rc=$?

exit $rc

echo "-- DONE --"
