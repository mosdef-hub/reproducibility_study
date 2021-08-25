#! /bin/bash
#SBATCH --job-name="rpd-uamethane"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=11:59:59
#SBATCH --export=ALL
#SBATCH --mail-user=zijiewu@udel.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --requeue
# SBATCH --partition=_workgroup_
#SBATCH --mem=8G

export VALET_PATH=$VALET_PATH:$WORKDIR/udsw/valet/etc
vpkg_require jlab-lammps

OPENMPI_VERSION='unknown'
OMPI_INFO="$(which ompi_info)"
if [ $? -eq 0 -a -n "$OMPI_INFO" ]; then
  OPENMPI_VERSION="$(${OMPI_INFO} --version | egrep 'v[0-9]+' | sed 's/^.*v//')"
fi
OPENMPI_FLAGS="--display-map --mca btl ^tcp"
if [ "${WANT_CPU_AFFINITY:-NO}" = "YES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --bind-to core"
fi
if [ "${WANT_NPROC:-0}" -gt 0 ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --np ${WANT_NPROC} --map-by node"
fi
if [ "${SHOW_MPI_DEBUGGING:-NO}" = "YES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --debug-devel --debug-daemons --display-devel-map --display-devel-allocation --mca mca_verbose 1 --mca coll_base_verbose 1 --mca ras_base_verbose 1 --mca ras_gridengine_debug 1 --mca ras_gridengine_verbose 1 --mca btl_base_verbose 1 --mca mtl_base_verbose 1 --mca plm_base_verbose 1 --mca pls_rsh_debug 1"
  if [ "${WANT_CPU_AFFINITY:-NO}" = "YES" ]; then
    OPENMPI_FLAGS="${OPENMPI_FLAGS} --report-bindings"
  fi
fi

if [ $((1)) ]
    then
    echo "GridEngine parameters:"
    echo "  mpirun        = "`which mpirun`
    echo "  nhosts        = $NHOSTS"
    echo "  nproc         = $NSLOTS"
    echo "  executable    = $MY_EXE"
    echo "  Open MPI vers = $OPENMPI_VERSION"
    echo "  MPI flags     = $OPENMPI_FLAGS"
    echo "-- begin OPENMPI run --"
    mpirun ${OPENMPI_FLAGS} lmp -in in.minimize
    rc=$?
    echo "-- end OPENMPI run --"
fi
