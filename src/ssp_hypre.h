#include "mpi.h"

#include "HYPRE_struct_mv.h"
#include "HYPRE_struct_ls.h"

struct ssp_hypre_struct{
  MPI_Comm            comm;
  int                 lb[ 3 ], ub[ 3 ];
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  A;
};
