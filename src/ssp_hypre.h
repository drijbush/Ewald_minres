#include "mpi.h"

#include "HYPRE_struct_mv.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_mv.h"
#include "HYPRE_sstruct_ls.h"

/* Struct to hold hypre stuff so can return as a pointer to Fortran */
struct ssp_hypre_struct{
  MPI_Comm            comm;
  int                 lb[ 3 ], ub[ 3 ];
  int                 *stencil_element_list; /* for convenience */
  HYPRE_StructGrid    grid;
  HYPRE_StructStencil stencil;
  HYPRE_StructMatrix  A;
};

struct ssp_hypre_semi_struct{
  MPI_Comm            comm;
  int                 lb[ 3 ], ub[ 3 ];
  int                 *stencil_element_list; /* for convenience */
  HYPRE_SStructGrid    grid;
  HYPRE_SStructStencil stencil;
  HYPRE_SStructGraph   graph;
  HYPRE_SStructMatrix  A;
  HYPRE_ParCSRMatrix   Aij;
};

struct ssp_hypre_ij{
  MPI_Comm             comm;
  int                  ilower, iupper;
  int                  jlower, jupper;
  HYPRE_SStructMatrix  A;
  HYPRE_ParCSRMatrix   Acsr;
};

/* Function Prototypes */

/* Struct interface */
struct ssp_hypre_struct *ssp_hypre_struct_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
						 int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] );
void ssp_hypre_struct_pfmg_solve( struct ssp_hypre_struct *data_for_hypre_struct, int n1, int n2, int n3, double b[ n3 ][ n2 ][ n1 ],
				  double x[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info );
void ssp_hypre_struct_free( struct ssp_hypre_struct *data_for_hypre_struct );

/* Semi struct interface */
struct ssp_hypre_semi_struct *ssp_hypre_semi_struct_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
							   int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] );
void ssp_hypre_semi_struct_pfmg_solve( struct ssp_hypre_semi_struct *data_for_hypre_semi_struct, int n1, int n2, int n3, double b[ n3 ][ n2 ][ n1 ],
				  double x[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info );
void ssp_hypre_semi_struct_free( struct ssp_hypre_semi_struct *data_for_hypre_semi_struct );

/* IJ Matrix interface */
struct ssp_hypre_struct *ssp_hypre_ij_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
						 int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] );
void ssp_hypre_ij_pfmg_solve( struct ssp_hypre_struct *data_for_hypre_struct, int n1, int n2, int n3, double b[ n3 ][ n2 ][ n1 ],
				  double x[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info );
void ssp_hypre_ij_free( struct ssp_hypre_struct *data_for_hypre_struct );

