#include <stdlib.h>

#include "mpi.h"

#include "HYPRE_struct_mv.h"
#include "HYPRE_struct_ls.h"

#include "ssp_hypre.h"

struct ssp_hypre_struct *ssp_hypre_struct_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
				   int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] ) {

  /* Set up the grid, stencil and matrix */

  struct ssp_hypre_struct *data_for_hypre_struct;

  int ele[ 3 ];

  int i;

  data_for_hypre_struct = malloc( sizeof( *data_for_hypre_struct ) );

  data_for_hypre_struct -> comm = MPI_Comm_f2c( comm );

  /* Note reverse order of dimensions to cope with C <-> Fortran mappings */
  for( i = 0; i < 3; i++ ) {
    data_for_hypre_struct -> lb[ i ] = lb[ 2 - i ];
    data_for_hypre_struct -> ub[ i ] = ub[ 2 - i ];
  }

  HYPRE_StructGridCreate( data_for_hypre_struct -> comm, 3,
			  &( data_for_hypre_struct -> grid ) );
  HYPRE_StructGridSetPeriodic( data_for_hypre_struct -> grid, n );
  HYPRE_StructGridSetExtents( data_for_hypre_struct -> grid,
			      data_for_hypre_struct -> lb, data_for_hypre_struct -> ub );
  HYPRE_StructGridAssemble( data_for_hypre_struct -> grid );

  /* The stencil */
  HYPRE_StructStencilCreate( 3, n_stencil, &( data_for_hypre_struct -> stencil ) );
  data_for_hypre_struct -> stencil_element_list =
    malloc( n_stencil * sizeof( *( data_for_hypre_struct -> stencil_element_list ) ) );
  for( i = 0; i < n_stencil; i++ ) {
    /* Note reverse order of dimensions to cope with C <-> Fortran mappings */
    ele[ 0 ] = stencil_elements[ i ][ 2 ];
    ele[ 1 ] = stencil_elements[ i ][ 1 ];
    ele[ 2 ] = stencil_elements[ i ][ 0 ];
    printf( "%d %d %d %d %24.12f\n", i, ele[ 0 ], ele[ 1 ], ele[ 2 ], stencil_values[ i ] );
    data_for_hypre_struct -> stencil_element_list[ i ]= i;
    HYPRE_StructStencilSetElement( data_for_hypre_struct -> stencil, i, ele );
  }

  printf( "%24.12f %24.12f %24.12f\n", stencil_values[ 0 ], stencil_values[ 1 ], stencil_values[ 234 ] );
  /* The matrix */
  HYPRE_StructMatrixCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, data_for_hypre_struct -> stencil,
			    &( data_for_hypre_struct -> A ) );
  HYPRE_StructMatrixSetSymmetric( data_for_hypre_struct -> A, 1 );
  HYPRE_StructMatrixInitialize  ( data_for_hypre_struct -> A );

  /* Set the values for the matrix */
  /* The stencil is constant across the whole grid */
  HYPRE_StructMatrixSetConstantEntries( data_for_hypre_struct -> A, n_stencil, data_for_hypre_struct -> stencil_element_list );
  HYPRE_StructMatrixSetConstantValues ( data_for_hypre_struct -> A, n_stencil, data_for_hypre_struct -> stencil_element_list, stencil_values );

  /* Assemble the matrix - blocking */
  HYPRE_StructMatrixAssemble( data_for_hypre_struct -> A );
  HYPRE_StructMatrixPrint( "matrix.dat", data_for_hypre_struct -> A, 1 );

  return data_for_hypre_struct;
  
}

void ssp_hypre_struct_pfmg_solve( struct ssp_hypre_struct *data_for_hypre_struct, int n1, int n2, int n3, double rhs[ n3 ][ n2 ][ n1 ],
				  double soln[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info ) {

  /* Solve the equations using the PFMG solver */
  HYPRE_StructVector b;
  HYPRE_StructVector x;
  HYPRE_StructSolver solver;
  
  /* Assume everything worked for the moment */
  *info = 0;

  /* Create the RHS */
  HYPRE_StructVectorCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, &b );
  /* Put the RHS from the Fortran array into the HYPRE object */
  HYPRE_StructVectorInitialize( b );
  HYPRE_StructVectorSetBoxValues( b, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( rhs[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the RHS */
  HYPRE_StructVectorAssemble( b );
  HYPRE_StructVectorPrint( "rhs.dat", b, 1 );
  
  /* Create the initial guess at solution */
  HYPRE_StructVectorCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, &x );
  /* Put the initial guess from the Fortran array into the HYPRE object */
  HYPRE_StructVectorInitialize( x );
  HYPRE_StructVectorSetBoxValues( x, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( soln[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the initial guess */
  HYPRE_StructVectorAssemble( x );
  HYPRE_StructVectorPrint( "guess.dat", x, 1 );

  /* Create the solver */
  HYPRE_StructPFMGCreate( data_for_hypre_struct -> comm, &solver );

  /* Set up the solver */
  HYPRE_StructPFMGSetup( solver, data_for_hypre_struct -> A, b, x );

  /* Solve the equations */
  HYPRE_StructPFMGSolve( solver, data_for_hypre_struct -> A, b, x );

  /* Get the solution from the HYPRE object into the Fortran array */
  HYPRE_StructVectorGetBoxValues( x, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( soln[ 0 ][ 0 ][ 0 ] ) );

  /* Get some interesting data */
  HYPRE_StructPFMGGetNumIterations( solver, n_iter );
  HYPRE_StructPFMGGetFinalRelativeResidualNorm( solver, residual );

  /* Destroy the solver */
  HYPRE_StructPFMGDestroy( solver );
  
  /* Destroy the vectors */
  HYPRE_StructVectorDestroy( x );
  HYPRE_StructVectorDestroy( b );
  
}

void ssp_hypre_struct_free( struct ssp_hypre_struct *data_for_hypre_struct ){

  /* Free up the memory used for the grid, stencil and matrix */

  HYPRE_StructMatrixDestroy ( data_for_hypre_struct -> A       );
  free                      ( data_for_hypre_struct -> stencil_element_list );
  HYPRE_StructStencilDestroy( data_for_hypre_struct -> stencil );
  HYPRE_StructGridDestroy   ( data_for_hypre_struct -> grid    );
    
  free( data_for_hypre_struct );
  
  return;
}

