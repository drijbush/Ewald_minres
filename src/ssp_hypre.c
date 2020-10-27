#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#include "HYPRE_struct_mv.h"
#include "HYPRE_struct_ls.h"

#include "ssp_hypre.h"

struct ssp_hypre_struct *ssp_hypre_struct_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
				   int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] ) {

  /* Set up the grid, stencil and matrix */

  struct ssp_hypre_struct *data_for_hypre_struct;
  
  int num_ghost[ 6 ];
  int ele[ 3 ];

  int max_stencil;
  int i, j;

  /* Find the halo size */
  max_stencil = -1;
  for( i = 0; i < n_stencil; i++ ) {
    for( j = 0; j < 3; j++ ) {
      if( abs( stencil_elements[ i ][ j ] ) > max_stencil ) {
	max_stencil = abs( stencil_elements[ i ][ j ] );
      }
    }
  }
  for( i = 0; i < 6; i++ ) {
    num_ghost[ i ] = max_stencil;
  }
  
  data_for_hypre_struct = malloc( sizeof( *data_for_hypre_struct ) );

  data_for_hypre_struct -> comm = MPI_Comm_f2c( comm );

  /* Note reverse order of dimensions to cope with C <-> Fortran mappings */
  for( i = 0; i < 3; i++ ) {
    /*        data_for_hypre_struct -> lb[ i ] = lb[ 2 - i ];
	      data_for_hypre_struct -> ub[ i ] = ub[ 2 - i ];  */
    data_for_hypre_struct -> lb[ i ] = lb[ i ];
    data_for_hypre_struct -> ub[ i ] = ub[ i ]; 
  }

  HYPRE_StructGridCreate( data_for_hypre_struct -> comm, 3,
			  &( data_for_hypre_struct -> grid ) );
  HYPRE_StructGridSetPeriodic( data_for_hypre_struct -> grid, n );
  HYPRE_StructGridSetExtents( data_for_hypre_struct -> grid,
			      data_for_hypre_struct -> lb, data_for_hypre_struct -> ub );
  HYPRE_StructGridSetNumGhost( data_for_hypre_struct -> grid, num_ghost );
  HYPRE_StructGridAssemble( data_for_hypre_struct -> grid );

  /* The stencil */
  HYPRE_StructStencilCreate( 3, n_stencil, &( data_for_hypre_struct -> stencil ) );
  data_for_hypre_struct -> stencil_element_list =
      malloc( n_stencil * sizeof( *( data_for_hypre_struct -> stencil_element_list ) ) );
  for( i = 0; i < n_stencil; i++ ) {
    ele[ 0 ] = stencil_elements[ i ][ 0 ];
    ele[ 1 ] = stencil_elements[ i ][ 1 ];
    ele[ 2 ] = stencil_elements[ i ][ 2 ]; 
    data_for_hypre_struct -> stencil_element_list[ i ]= i;
    HYPRE_StructStencilSetElement( data_for_hypre_struct -> stencil, i, ele );
    printf( "ele %d %d %d %d %20.12f\n", i, ele[ 0 ], ele[ 1 ], ele[ 2 ], stencil_values[ i ] );
  }

  /* The matrix */
  HYPRE_StructMatrixCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, data_for_hypre_struct -> stencil,
			    &( data_for_hypre_struct -> A ) );
  HYPRE_StructMatrixSetSymmetric( data_for_hypre_struct -> A, 1 );
  HYPRE_StructMatrixSetNumGhost ( data_for_hypre_struct -> A, num_ghost );
  HYPRE_StructMatrixInitialize  ( data_for_hypre_struct -> A );

  /* Set the values for the matrix */
  /* The stencil is constant across the whole grid */
  
  HYPRE_StructMatrixSetConstantEntries( data_for_hypre_struct -> A, n_stencil, data_for_hypre_struct -> stencil_element_list ); 
  HYPRE_StructMatrixAddToConstantValues ( data_for_hypre_struct -> A, n_stencil, data_for_hypre_struct -> stencil_element_list, stencil_values ); 


  /* Assemble the matrix - blocking */
  HYPRE_StructMatrixAssemble( data_for_hypre_struct -> A );
  /* HYPRE_StructMatrixPrint( "matrix.dat", data_for_hypre_struct -> A, 1 ); */

  return data_for_hypre_struct;
  
}

void ssp_hypre_struct_pfmg_solve( struct ssp_hypre_struct *data_for_hypre_struct, int n1, int n2, int n3, double rhs[ n3 ][ n2 ][ n1 ],
				  double soln[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info ) {

  /* Solve the equations using the PFMG solver */
  HYPRE_StructVector b;
  HYPRE_StructVector x;
  HYPRE_StructSolver solver;

  int retval;
  
  /* Assume everything worked for the moment */
  *info = 0;

  /* Create the RHS */
  HYPRE_StructVectorCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, &b );
  /* Put the RHS from the Fortran array into the HYPRE object */
  HYPRE_StructVectorInitialize( b );
  HYPRE_StructVectorSetBoxValues( b, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( rhs[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the RHS */
  HYPRE_StructVectorAssemble( b );
  /* HYPRE_StructVectorPrint( "rhs.dat", b, 1 ); */
  
  /* Create the initial guess at solution */
  HYPRE_StructVectorCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, &x );
  /* Put the initial guess from the Fortran array into the HYPRE object */
  printf( "!!!%f\n", soln[ 0 ][ 0 ][ 0 ] );
  HYPRE_StructVectorInitialize( x );
  HYPRE_StructVectorSetBoxValues( x, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( soln[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the initial guess */
  HYPRE_StructVectorAssemble( x );
  /* HYPRE_StructVectorPrint( "guess.dat", x, 1 ); */

  /* Create the solver */
  HYPRE_StructPFMGCreate( data_for_hypre_struct -> comm, &solver );

  /* Set up the solver */
  /* HYPRE_StructPFMGSetTol( solver, 5.0e-8 ); */
  HYPRE_StructPFMGSetLogging( solver, 100 );
  HYPRE_StructPFMGSetPrintLevel(solver, 1 );
  HYPRE_StructPFMGSetMaxIter( solver, 100 );
  HYPRE_StructPFMGSetRAPType( solver,  1 );
  /* HYPRE_StructPFMGSetMaxLevels( solver,  9 ); */
  /* HYPRE_StructPFMGSetRelaxType( solver,  2 ); */
  printf( "setup\n" );
  retval = HYPRE_StructPFMGSetup( solver, data_for_hypre_struct -> A, b, x );
  printf( "solve 1 %d\n", retval );

  /* Solve the equations */
  retval = HYPRE_StructPFMGSolve( solver, data_for_hypre_struct -> A, b, x );
  printf( "solve 2 %d\n", retval );

  /* Get the solution from the HYPRE object into the Fortran array */
  retval = HYPRE_StructVectorGetBoxValues( x, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( soln[ 0 ][ 0 ][ 0 ] ) );
  printf( "extract soln %d\n", retval );

  /* Get some interesting data */
  HYPRE_StructPFMGGetNumIterations( solver, n_iter );
  retval = HYPRE_StructPFMGGetFinalRelativeResidualNorm( solver, residual );
  printf( "Residual %d %f\n", retval, *residual );

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

struct ssp_hypre_semi_struct *ssp_hypre_semi_struct_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
				   int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] ) {

  /* Set up the grid, stencil and matrix */

  struct ssp_hypre_semi_struct *data_for_hypre_semi_struct;

  double *matrix_entries;
  
  int num_ghost[ 6 ];
  int ele[ 3 ];

  int max_stencil;
  int var_type;
  int n_mat;
  int i, j, k, m, i_stencil;

  /* Find the halo size */
  max_stencil = -1;
  for( i = 0; i < n_stencil; i++ ) {
    for( j = 0; j < 3; j++ ) {
      if( abs( stencil_elements[ i ][ j ] ) > max_stencil ) {
	max_stencil = abs( stencil_elements[ i ][ j ] );
      }
    }
  }
  for( i = 0; i < 6; i++ ) {
    num_ghost[ i ] = max_stencil;
  }
  
  data_for_hypre_semi_struct = malloc( sizeof( *data_for_hypre_semi_struct ) );

  data_for_hypre_semi_struct -> comm = MPI_Comm_f2c( comm );

  for( i = 0; i < 3; i++ ) {
    data_for_hypre_semi_struct -> lb[ i ] = lb[ i ];
    data_for_hypre_semi_struct -> ub[ i ] = ub[ i ]; 
  }

  HYPRE_SStructGridCreate( data_for_hypre_semi_struct -> comm, 3,
			   1, &( data_for_hypre_semi_struct -> grid ) );
  HYPRE_SStructGridSetPeriodic( data_for_hypre_semi_struct -> grid, 0, n );
  HYPRE_SStructGridSetExtents( data_for_hypre_semi_struct -> grid, 0, 
			      data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub );
  var_type = HYPRE_SSTRUCT_VARIABLE_CELL;
  HYPRE_SStructGridSetVariables( data_for_hypre_semi_struct -> grid, 0, 1, &var_type );
  HYPRE_SStructGridSetNumGhost( data_for_hypre_semi_struct -> grid, num_ghost );
  HYPRE_SStructGridAssemble( data_for_hypre_semi_struct -> grid );

  HYPRE_SStructStencilCreate( 3, n_stencil, &( data_for_hypre_semi_struct -> stencil ) );
  data_for_hypre_semi_struct -> stencil_element_list =
      malloc( n_stencil * sizeof( *( data_for_hypre_semi_struct -> stencil_element_list ) ) );
  for( i = 0; i < n_stencil; i++ ) {
    ele[ 0 ] = stencil_elements[ i ][ 0 ];
    ele[ 1 ] = stencil_elements[ i ][ 1 ];
    ele[ 2 ] = stencil_elements[ i ][ 2 ]; 
    data_for_hypre_semi_struct -> stencil_element_list[ i ]= i;
    HYPRE_SStructStencilSetEntry( data_for_hypre_semi_struct -> stencil, i, ele, 0 );
    printf( "ele %d %d %d %d %20.12f\n", i, ele[ 0 ], ele[ 1 ], ele[ 2 ], stencil_values[ i ] );
  }

  HYPRE_SStructGraphCreate( data_for_hypre_semi_struct -> comm, data_for_hypre_semi_struct -> grid,
			    &(  data_for_hypre_semi_struct -> graph ) );
  HYPRE_SStructGraphSetStencil( data_for_hypre_semi_struct -> graph, 0, 0,
				data_for_hypre_semi_struct -> stencil );
  HYPRE_SStructGraphAssemble( data_for_hypre_semi_struct -> graph );
  
  HYPRE_SStructMatrixCreate( data_for_hypre_semi_struct -> comm, data_for_hypre_semi_struct -> graph, &( data_for_hypre_semi_struct -> A ) );
  HYPRE_SStructMatrixSetSymmetric( data_for_hypre_semi_struct -> A, 0, 0, 0, 1 );
  /* HYPRE_SStructMatrixSetObjectType( data_for_hypre_semi_struct -> A, HYPRE_PARCSR ); */
  HYPRE_SStructMatrixSetObjectType( data_for_hypre_semi_struct -> A, HYPRE_SSTRUCT );
  HYPRE_SStructMatrixInitialize  ( data_for_hypre_semi_struct -> A );

  
  /*  Hypre_Semi_StructMatrixSetConstantEntries( data_for_hypre_semi_struct -> A, n_stencil, data_for_hypre_semi_struct -> stencil_element_list ); 
  Hypre_Semi_StructMatrixSetConstantValues ( data_for_hypre_semi_struct -> A, n_stencil, data_for_hypre_semi_struct -> stencil_element_list, stencil_values ); 
  */

  /* Can't use constant values for semi struct grids */
  n_mat  = ( ub[ 2 ] - lb[ 2 ] + 1 );
  n_mat *= ( ub[ 1 ] - lb[ 1 ] + 1 );
  n_mat *= ( ub[ 0 ] - lb[ 0 ] + 1 );
  matrix_entries = malloc( n_stencil * n_mat * sizeof( *matrix_entries ) );
  m = 0;
  for( k = lb[ 2 ]; k <= ub[ 2 ]; k++ ) {
    for( j = lb[ 1 ]; j <= ub[ 1 ]; j++ ) {
      for( i = lb[ 0 ]; i <= ub[ 0 ]; i++ ) {
	for( i_stencil = 0; i_stencil < n_stencil; i_stencil++ ) {
	  matrix_entries[ m ] = stencil_values[ i_stencil ];
	  m++;
	}
      }
    }
  }
  HYPRE_SStructMatrixSetBoxValues( data_for_hypre_semi_struct -> A, 0, lb, ub, 0, n_stencil,
				   data_for_hypre_semi_struct -> stencil_element_list, matrix_entries );
  free( matrix_entries );

  HYPRE_SStructMatrixAssemble( data_for_hypre_semi_struct -> A );
  HYPRE_SStructMatrixPrint( "matrix.dat", data_for_hypre_semi_struct -> A, 1 );

  HYPRE_SStructMatrixGetObject( data_for_hypre_semi_struct -> A,
				    (void **) &( data_for_hypre_semi_struct -> Aij ) );

  return data_for_hypre_semi_struct;
  
}

void ssp_hypre_semi_struct_pfmg_solve_gmres( struct ssp_hypre_semi_struct *data_for_hypre_semi_struct, int n1, int n2, int n3, double rhs[ n3 ][ n2 ][ n1 ],
				  double soln[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info ) {

  /* Solve the equations using the SYSPFMG solver */
  HYPRE_SStructVector b;
  HYPRE_SStructVector x;
  HYPRE_SStructSolver solver;

  int retval;
  
  /* Assume everything worked for the moment */
  *info = 0;

  /* Create the RHS */
  HYPRE_SStructVectorCreate( data_for_hypre_semi_struct -> comm, data_for_hypre_semi_struct -> grid, &b );
  /* Put the RHS from the Fortran array into the HYPRE object */
  HYPRE_SStructVectorSetObjectType( b, HYPRE_PARCSR );
  HYPRE_SStructVectorInitialize( b );
  HYPRE_SStructVectorSetBoxValues( b, 0, data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub, 0,
				   &( rhs[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the RHS */
  HYPRE_SStructVectorAssemble( b );
  HYPRE_SStructVectorPrint( "rhs.dat", b, 1 );
  
  /* Create the initial guess at solution */
  HYPRE_SStructVectorCreate( data_for_hypre_semi_struct -> comm, data_for_hypre_semi_struct -> grid, &x );
  /* Put the initial guess from the Fortran array into the HYPRE object */
  printf( "!!!%f\n", soln[ 0 ][ 0 ][ 0 ] );
  HYPRE_SStructVectorSetObjectType( x, HYPRE_PARCSR );
  HYPRE_SStructVectorInitialize( x );
  HYPRE_SStructVectorSetBoxValues( x, 0, data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub, 0,
				   &( soln[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the initial guess */
  HYPRE_SStructVectorAssemble( x );
  HYPRE_SStructVectorPrint( "guess.dat", x, 1 );

  /* Create the solver */
  HYPRE_SStructGMRESCreate( data_for_hypre_semi_struct -> comm, &solver );

  /* Set up the solver */
  /* HYPRE_SStructGMRESSetTol( solver, 5.0e-8 ); */
  HYPRE_SStructGMRESSetLogging( solver, 100 );
  HYPRE_SStructGMRESSetPrintLevel(solver, 100 );
  HYPRE_SStructGMRESSetMaxIter( solver, 1000 );
  /* HYPRE_SStructGMRESSetRAPType( solver,  1 ); */
  /* HYPRE_SStructGMRESSetMaxLevels( solver,  9 ); */
  /*  HYPRE_SStructGMRESSetRelaxType( solver,  2 ); */
  printf( "setup\n" );
  retval = HYPRE_SStructGMRESSetup( solver, data_for_hypre_semi_struct -> A, b, x );
  printf( "solve 1 %d\n", retval );

  /* Solve the equations */
  retval = HYPRE_SStructGMRESSolve( solver, data_for_hypre_semi_struct -> A, b, x );
  printf( "solve 2 %d\n", retval );
  HYPRE_SStructVectorPrint( "soln.dat", x, 1 );

  /* Get the solution from the HYPRE object into the Fortran array */
  retval = HYPRE_SStructVectorGetBoxValues( x, 0, data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub, 0,
					    &( soln[ 0 ][ 0 ][ 0 ] ) );
  printf( "extract soln %d\n", retval );

  /* Get some interesting data */
  HYPRE_SStructGMRESGetNumIterations( solver, n_iter );
  retval = HYPRE_SStructGMRESGetFinalRelativeResidualNorm( solver, residual );
  printf( "Residual %d %f\n", retval, *residual );

  /* Destroy the solver */
  HYPRE_SStructGMRESDestroy( solver );
  
  /* Destroy the vectors */
  HYPRE_SStructVectorDestroy( x );
  HYPRE_SStructVectorDestroy( b );
  
}

void ssp_hypre_semi_struct_free( struct ssp_hypre_semi_struct *data_for_hypre_semi_struct ){

  /* Free up the memory used for the grid, stencil and matrix */

  HYPRE_SStructMatrixDestroy ( data_for_hypre_semi_struct -> A       );
  HYPRE_SStructGraphDestroy  ( data_for_hypre_semi_struct -> graph       );
  free                       ( data_for_hypre_semi_struct -> stencil_element_list );
  HYPRE_SStructStencilDestroy( data_for_hypre_semi_struct -> stencil );
  HYPRE_SStructGridDestroy   ( data_for_hypre_semi_struct -> grid    );
    
  free( data_for_hypre_semi_struct );
  
  return;
}

void ssp_hypre_semi_struct_pfmg_solve( struct ssp_hypre_semi_struct *data_for_hypre_semi_struct, int n1, int n2, int n3, double rhs[ n3 ][ n2 ][ n1 ],
				  double soln[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info ) {

  /* Solve the equations using the SYSPFMG solver */
  HYPRE_SStructVector b;
  HYPRE_SStructVector x;
  HYPRE_ParVector bij;
  HYPRE_ParVector xij;
  HYPRE_SStructSolver solver;

  int retval;
  
  /* Assume everything worked for the moment */
  *info = 0;

  /* Create the RHS */
  HYPRE_SStructVectorCreate( data_for_hypre_semi_struct -> comm, data_for_hypre_semi_struct -> grid, &b );
  /* Put the RHS from the Fortran array into the HYPRE object */
  HYPRE_SStructVectorInitialize( b );
  HYPRE_SStructVectorSetBoxValues( b, 0, data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub, 0,
				   &( rhs[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the RHS */
  HYPRE_SStructVectorAssemble( b );
  HYPRE_SStructVectorPrint( "rhs.dat", b, 1 );
  
  /* Create the initial guess at solution */
  HYPRE_SStructVectorCreate( data_for_hypre_semi_struct -> comm, data_for_hypre_semi_struct -> grid, &x );
  /* Put the initial guess from the Fortran array into the HYPRE object */
  printf( "!!!%f\n", soln[ 0 ][ 0 ][ 0 ] );
  HYPRE_SStructVectorInitialize( x );
  HYPRE_SStructVectorSetBoxValues( x, 0, data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub, 0,
				   &( soln[ 0 ][ 0 ][ 0 ] ) );
  /* Assemble the initial guess */
  HYPRE_SStructVectorAssemble( x );
  HYPRE_SStructVectorPrint( "guess.dat", x, 1 );

  /* Create the solver */
  HYPRE_SStructGMRESCreate( data_for_hypre_semi_struct -> comm, &solver );

  /* Set up the solver */
  /* HYPRE_SStructGMRESSetTol( solver, 5.0e-8 ); */
  HYPRE_SStructGMRESSetLogging( solver, 100 );
  HYPRE_SStructGMRESSetPrintLevel(solver, 100 );
  HYPRE_SStructGMRESSetMaxIter( solver, 1000 );
  /* HYPRE_SStructGMRESSetRAPType( solver,  1 ); */
  /* HYPRE_SStructGMRESSetMaxLevels( solver,  9 ); */
  /*  HYPRE_SStructGMRESSetRelaxType( solver,  2 ); */
  printf( "setup\n" );
  retval = HYPRE_SStructGMRESSetup( solver, data_for_hypre_semi_struct -> A, b, x );
  printf( "solve 1 %d\n", retval );

  /* Solve the equations */
  retval = HYPRE_SStructGMRESSolve( solver, data_for_hypre_semi_struct -> A, b, x );
  printf( "solve 2 %d\n", retval );
  HYPRE_SStructVectorPrint( "soln.dat", x, 1 );

  /* Get the solution from the HYPRE object into the Fortran array */
  retval = HYPRE_SStructVectorGetBoxValues( x, 0, data_for_hypre_semi_struct -> lb, data_for_hypre_semi_struct -> ub, 0,
					    &( soln[ 0 ][ 0 ][ 0 ] ) );
  printf( "extract soln %d\n", retval );

  /* Get some interesting data */
  HYPRE_SStructGMRESGetNumIterations( solver, n_iter );
  retval = HYPRE_SStructGMRESGetFinalRelativeResidualNorm( solver, residual );
  printf( "Residual %d %f\n", retval, *residual );

  /* Destroy the solver */
  HYPRE_SStructGMRESDestroy( solver );
  
  /* Destroy the vectors */
  HYPRE_SStructVectorDestroy( x );
  HYPRE_SStructVectorDestroy( b );
  
}
