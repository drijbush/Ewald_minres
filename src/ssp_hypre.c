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

void ssp_hypre_struct_pfmg_solve( struct ssp_hypre_struct *data_for_hypre_struct, int n1, int n2, int n3,
				  double rtol, double rhs[ n3 ][ n2 ][ n1 ],
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
  /* HYPRE_StructVectorPrint( "rhs.dat", b, 1 );  */
  
  /* Create the initial guess at solution */
  HYPRE_StructVectorCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, &x );
  /* Put the initial guess from the Fortran array into the HYPRE object */
  HYPRE_StructVectorInitialize( x );
  HYPRE_StructVectorSetBoxValues( x, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( soln[ 0 ][ 0 ][ 0 ] ) ); 
  /* Assemble the initial guess */
  HYPRE_StructVectorAssemble( x ); 
  /* HYPRE_StructVectorPrint( "guess.dat", x, 1 ); */

  /* Create the solver */
  HYPRE_StructPFMGCreate( data_for_hypre_struct -> comm, &solver );

  /* Set up the solver */
  /* HYPRE_StructPFMGSetZeroGuess( solver ); */
  HYPRE_StructPFMGSetTol( solver, rtol ); 
  HYPRE_StructPFMGSetLogging( solver, 1 );
  /* HYPRE_StructPFMGSetPrintLevel( solver, 1 ); */
  /* HYPRE_StructPFMGSetMaxIter( solver, 100 ); */
  /* HYPRE_StructPFMGSetRAPType( solver,  1 ); */
  /* HYPRE_StructPFMGSetMaxLevels( solver,  9 ); */
  /* HYPRE_StructPFMGSetRelaxType( solver,  2 ); */
  HYPRE_StructPFMGSetup( solver, data_for_hypre_struct -> A, b, x );

  /* Solve the equations */
  HYPRE_StructPFMGSolve( solver, data_for_hypre_struct -> A, b, x );

  /* Get the solution from the HYPRE object into the Fortran array */
  HYPRE_StructVectorGetBoxValues( x, data_for_hypre_struct -> lb, data_for_hypre_struct -> ub, &( soln[ 0 ][ 0 ][ 0 ] ) );
  /* HYPRE_StructVectorPrint( "soln.dat", x, 1 ); */

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
  // HYPRE_ParVector bij;
  // HYPRE_ParVector xij;
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






















struct ssp_hypre_ij *ssp_hypre_ij_setup( int comm, int n[ 3 ], int lb[ 3 ], int ub[ 3 ],
					 int n_stencil, int stencil_elements[ n_stencil ][ 3 ], double stencil_values[ n_stencil ] ) {

  /* Set up the grid, stencil and matrix */
  struct ssp_hypre_ij *data_for_hypre_ij;

  int n_loc[ 3 ];

  int *row_sizes;
  int *col_indices;
  int *col_base_indices;

  int n_tot;
  int n_my_rows;
  int n_init;
  int my_first_row, my_last_row;
  int row, col;
  int i_stencil;

  /* Order of matrix */
  n_tot = n[ 0 ] * n[ 1 ] * n[ 2 ];

  /* Local dimensions */
  n_loc[ 0 ] = ub[ 0 ] - lb[ 0 ] + 1;
  n_loc[ 1 ] = ub[ 1 ] - lb[ 1 ] + 1;
  n_loc[ 2 ] = ub[ 2 ] - lb[ 2 ] + 1;
  
  data_for_hypre_ij = malloc( sizeof( *data_for_hypre_ij ) );

  data_for_hypre_ij -> comm = MPI_Comm_f2c( comm );

  /* Currently distribute the matrix by rows */
  /* For this need to know what my first row is */
  n_my_rows  = n_loc[ 0 ] * n_loc[ 1 ] * n_loc[ 2 ];
  /* Use MPI_Scan to find the global index of my first row */
  MPI_Scan( &n_my_rows, &my_first_row, 1, MPI_INT, MPI_SUM, data_for_hypre_ij -> comm );
  /* Correct inclusive -> exclusive scan */
  my_first_row = my_first_row - n_my_rows;
  my_last_row  = my_first_row + n_my_rows - 1;
  printf( "First last rows %d %d %d %d\n", my_first_row, my_last_row, n_tot, n_stencil );
  
  /* Create the matrix - Distribute currently by rows */
  data_for_hypre_ij -> ilower = my_first_row;
  data_for_hypre_ij -> iupper = my_last_row;
  data_for_hypre_ij -> jlower = 0;
  data_for_hypre_ij -> jupper = n_tot - 1;
  HYPRE_IJMatrixCreate( data_for_hypre_ij -> comm,
			data_for_hypre_ij -> ilower, data_for_hypre_ij -> iupper,
			data_for_hypre_ij -> jlower, data_for_hypre_ij -> jupper,
			&( data_for_hypre_ij -> A ) );
  /* local # rows or global # rows? Doc not clear */
  n_init = n_tot; /* CHECK THIS */
  row_sizes = malloc( n_init * sizeof( *row_sizes ) );
  for( row = 0; row < n_init; row ++ )
    row_sizes[ row ] = n_stencil;
  HYPRE_IJMatrixSetRowSizes  ( data_for_hypre_ij -> A, row_sizes );
  free( row_sizes );

  /* The object type */
  HYPRE_IJMatrixSetObjectType( data_for_hypre_ij -> A, HYPRE_PARCSR );
  printf( "!! Set type\n" );

  /* Prepare matrix for setting */
  HYPRE_IJMatrixInitialize( data_for_hypre_ij -> A );
  printf( "!! Init\n" );

  /* Set matrix a row at a time */
  /* set up the base column indices */  
  col_base_indices = malloc( n_stencil * sizeof( *col_base_indices ) );
  for( i_stencil = 0; i_stencil < n_stencil; i_stencil++ )
    col_base_indices[ i_stencil ] =
      stencil_elements[ i_stencil ][ 2 ]                           +
      stencil_elements[ i_stencil ][ 1 ] * n_loc[ 2 ]              +
      stencil_elements[ i_stencil ][ 0 ] * n_loc[ 2 ] * n_loc[ 1 ];
  printf( "!! Base index\n" );
  
  col_indices = malloc( n_stencil * sizeof( *col_indices ) );
  for( row = my_first_row; row <= my_last_row; row++ ){
    /* For this row work out what the actual indices are */
    for( i_stencil = 0; i_stencil < n_stencil; i_stencil++ ){
      /* Because C is a stupid language the % operator is not well defined for negative values
	 Hence do it this way */
      col = row + col_base_indices[ i_stencil ];
      if( col > n_tot - 1 )
	col = col - n_tot;
      else if( col < 0 )
	col = col + n_tot;
      if( col < 0 || col > n_tot - 1 ) printf( "Wibble!!\n" );
      /* And store the value */
      col_indices[ i_stencil ] = col;
    }
    /* Set this row. Use AddTo to copy with case where stencil wraps around to add to an already
       intialised location - as can happen for small grids compared to the stencil */
    HYPRE_IJMatrixAddToValues( data_for_hypre_ij -> A, 1, &n_stencil, &row, col_indices, stencil_values );
  }
  free( col_indices );
  free( col_base_indices );
  printf( "!! Add loop\n" );

  /* Assemble the matrix and get the ParCSR handle */
  HYPRE_IJMatrixAssemble ( data_for_hypre_ij -> A );
  printf( "!! assemble\n" );

  HYPRE_IJMatrixPrint( data_for_hypre_ij -> A, "ij_matrix.dat" );

  
  HYPRE_IJMatrixGetObject( data_for_hypre_ij -> A, (void **) & (data_for_hypre_ij -> Acsr )  );

  printf( "!! exit set up\n" );
  
  return data_for_hypre_ij;
  
}

void ssp_hypre_ij_solve( struct ssp_hypre_ij *data_for_hypre_ij, int n1, int n2, int n3, double b[ n3 ][ n2 ][ n1 ],
				  double x[ n3 ][ n2 ][ n1 ], int *n_iter, double *residual, int *info )
{

  HYPRE_IJVector  bij;
  HYPRE_ParVector bij_csr;
  HYPRE_IJVector  xij;
  HYPRE_ParVector xij_csr;

  HYPRE_Solver solver;
  
  int *indices;

  int nvalues;
  int i; 
  
  /* Assume everything worked for the moment */
  *info = 0;      

  /* Set up the indices list */
  nvalues = data_for_hypre_ij -> iupper - data_for_hypre_ij -> ilower + 1;
  indices = malloc( nvalues * sizeof( *indices ) );
  for( i = 0; i < nvalues; i++ )
    indices[ i ] = data_for_hypre_ij -> ilower + i; 
    /*   indices[ i ] = i; */


  /* Create the RHS */
  HYPRE_IJVectorCreate( data_for_hypre_ij -> comm, data_for_hypre_ij -> ilower, data_for_hypre_ij -> iupper,
			&bij );
  HYPRE_IJVectorSetObjectType( bij, HYPRE_PARCSR);

  /* Get ready to set up RHS vector */
  HYPRE_IJVectorInitialize( bij );

  /* Set the RHS */
  HYPRE_IJVectorSetValues( bij, nvalues, indices, &( b[ 0 ][ 0 ][ 0 ] ) );

  /* Assemble the RHS and get the object handle */
  HYPRE_IJVectorAssemble ( bij );
  HYPRE_IJVectorGetObject( bij, ( void **) & bij_csr );

  HYPRE_IJVectorPrint( bij, "rhs_before.dat" );


  /* Create the Solution */
  HYPRE_IJVectorCreate( data_for_hypre_ij -> comm, data_for_hypre_ij -> ilower, data_for_hypre_ij -> iupper,
			&xij );
  HYPRE_IJVectorSetObjectType( xij, HYPRE_PARCSR);

  /* Get ready to set up Solution vector */
  HYPRE_IJVectorInitialize( xij );

  /* Set the intial guess at the Solution */
  HYPRE_IJVectorSetValues( xij, nvalues, indices, &( x[ 0 ][ 0 ][ 0 ] ) );

  /* Assemble the initial guess and get the object handle */
  HYPRE_IJVectorAssemble ( xij );
  HYPRE_IJVectorGetObject( xij, ( void **) & xij_csr );

  HYPRE_IJVectorPrint( xij, "soln_before.dat" );

  /* Solve the equations */
  HYPRE_BoomerAMGCreate( &solver );
  HYPRE_BoomerAMGSetTol( solver, 1e-12 );
  HYPRE_BoomerAMGSetup( solver, data_for_hypre_ij -> Acsr, bij_csr, xij_csr );
  HYPRE_BoomerAMGSolve( solver, data_for_hypre_ij -> Acsr, bij_csr, xij_csr );
  HYPRE_BoomerAMGGetFinalRelativeResidualNorm( solver, residual );
  HYPRE_BoomerAMGGetNumIterations( solver, n_iter );
  HYPRE_BoomerAMGDestroy( solver );

  /* Get back the solution */
  HYPRE_IJVectorGetValues( xij, nvalues, indices, &( x[ 0 ][ 0 ][ 0 ] ) );  

  HYPRE_IJVectorPrint( xij, "soln_after.dat" );

  /* Tidy up */
  HYPRE_IJVectorDestroy( xij );
  HYPRE_IJVectorDestroy( bij );
  free( indices );
  
}

void ssp_hypre_ij_free( struct ssp_hypre_ij *data_for_hypre_ij ){

  /* Free up the memory used for the grid, stencil and matrix */

  HYPRE_IJMatrixDestroy ( data_for_hypre_ij -> A );
    
  free( data_for_hypre_ij );
  
  return;
}
