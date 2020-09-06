#include <stdlib.h>

#include "mpi.h"

#include "HYPRE_struct_mv.h"
#include "HYPRE_struct_ls.h"

#include "ssp_hypre.h"

struct ssp_hypre_struct *ssp_hypre_struct_setup( int comm, int lb[ 3 ], int ub[ 3 ],
				   int n_stencil, int stencil_elements[][3], double stencil_values[] ) {

  struct ssp_hypre_struct *data_for_hypre_struct;

  int *ele_list;
  
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
  HYPRE_StructGridSetExtents( data_for_hypre_struct -> grid,
			      data_for_hypre_struct -> lb, data_for_hypre_struct -> ub );
  HYPRE_StructGridAssemble( data_for_hypre_struct -> grid );

  /* The stencil */
  HYPRE_StructStencilCreate( 3, n_stencil, &( data_for_hypre_struct -> stencil ) );
  /* Note reverse order of dimensions to cope with C <-> Fortran mappings */
  ele_list = malloc( n_stencil * sizeof( *ele_list ) );
  for( i = 0; i < n_stencil; i++ ) {
    ele[ 0 ] = stencil_elements[ i ][ 2 ];
    ele[ 1 ] = stencil_elements[ i ][ 1 ];
    ele[ 2 ] = stencil_elements[ i ][ 0 ];
    ele_list[ i ] = i;
    HYPRE_StructStencilSetElement( data_for_hypre_struct -> stencil, i, ele );
  }

  /* The matrix */
  HYPRE_StructMatrixCreate( data_for_hypre_struct -> comm, data_for_hypre_struct -> grid, data_for_hypre_struct -> stencil,
			    &( data_for_hypre_struct -> A ) );
  HYPRE_StructMatrixSetSymmetric( data_for_hypre_struct -> A, 1 );
  HYPRE_StructMatrixInitialize( data_for_hypre_struct -> A );
  HYPRE_StructMatrixSetConstantValues( data_for_hypre_struct -> A, n_stencil, ele_list, stencil_values );
  HYPRE_StructMatrixAssemble( data_for_hypre_struct -> A );

  free( ele_list );
  
  return data_for_hypre_struct;
  
}

void ssp_hypre_struct_free( struct ssp_hypre_struct *data_for_hypre_struct ){

  HYPRE_StructMatrixDestroy ( data_for_hypre_struct -> A       );
  HYPRE_StructStencilDestroy( data_for_hypre_struct -> stencil );
  HYPRE_StructGridDestroy   ( data_for_hypre_struct -> grid    );
    
  free( data_for_hypre_struct );
  
  return;
}
