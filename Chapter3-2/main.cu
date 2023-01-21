#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "param.h"
#include "go.cuh"
#include "go_solve.cuh"
#include "go_ion.cuh"
#include "gr.cuh"
#include "gr_solve.cuh"
#include "gr_ion.cuh"
#include "pkj.cuh"
#include "pkj_solve.cuh"
#include "pkj_ion.cuh"
#include "io.cuh"
#include "io_solve.cuh"
#include "io_ion.cuh"
#include "solve_bem.cuh"
#include "solve_rkc.cuh"
#include "syn.cuh"
#include "gap.cuh"
#include "output.cuh"
#include <cuda_profiler_api.h>


// set solver name in argv [ 1 ] from [ BE or CN or RKC ] for GR
// set GR_X and GR_Y in argv [ 2 ] and [ 3 ] 

// set solver name in argv [ 4 ] from [ BE or CN or RKC ] for GO
// set GO_X and GO_Y in argv [ 5 ] and [ 6 ] 

// set solver name in argv [ 7 ] from [ BE or CN or RKC ] for PKJ
// set PKJ_X and PKJ_Y in argv [ 8 ] and [ 9 ] 

int main ( int argc, char *argv [] ) 
{
  //srand ( ( int ) time ( NULL ) );
  srand ( 0 );
  bool network_sim = false;
  // Debug
  printf ( "argc = %d\n", argc );
  for ( int i = 0; i < argc; i++ ) {
    printf( "argv[%d] = %s\n", i, argv [ i ] );
  }

  // set mh_tau_rate
  double mh_tau_rate = atof ( argv [ 13 ] );
  printf ( "%lf\n", mh_tau_rate );

  // set solver name 
  char solver_name_gr [ 128 ], solver_name_go [ 128 ], solver_name_pkj [ 128 ], solver_name_io [ 128 ];
  strcpy ( solver_name_gr,  argv [ 1 ] );
  strcpy ( solver_name_go,  argv [ 4 ] );
  strcpy ( solver_name_pkj, argv [ 7 ] );
  strcpy ( solver_name_io,  argv [ 10 ] );
  printf ( "GR  solver -> " ); puts ( solver_name_gr );  
  printf ( "GO  solver -> " ); puts ( solver_name_go );  
  printf ( "PKJ solver -> " ); puts ( solver_name_pkj ); 
  printf ( "IO  solver -> " ); puts ( solver_name_io );  printf ("\n");

  FILE *f_gr_v; FILE *f_go_v; FILE *f_pkj_v; FILE *f_io_v; 
  FILE *f_gr_raster; FILE *f_go_raster; FILE *f_pkj_raster; FILE *f_io_raster; 

  if ( 0 == strncmp ( solver_name_gr, "BE",  2 ) ) {
    f_gr_v = fopen  ( "gr_v_by_BE.csv",  "w" );  
    f_gr_raster = fopen  ( "gr_raster_by_BE.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_gr, "CN",  2 ) ) {
    f_gr_v = fopen  ( "gr_v_by_CN.csv",  "w" );  
    f_gr_raster = fopen  ( "gr_raster_by_CN.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_gr, "RKC",  3 ) ) {
    f_gr_v = fopen  ( "gr_v_by_RKC.csv",  "w" );  
    f_gr_raster = fopen  ( "gr_raster_by_RKC.csv",  "w" );  
  } else { printf ( "Error!! -> set solver name in argv [ 1 ]\n" ); return 0; }

  if ( 0 == strncmp ( solver_name_go, "BE",  2 ) ) {
    f_go_v = fopen  ( "go_v_by_BE.csv",  "w" );  
    f_go_raster = fopen  ( "go_raster_by_BE.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_go, "CN",  2 ) ) {
    f_go_v = fopen  ( "go_v_by_CN.csv",  "w" );  
    f_go_raster = fopen  ( "go_raster_by_CN.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_go, "RKC",  3 ) ) {
    f_go_v = fopen  ( "go_v_by_RKC.csv",  "w" );  
    f_go_raster = fopen  ( "go_raster_by_RKC.csv",  "w" );  
  } else { printf ( "Error!! -> set solver name in argv [ 4 ]\n" ); return 0; }

  if ( 0 == strncmp ( solver_name_pkj, "BE",  2 ) ) {
    f_pkj_v = fopen  ( "pkj_v_by_BE.csv",  "w" );  
    f_pkj_raster = fopen  ( "pkj_raster_by_BE.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_pkj, "CN",  2 ) ) {
    f_pkj_v = fopen  ( "pkj_v_by_CN.csv",  "w" );  
    f_pkj_raster = fopen  ( "pkj_raster_by_CN.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_pkj, "RKC",  3 ) ) {
    f_pkj_v = fopen  ( "pkj_v_by_RKC.csv",  "w" );  
    f_pkj_raster = fopen  ( "pkj_raster_by_RKC.csv",  "w" );  
  } else { printf ( "Error!! -> set solver name in argv [ 7 ]\n" ); return 0; }
  
  if ( 0 == strncmp ( solver_name_io, "BE",  2 ) ) {
    f_io_v = fopen  ( "io_v_by_BE.csv",  "w" );  
    f_io_raster = fopen  ( "io_raster_by_BE.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_io, "CN",  2 ) ) {
    f_io_v = fopen  ( "io_v_by_CN.csv",  "w" );  
    f_io_raster = fopen  ( "io_raster_by_CN.csv",  "w" );  
  } else if ( 0 == strncmp ( solver_name_io, "RKC",  3 ) ) {
    f_io_v = fopen  ( "io_v_by_RKC.csv",  "w" );  
    f_io_raster = fopen  ( "io_raster_by_RKC.csv",  "w" );  
  } else { printf ( "Error!! -> set solver name in argv [ 10 ]\n" ); return 0; }
    
  // set gr
  int n_gr = atoi ( argv [ 2 ] ) * atoi ( argv [ 3 ] );
  //printf ( "n_gr -> %d\n", n_gr );
  neuron_t        *p_gr       = ( neuron_t * ) malloc ( sizeof ( neuron_t ) );
  neuron_solve_t  *p_gr_solve = ( neuron_solve_t * ) malloc ( sizeof ( neuron_solve_t ) );
  neuron_t        *d_gr       = gr_initialize ( atoi ( argv [ 2 ] ), atoi ( argv [ 3 ] ), solver_name_gr, p_gr );  
  neuron_solve_t  *d_gr_solve = gr_solve_initialize ( p_gr_solve, solver_name_gr, p_gr, d_gr );

  // set go
  int n_go = atoi ( argv [ 5 ] ) * atoi ( argv [ 6 ] );
  //printf ( "n_go -> %d\n", n_go );
  neuron_t        *p_go       = ( neuron_t * ) malloc ( sizeof ( neuron_t ) );
  neuron_solve_t  *p_go_solve = ( neuron_solve_t * ) malloc ( sizeof ( neuron_solve_t ) );
  neuron_t        *d_go       = go_initialize ( atoi ( argv [ 5 ] ), atoi ( argv [ 6 ] ), solver_name_go, p_go );  
  neuron_solve_t  *d_go_solve = go_solve_initialize ( p_go_solve, solver_name_go, p_go, d_go );
  
  // set pkj
  int n_pkj = atoi ( argv [ 8 ] ) * atoi ( argv [ 9 ] );
  //printf ( "n_pkj -> %d\n", n_pkj );
  neuron_t        *p_pkj       = ( neuron_t * ) malloc ( sizeof ( neuron_t ) );
  neuron_solve_t  *p_pkj_solve = ( neuron_solve_t * ) malloc ( sizeof ( neuron_solve_t ) );
  neuron_t        *d_pkj       = pkj_initialize ( atoi ( argv [ 8 ] ), atoi ( argv [ 9 ] ), solver_name_pkj, p_pkj );  
  neuron_solve_t  *d_pkj_solve = pkj_solve_initialize ( p_pkj_solve, solver_name_pkj, p_pkj, d_pkj );
  
  // set io
  int n_io = atoi ( argv [ 11 ] ) * atoi ( argv [ 12 ] );
  //printf ( "n_io -> %d\n", n_io );
  neuron_t        *p_io       = ( neuron_t * ) malloc ( sizeof ( neuron_t ) );
  neuron_solve_t  *p_io_solve = ( neuron_solve_t * ) malloc ( sizeof ( neuron_solve_t ) );
  neuron_t        *d_io       = io_initialize ( atoi ( argv [ 11 ] ), atoi ( argv [ 12 ] ), solver_name_io, p_io );  
  neuron_solve_t  *d_io_solve = io_solve_initialize ( p_io_solve, solver_name_io, p_io, d_io );
  
  printf ( "\n" );

  // create synapse
  synapse_t *d_mfgr   = mfgr_create  ( n_gr );
  synapse_t *d_gogr   = gogr_create  ( p_go -> nx, p_go -> ny, p_gr -> nx, p_gr -> ny );
  synapse_t *d_grgo   = grgo_create  ( p_gr -> nx, p_gr -> ny, p_go -> nx, p_go -> ny );
  synapse_t *d_grpkj  = grpkj_create ( p_gr -> nx, p_gr -> ny, p_pkj -> nx, p_pkj -> ny );
  synapse_t *d_mlipkj = mlipkj_create ( n_pkj, n_gr );
  
  // create gap_junction ( GJ )
  gap_t *d_io_gap = io_gap_create ( atoi ( argv [ 11 ] ), atoi ( argv [ 12 ] ) );

  
  double *memv_test_gr = ( double * ) malloc ( p_gr -> nc * sizeof ( double ) );
  if ( p_gr -> nc > 0 )
    cudaMemcpy ( memv_test_gr, p_gr -> elem [ v ], p_gr -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
    double *memv_test_go = ( double * ) malloc ( p_go -> nc * sizeof ( double ) );
  if ( p_go -> nc > 0 )
    cudaMemcpy ( memv_test_go, p_go -> elem [ v ], p_go -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
    double *memv_test_pkj = ( double * ) malloc ( p_pkj -> nc * sizeof ( double ) );
  if ( p_pkj -> nc > 0 )
    cudaMemcpy ( memv_test_pkj, p_pkj -> elem [ v ], p_pkj -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
    double *memv_test_io = ( double * ) malloc ( p_io -> nc * sizeof ( double ) );
  if ( p_io -> nc > 0 )
    cudaMemcpy ( memv_test_io, p_io -> elem [ v ], p_io -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );

  double *memv_raster_gr = ( double * ) malloc ( p_gr -> n * sizeof ( double ) );
  for ( int i = 0; i < p_gr -> n; i++ ) memv_raster_gr [ i ] = 100;
  double *memv_raster_go = ( double * ) malloc ( p_go -> n * sizeof ( double ) );
  for ( int i = 0; i < p_go -> n; i++ ) memv_raster_go [ i ] = 100;
  double *memv_raster_pkj = ( double * ) malloc ( p_pkj -> n * sizeof ( double ) );
  for ( int i = 0; i < p_pkj -> n; i++ ) memv_raster_pkj [ i ] = 100;
  double *memv_raster_io = ( double * ) malloc ( p_io -> n * sizeof ( double ) );
  for ( int i = 0; i < p_io -> n; i++ ) memv_raster_io [ i ] = 100;

      // open input file
  int i;
      FILE *input = fopen("input_pos.txt","r");
      if( input == NULL ){
	printf("input file open error\n");
	exit(1);
      }
      int pulse_num, pat_num;
      int cnt = 0;
      fscanf(input, "%d", &pulse_num);
      fscanf(input, "%d", &pat_num);
      
      int *pos = (int *)malloc( sizeof(int) * pulse_num * pat_num); 
      while(fscanf(input, "%d", &pos[cnt]) != EOF ){cnt++;};
      fclose(input);
      printf("%d",cnt);
      for( i = 0; i<pulse_num*pat_num;i++){
        printf("%d, ", pos[i]);
      }
      printf("\n");

      input = fopen("input_time.txt","r");
      if( input == NULL ){
	printf("input file open error\n");
	exit(1);
      }
      cnt = 0;      
      double *time = (double *)malloc( sizeof(double) * pulse_num * pat_num); 
      while(fscanf(input, "%lf", &time[cnt]) != EOF ){cnt++;};
      fclose(input);
      printf("%d",cnt);
      for( i = 0; i<pulse_num*pat_num;i++){
        printf("%lf, ", time[i]);
      }
      printf("\n");


      input = fopen("input_size.txt","r");
      if( input == NULL ){
	printf("input file open error\n");
	exit(1);
      }
      cnt = 0;      
      double *size = (double *)malloc( sizeof(double) * pulse_num * pat_num); 
      while(fscanf(input, "%lf", &size[cnt]) != EOF ){cnt++;};
      fclose(input);
      printf("%d",cnt);
      for( i = 0; i<pulse_num*pat_num;i++){
        printf("%lf, ", size[i]);
      }
      printf("\n");



      int *d_pos;
      double *d_time, *d_size;
      cudaMalloc( (void **) &d_pos, pulse_num * pat_num * sizeof( int ));
      cudaMalloc( (void **) &d_time, pulse_num * pat_num * sizeof( double ));
      cudaMalloc( (void **) &d_size, pulse_num * pat_num * sizeof( double ));
      cudaMemcpy ( d_pos, pos, pulse_num * pat_num * sizeof ( int ), cudaMemcpyHostToDevice );
      cudaMemcpy ( d_time, time, pulse_num * pat_num * sizeof ( double ), cudaMemcpyHostToDevice );
      cudaMemcpy ( d_size, size, pulse_num * pat_num * sizeof ( double ), cudaMemcpyHostToDevice );
      
      free(pos);
      free(time);
      free(size);

      

  // Debug
  FILE *debug_file;
  debug_file = fopen ( "debug_file.csv", "w" );
  fclose ( debug_file );

  double t = 0.0;
  int t_count = 0;
  clock_t start, half, end;
  start = clock();

  cudaProfilerStart();
  while ( t < TS - 0.0001 ) 
  {
    double l_dt = 0.0;
    if ( t_count % 8 == 0 ) //40 
    {   
      printf ( "time:%f\n", t ); 
      //printf(  "DT = %f \n ", p_go -> DT ); 
    }
 
    // update gr
    if ( n_gr > 0 ) 
    {
      l_dt = 0.0;
      while ( l_dt < 0.125 ) 
      {
        if ( network_sim == false ) {
          gr_set_current <<< ( int ) ( p_gr -> n / 128 ) + 1, 128 >>> ( d_gr, t + l_dt );
        }
        gr_solve_update_v ( d_gr, d_gr_solve, p_gr, p_gr_solve, d_mfgr, d_gogr );      
        if ( isnan ( memv_test_gr [ 1 ] ) ) { printf ( "isnan in gr, time=%f\n", t ); exit ( 1 ); }
        l_dt += p_gr -> DT;
      }
      switch ( p_gr -> n_solver ) 
      {
        case kBE:
        case kCN:
          cudaMemcpy ( memv_test_gr, p_gr -> elem [ v ], p_gr -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
          break;
        case kRKC:
          cudaMemcpy ( memv_test_gr, p_gr_solve -> vec [ v_new ], p_gr -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
      }     
      gr_output_file ( p_gr, memv_test_gr, memv_raster_gr, t + l_dt, f_gr_v, f_gr_raster );
      //mf_output_file ( d_mfgr, t + l_dt, p_gr );
    }
    // update go
    if ( n_go > 0 ) 
    {
      l_dt = 0.0;
      while ( l_dt < 0.125 ) 
      {
        if ( network_sim == false ) {
          go_set_current <<< ( int ) ( p_go -> n / 128 ) + 1, 128 >>> ( d_go, t + l_dt );
        }
        go_solve_update_v ( d_go, d_go_solve, p_go, p_go_solve, d_grgo );
        if ( isnan ( memv_test_go [ 1 ] ) ) { printf ( "isnan in go, time=%f\n", t ); exit ( 1 ); }
        l_dt += p_go -> DT;
      }
      switch ( p_go -> n_solver ) 
      {
        case kBE:
        case kCN:
          cudaMemcpy ( memv_test_go, p_go -> elem [ v ], p_go -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
          break;
        case kRKC:
          cudaMemcpy ( memv_test_go, p_go_solve -> vec [ v_new ], p_go -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
      } 
      go_output_file ( p_go, memv_test_go, memv_raster_go, t + l_dt, f_go_v, f_go_raster );
    }

    // update pkj
    if ( n_pkj > 0 ) 
    {
      l_dt = 0.0;

      // poisson input
      //pkj_set_current_Poisson ( p_pkj -> n, d_pkj, t );

      while ( l_dt < 0.125 ) 
      {
        // constant input
        if ( network_sim == false ) {
          pkj_set_current <<< ( int ) ( p_pkj -> n / 128 ) + 1, 128 >>> ( d_pkj, t + l_dt, pulse_num, pat_num, d_pos, d_time, d_size);
        }
        // update
        pkj_solve_update_v ( d_pkj, d_pkj_solve, p_pkj, p_pkj_solve, d_grpkj, d_mlipkj, mh_tau_rate );
        if ( isnan ( memv_test_pkj [ 1 ] ) ) { printf ( "isnan in pkj, time=%f\n", t ); exit ( 1 ); }
        l_dt += p_pkj -> DT;
      }
      switch ( p_pkj -> n_solver ) 
      {
        case kBE:
        case kCN:
          cudaMemcpy ( memv_test_pkj, p_pkj -> elem [ v ], p_pkj -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
          break;
        case kRKC:
          cudaMemcpy ( memv_test_pkj, p_pkj_solve -> vec [ v_new ], p_pkj -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
      }
      pkj_output_file ( p_pkj, memv_test_pkj, memv_raster_pkj, t + l_dt, f_pkj_v, f_pkj_raster );
      //mli_output_file ( d_mlipkj, t + l_dt, p_pkj );

    }
    
    // update io
    if ( n_io > 0 ) 
    {
      l_dt = 0.0;
      while ( l_dt < 0.125 ) 
      {
        if ( network_sim == false ){
          io_set_current <<< ( int ) ( p_io -> n / 128 ) + 1, 128 >>> ( d_io, t + l_dt );
        }
        io_solve_update_v ( d_io, d_io_solve, p_io, p_io_solve, d_io_gap );
        if ( isnan ( memv_test_io [ 1 ] ) ) { printf ( "isnan,time=%f\n", t ); exit ( 1 ); }
        l_dt += p_io -> DT;
      }      
      switch ( p_io -> n_solver ) 
      {
        case kBE:
        case kCN:
          cudaMemcpy ( memv_test_io, p_io -> elem [ v ], p_io -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
          break;
        case kRKC:
          cudaMemcpy ( memv_test_io, p_io_solve -> vec [ v_new ], p_io -> nc * sizeof ( double ), cudaMemcpyDeviceToHost );
      }
      io_output_file ( p_io, memv_test_io, memv_raster_io, t + l_dt, f_io_v, f_io_raster );
    }
    
    // Update synapses
    if ( network_sim == true ) {
      gr_synapse_update ( t, 0.125, d_mfgr, d_gogr, d_go, p_gr_solve );
      go_synapse_update ( t, 0.125, d_grgo,  d_gr, p_go_solve );
      pkj_synapse_update ( t, 0.125, d_grpkj, d_mlipkj, d_gr, p_pkj_solve );
    }
    t_count++; t += 0.125; 
    // record 1000msec time
    if ( t >= 500.0 && 500.0 + 0.125 > t ) { half = clock(); }
  }
  cudaProfilerStop();

  // Debug
  printf ( "time:%f \n", t );  
  end = clock();
  printf ( "computational time = %.2f \n", ( double ) ( end - start ) / CLOCKS_PER_SEC );

  // output computational time
  if ( TS <= 500.0 )
    output_time ( n_gr, solver_name_gr, n_go, solver_name_go, n_pkj, solver_name_pkj, n_io, solver_name_io, ( double ) ( end - start ) / CLOCKS_PER_SEC, TS );
  else
    output_time2 ( n_gr, solver_name_gr, n_go, solver_name_go, n_pkj, solver_name_pkj, n_io, solver_name_io, ( double ) ( end - start ) / CLOCKS_PER_SEC, ( double ) ( end - half ) / CLOCKS_PER_SEC, TS );
  //raster_plot ( n_go, solver_name );
  //fflush( f_out );

  // synapse_finalize
  mfgr_finalize ( d_mfgr, n_gr );
  gogr_finalize ( d_gogr, n_go * n_gr );
  grgo_finalize ( d_grgo, n_go * n_gr );
  grpkj_finalize ( d_grpkj, n_gr * n_pkj );
  mlipkj_finalize ( d_mlipkj, n_gr * n_pkj );
  io_gap_finalize ( d_io_gap );

  // gr_finalize
  gr_finalize ( p_gr, d_gr, f_gr_v, f_gr_raster );
  gr_solve_finalize ( n_gr, d_gr_solve, p_gr_solve );

  // go_finalize
  go_finalize ( p_go, d_go, f_go_v, f_go_raster );
  go_solve_finalize ( n_go, d_go_solve, p_go_solve );

  // pkj_finalize
  pkj_finalize ( p_pkj, d_pkj, f_pkj_v, f_pkj_raster );
  pkj_solve_finalize ( n_pkj, d_pkj_solve, p_pkj_solve );

  // io_finalize
  io_finalize ( p_io, d_io, f_io_v, f_io_raster );
  io_solve_finalize ( n_io, d_io_solve, p_io_solve );

  free  ( memv_test_gr );   free  ( memv_test_go );   free  ( memv_test_pkj );   free  ( memv_test_io );
  free  ( memv_raster_gr ); free  ( memv_raster_go ); free  ( memv_raster_pkj ); free  ( memv_raster_io );
  printf("%d,%d\n",pulse_num,pat_num);
  //cudaFree  ( d_pos );            cudaFree  ( d_time );           cudaFree  ( d_size );

  return 0;
}
