#include "pkj.cuh"

__host__ static void initialize_host_pkj ( neuron_t *pkj )
{ 
  int nc = pkj -> nc; // # of all compartments
  pkj -> elem = ( double ** ) malloc ( n_elem * sizeof ( double *) );
  for ( int i = 0; i < n_elem; i++ ) {
    pkj -> elem [ i ] = ( double * ) malloc ( nc * sizeof ( double ) );
  }

  pkj -> cond = ( double ** ) malloc ( pkj_n_cond * sizeof ( double *) );
  for ( int i = 0; i < pkj_n_cond; i++ ) {
    pkj -> cond [ i ] = ( double * ) malloc ( nc * sizeof ( double ) );
  }

  pkj -> ion = ( double ** ) malloc ( pkj_n_ion * sizeof ( double *) );
  for ( int i = 0; i < pkj_n_ion; i++ ) {
    pkj -> ion [ i ] = ( double * ) malloc ( nc * sizeof ( double ) );
  }
  
  pkj -> shell = ( double * ) malloc ( nc * sizeof ( double ) );
  //pkj -> ca2 = ( double * ) malloc ( nc * sizeof ( double ) );
  pkj -> rev_ca2 = ( double * ) malloc ( nc * sizeof ( double ) );

  double rad [ PKJ_COMP ], len [ PKJ_COMP ];//, Ra [ PKJ_COMP ];

  FILE *file = fopen ( PARAM_FILE_PKJ, "r" );
  if ( ! file ) { fprintf ( stderr, "no such file %s\n", PARAM_FILE_PKJ ); exit ( 1 ); }

  for ( int i = 0; i < PKJ_COMP; i++ ) {
    int i1, i2, i3;
    double x, y, z, d, area_val;

    if ( fscanf ( file, "%d %d %lf %lf %lf %lf %d ", &i1, &i2, &x, &y, &z, &d, &i3 ) == ( EOF ) )
    {
        printf ( "PKJ_PARAM_FILE_READING_ERROR\n" );
        exit ( 1 );
    }

    if ( i3 == PKJ_soma ) 
    {
      rad [ i1 ] = 0.5 * d * 1e-4; // [mum -> cm]
      len [ i1 ] = d * 1e-4; // [mum -> cm]
      //Ra  [ i1 ] = PKJ_RI; // [kohm-cm]
  
      pkj -> elem [ connect ] [ i1 ] = ( double ) i2;
      pkj -> elem [ compart ] [ i1 ] = ( double ) i3;
      pkj -> elem [ area    ] [ i1 ] = 4.0 * M_PI * rad [ i1 ] * rad [ i1 ]; // [cm^2]
      //pkj -> elem [ area    ] [ i1 ] = 2.0 * M_PI * rad [ i1 ] * len [ i1 ]; // [cm^2]

      area_val = pkj -> elem [ area ] [ i1 ];
      pkj -> elem [ Cm      ] [ i1 ] = PKJ_CM * area_val; // [muF]
      pkj -> elem [ i_ext   ] [ i1 ] = 0.0;
      pkj -> shell [ i1 ] = ( PKJ_SHELLD - 2 * pow ( PKJ_SHELLD, 2 ) / ( d * 1e-4 )
                                           + 4 * pow ( PKJ_SHELLD, 3 ) / ( 3 * pow ( d * 1e-4, 2 ) ) ) * pkj -> elem [ area    ] [ i1 ];
      pkj -> rev_ca2 [ i1 ] = 12.5 * log ( Ca1OUT_PKJ / Ca1_0_PKJ );        
    } 
    else 
    { // pkj_dend // main, spiny, smooth
      rad [ i1 ] = 0.5 * d * 1e-4; // [mum -> cm]
      len [ i1 ] = sqrt(pow(x,2)+pow(y,2)+pow(z,2))*1e-4; // [mum -> cm]
      //Ra  [ i1 ] = PKJ_RI; // [kohm-cm]
  
      pkj -> elem [ connect ] [ i1 ] = ( double ) i2;
      pkj -> elem [ compart ] [ i1 ] = ( double ) i3;
      pkj -> elem [ area    ] [ i1 ] = 2.0 * M_PI * rad [ i1 ] * len [ i1 ]; // [cm^2]
      area_val = pkj -> elem [ area ] [ i1 ];
      pkj -> shell [ i1 ] = ( PKJ_SHELLD - pow ( PKJ_SHELLD, 2 ) / ( d * 1e-4 ) ) * pkj -> elem [ area    ] [ i1 ];
      pkj -> rev_ca2 [ i1 ] = 13.361624877 * log ( Ca1OUT_PKJ / Ca1_0_PKJ );
    
      if ( i3 == PKJ_spiny ) 
      {
        double other_spines = area_val + len [ i1 ] * PKJ_SPINE_AREA * PKJ_SPINE_DENS;
        pkj -> elem [ Cm      ] [ i1 ] = PKJ_CM * other_spines; // [muF]
      }
      else 
      {   
        pkj -> elem [ Cm      ] [ i1 ] = PKJ_CM * area_val; // [muF]
      }
        pkj -> elem [ i_ext   ] [ i1 ] = 0.0;
    }

    if ( i3 == PKJ_soma )
    {
      pkj -> cond [ g_leak_pkj ] [ i1 ] = 1.0 / PKJ_RM_S * area_val;
      pkj -> cond [ g_NaF_pkj ]  [ i1 ] = PKJ_S_G_NaF * area_val;
      pkj -> cond [ g_NaP_pkj ]  [ i1 ] = PKJ_S_G_NaP * area_val;
      pkj -> cond [ g_CaP_pkj ]  [ i1 ] = PKJ_S_G_CaP * area_val;
      pkj -> cond [ g_CaT_pkj ]  [ i1 ] = PKJ_S_G_CaT * area_val;
      pkj -> cond [ g_Kh_pkj ]   [ i1 ] = PKJ_S_G_KH * area_val;
      pkj -> cond [ g_Kdr_pkj ]  [ i1 ] = PKJ_S_G_KDR * area_val;
      pkj -> cond [ g_KM_pkj ]   [ i1 ] = PKJ_S_G_KM * area_val;
      pkj -> cond [ g_KA_pkj ]   [ i1 ] = PKJ_S_G_KA * area_val;
      pkj -> cond [ g_KC_pkj ]   [ i1 ] = PKJ_S_G_KC * area_val;
      pkj -> cond [ g_K2_pkj ]   [ i1 ] = PKJ_S_G_K2 * area_val;
    } 
    else if ( i3 == PKJ_main ) 
    {
      pkj -> cond [ g_leak_pkj ] [ i1 ] = 1.0 / PKJ_RM_D * area_val;
      pkj -> cond [ g_NaF_pkj ]  [ i1 ] = PKJ_MD_G_NaF * area_val;
      pkj -> cond [ g_NaP_pkj ]  [ i1 ] = PKJ_MD_G_NaP * area_val;
      pkj -> cond [ g_CaP_pkj ]  [ i1 ] = PKJ_MD_G_CaP * area_val;
      pkj -> cond [ g_CaT_pkj ]  [ i1 ] = PKJ_MD_G_CaT * area_val;
      pkj -> cond [ g_Kh_pkj ]   [ i1 ] = PKJ_MD_G_KH * area_val;
      pkj -> cond [ g_Kdr_pkj ]  [ i1 ] = PKJ_MD_G_KDR * area_val;
      pkj -> cond [ g_KM_pkj ]   [ i1 ] = PKJ_MD_G_KM * area_val;
      pkj -> cond [ g_KA_pkj ]   [ i1 ] = PKJ_MD_G_KA * area_val;
      pkj -> cond [ g_KC_pkj ]   [ i1 ] = PKJ_MD_G_KC * area_val;
      pkj -> cond [ g_K2_pkj ]   [ i1 ] = PKJ_MD_G_K2 * area_val;
    } 
    else if ( i3 == PKJ_smooth )
    {
      pkj -> cond [ g_leak_pkj ] [ i1 ] = 1.0 / PKJ_RM_D * area_val;
      pkj -> cond [ g_NaF_pkj ]  [ i1 ] = PKJ_RD_G_NaF * area_val;
      pkj -> cond [ g_NaP_pkj ]  [ i1 ] = PKJ_RD_G_NaP * area_val;
      pkj -> cond [ g_CaP_pkj ]  [ i1 ] = PKJ_RD_G_CaP * area_val;
      pkj -> cond [ g_CaT_pkj ]  [ i1 ] = PKJ_RD_G_CaT * area_val;
      pkj -> cond [ g_Kh_pkj ]   [ i1 ] = PKJ_RD_G_KH * area_val;
      pkj -> cond [ g_Kdr_pkj ]  [ i1 ] = PKJ_RD_G_KDR * area_val;
      pkj -> cond [ g_KM_pkj ]   [ i1 ] = PKJ_RD_G_KM * area_val;
      pkj -> cond [ g_KA_pkj ]   [ i1 ] = PKJ_RD_G_KA * area_val;
      pkj -> cond [ g_KC_pkj ]   [ i1 ] = PKJ_RD_G_KC * area_val;
      pkj -> cond [ g_K2_pkj ]   [ i1 ] = PKJ_RD_G_K2 * area_val;
    } 
    else
    { // spiny
      double other_spines = area_val + len [ i1 ] * PKJ_SPINE_AREA * PKJ_SPINE_DENS;
      pkj -> cond [ g_leak_pkj ] [ i1 ] = 1.0 / PKJ_RM_D * other_spines;
      pkj -> cond [ g_NaF_pkj ]  [ i1 ] = PKJ_RD_G_NaF * area_val;
      pkj -> cond [ g_NaP_pkj ]  [ i1 ] = PKJ_RD_G_NaP * area_val;
      pkj -> cond [ g_CaP_pkj ]  [ i1 ] = PKJ_RD_G_CaP * area_val;
      pkj -> cond [ g_CaT_pkj ]  [ i1 ] = PKJ_RD_G_CaT * area_val;
      pkj -> cond [ g_Kh_pkj ]   [ i1 ] = PKJ_RD_G_KH * area_val;
      pkj -> cond [ g_Kdr_pkj ]  [ i1 ] = PKJ_RD_G_KDR * area_val;
      pkj -> cond [ g_KM_pkj ]   [ i1 ] = PKJ_RD_G_KM * area_val;
      pkj -> cond [ g_KA_pkj ]   [ i1 ] = PKJ_RD_G_KA * area_val;
      pkj -> cond [ g_KC_pkj ]   [ i1 ] = PKJ_RD_G_KC * area_val;
      pkj -> cond [ g_K2_pkj ]   [ i1 ] = PKJ_RD_G_K2 * area_val;
    }      
  }
  fclose ( file );

  for ( int i = 1; i < pkj -> n; i++ ) 
  {
    for ( int j = 0; j < PKJ_COMP; j++ ) 
    {
      pkj -> elem [ connect ] [ j + PKJ_COMP * i ] = ( pkj -> elem [ connect] [ j + PKJ_COMP * i ] > 0 ?
						    pkj -> elem [ connect ] [ j ] + PKJ_COMP * i : -1 );
      pkj -> elem [ compart ] [ j + PKJ_COMP * i ] = pkj -> elem [ compart ] [ j ];
      pkj -> elem [ Cm      ] [ j + PKJ_COMP * i ] = pkj -> elem [ Cm      ] [ j ];
      pkj -> elem [ area    ] [ j + PKJ_COMP * i ] = pkj -> elem [ area    ] [ j ];
      pkj -> elem [ i_ext   ] [ j + PKJ_COMP * i ] = pkj -> elem [ i_ext   ] [ j ];
      pkj -> shell   [ j + PKJ_COMP * i ] = pkj -> shell   [ j ];
      pkj -> rev_ca2 [ j + PKJ_COMP * i ] = pkj -> rev_ca2 [ j ];

      pkj -> cond [ g_leak_pkj ] [ j + PKJ_COMP * i ] = pkj -> cond [ g_leak_pkj ] [ j ];
      pkj -> cond [ g_NaF_pkj ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_NaF_pkj ] [ j ];
      pkj -> cond [ g_NaP_pkj ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_NaP_pkj ] [ j ];
      pkj -> cond [ g_CaP_pkj  ] [ j + PKJ_COMP * i ] = pkj -> cond [ g_CaP_pkj  ] [ j ];
      pkj -> cond [ g_CaT_pkj  ] [ j + PKJ_COMP * i ] = pkj -> cond [ g_CaT_pkj  ] [ j ];
      pkj -> cond [ g_Kh_pkj  ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_Kh_pkj  ] [ j ];
      pkj -> cond [ g_Kdr_pkj ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_Kdr_pkj ] [ j ];
      pkj -> cond [ g_KM_pkj  ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_KM_pkj  ] [ j ];
      pkj -> cond [ g_KA_pkj ] [ j + PKJ_COMP * i ]   = pkj -> cond [ g_KA_pkj ] [ j ];
      pkj -> cond [ g_KC_pkj  ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_KC_pkj  ] [ j ];
      pkj -> cond [ g_K2_pkj  ] [ j + PKJ_COMP * i ]  = pkj -> cond [ g_K2_pkj  ] [ j ];
    }
  }
}

__host__
static void finalize_host_pkj ( neuron_t *pkj )
{  
  if ( ( pkj -> n ) > 0 )
  {
    for ( int i = 0; i < n_elem; i++ ) { free ( pkj -> elem [ i ] ); }
    free ( pkj -> elem );

    for ( int i = 0; i < pkj_n_cond; i++ ) { free ( pkj -> cond [ i ] ); }
    free ( pkj -> cond );

    for ( int i = 0; i < pkj_n_ion; i++ ) { free ( pkj -> ion [ i ] ); }
    free ( pkj -> ion );

    free ( pkj -> shell );
    //free ( pkj -> ca2 );
    free ( pkj -> rev_ca2 );
  }
  free ( pkj );
}

__global__
static void device_mem_allocation ( const double DT, const int nx, const int ny, 
                                    neuron_t* d_pkj, double **d_elem, double **d_cond, double **d_ion,
                                    double *d_shell, double *d_rev ) 
{
  d_pkj -> elem    = d_elem;
  d_pkj -> cond    = d_cond;
  d_pkj -> ion     = d_ion;
  d_pkj -> shell     = d_shell;
  d_pkj -> rev_ca2 = d_rev;
  d_pkj -> nx = nx; // PKJ_X 
  d_pkj -> ny = ny; // PKJ_Y
  int n = nx * ny; int nc = n * PKJ_COMP;
  d_pkj -> n = n;   // # of neurons
  d_pkj -> nc = nc; // # of all compartments
  d_pkj -> DT = DT;
  //Debug
  printf ( "pkj -> n = %d, nc = %d\n", d_pkj -> n, d_pkj -> nc );
}

__global__
static void device_mem_allocation2 (const int n, double ** dev, double *ptr )
{
  dev [ n ] = ptr;
}

__host__ static void memcpy_neuron ( neuron_t *p_pkj, neuron_t *h_pkj )
{
  int nc = h_pkj -> nc;
  cudaMemcpy ( p_pkj -> elem [ v     ], h_pkj -> elem [ v     ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> elem [ Ca    ], h_pkj -> elem [ Ca    ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> elem [ Cm    ], h_pkj -> elem [ Cm    ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> elem [ area  ], h_pkj -> elem [ area  ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> elem [ i_ext ], h_pkj -> elem [ i_ext ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> elem [ connect ], h_pkj -> elem [ connect ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> elem [ compart ], h_pkj -> elem [ compart ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  
  for ( int i = 0; i < pkj_n_ion; i++ )
    cudaMemcpy ( p_pkj -> ion [ i ], h_pkj -> ion [ i ], nc * sizeof ( double ), cudaMemcpyHostToDevice );
  for ( int i = 0; i < pkj_n_cond; i++ )
    cudaMemcpy ( p_pkj -> cond [ i ], h_pkj -> cond [ i ], nc * sizeof ( double ), cudaMemcpyHostToDevice );

  cudaMemcpy ( p_pkj -> rev_ca2, h_pkj -> rev_ca2, nc * sizeof ( double ), cudaMemcpyHostToDevice );
  cudaMemcpy ( p_pkj -> shell,   h_pkj -> shell,   nc * sizeof ( double ), cudaMemcpyHostToDevice );

  cudaDeviceSynchronize ();
} 

__host__
neuron_t *pkj_initialize ( const int nx, const int ny, const char *type, neuron_t *p_pkj )
{
  p_pkj -> nx = nx; // PKJ_X 
  p_pkj -> ny = ny; // PKJ_Y
  int n = nx * ny; int nc = n * PKJ_COMP;
  p_pkj -> n = n; // # of neurons
  p_pkj -> nc = nc; // # of all compartments
  p_pkj -> neuron_type = N_TYPE_PKJ;  
  // set DT
  if      ( 0 == strncmp ( type, "BE", 2 ) ) { p_pkj -> DT = BE_DT;  p_pkj -> n_solver = kBE; }
  else if ( 0 == strncmp ( type, "CN", 2 ) ) { p_pkj -> DT = CN_DT;  p_pkj -> n_solver = kCN; }
  else if ( 0 == strncmp ( type, "RKC", 3 ) ) { p_pkj -> DT = 0.125; p_pkj -> n_solver = kRKC; }
  else { printf ("error in pkj_initialize\n"); exit ( 1 ); }

  neuron_t *d_pkj;  
  cudaMalloc ( ( neuron_t ** ) &d_pkj, sizeof ( neuron_t ) );

  if ( n == 0 ) { printf ( "pkj -> n = 0\n" ); return d_pkj; }
  
  /**/
  double **d_elem;
  cudaMalloc ( ( double *** ) &d_elem, n_elem * sizeof ( double * ) );
  double **d_cond;
  cudaMalloc ( ( double *** ) &d_cond, pkj_n_cond * sizeof ( double * ) );
  double **d_ion;
  cudaMalloc ( ( double *** ) &d_ion , pkj_n_ion  * sizeof ( double * ) );
  double *d_shell;
  cudaMalloc ( ( double **  ) &d_shell , nc * sizeof ( double ) );
  double *d_rev;
  cudaMalloc ( ( double **  ) &d_rev , nc * sizeof ( double ) );  

  p_pkj -> shell  = d_shell;
  p_pkj -> rev_ca2 = d_rev;

  p_pkj -> elem = ( double ** ) malloc ( n_elem * sizeof ( double * ) );
  for ( int i = 0; i < n_elem; i++ ) {
    cudaMalloc ( ( double ** ) ( & ( p_pkj -> elem [ i ] ) ), nc * sizeof ( double ) );
    device_mem_allocation2 <<< 1, 1 >>> ( i, d_elem, p_pkj -> elem [ i ] );
  }
  p_pkj -> cond = ( double ** ) malloc ( pkj_n_cond * sizeof ( double * ) );
  for ( int i = 0; i < pkj_n_cond; i++ ) {
    cudaMalloc ( ( double ** ) ( & ( p_pkj -> cond [ i ] ) ), nc * sizeof ( double ) );
    device_mem_allocation2 <<< 1, 1 >>> ( i, d_cond, p_pkj -> cond [ i ] );
  }
  p_pkj -> ion = ( double ** ) malloc ( pkj_n_ion * sizeof ( double * ) );
  for ( int i = 0; i < pkj_n_ion; i++ ) {
    cudaMalloc ( ( double ** ) ( & ( p_pkj -> ion [ i ] ) ), nc * sizeof ( double ) );
    device_mem_allocation2 <<< 1, 1 >>> ( i, d_ion, p_pkj -> ion [ i ] );
  }

  device_mem_allocation <<< 1, 1 >>> ( p_pkj -> DT, nx, ny, d_pkj, d_elem, d_cond, d_ion, d_shell, d_rev );
  cudaDeviceSynchronize();
  
  // set temporary host pkj
  neuron_t *h_pkj = ( neuron_t * ) malloc ( sizeof ( neuron_t ) );
  h_pkj -> nx = nx; h_pkj -> ny = ny;
  h_pkj -> n = n; h_pkj -> nc = nc;
  initialize_host_pkj ( h_pkj );
  pkj_initialize_ion  ( h_pkj );
  
  // copy host pkj -> device pkj
  memcpy_neuron ( p_pkj, h_pkj );
  finalize_host_pkj ( h_pkj );
 
  return d_pkj;
}

__host__
void pkj_finalize ( neuron_t *p_pkj, neuron_t *d_pkj, FILE *f_out, FILE *f_out_raster )
{
  if ( p_pkj -> n > 0 )
  {
    for ( int i = 0; i < n_elem; i++ ) { cudaFree ( p_pkj -> elem [ i ] ); }
    for ( int i = 0; i < pkj_n_cond; i++ ) { cudaFree ( p_pkj -> cond [ i ] ); }
    for ( int i = 0; i < pkj_n_ion;  i++ ) { cudaFree ( p_pkj -> ion [ i ] ); }
    free ( p_pkj -> elem  );
    free ( p_pkj -> cond  );
    free ( p_pkj -> ion  );
    cudaFree ( p_pkj -> shell );
    cudaFree ( p_pkj -> rev_ca2 );

    cudaMemcpy ( p_pkj, d_pkj, sizeof ( neuron_t ), cudaMemcpyDeviceToHost );
    cudaFree ( p_pkj -> elem );
    cudaFree ( p_pkj -> cond );
    cudaFree ( p_pkj -> ion  );
  }
  fclose ( f_out );
  fclose ( f_out_raster );
  free ( p_pkj );
  cudaFree ( d_pkj );
}

__global__ void pkj_set_current ( neuron_t *d_pkj, const double t, const int pulse_num, const int pat_num, int *pos, double *time, double *size )
{
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if ( id < d_pkj -> n && id < pat_num){
    double *i_ext_val;
    for ( int i = 0; i<pulse_num; i++ ){
      i_ext_val = & d_pkj -> elem [ i_ext ] [ pos[ id * pulse_num + i ] + id * PKJ_COMP ];
      *i_ext_val = 0;
      if( time[id * pulse_num + i] <= t && t < time[ id * pulse_num + i] + 1){
	*i_ext_val = size[ id * pulse_num + i ] *  ( 1.0 - exp ( - ( t - time[id * pulse_num + i] ) / 50 )) ;
      }
     }
  }
   
  //if ( id < d_pkj -> n )
  //{
  /*
    if(id == 0){
                        
      //d-p sitimulus
    double *i_ext_val = & d_pkj -> elem [ i_ext ] [ 902 + id * PKJ_COMP ];
    *i_ext_val = 0;
    if ( (3000 <= t && t < 3001 )  ) { *i_ext_val = 0.335* ( 1.0 - exp ( - ( t - 3000 ) / 50 )    ); }
    double *i_ext_val1 = & d_pkj -> elem [ i_ext ] [ 1007 + id * PKJ_COMP ];
    *i_ext_val1 = 0;
    if ( (3100 <= t && t < 3101 )  ) { *i_ext_val1 = 0.335 * ( 1.0 -  exp (  - ( t - 3100 ) / 50 )    ); }
    double *i_ext_val2 = & d_pkj -> elem [ i_ext ] [ 1061 + id * PKJ_COMP ];
    *i_ext_val2 = 0;
    if ( (3200 <= t && t < 3201 )  ) { *i_ext_val2 = 0.335 * ( 1.0 - exp ( - ( t - 3200 ) / 50 )    ); }
    double *i_ext_val3 = & d_pkj -> elem [ i_ext ] [ 1106 + id * PKJ_COMP ];
    *i_ext_val3 = 0;
    if ( (3300 <= t && t < 3301 )  ) { *i_ext_val3 = 0.335 * (1.0 - exp (  - ( t - 3300 ) / 50 )    ); }
    double *i_ext_val4 = & d_pkj -> elem [ i_ext ] [ 1210 + id * PKJ_COMP ];
    *i_ext_val4 = 0;
    if ( (3400 <= t && t < 3401 )  ) { *i_ext_val4 = 0.335 * ( 1.0 - exp (  - ( t - 3400 ) / 50 )    ); }
    double *i_ext_val5 = & d_pkj -> elem [ i_ext ] [ 1258 + id * PKJ_COMP ];
    *i_ext_val5 = 0;
    if ( (3500 <= t && t < 3501 )  ) { *i_ext_val5 = 0.335 * ( 1.0 - exp ( - ( t - 3500 ) / 50 )    ); }	
      
    }
   
   */
    /*              
    //p-d stimulus
    if(id <= 1){
    double *i_ext_val = & d_pkj -> elem [ i_ext ] [ 1258 + id * PKJ_COMP ];
    *i_ext_val = 0;
    if ( (3000 <= t && t < 3001 )  ) { *i_ext_val = 0.335 * ( 1.0 - exp ( - ( t - 3000 ) / 50 )    ); }
    double *i_ext_val1 = & d_pkj -> elem [ i_ext ] [ 1210 + id * PKJ_COMP ];
    *i_ext_val1 = 0;
    if ( (3100 <= t && t < 3101 )  ) { *i_ext_val1 = 0.335 * ( 1.0 -  exp (  - ( t - 3100 ) / 50 )    ); }
    double *i_ext_val2 = & d_pkj -> elem [ i_ext ] [ 1106 + id * PKJ_COMP ];
    *i_ext_val2 = 0;
    if ( (3200 <= t && t < 3201 )  ) { *i_ext_val2 = 0.335 * ( 1.0 - exp ( - ( t - 3200 ) / 50 )    ); }
    double *i_ext_val3 = & d_pkj -> elem [ i_ext ] [ 1061 + id * PKJ_COMP ];
    *i_ext_val3 = 0;
    if ( (3300 <= t && t < 3301 )  ) { *i_ext_val3 = 0.335 * (1.0 - exp (  - ( t - 3300 ) / 50 )    ); }
    double *i_ext_val4 = & d_pkj -> elem [ i_ext ] [ 1007 + id * PKJ_COMP ];
    *i_ext_val4 = 0;
    if ( (3400 <= t && t < 3401 )  ) { *i_ext_val4 = 0.335 * ( 1.0 - exp (  - ( t - 3400 ) / 50 )    ); }
    double *i_ext_val5 = & d_pkj -> elem [ i_ext ] [ 902 + id * PKJ_COMP ];
    *i_ext_val5 = 0;
    if ( (3500 <= t && t < 3501 )  ) { *i_ext_val5 = 0.335 * ( 1.0 - exp ( - ( t - 3500 ) / 50 )    ); }	
    }
    
    */
    //double curr = 0.215;
    //random sequence
    /*    
    if(id == 0){
    double *i_ext_val = & d_pkj -> elem [ i_ext ] [ 1305 + id * PKJ_COMP ];
    *i_ext_val = 0;
    if ( (3000 <= t && t < 3001 )  ) { *i_ext_val = curr * ( 1.0 - exp ( - ( t - 3000 ) / 50 )    ); }
    double *i_ext_val1 = & d_pkj -> elem [ i_ext ] [ 5 + id * PKJ_COMP ];
    *i_ext_val1 = 0;
    if ( (3040 <= t && t < 3041 )  ) { *i_ext_val1 = curr * ( 1.0 -  exp (  - ( t - 3040 ) / 50 )    ); }
    double *i_ext_val2 = & d_pkj -> elem [ i_ext ] [ 1057 + id * PKJ_COMP ];
    *i_ext_val2 = 0;
    if ( (3080 <= t && t < 3081 )  ) { *i_ext_val2 = curr * ( 1.0 - exp ( - ( t - 3080 ) / 50 )    ); }
    double *i_ext_val3 = & d_pkj -> elem [ i_ext ] [ 10 + id * PKJ_COMP ];
    *i_ext_val3 = 0;
    if ( (3120 <= t && t < 3121 )  ) { *i_ext_val3 = curr * (1.0 - exp (  - ( t - 3120 ) / 50 )    ); }
    double *i_ext_val4 = & d_pkj -> elem [ i_ext ] [ 255 + id * PKJ_COMP ];
    *i_ext_val4 = 0;
    if ( (3160 <= t && t < 3161 )  ) { *i_ext_val4 = curr * ( 1.0 - exp (  - ( t - 3160 ) / 50 )    ); }
    double *i_ext_val5 = & d_pkj -> elem [ i_ext ] [ 1092 + id * PKJ_COMP ];
    *i_ext_val5 = 0;
    if ( (3200 <= t && t < 3201 )  ) { *i_ext_val5 = curr * ( 1.0 - exp ( - ( t - 3200 ) / 50 )    ); }	
    double *i_ext_val6 = & d_pkj -> elem [ i_ext ] [ 679 + id * PKJ_COMP ];
    *i_ext_val6 = 0;
    if ( (3240 <= t && t < 3241 )  ) { *i_ext_val6 = curr * ( 1.0 - exp ( - ( t - 3240 ) / 50 )    ); }
    double *i_ext_val7 = & d_pkj -> elem [ i_ext ] [ 1356 + id * PKJ_COMP ];
    *i_ext_val7 = 0;
    if ( (3280 <= t && t < 3281 )  ) { *i_ext_val7 = curr * ( 1.0 -  exp (  - ( t - 3280 ) / 50 )    ); }
    double *i_ext_val8 = & d_pkj -> elem [ i_ext ] [ 455 + id * PKJ_COMP ];
    *i_ext_val8 = 0;
    if ( (3320 <= t && t < 3321 )  ) { *i_ext_val8 = curr * ( 1.0 - exp ( - ( t - 3320 ) / 50 )    ); }
    double *i_ext_val9 = & d_pkj -> elem [ i_ext ] [379 + id * PKJ_COMP ];
    *i_ext_val9 = 0;
    if ( (3360 <= t && t < 3361 )  ) { *i_ext_val9 = curr * ( 1.0 - exp ( - ( t - 3360 ) / 50 )    ); }

    }
    */
    /*
    if(id == 0){
    double *i_ext_val = & d_pkj -> elem [ i_ext ] [ 379 + id * PKJ_COMP ];
    *i_ext_val = 0;
    if ( (3000 <= t && t < 3001 )  ) { *i_ext_val = curr * ( 1.0 - exp ( - ( t - 3000 ) / 50 )    ); }
    double *i_ext_val1 = & d_pkj -> elem [ i_ext ] [ 455 + id * PKJ_COMP ];
    *i_ext_val1 = 0;
    if ( (3040 <= t && t < 3041 )  ) { *i_ext_val1 = curr * ( 1.0 -  exp (  - ( t - 3040 ) / 50 )    ); }
    double *i_ext_val2 = & d_pkj -> elem [ i_ext ] [ 1356 + id * PKJ_COMP ];
    *i_ext_val2 = 0;
    if ( (3080 <= t && t < 3081 )  ) { *i_ext_val2 = curr * ( 1.0 - exp ( - ( t - 3080 ) / 50 )    ); }
    double *i_ext_val3 = & d_pkj -> elem [ i_ext ] [ 679 + id * PKJ_COMP ];
    *i_ext_val3 = 0;
    if ( (3120 <= t && t < 3121 )  ) { *i_ext_val3 = curr * (1.0 - exp (  - ( t - 3120 ) / 50 )    ); }
    double *i_ext_val4 = & d_pkj -> elem [ i_ext ] [ 1092 + id * PKJ_COMP ];
    *i_ext_val4 = 0;
    if ( (3160 <= t && t < 3161 )  ) { *i_ext_val4 = curr * ( 1.0 - exp (  - ( t - 3160 ) / 50 )    ); }
    double *i_ext_val5 = & d_pkj -> elem [ i_ext ] [ 255 + id * PKJ_COMP ];
    *i_ext_val5 = 0;
    if ( (3200 <= t && t < 3201 )  ) { *i_ext_val5 = curr * ( 1.0 - exp ( - ( t - 3200 ) / 50 )    ); }	
    double *i_ext_val6 = & d_pkj -> elem [ i_ext ] [ 10 + id * PKJ_COMP ];
    *i_ext_val6 = 0;
    if ( (3240 <= t && t < 3241 )  ) { *i_ext_val6 = curr * ( 1.0 - exp ( - ( t - 3240 ) / 50 )    ); }
    double *i_ext_val7 = & d_pkj -> elem [ i_ext ] [ 1057 + id * PKJ_COMP ];
    *i_ext_val7 = 0;
    if ( (3280 <= t && t < 3281 )  ) { *i_ext_val7 = curr * ( 1.0 -  exp (  - ( t - 3280 ) / 50 )    ); }
    double *i_ext_val8 = & d_pkj -> elem [ i_ext ] [ 5 + id * PKJ_COMP ];
    *i_ext_val8 = 0;
    if ( (3320 <= t && t < 3321 )  ) { *i_ext_val8 = curr * ( 1.0 - exp ( - ( t - 3320 ) / 50 )    ); }
    double *i_ext_val9 = & d_pkj -> elem [ i_ext ] [ 1305 + id * PKJ_COMP ];
    *i_ext_val9 = 0;
    if ( (3360 <= t && t < 3361 )  ) { *i_ext_val9 = curr * ( 1.0 - exp ( - ( t - 3360 ) / 50 )    ); }
     
    }
    */
    
    /*    
    //fig3a~e
   
    if ( ( 550 <= t && t < 1300 ) && ( id == 0 ) ) { *i_ext_val = 0.1e-3; }
    if ( ( 550 <= t && t < 1300 ) && ( id == 1 ) ) { *i_ext_val = 0.5e-3; } 
    if ( ( 550 <= t && t < 1300 ) && ( id == 2 ) ) { *i_ext_val = 1.0e-3; } 
    if ( ( 550 <= t && t < 1300 ) && ( id == 3 ) ) { *i_ext_val = 2.0e-3; } 
    if ( ( 550 <= t && t < 1300 ) && ( id == 4 ) ) { *i_ext_val = 3.0e-3; } 
    */

    //if ( ( 500 <= t && t < 800 )  ) { *i_ext_val = 3.0e-3; }
    //if ( ( 800 <= t && t < 900 )  ) { *i_ext_val = -1.0e-3; }

    // fig3f
    //if ( ( 550 <= t && t < 1300 )  ) { *i_ext_val = PKJ_CURRENT; }
    //    if ( (200 <= t && t < 201 )  ) { *i_ext_val = 3.0e-3; }
    //
    //if ( ( 550 <= t && t < 1300 )  ) { *i_ext_val = 3.0e-3 * ( 1.0 - exp ( - ( t - 550 ) / 50 )    ); }
    //}
}

__global__ static void set_current_poisson ( neuron_t *d_pkj, const double val_current )
{  
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  if ( id < d_pkj -> n )
  {     
    double *i_ext_val = & d_pkj -> elem [ i_ext ] [ 1599 + id * PKJ_COMP ];
    *i_ext_val = 0.1e-3 * val_current;
  }
}

__host__ void pkj_set_current_Poisson ( const int num, neuron_t *d_pkj, const double t )
{
  // Poislon
  int firing_flag;
  double fr;
  static double rate = 25.0;
  static double val_current = 0.0;
  ( 350 <= t && t < 1000 )? fr = rate : fr = 0.0;
  double f = ( ( double ) rand ( ) + 1.0 ) / ( ( double ) RAND_MAX + 2.0 );
  ( fr * 0.125 * 0.001 > f )? firing_flag = 1 : firing_flag = 0;
  val_current = val_current * 0.81193634615 + firing_flag;  // 0.81193634615 = exp ( -0.125 / tau ), tau = 0.6ms
  set_current_poisson <<< ( num + 127 ) / 128, 128 >>> ( d_pkj, val_current ); 
} 

__host__
void pkj_output_file ( neuron_t *pkj, double *memv_test, double *memv_raster, const double t, FILE *f_out, FILE *f_out_raster )
{
  
  FILE *f = f_out;
  FILE *f_raster = f_out_raster;  

  //specific comps
  ///*
  fprintf ( f, "%lf,", t );    
  for ( int j = 0; j < pkj -> n; j++ ) 
  {
    fprintf ( f, "%lf,", memv_test [ 1599 + j * PKJ_COMP ] );
    //fprintf ( f, "%lf,", memv_test [ 1599 + j * PKJ_COMP ] );
    //fprintf ( f, "%lf,", memv_test [ 1599 + j * PKJ_COMP ] );
    //fprintf ( f, "%lf,", memv_test [ 1599 + j * PKJ_COMP ] );
  }
  fprintf ( f, "\n" );
  //*/
  // all comps
  /*  
  for ( int j = 0; j < pkj -> n; j++ ) 
]  {
    fprintf ( f, "%lf,", t );    
    for ( int i = 0; i < PKJ_COMP; i++ ) { 
      // V
      fprintf ( f, "%lf,", memv_test [ i + j * PKJ_COMP ] );
      // Ca
      //fprintf ( f, "%lf,", 1000.0 * memv_test [ i + j * PKJ_COMP ] );
    }
    fprintf ( f, "," );
  }
  fprintf ( f, "\n" );
  */ 
 
  //fprintf ( f_raster, "%lf,", t );
  for ( int j = 0; j < pkj -> n; j++ ) {
    if ( memv_test [ 1599 + j * PKJ_COMP ] > -10.0 && memv_raster [ j ] <= -10.0 )
      fprintf ( f_raster, "%lf,%d\n", t, j );
    //else 
      //fprintf ( f_raster, "," );
    memv_raster [ j ] = memv_test [ 1599 + j * PKJ_COMP ];
  }
  //fprintf ( f_raster, "\n" );

}
