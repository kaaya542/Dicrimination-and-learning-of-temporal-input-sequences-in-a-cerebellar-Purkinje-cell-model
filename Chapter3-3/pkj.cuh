#ifndef __PKJ_CUH__
#define __PKJ_CUH__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "param.h"
#include "pkj_ion.cuh"
//#include <cuda_runtime.h>

#define PKJ_COMP ( 1600 ) // if change PKJ_COMP -> change pkj_ion soma num
#define PARAM_FILE_PKJ "./Pkj.txt"
//#define PARAM_FILE_PKJ "./Pkj_test.txt"
#define GJ_EACH_PKJ 26
#define PKJ_CURRENT ( 2.0e-4 ) // ]nA]

typedef enum { PKJ_soma, PKJ_main, PKJ_smooth, PKJ_spiny, pkj_n_comptype } pkj_comptype_t;

#define V_INIT_PKJ   ( -68.0 )
#define V_LEAK_PKJ   ( -80.0 )//Tsuyuki -55 // De schtter -80
#define V_Na_PKJ     ( 45.0 )
#define V_Ca_PKJ     ( 135.0 )
#define V_KH_PKJ     ( -30.0 )
#define V_K_PKJ      ( -85.0 )
#define V_EX_PKJ     ( 0.0)
#define G_GJ_PKJ    ( 0.0 ) //( 2.5e-6 ) // Gap Junctional Conductance [mS]

// g: Conductance [mS/cm^2] // Noter : De schtter [ S / m^2 ]
// Soma
#define PKJ_S_G_NaF   (7500.0)
#define PKJ_S_G_NaP   (1.0)
#define PKJ_S_G_CaP   (0.0)
#define PKJ_S_G_CaT   (0.5)//0.5
#define PKJ_S_G_KH    (0.3)
#define PKJ_S_G_KDR   900.0//600.0//900.0//(600.0)
#define PKJ_S_G_KM    0.14//0.040//0.14//(0.040)
#define PKJ_S_G_KA    (15.0)
#define PKJ_S_G_KC    (0.0)
#define PKJ_S_G_K2    (0.0)

// Main Dendrite
#define PKJ_MD_G_NaF  (0.0)
#define PKJ_MD_G_NaP  (0.0)
#define PKJ_MD_G_CaP  (4.0)//4.0//4.5
#define PKJ_MD_G_CaT  (0.5)//0.5
#define PKJ_MD_G_KH   (0.0)
#define PKJ_MD_G_KDR  90.0//60.0//90.0//(60.0)
#define PKJ_MD_G_KM   0.04//0.01//0.04//(0.010)
#define PKJ_MD_G_KA   (2.0)
#define PKJ_MD_G_KC   (80.0)//80.0
#define PKJ_MD_G_K2   (0.39)

// Rest of Dendrite
#define PKJ_RD_G_NaF  (0.0)
#define PKJ_RD_G_NaP  (0.0)
#define PKJ_RD_G_CaP  (4.5)//4.5
#define PKJ_RD_G_CaT  (0.5)//0.5
#define PKJ_RD_G_KH   (0.0)
#define PKJ_RD_G_KDR  (0.0)
#define PKJ_RD_G_KM   0.013//0.10//(0.013)
#define PKJ_RD_G_KA   (0.0)
#define PKJ_RD_G_KC   (80.0)//80.0
#define PKJ_RD_G_K2   (0.39)
#define PKJ_G_EX      (0.7e-6) // [mS]

// spine
#define PKJ_SPINE_DENS (13.0e4) // [1/cm]
#define PKJ_SPINE_AREA (1.33e-8) // [cm^2]

// Membrane
#define PKJ_CM        (1.64) // [muF/cm^2]
#define PKJ_RM_S      (10.0) // Soma [kohmcm^2]
#define PKJ_RM_D      (30.0) // Dendrite [kohmcm^2]
#define PKJ_RI        (0.25) // [kohmcm]

// Calcium
#define B_Ca1_PKJ      (10.0  )//(1.0/2.0/96494.0) //(5.18e-6)//(0.0518)//52
//#define PKJ_TAU_Ca    (0.1   ) // [ms]
#define PKJ_SHELLD    (0.2e-4) // [cm]
#define Ca1OUT_PKJ    (2.4   ) // [mM]
#define Ca1_0_PKJ     (4.0e-5) // [muM]
#define F_PKJ         ( 9.6485e4 ) // Faraday constant [s*A/mol]


#define pkj_n_ion 18
typedef enum { m_NaF_pkj, h_NaF_pkj, m_NaP_pkj, m_CaP_pkj, h_CaP_pkj, m_CaT_pkj, h_CaT_pkj, m_KA_pkj, h_KA_pkj, 
              m_KC_pkj, m_Kh1_pkj, m_Kh2_pkj, m_Kdr_pkj, h_Kdr_pkj, m_KM_pkj, z_KC_pkj, m_K2_pkj, z_K2_pkj } pkj_ion_t;

#define pkj_n_cond 11
typedef enum { g_leak_pkj, g_NaF_pkj, g_NaP_pkj, g_CaP_pkj, g_CaT_pkj, 
            g_Kh_pkj, g_Kdr_pkj, g_KM_pkj, g_KA_pkj, g_KC_pkj, g_K2_pkj } pkj_cond_t;

__host__ neuron_t *pkj_initialize ( const int, const int, const char *, neuron_t * );
__host__ void pkj_finalize ( neuron_t *, neuron_t *, FILE *, FILE * );
__global__ void pkj_set_current ( neuron_t *, const double, const int, const int, int *, double *, double * );
__host__ void pkj_set_current_Poisson ( const int, neuron_t *, const double );
__host__ void pkj_output_file  ( neuron_t *, double *, double *, const double, FILE *, FILE * );

#endif // __PKJ_CUH__
