#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<stdbool.h>
//#define CN 2000     //cell num
//#define PN 6 // pulse num
#define TM 5000 // simulation time
#define CF_IN 2.0
#define DT 1.
#define START 3000
#define TAU_PRE 600.
#define TAU_POST 100.
#define MIN_PULSE 2
#define MAX_PULSE 20
#define A1 0.9
#define A2 0.1
#define A3 0.0001


void updateTrace( double *trace_pre, double *trace_post, bool *spike_pre, bool *spike_post, int pulse_num, int seq_num){
  int i;
  for( i = 0; i < pulse_num * seq_num; i++){
    trace_pre[i] = exp ( -DT / TAU_PRE ) * trace_pre[i] + spike_pre[i];
    trace_post[i] = exp ( -DT / TAU_POST ) * trace_post[i] + spike_post[i];
  }  
}


bool check_spike( int cell, int num, double t, double *time_pre, double *time_post, int pulse_num, int seq_num ){
  int i;
  if( cell == 0 ){ 
    if( t != 0 && t == time_pre[num] ){
      return true;
    }
  }else{ 
    if( t == time_post[num] ){
      return true;
    }
  }
  return false;
}

int main(){


  FILE *pos_file;
  FILE *size_file;
  FILE *time_file;
  FILE *output_file;

  int pulse_num;
  int seq_num;
  int seq_max;
  int count;
  int i,j;
  int n,dummy;

  double t;
  double dw;

  char pos_txt[] = "input_pos_";
  char size_txt[] = "input_size_";
  char time_txt[]  = "input_time_";
  char pos_path[20];
  char size_path[20];
  char time_path[20];
  char output_path[20];
  char num_txt[3];

  double *syn;
  double *trace_pre, *trace_post, *time_pre, *time_post; 
  bool *spike_pre, *spike_post;

  for( pulse_num = MIN_PULSE-1; pulse_num <= MAX_PULSE-1; pulse_num++){

    /* initialize */
    t = START;
    count = 0;
    
    /*make path name*/
    sprintf( num_txt, "%d", pulse_num );
    sprintf( pos_path, "%s%s.txt", pos_txt, num_txt );
    sprintf( size_path, "%s%s.txt", size_txt, num_txt );
    sprintf( time_path, "%s%s.txt", time_txt, num_txt );

    printf ("pulse_num = %d, processing %s\n",pulse_num, size_path); 

    /*file open*/
    pos_file = fopen ( pos_path, "r" );
    size_file = fopen ( size_path, "r" );
    time_file = fopen ( time_path, "r");

    if( pos_file == NULL ){
      printf("pos file open error");
      exit(1);
    }
    if( size_file == NULL ){
      printf("size file open error");
      exit(1);
    }
    if( time_file == NULL ){
      printf("time file open error");
      exit(1);
    }

    fscanf( pos_file, "%d %d", &dummy, &seq_num );
    seq_max = seq_num;
    printf("dummy = %d, seq_num = %d\n", dummy, seq_num);
    fclose(pos_file);


    syn = (double *)malloc( dummy * seq_num * sizeof(double));

    while( fscanf(size_file, "%lf ", &syn[count]) != EOF){
      //printf("%lf,", syn[count]);
      count++;
    }
    fclose(size_file);

    
    
    trace_pre = ( double *) malloc ( pulse_num * seq_num * sizeof( double ));
    trace_post = ( double *) malloc ( pulse_num * seq_num *sizeof( double ));
    spike_pre = ( bool *) malloc ( pulse_num * seq_num *sizeof( bool ));
    spike_post = ( bool *) malloc ( pulse_num * seq_num *sizeof( bool ));
    time_pre = ( double *) malloc (  pulse_num *  seq_num * sizeof( double ));
    time_post = ( double *) malloc (  pulse_num *  seq_num * sizeof( double ));
    
    count = 0;
    while( fscanf(time_file, "%lf", &time_pre[count]) != EOF){
      count++;
    }
    fclose( time_file );

    for( i = 0; i < seq_num; i++){
      for ( j = 0; j < pulse_num; j++){
	time_post[pulse_num * i + j] = time_pre[pulse_num * i + pulse_num - 1] + CF_IN;
      }
    }
      
 
    for( i = 0; i < pulse_num * seq_num; i++ ){
      trace_pre[i] = 0;
      trace_post[i] = 0;
      spike_pre[i] = false;
      spike_post[i] = false;
    }
    

    /*calculate synaptic weights*/

    

    dw =0;
    while( t < TM ){      
      for( i = 0; i < pulse_num * seq_num; i++){
	spike_pre[i] = check_spike( 0, i, t, time_pre, time_post, pulse_num, seq_num);
	spike_post[i] = check_spike( 1, i, t, time_pre, time_post, pulse_num, seq_num);
	//if(spike_post[i] == true || spike_pre[i] == true) printf("%lf, %d, %d\n", t, spike_pre[i], spike_post[i]);
      }
      updateTrace( trace_pre, trace_post, spike_pre, spike_post,  pulse_num, seq_num);
      for( i = 0; i < pulse_num * seq_num; i++){
	dw =  -A1 * trace_pre[i] * spike_post[i] + A2 * trace_post[i] * spike_pre[i] + A3;
	//if(dw != 0. && t<3700)printf("%d, %lf, %lf,%d, %d\n ", i, t, dw, spike_pre[i], spike_post[i]);
	syn[i] = syn[i] + dw;
	//syn[PN*2*n+i] = ( syn[PN*2*n+i] <0) ? 0. : syn[PN*2*n+i];
      }      
      t += DT;
    }
    
    for(j=0;j<pulse_num * seq_num;j++){
      syn[j] = ( syn[j] < 0. ) ? 0. : syn[j];
      syn[j] = ( syn[j] > 0.65 ) ? 0.65 : syn[j];
      //printf("%lf,", syn[PN*2*n+j]);
    }
    

    /*output file*/
    

    sprintf( output_path, "%s%s_l.txt", size_txt, num_txt );
    output_file = fopen( output_path , "w" );
    for( i = 0; i < pulse_num * seq_num ;i++){
      fprintf(output_file, "%lf ", syn[i]);
    }
    fclose(output_file); 
    
  }

}
