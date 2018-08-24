//
//  dis_hist.c
//  
//
//  Created by 高山雄揮 on 2018/01/20.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <time.h>
#include <string.h>
#include "MT.h"

#define DIMENSION (3)
#define NUMBER (6193)      //粒子数
#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PI (M_PI)
#define MEMBRAIN_RADIUS_MINIMUM  (32.0) //膜の半径の最小値
#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )
#define DIV_NUMBER (100)
#define DIV_DELTA (23.7)
#define RANGE_MIN (5880) //グラフに表示させる最小値

#define TOP_NUMBER (100)
#define BOTTOM_NUMBER (100)
#define RAN_NUMBER (100)

typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

enum label{ X, Y, Z};

typedef struct particle {           //構造体の型宣言
    CHAIN chr_no;
    TYPE particle_type;
    double position[DIMENSION];
    //int ex_state;
    
} Particle;

Particle *ptr;

Particle part[NUMBER];

void read_coordinate( int time ){       //初期値設定
    
    int i, i_dummy;
    
    double d_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    sprintf (filename, "fission_result_%d.txt", time);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++){
        
        fscanf(fpr, "%d %d %d %lf %lf %lf %lf %lf %lf %lf \n", &i_dummy, &part[i].chr_no, &part[i].particle_type,
               &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
               &d_dummy, &d_dummy, &d_dummy, &d_dummy);
        //fgets(dummy, 128, fpr);
        
    }
    
    fclose (fpr);
    
}

void read_expression_data(unsigned int top_list[TOP_NUMBER], unsigned int bottom_list[BOTTOM_NUMBER]) {
    
    FILE *fpr;
    
    int i, i_dummy, state, top_count = 0, bottom_count = 0;
    
    char filename[128], dummy[128];
    
    sprintf (filename, "gene_best100.txt");

    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("\n error \n");
        
        exit (1);
    }
    
    for (i=0; i<TOP_NUMBER; i++) {
        
        fscanf (fpr,"%d %d\n", &i_dummy, &top_list[i]);
    }
    fclose (fpr);
    
    
    sprintf (filename, "gene_worst100.txt");
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("\n error \n");
        
        exit (1);
    }
    
    for (i=0; i<TOP_NUMBER; i++) {
        
        fscanf (fpr,"%d %d\n", &i_dummy, &bottom_list[i]);
    }
    fclose (fpr);
    
    
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

void write_data (const unsigned int top_ppv_hist[DIV_NUMBER], const unsigned int bottom_ppv_hist[DIV_NUMBER],
                 const unsigned int ran_ppv_hist[DIV_NUMBER], const char *output_file) {
    
    unsigned int i, j;
    
    FILE *fpw;

    char result[128], str[128];
    
    sprintf (result, "%s", output_file);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    for (i=0; i<DIV_NUMBER; i++) {
        
        fprintf (fpw, "%d %d %d %d\n", i*DIV_DELTA + RANGE_MIN, top_ppv_hist[i], bottom_ppv_hist[i],
                 ran_ppv_hist[i]);
    }
    
}

int main ( int argc, char **argv) {
    
    unsigned int i, j, t, counter;
    
    int start_number = atoi (argv[1]);
    int calculate_number = atoi (argv[2]);
    int width = atoi (argv[3]);
    
    Particle *top[TOP_NUMBER], *bottom[BOTTOM_NUMBER], *ran[RAN_NUMBER];
    unsigned int top_list[TOP_NUMBER], bottom_list[BOTTOM_NUMBER], ran_list[RAN_NUMBER], num_list[NUMBER];
    unsigned int top_ppv_hist[DIV_NUMBER], bottom_ppv_hist[DIV_NUMBER], ran_ppv_hist[DIV_NUMBER];
    double top_ppv, bottom_ppv, ran_ppv;
    
    double top_dist, bottom_dist, ran_dist;
    
    unsigned int no_counter;
    
    //printf ("%d %d\n", start_number, calculate_number);
    
    ptr = (Particle *)malloc(NUMBER * sizeof(Particle));
    
    if (ptr == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    init_genrand((unsigned)time(NULL));
    
    //ランダムにRAN_NUMBER粒子取り出す
    for (i=0; i<NUMBER; i++){
        
        num_list[i] = i;
    }
    init_genrand (10);
    ran_list[0] = genrand_int32() % NUMBER;
    for (i=1; i<RAN_NUMBER; i++){
        
        for (j=ran_list[i-1]; j<NUMBER-i; j++){
            
            num_list[j] += 1;
        }
        ran_list[i] = genrand_int32() % (NUMBER - i);
        
        //printf ("ran_%d = %d\n", i, ran_list[i]);
    }
    
    read_expression_data (top_list, bottom_list);
    
    for (i=0; i<DIV_NUMBER; i++){
        
        top_ppv_hist[i] = 0;
        bottom_ppv_hist[i] = 0;
        ran_ppv_hist[i] = 0;
    }
    
    for ( t=0; t < calculate_number; t+= width) {
        
        //printf ("t = %d\r", t);
        
        //座標の読み込み
        read_coordinate ( t + start_number);
        
        //上位下位ランダムの粒子の座標読み込み
        for (i=0; i<TOP_NUMBER; i++) {
            
            top[i] = &part[top_list[i]];
        }
        for (i=0; i<BOTTOM_NUMBER; i++) {
            
            bottom[i] = &part[bottom_list[i]];
        }
        for (i=0; i<RAN_NUMBER; i++) {
            
            ran[i] = &part[ran_list[i]];
        }
        
        top_ppv = 0.0;
        bottom_ppv = 0.0;
        ran_ppv = 0.0;
        
        no_counter = 0;
        
        for (i=0; i<TOP_NUMBER; i++) {
            
            for (j=i+1; j<TOP_NUMBER; j++) {
                
                top_dist = Euclid_norm (top[i]->position, top[j]->position) * 0.04 * 0.71/ 1.28;
                
                if (top_dist > 0.1713)  top_ppv += -log((top_dist - 0.1713)/1.0965) / 0.6865;
                else no_counter++;

            }
        }
        
        //top_ppv /= TOP_NUMBER*(TOP_NUMBER-1)/2.0 - no_counter;
        top_ppv_hist[(int)((top_ppv - RANGE_MIN)/ DIV_DELTA)] += 1;
        
        no_counter = 0;
        
        for (i=0; i<BOTTOM_NUMBER; i++) {
            
            for (j=i+1; j<BOTTOM_NUMBER; j++) {
                
                bottom_dist = Euclid_norm (bottom[i]->position, bottom[j]->position) * 0.04 * 0.71/ 1.28;
                
                if (bottom_dist > 0.1713) bottom_ppv += -log((bottom_dist - 0.1713)/1.0965) / 0.6865;
                else no_counter++;
            }
        }
        
        //bottom_ppv /= BOTTOM_NUMBER*(BOTTOM_NUMBER-1)/2.0 - no_counter;
        bottom_ppv_hist[(int)((bottom_ppv - RANGE_MIN)/ DIV_DELTA)] += 1;
        
        no_counter = 0;
        
        for (i=0; i<RAN_NUMBER; i++) {
            
            for (j=i+1; j<RAN_NUMBER; j++) {
                
                ran_dist = Euclid_norm (ran[i]->position, ran[j]->position) * 0.04 * 0.71/ 1.28;
                
                if (ran_dist > 0.1713) ran_ppv += -log((ran_dist - 0.1713)/1.0965) / 0.6865;
                else no_counter++;
            }
        }
        
        //ran_ppv /= RAN_NUMBER * (RAN_NUMBER-1)/2.0 - no_counter;
        ran_ppv_hist[(int)((ran_ppv - RANGE_MIN)/ DIV_DELTA)] += 1;
        
    }
    
    write_data (top_ppv_hist, bottom_ppv_hist, ran_ppv_hist, argv[4]);
    
    
    return (0);
}
