//
//  common.h
//  
//
//  Created by tkym on 2018/06/06.
//

#ifndef common_h
#define common_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "dSFMT/dSFMT.h"
#include <time.h>
#indef _OPENMP
#include <omp.h>
#endif

#define DIMENSION (3) //次元
#define NANO (4.0e-8)   // 長さの単位
#define KBT (1.38064852e-23 / NANO / NANO) //ボルツマン
#define TEMPARTURE (300)
#define M_A (1.82527596e+6)     //分子量
#define N_A (6.022140857e+23)   　//アボガドロ

#define NUMBER (6193)      //粒子数
#define DELTA (1.0e-11)  //刻み幅
#define PI (M_PI)

#define K_EXCLUDE (1.0)    //粒子間の排除体積効果
#define K_BOND (1.0)    //１つ隣ばね
#define K_BOND_2 (1.0e-4)  //２つ隣ばね
#define K_MEMBRAIN_EXCLUDE (1.0)     //膜との排除体積効果

#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PARTICLE_MASS ( M_A / N_A / 1000.0)     //染色体粒子の質量
#define PARTICLE_MYU (2.0 * DIMENSION * PI * PARTICLE_RADIUS * NANO * 0.000890 / 100 )    //粘性抵抗の強さ

#define SPB_RADIUS ( 3.0 )      //SPBの半径
#define SPB_MASS ( 9.0 * PARTICLE_MASS)      //SPBの質量
#define SPB_MYU (2.0 * DIMENSION * PI * SPB_RADIUS * NANO * 0.000890 / 100 )  //SPBの粘性

#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )

unsigned int start_no, save_max, loop;
double membrain_radius, nucleolus_myu, nucleolus_mass, k_expression;

char result_name[256];

typedef enum coordinate {
    
    X, Y, Z

}COORDINATE;

typedef enum type{
    
    Normal, Centromere, Telomere, rDNA

}TYPE;

typedef enum chromosome_type{
    
    A, B, C

}CHROMOSOME_TYPE;

typedef struct particle {
    
    CHROMOSOME_TYPE chr_no;
    TYPE particle_type;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_old[DIMENSION];
    double velocity[DIMENSION];
    double velocity_2[DIMENSION];
    double velocity_new[DIMENSION];
    int list_no;
    int *list;
    
} Particle;

Particle *ptr;

Particle part[NUMBER];

Particle spb;

Particle nucleolus;


double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

void Set_result ( int argc, char **argv ) {
    
    if (argc < 1 ) {
        
        printf ("\n\n   error : input result_name \n\n");
    }
    
    sprintf (result_name, "result/%s", argv[1]);
}

/// Input setting information and make setting file ////
void Load_setting ( int argc, char **argv ) {
    
    Set_result ( argc, argv);
    
    char filename[128], dummy[256];
    FILE *fpw;
    
    printf ("Input : Start, Loop \n");
    
    scanf ("%d %d %d ", &start_no, &save_max, &loop);
    
    if ( (fpw = fopen(filename, "w")) == NULL ) {
        
        printf ("\n\n   error : can not make setting file \n\n");
    }
    
    fprintf (fpw, "%d %d %d", start_no, save_max, loop);
    fprintf (fpw, "\n\n #0 : start number \n #1 : save times \n #2 : loop ");
    
    fclose (fpw);
    
    printf ("   MASS = %lf \n", PARTICLE_MASS);
    
}


#endif /* common_h */
