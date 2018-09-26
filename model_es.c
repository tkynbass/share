//
//  model_es.c
//  
//
//  Created by tkym on 2018/09/26.
//

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "dSFMT/dSFMT.h"
//#include "MT.h"
#include <time.h>
#include <omp.h>

#define DIMENSION (3) //次元
#define NANO (4.0e-8)   // 長さの単位
#define KBT (1.38064852e-23 / NANO / NANO) //ボルツマン
#define TEMPARTURE (300)
#define M_A (1.82527596e+6)
#define N_A (6.022140857e+23)

#define NUMBER (6193)      //粒子数
#define PARTICLE_MASS (M_A / N_A / 1000)      //染色体粒子の質量
#define SPB_MASS ( 9.0 * PARTICLE_MASS)      //SPBの質量
#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PI (M_PI)
#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )

#define K_EXCLUDE (1.0)    //排除体積効果の強さ
#define K_BOND (1.0)    //ばね定数
//#define K_BOND_2 (1.0e-4)  //ひもの硬さ
#define DELTA (1.0e-11)  //刻み幅
#define PARTICLE_MYU (2.0 * DIMENSION * PI * PARTICLE_RADIUS * NANO * 0.000890 / 100 )    //粘性抵抗の強さ
#define MEMBRAIN_EXCLUDE (1.0)     //膜との衝突
#define MEMBRAIN_EXCLUDE_SPB (1.0) //SPBとの衝突

#define SPB_RADIUS ( 3.0 )      //SPBの半径
#define SPB_MYU (2.0 * DIMENSION * PI * SPB_RADIUS * NANO * 0.000890 / 100 )  //SPBの粘性

//k_bond2 k_expression

typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

typedef struct particle {           //構造体の型宣言
    CHAIN chr_no;
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

void read_coordinate_init ( const char *str,  int start ){       //初期値設定
    
    int i, i_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    sprintf (filename, "fission_result_%d.txt", start);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++){
        
        fscanf(fpr, "%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", &i_dummy, &part[i].chr_no, &part[i].particle_type,
               &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
               &part[i].velocity[X], &part[i].velocity[Y], &part[i].velocity[Z],
               &membrain_radius);
        //fgets(dummy, 128, fpr);
        
        //printf ("%d %lf %lf %lf \n", i, part[i].position[X], part[i].position[Y], part[i].position[Z]);
        
    }
    
    fscanf (fpr, "%s %s %lf %lf %lf %lf %lf %lf", dummy, dummy, &spb.position[X], &spb.position[Y], &spb.position[Z],
            &spb.velocity[X], &spb.velocity[Y], &spb.velocity[Z]);
    fgets(dummy, 128, fpr);
    
    //fscanf (fpr, "%s %s %lf %lf %lf", dummy, dummy, &Nucleolus_circle_center[X], &Nucleolus_circle_center[Y], &Nucleolus_circle_center[Z]);
    
    fclose (fpr);
    
}

