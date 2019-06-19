// 粒子をランダム配置　→　楕円体近似の核内に緩和

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <dSFMT/dSFMT.h>
#include <time.h>

#define DIMENSION ( 3 ) //次元
#define LENGTH ( 7.0e-8 )   //長さの単位 (粒子径)

#define M_A ( 1.85131596e+6 )   // g/mol/particle
#define N_A ( 6.022140857e+23 )
#define PASTIS_SCALING ( 1.8e-6 / 75 / LENGTH)

#define NUMBER_MAX ( 2516 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )

#define K_BOND ( 1.0e+0 )       //1つ隣　ばね定数
#define K_BOND_2 ( 1.0e+0 )     //2つ隣
#define K_BOND_3 ( 1.0e+0 )     //3つ隣
#define K_EXCLUDE ( 1.0e+1 )

#define DELTA ( 1.0e-3 )  //刻み幅
#define MITIGATION_INTERVAL (1.0e+3)
#define LIST_INTERVAL ( 50 )   // リスト化の間隔
#define LIST_RADIUS ( 10.0 * PARTICLE_RADIUS)

#define MEMBRANE_EXCLUDE ( 1.0 )     //膜との衝突
#define MEMBRANE_EXCLUDE_SPB ( 1.0 ) //SPBとの衝突

#define BOND_DISTANCE ( 2.0 * PARTICLE_RADIUS * 0.8 )   // １個隣ばねの自然長

#define SPB_RADIUS (  3.0  )      //SPBの半径

// Ellipsoid axes parameter of nucleus & nucleolus //

#define MEMBRANE_AXIS_1 ( 1.889011e-6 / LENGTH )
#define MEMBRANE_AXIS_2 ( 0.85 * MEMBRANE_AXIS_1 )
#define MEMBRANE_AXIS_3 ( 0.75 * MEMBRANE_AXIS_1 )

#define NUCLEOLUS_AXIS_1 ( 1.1426593e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.9 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.8 * NUCLEOLUS_AXIS_1 )

const unsigned int CENT_LIST[] = { 754, 1440, 2244 };
const unsigned int TELO_LIST[] = { 0, 1115, 1116, 2023};
const unsigned int rDNA_LIST[] = { 2024, 2515};
const double ORIGIN[] = { 0.0, 0.0, 0.0};
const double SPB_POS[] = { 0.9e-6 / LENGTH, 0.85e-6 / LENGTH, 1.0e-6 / LENGTH };
const double NUCLEOLUS_POS[] = { -0.25e-6 / LENGTH, -0.365e-6 / LENGTH, 0.3e-6 / LENGTH};

//#define POTENTIAL_DELTA (1.0e-7)

typedef enum chain {
    A, B, C
}CHAIN;

typedef enum type {
    Normal, Centromere, Telomere, rDNA
}TYPE;

// 染色体末端粒子のリスト UP:上流 DOWN:下流 //
typedef enum chr_end {
    
    TELO1_UP = 0, TELO1_DOWN = 1115,
    TELO2_UP = 1116, TELO2_DOWN = 2023,
    rDNA_UP = 2024, rDNA_DOWN = 2515
    
}CHR_END;

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    int pastis_no;
    unsigned int chr_no;
    unsigned int particle_type;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_init[DIMENSION];
    double velocity[DIMENSION];
    //double velocity_h[DIMENSION];
    double force[DIMENSION];
    unsigned int list_no;
    unsigned int *list;
    
} Particle;

enum label{ X, Y, Z};

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

//  内積計算    //
double Inner_product (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    return ( pos_1[X] * pos_2[X] + pos_1[Y] * pos_2[Y] + pos_1[Z] * pos_2[Z]);
}

void secure_main_memory (Particle **part) {   // メモリ確保 //
    
    if ( (*part = (Particle *)malloc(NUMBER_MAX * sizeof(Particle))) == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    for (unsigned int loop = 0; loop < NUMBER_MAX; loop++) {
        
        ( *part )[loop].list = (unsigned int *) malloc (NUMBER_MAX * sizeof (unsigned int));
        
        if ( ( *part )[loop].list == NULL) {
            
            printf ("\n error : can not secure the memory part.list \n");
            exit (1);
        }
    }
    
}

void random_init_pos (Particle *part, dsfmt_t *dsfmt) {
    
    const double ratio = 1.5;
    const double X_max = ratio * MEMBRANE_AXIS_1, Y_max = ratio * MEMBRANE_AXIS_2, Z_max = ratio * MEMBRANE_AXIS_3;
    Particle *part_1;
    
    // 初期値位置の設定 粒子間距離が一定値未満ならやり直す //
    for (count=0; count < NUMBER_MAX; count++) {
        
        part_1 = &part [count];
        
        part_1[X] = -X_max + 2.0 * X_max * dsfmt_genrand_open_open (dsfmt);
        part_1[Y] = -Y_max + 2.0 * Y_max * dsfmt_genrand_open_open (dsfmt);
        part_1[Z] = -Z_max + 2.0 * Z_max * dsfmt_genrand_open_open (dsfmt);
        
        for (count2 = 0; count2 < count; count2++) {
            
            while (Euclid_norm ( part_1->position, part[count2].position ) > 1.0e-2) {
                
                part_1[X] = -X_max + 2.0 * X_max * dsfmt_genrand_open_open (dsfmt);
                part_1[Y] = -Y_max + 2.0 * Y_max * dsfmt_genrand_open_open (dsfmt);
                part_1[Z] = -Z_max + 2.0 * Z_max * dsfmt_genrand_open_open (dsfmt);
            }
        }
        
        
    }
    
}

int main ( int argc, char **argv) {
    
    Particle *part;
    
    
    
    // dSFMT
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand (&dsfmt, (unsigned)time (NULL));
    
}
