// 粒子をランダム配置　→　楕円体近似の核内に緩和

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "dSFMT/dSFMT.h"
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
#define K_BOND_2 ( 1.0e-2 )     //2つ隣
#define K_BOND_3 ( 1.0e-2 )     //3つ隣
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

void labeling (Particle *part) {
    
    unsigned int loop;
    
    for ( loop = 0; loop < NUMBER_MAX; loop++)  {
        
        part[loop].pastis_no = -1;
        
        if ( loop <= 1115 ) part[loop].chr_no = A;
        else if ( loop <= 2023) part[loop].chr_no = B;
        else part[loop].chr_no = C;
        
    }
    
    for ( loop = 0; loop < NUMBER_MAX; loop++ ) part [loop].particle_type = Normal;
    
    for ( loop = 0; loop < sizeof (CENT_LIST) / sizeof (CENT_LIST[0]); loop++ ) part [CENT_LIST [loop]].particle_type = Centromere;
    for ( loop = 0; loop < sizeof (TELO_LIST) / sizeof (TELO_LIST[0]); loop++ ) part [TELO_LIST [loop]].particle_type = Telomere;
    for ( loop = 0; loop < sizeof (rDNA_LIST) / sizeof (rDNA_LIST[0]); loop++ ) part [rDNA_LIST [loop]].particle_type = rDNA;
}


void random_init_pos (Particle *part, dsfmt_t *dsfmt) {
    
    const double ratio = 1.5;
    const double X_max = ratio * MEMBRANE_AXIS_1, Y_max = ratio * MEMBRANE_AXIS_2, Z_max = ratio * MEMBRANE_AXIS_3;
    Particle *part_1;
    
    // 初期値位置の設定 粒子間距離が一定値未満ならやり直す //
    for (unsigned int count=0; count < NUMBER_MAX; count++) {
        
        part_1 = &part [count];
        
        part_1->position[X] = -X_max + 2.0 * X_max * dsfmt_genrand_open_open (dsfmt);
        part_1->position[Y] = -Y_max + 2.0 * Y_max * dsfmt_genrand_open_open (dsfmt);
        part_1->position[Z] = -Z_max + 2.0 * Z_max * dsfmt_genrand_open_open (dsfmt);
        
        for (unsigned int count2 = 0; count2 < count; count2++) {
            
            while (Euclid_norm ( part_1->position, part[count2].position ) < 1.0e-2) {
                
                part_1->position[X] = -X_max + 2.0 * X_max * dsfmt_genrand_open_open (dsfmt);
                part_1->position[Y] = -Y_max + 2.0 * Y_max * dsfmt_genrand_open_open (dsfmt);
                part_1->position[Z] = -Z_max + 2.0 * Z_max * dsfmt_genrand_open_open (dsfmt);
            }
        }
        
        
    }
    
}

//　ばねによる力 part_1側の力計算//
void spring (Particle *part_1, const Particle *part_2, unsigned int interval) {
    
    // 線形バネの強さ　0:spb-centromere, 1,2,3: n個隣 //
    const double bonding_power[] = { K_BOND, K_BOND, K_BOND_2, K_BOND_3 };
    double dist, f, dist_0;
    
    switch (interval) {
        case 0:
        
            dist_0 =  PARTICLE_RADIUS + SPB_RADIUS;
        
            dist = Euclid_norm (part_1->position, SPB_POS);
        
            f = bonding_power[interval] * (dist_0 - dist) / dist;
        
            part_1->force[X] += f * (part_1->position[X] - SPB_POS[X]);
            part_1->force[Y] += f * (part_1->position[Y] - SPB_POS[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - SPB_POS[Z]);
            break;
        
//        case 1:
//
//            dist_0 = BOND_DISTANCE * interval;
//            dist = Euclid_norm (part_1->position, part_2->position);
//
//            f = bonding_power[interval] * (dist_0 - dist) / dist;
//
//            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
//            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
//            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
//            break;
//
        default:
                
            dist_0 = BOND_DISTANCE * interval;
            dist = Euclid_norm (part_1->position, part_2->position);
        
            f = bonding_power[interval] * (dist_0 - dist) / dist;
        
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                
            break;
    }
}

void membrane_interaction_fix ( Particle *part_1, char interaction_type /* F: fix, E: exclude */) {
    
    double dist = Euclid_norm (part_1->position, ORIGIN);
    
    double ellipsoid_dist = part_1->position[X] * part_1->position[X] / ( MEMBRANE_AXIS_1 * MEMBRANE_AXIS_1 )
    + part_1->position[Y] * part_1->position[Y] / ( MEMBRANE_AXIS_2 * MEMBRANE_AXIS_2 )
    + part_1->position[Z] * part_1->position[Z] / ( MEMBRANE_AXIS_3 * MEMBRANE_AXIS_3 );
    
    if ( interaction_type == 'F' || ellipsoid_dist - 1 > 0 ) {
        
        // 法線ベクトル
        double normal_vector[] = { 2.0 * part_1->position[X] / ( MEMBRANE_AXIS_1 * MEMBRANE_AXIS_1),
            2.0 * part_1->position[Y] / ( MEMBRANE_AXIS_2 * MEMBRANE_AXIS_2),
            2.0 * part_1->position[Z] / ( MEMBRANE_AXIS_3 * MEMBRANE_AXIS_3) };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * ( part_1->position[X] * normal_vector[X]
                                                                  + part_1->position[Y] * normal_vector[Y]
                                                                  + part_1->position[Z] * normal_vector[Z]);
        
        part_1->force[X] += f * normal_vector[X] / normal_vector_norm;
        part_1->force[Y] += f * normal_vector[Y] / normal_vector_norm;
        part_1->force[Z] += f * normal_vector[Z] / normal_vector_norm;
    }
}

// nucleolus interaction //
void nucleolus_interaction ( Particle *part_1, const char interaction_type ) {
    
    //核小体中心から粒子へのベクトル
    double nuc_to_pos[DIMENSION] = { part_1->position[X] - NUCLEOLUS_POS[X],
        part_1->position[Y] - NUCLEOLUS_POS[Y],
        part_1->position[Z] - NUCLEOLUS_POS[Z]};
    
    //位置座標をz軸まわりに-10度回転
    rotate_about_z (nuc_to_pos, - PI / 18);
    
    double ellipsoid_dist =  nuc_to_pos[X] * nuc_to_pos[X] / ( NUCLEOLUS_AXIS_1 * NUCLEOLUS_AXIS_1 )
    + nuc_to_pos[Y] * nuc_to_pos[Y] / ( NUCLEOLUS_AXIS_3 * NUCLEOLUS_AXIS_3 )
    + nuc_to_pos[Z] * nuc_to_pos[Z] / ( NUCLEOLUS_AXIS_2 * NUCLEOLUS_AXIS_2 );
    
    if ( interaction_type == 'F' || ellipsoid_dist < 1.0 ) {
        
        // 法線ベクトル
        double normal_vector[] = { 2.0 * nuc_to_pos[X] / ( NUCLEOLUS_AXIS_1 * NUCLEOLUS_AXIS_1),
            2.0 * nuc_to_pos[Y] / ( NUCLEOLUS_AXIS_3 * NUCLEOLUS_AXIS_3),
            2.0 * nuc_to_pos[Z] / ( NUCLEOLUS_AXIS_2 * NUCLEOLUS_AXIS_2) };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * Inner_product (nuc_to_pos, normal_vector);
        
        rotate_about_z (normal_vector, - PI / 18.0);
        
        part_1->force[X] += f * normal_vector[X] / normal_vector_norm;
        part_1->force[Y] += f * normal_vector[Y] / normal_vector_norm;
        part_1->force[Z] += f * normal_vector[Z] / normal_vector_norm;
    }
}

void spb_exclusion (Particle *part_1) {
    
    // spb_exclusion //
    double dist = Euclid_norm (part_1->position, SPB_POS);
    double f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
    
    if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
        
        part_1->force[X] += f * (part_1->position[X] -  SPB_POS[X]);
        part_1->force[Y] += f * (part_1->position[Y] -  SPB_POS[Y]);
        part_1->force[Z] += f * (part_1->position[Z] -  SPB_POS[Z]);
    }
}

//　ひも同士のリスト化 //
void make_ve_list (Particle *part, Particle *part_1, const unsigned int target) {
    
    unsigned int loop, list_count = 0;
    double dist;
    Particle *part_2;
    
    part_1->list_no = 0;
    
    for ( loop = 0; loop < NUMBER_MAX; loop++ ){
        
        part_2 = &part[loop];
        dist = Euclid_norm ( part_1->position, part_2->position);
        
        // 1個隣 そうでなくてもテロメア同士 (番号は隣だが染色体No.が異なる) //
        if ( dist < LIST_RADIUS && (abs (target - loop) > 1 || part_1->chr_no != part_2->chr_no ) ) {
            
            list_count++;
            part_1->list_no = list_count;
            part_1->list [list_count] = loop;
        }
    }
}

// Volume exclusion between particles //
void particle_exclusion (Particle *part, Particle *part_1) {
    
    unsigned int loop;
    double dist, f;
    Particle *part_2;
    
    for ( loop = 1; loop <= part_1->list_no; loop++) {
        
        part_2 = &part [part_1->list [loop]];
        dist = Euclid_norm (part_1->position, part_2->position);
        
        if ( dist < 2.0 * PARTICLE_RADIUS ){
            
            f = K_EXCLUDE * (2.0 * PARTICLE_RADIUS - dist) / dist;
            
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
        }
    }
}

// 各stepごとの座標計算 //
void calculation (Particle *part, const unsigned int mitigation, double mem_al[3] ) {
    
    unsigned int loop;
    Particle *part_1, *part_2, *part_3;
    
    // 位置の計算 & 力の初期化 //
    for ( loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        /*
         part_1->position[X] += DELTA * part_1->velocity_h[X];
         part_1->position[Y] += DELTA * part_1->velocity_h[Y];
         part_1->position[Z] += DELTA * part_1->velocity_h[Z];
         */
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
    }
    
#pragma omp parallel for private (part_1) num_threads (8)
    for ( loop = 0; loop < NUMBER_MAX; loop++ ){
        
        part_1 = &part [loop];
        
        switch (part_1->particle_type) {
            case Normal:
            
                // spring 1 //
                spring (part_1, &part [loop + 1], 1);
                spring (part_1, &part [loop - 1], 1);
            
                // spring 2 & 3 //
                switch (loop)  {
                    
                    case TELO1_UP + 1:
                    case TELO2_UP + 1:
                    case rDNA_UP + 1:
                    
                        spring (part_1, &part[loop + 2], 2);
                    
                        spring (part_1, &part[loop + 3], 3);
                        break;
                    
                    case TELO1_UP + 2:
                    case TELO2_UP + 2:
                    case rDNA_UP + 2:
                    
                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);
                    
                        spring (part_1, &part[loop + 3], 3);
                        break;
                    
                    case TELO1_DOWN - 1:
                    case TELO2_DOWN - 1:
                    case rDNA_DOWN - 1:
                    
                        spring (part_1, &part[loop - 2], 2);
                    
                        spring (part_1, &part[loop - 3], 3);
                        break;
                    
                    case TELO1_DOWN - 2:
                    case TELO2_DOWN - 2:
                    case rDNA_DOWN - 2:
                    
                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);
                    
                        spring (part_1, &part[loop - 3], 3);
                        break;
                    
                    default:
                    
                    spring (part_1, &part[loop + 2], 2);
                    spring (part_1, &part[loop - 2], 2);
                    
                    spring (part_1, &part[loop + 3], 3);
                    spring (part_1, &part[loop - 3], 3);
                    break;
                }
            
                spb_exclusion (part_1);
                nucleolus_interaction (part_1, 'E');
                membrane_interaction_change (part_1, 'E', mem_al);
            
                break;
            
            case Centromere:
            
                spring (part_1, &part[loop + 1], 1);
                spring (part_1, &part[loop - 1], 1);
            
                spring (part_1, &part[loop + 2], 2);
                spring (part_1, &part[loop - 2], 2);
            
                spring (part_1, &part[loop + 3], 3);
                spring (part_1, &part[loop - 3], 3);
            
                nucleolus_interaction (part_1, 'E');
                membrane_interaction_change (part_1, 'E', mem_al);
                spring (part_1, NULL, 0);
            
                break;
            
            case Telomere:
            
                switch (loop) {
                    case TELO1_UP:
                    case TELO2_UP:
                    
                        spring (part_1, &part[loop + 1], 1);
                    
                        spring (part_1, &part[loop + 2], 2);
                    
                        spring (part_1, &part[loop + 3], 3);
                    
                        break;
                    
                    case TELO1_DOWN:
                    case TELO2_DOWN:
                    
                        spring (part_1, &part[loop - 1], 1);
                    
                        spring (part_1, &part[loop - 2], 2);
                    
                        spring (part_1, &part[loop - 3], 3);
                    
                        break;
                }
            
                spb_exclusion (part_1);
                membrane_interaction_change (part_1, 'F', mem_al);
                nucleolus_interaction (part_1, 'E');
            
                break;
            
            case rDNA:
            
                switch (loop) {
                    case rDNA_UP:
                    
                        spring (part_1, &part[loop + 1], 1);
                    
                        spring (part_1, &part[loop + 2], 2);
                    
                        spring (part_1, &part[loop + 3], 3);
                    
                        break;
                    
                    case rDNA_DOWN:
                    
                        spring (part_1, &part[loop - 1], 1);
                    
                        spring (part_1, &part[loop - 2], 2);
                    
                        spring (part_1, &part[loop - 3], 3);
                    
                        break;
                }
            
                spb_exclusion (part_1);
                membrane_interaction_change (part_1, 'E', mem_al);
                nucleolus_interaction (part_1, 'F');
            
                break;
            
            default:
                printf ("\t Labeling error occured. \n");
                exit(1);
        }
        
        if ( mitigation % LIST_INTERVAL == 0 ) make_ve_list (part, part_1, loop);
        particle_exclusion (part, part_1);
        
        /*
         part_1->velocity[X] = ( 2.0 * PARTICLE_MASS * part_1->velocity_h[X] + DELTA * part_1->force[X] ) / ( 2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
         part_1->velocity[Y] = ( 2.0 * PARTICLE_MASS * part_1->velocity_h[Y] + DELTA * part_1->force[Y] ) / ( 2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
         part_1->velocity[Z] = ( 2.0 * PARTICLE_MASS * part_1->velocity_h[Z] + DELTA * part_1->force[Z] ) / ( 2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
         
         part_1->velocity_h[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - PARTICLE_MYU * part_1->velocity[X] ) / (2.0 * PARTICLE_MASS);
         part_1->velocity_h[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / (2.0 * PARTICLE_MASS);
         part_1->velocity_h[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / (2.0 * PARTICLE_MASS);
         */
    }
    
    for ( loop = 0; loop < NUMBER_MAX; loop++ ) {
        
        part_1 = &part[loop];
        
        part_1->position[X] += DELTA * part_1->force[X];
        part_1->position[Y] += DELTA * part_1->force[Y];
        part_1->position[Z] += DELTA * part_1->force[Z];
    }
}


void write_coordinate (Particle *part, const unsigned int time) {

    unsigned int loop, mitigation, start, calculation_max;
    char output_file[256];

    Particle *part_1;

    FILE *fpw;

    char result[128], str[128];

    sprintf (result, "result_%d.txt", time);

    if ((fpw = fopen (result, "w")) == NULL) {

        printf (" \t error : cannot write coordinate. \n");

        exit (1);
    }

    for (loop = 0; loop < NUMBER_MAX; loop++) {

        part_1 = &part[loop];
        fprintf (fpw, "%d 0 %d %d %lf %lf %lf\n", loop, part_1->chr_no, part_1->particle_type,
                 part_1->position[X],part_1->position[Y], part_1->position[Z]);
    }

    fprintf (fpw, "Radius %1.1e\n", LENGTH);
    fclose (fpw);
}

void write_coordinate (Particle *part) {

    unsigned int loop, mitigation, start, calculation_max;
    char output_file[256];
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "result_test.txt");
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \t error : cannot write coordinate. \n");
        
        exit (1);
    }
    
    for (loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        fprintf (fpw, "%d 0 %d %d %lf %lf %lf\n", loop, part_1->chr_no, part_1->particle_type,
                 part_1->position[X],part_1->position[Y], part_1->position[Z]);
    }
    
    fprintf (fpw, "Radius %1.1e\n", LENGTH);
    fclose (fpw);
}

int main ( int argc, char **argv) {
    
    unsigned int loop;
    Particle *part, *part_1;
    
    // dSFMT
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand (&dsfmt, (unsigned)time (NULL));
    
    
    secure_main_memory (&part);
    
    labeling (part);
    
    random_init_pos (part, &dsfmt);
    
    write_coordinate (part, 0)
    
    printf ("\n\t PARTICLE_RADIUS = %1.1e, DELTA = %1.1e, mitigation = %1.1e \n\n", LENGTH, DELTA, MITIGATION_INTERVAL);
    
    for ( unsigned int time = 1; time <= calculation_max; time++) {

        printf ("\t Now calculating...  time = %d \r", time);
        fflush (stdout);
        
        for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){
            
            calculation (part, mitigation, mem_al);
        }
        
        write_coordinate (part, time);
    }
    
    // メモリ解放 //
    for ( loop = 0 ; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part [loop];
        free (part_1->list);
    }
    free (part);
    
    return ( 0 );
    
}
