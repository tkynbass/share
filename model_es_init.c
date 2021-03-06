//
//  model_es.c
//  
//
//  Created by tkym on 2018/09/26.
//

////noma 3C data　から　核膜を楕円体、粒子径3.0e-8で初期構造を作る

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "dSFMT/dSFMT.h"
#include <time.h>
#include <omp.h>

#define DIMENSION (3) //次元
#define LENGTH (4.0e-8)   // 長さの単位
#define KBT (1.38064852e-23 / LENGTH / LENGTH) //ボルツマン
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
#define K_BOND_2 (1.0e-4)  //ひもの硬さ
#define DELTA (1.0e-11)  //刻み幅
#define PARTICLE_MYU (2.0 * DIMENSION * PI * PARTICLE_RADIUS * LENGTH * 0.000890 / 100 )    //粘性抵抗の強さ
#define MEMBRANE_EXCLUDE (1.0)     //膜との衝突
#define MEMBRANE_EXCLUDE_SPB (1.0) //SPBとの衝突

#define SPB_RADIUS ( 3.0 )      //SPBの半径
#define SPB_MYU (2.0 * DIMENSION * PI * SPB_RADIUS * LENGTH * 0.000890 / 100 )  //SPBの粘性


#define MEMBRANE_INIT_RADIUS (31.8)
// Ellipsoid axes parameter of nucleus & nucleolus //
#define MEMBRANE_AXIS_1 ( 1.64032e-6 / LENGTH )
#define MEMBRANE_AXIS_2 ( 0.85 * MEMBRANE_AXIS_1 )
#define MEMBRANE_AXIS_3 ( 0.7 * MEMBRANE_AXIS_1 )

#define NUCLEOLUS_AXIS_1 ( 0.98325e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.85 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.65 * NUCLEOLUS_AXIS_1 )


//k_bond2 k_expression
double set_al_1;

typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

const double origin[] = { 0.0, 0.0, 0.0};

typedef struct particle {           //構造体の型宣言
    CHAIN chr_no;
    TYPE particle_type;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_old[DIMENSION];
    double velocity[DIMENSION];
    double velocity_2[DIMENSION];
    double force[DIMENSION];
    int list_no;
    int *list;
    
} Particle;

typedef struct ellipsoid {      //楕円体主成分の構造体
    
    double al_1;
    double al_2;
    double al_3;
} Ellipsoid;

Particle *part;

Particle spb;

Ellipsoid mem;
Ellipsoid nuc;

enum label{ X, Y, Z};

void read_coordinate_init ( int start ){       //初期値設定
    
    int i, i_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    sprintf (filename, "fission_result_%d.dat", start);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++){
        
        fscanf(fpr, "%d %d %d %lf %lf %lf %lf %lf %lf", &i_dummy, &part[i].chr_no, &part[i].particle_type,
               &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
               &part[i].velocity[X], &part[i].velocity[Y], &part[i].velocity[Z]);
        fgets(dummy, 128, fpr);
        
        //printf ("%d %lf %lf %lf \n", i, part[i].position[X], part[i].position[Y], part[i].position[Z]);
        
    }
    
    fscanf (fpr, "%s %s %lf %lf %lf %lf %lf %lf", dummy, dummy, &spb.position[X], &spb.position[Y], &spb.position[Z],
            &spb.velocity[X], &spb.velocity[Y], &spb.velocity[Z]);
    
    fscanf (fpr, "%s %lf %lf %lf", dummy, &mem.al_1, &mem.al_2, &mem.al_3);
    
    fclose (fpr);
    
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

void spring (const Particle *part_1, const Particle *part_2, double force[DIMENSION]) {     //ばね
    
    double dist, dist_0;
    
    double f;
    
    Particle *part_3;
    
    part_3 = &spb;
    
    //dist_0 = Euclid_norm (part_1->position_init, part_2->position_init);
    
    if (part_1 == part_3) dist_0 = SPB_RADIUS + PARTICLE_RADIUS;
    else dist_0 = INIT_DISTANCE;
    
    dist = Euclid_norm (part_1->position, part_2->position);
    
    f = K_BOND * (dist_0 - dist) / dist;
    
    
    force[X] += f * (part_1->position[X] - part_2->position[X]);
    force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
    force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
}

void spb_list (Particle *part_1){        //リスト化
    
    int m = 0, j;
    double dist;
    
    Particle *part_2;
    
    for(j=0; j<NUMBER; j++){
        
        part_2 = &part[j];
        
        dist = Euclid_norm (part_1->position, part_2->position);
        
        if (dist < 5.0 * PARTICLE_RADIUS && part_2->particle_type != Centromere){
            
            m++;
            part_1->list_no = m;
            part_1->list[m] = j;
        }
    }
    
    if (m == 0){
        
        part_1->list_no = 0;
    }
}

void membrane_exclude ( Particle *part_1 ) {
    
    double dist = Euclid_norm (part_1->position, origin);
    
    double ellipsoid_dist = part_1->position[X] * part_1->position[X] / ( mem.al_1 * mem.al_1 )
    + part_1->position[Y] * part_1->position[Y] / ( mem.al_2 * mem.al_2 )
    + part_1->position[Z] * part_1->position[Z] / ( mem.al_3 * mem.al_3 );
    
    // 法線ベクトル
    double normal_vector[] = { 2.0 * part_1->position[X] / ( mem.al_1 * mem.al_1),
        2.0 * part_1->position[Y] / ( mem.al_2 * mem.al_2),
        2.0 * part_1->position[Z] / ( mem.al_3 * mem.al_3) };
    
    double normal_vector_norm = Euclid_norm (normal_vector, origin);
    
    double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * ( part_1->position[X] * normal_vector[X]
                                                              + part_1->position[Y] * normal_vector[Y]
                                                              + part_1->position[Z] * normal_vector[Z]);
    
    if ( ellipsoid_dist - 1 > 0 ) {
        
        part_1->force[X] += f * normal_vector[X] / normal_vector_norm;
        part_1->force[Y] += f * normal_vector[Y] / normal_vector_norm;
        part_1->force[Z] += f * normal_vector[Z] / normal_vector_norm;
    }
}

void membrane_fix ( Particle *part_1 ) {
    
    double dist = Euclid_norm (part_1->position, origin);
    
    double ellipsoid_dist = part_1->position[X] * part_1->position[X] / ( mem.al_1 * mem.al_1 )
    + part_1->position[Y] * part_1->position[Y] / ( mem.al_2 * mem.al_2 )
    + part_1->position[Z] * part_1->position[Z] / ( mem.al_3 * mem.al_3 );
    
    // 法線ベクトル
    double normal_vector[] = { 2.0 * part_1->position[X] / ( mem.al_1 * mem.al_1),
        2.0 * part_1->position[Y] / ( mem.al_2 * mem.al_2),
        2.0 * part_1->position[Z] / ( mem.al_3 * mem.al_3) };
    
    double normal_vector_norm = Euclid_norm (normal_vector, origin);
    
    double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * ( part_1->position[X] * normal_vector[X]
                                                              + part_1->position[Y] * normal_vector[Y]
                                                              + part_1->position[Z] * normal_vector[Z]);
    
    part_1->force[X] += f * normal_vector[X] / normal_vector_norm;
    part_1->force[Y] += f * normal_vector[Y] / normal_vector_norm;
    part_1->force[Z] += f * normal_vector[Z] / normal_vector_norm;
}

void rotate_position_z ( double pos[DIMENSION], const double theta ) {
    
    double pos_new[DIMENSION];
    
    pos_new[X] = cos (theta) * pos[X] - sin (theta) * pos[Z];
    pos_new[Z] = sin (theta) * pos[X] + cos (theta) * pos[Z];
    
    pos[X] = pos_new[X];
    pos[Z] = pos_new[Z];
}

void init_SPB_calculate ( dsfmt_t *dsfmt) {
    
    int k, j;
    double p1, p2, theta, psi, dist, f;
    
    Particle *part_2;
    
    p1 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
    
    spb.force[X] = p1 * sin(theta) / sqrt(DELTA);
    spb.force[Y] = p1 * cos(theta) / sqrt(DELTA);
    spb.force[Z] = p2 * sin(psi) / sqrt(DELTA);
    
    membrane_fix ( &spb );
    
    spring (&spb, &part[1880], spb.force);       //セントロメアとのバネによる力
    spring (&spb, &part[3561], spb.force);
    spring (&spb, &part[5542], spb.force);
    
    spb_list (&spb);
    
    //ひも粒子との排除体積
    if (spb.list_no != 0){
        
        for (j = 1; j <= spb.list_no; j++){
            
            k = spb.list[j];
            
            part_2 = &part[k];
            dist = Euclid_norm (spb.position, part_2->position);
            
            if (dist < PARTICLE_RADIUS + SPB_RADIUS){
                
                f = K_EXCLUDE * (PARTICLE_RADIUS + SPB_RADIUS - dist) / dist;
                
                spb.force[X] += f * (spb.position[X] - part_2->position[X]);
                spb.force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                spb.force[Z] += f * (spb.position[Z] - part_2->position[Z]);
            }
        }
    }
    
    //粘性抵抗
    spb.force[X] += - SPB_MYU * spb.velocity[X];
    spb.force[Y] += - SPB_MYU * spb.velocity[Y];
    spb.force[Z] += - SPB_MYU * spb.velocity[Z];
    
    spb.velocity_2[X] = spb.velocity[X] + DELTA * ( spb.force[X] - SPB_MYU * spb.velocity[X] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Y] = spb.velocity[Y] + DELTA * ( spb.force[Y] - SPB_MYU * spb.velocity[Y] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Z] = spb.velocity[Z] + DELTA * ( spb.force[Z] - SPB_MYU * spb.velocity[Z] ) / ( 2.0 * SPB_MASS );
    
    spb.position_old[X] = spb.position[X];
    spb.position_old[Y] = spb.position[Y];
    spb.position_old[Z] = spb.position[Z];
    
    spb.position[X] += DELTA * spb.velocity_2[X];
    spb.position[Y] += DELTA * spb.velocity_2[Y];
    spb.position[Z] += DELTA * spb.velocity_2[Z];
    
    
}

void SPB_calculate (const unsigned int l, dsfmt_t *dsfmt){
    
    int k, j;
    double p1, p2, theta, psi, dist, f;

    Particle *part_2;
    
    p1 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
    
    spb.force[X] = p1 * sin(theta) / sqrt(DELTA);
    spb.force[Y] = p1 * cos(theta) / sqrt(DELTA);
    spb.force[Z] = p2 * sin(psi) / sqrt(DELTA);
    
    membrane_fix ( &spb );
    
    spring (&spb, &part[1880], spb.force);       //セントロメアとのバネによる力
    spring (&spb, &part[3561], spb.force);
    spring (&spb, &part[5542], spb.force);
    
    if ( l%500 == 0) spb_list (&spb);
    
    //ひも粒子との排除体積
    if (spb.list_no != 0){
        
        for (j = 1; j <= spb.list_no; j++){
            
            k = spb.list[j];
            
            part_2 = &part[k];
            dist = Euclid_norm (spb.position, part_2->position);
            
            if (dist < PARTICLE_RADIUS + SPB_RADIUS){
                
                f = K_EXCLUDE * (PARTICLE_RADIUS + SPB_RADIUS - dist) / dist;
                
                spb.force[X] += f * (spb.position[X] - part_2->position[X]);
                spb.force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                spb.force[Z] += f * (spb.position[Z] - part_2->position[Z]);
            }
        }
    }
    
    
    spb.velocity[X] = (2.0 * SPB_MASS * spb.velocity_2[X] + DELTA * spb.force[X]) / (2.0 * SPB_MASS + SPB_MYU * DELTA);
    spb.velocity[Y] = (2.0 * SPB_MASS * spb.velocity_2[Y] + DELTA * spb.force[Y]) / (2.0 * SPB_MASS + SPB_MYU * DELTA);
    spb.velocity[Z] = (2.0 * SPB_MASS * spb.velocity_2[Z] + DELTA * spb.force[Z]) / (2.0 * SPB_MASS + SPB_MYU * DELTA);
    
    spb.velocity_2[X] = spb.velocity[X] + DELTA * ( spb.force[X] - SPB_MYU * spb.velocity[X] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Y] = spb.velocity[Y] + DELTA * ( spb.force[Y] - SPB_MYU * spb.velocity[Y] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Z] = spb.velocity[Z] + DELTA * ( spb.force[Z] - SPB_MYU * spb.velocity[Z] ) / ( 2.0 * SPB_MASS );
    
    spb.position_new[X] = spb.position[X] + DELTA * spb.velocity_2[X];
    spb.position_new[Y] = spb.position[Y] + DELTA * spb.velocity_2[Y];
    spb.position_new[Z] = spb.position[Z] + DELTA * spb.velocity_2[Z];
    
    
}

void init_particle_calculate( dsfmt_t *dsfmt/*, const unsigned int gene_list [CLUSTER_GENE_NUMBER] */){
    
    int i, k, j, m, gene_counter=0;
    
    double f, f_2, f_3, dist;
    double p1, p2, theta, psi;
    
    Particle *part_1, *part_2, *part_3;
    
    for (i = 0; i < NUMBER; i++){
        
        part_1 = &part[i];
        
        
        //noise
        p1 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
        p2 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
        theta = 2.0 * PI * dsfmt_genrand_open_close(dsfmt) ;
        psi = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
        
        part_1->force[X] = p1 * sin(theta) / sqrt(DELTA);
        part_1->force[Y] = p1 * cos(theta) / sqrt(DELTA);
        part_1->force[Z] = p2 * sin(psi) / sqrt(DELTA);
        
        part_1->force[X] += - PARTICLE_MYU * part_1->velocity[X];
        part_1->force[Y] += - PARTICLE_MYU * part_1->velocity[Y];
        part_1->force[Z] += - PARTICLE_MYU * part_1->velocity[Z];
        
        switch (part[i].particle_type) {
            case Normal:
                
                part_2 = &part[i-1];
                part_3 = &part[i+1];

                //spring
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //membrane_exclude
                membrane_exclude (part_1);
                
                //spb_exclude
                dist = Euclid_norm (part_1->position, spb.position);
                f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                
                if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
                    
                    part_1->force[X] += f * (part_1->position[X] - spb.position[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                }
                
                //spring2
                switch (i) {
                    case 1:
                    case 2772:
                    case 5013:
                        
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    case 2769:
                    case 5010:
                    case 6191:
                        
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-2];
                        part_3 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        dist = Euclid_norm (part_1->position, part_3->position);
                        f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        
                        break;
                }
                
                break;
                
            case Centromere:
                
                //centromere    spb とのバネ
                dist = Euclid_norm( part_1->position , spb.position);
                f = K_BOND * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                part_1->force[X] += f * (part_1->position[X] - spb.position[X]);
                part_1->force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                part_1->force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                
                //spring
                part_2 = &part[i-1];
                part_3 = &part[i+1];
                
                dist = Euclid_norm(part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //spring2
                part_2 = &part[i-2];
                part_3 = &part[i+2];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm (part_1->position, part_3->position);
                f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //membrane_exclude
                membrane_exclude ( part_1 );
                
                break;
                
            case Telomere:
                
                switch (i) {
                    case 0:
                    case 2771:
                    case 5012:
                        
                        part_2 = &part[i+1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 5012) { //telomere
                            
                            membrane_fix ( part_1 );
                            
                        }
                        else { //telomere_3
                            
                            membrane_fix ( part_1 );
                        }
                        
                        //spring2
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 6192) {//telomere
                            
                            membrane_fix ( part_1 );
                        }
                        else { //telomere_3
                            
                            membrane_fix ( part_1 );
                        }
                        
                        //spring2
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                }
                break;
        }
        
        //list
        
        m = 0;
        
        for(j=0; j<NUMBER; j++){
            
            part_2 = &part[j];
            
            dist = Euclid_norm (part_1->position, part_2->position);
            
            if (dist < 5.0 * PARTICLE_RADIUS && abs(i-j) > 1){
                
                m++;
                part_1->list_no = m;
                part_1->list[m] = j;
            }
        }
        if (m == 0){
            
            part_1->list_no = 0;
        }
        
        
        //particle_exclude
        if (part_1->list_no != 0){
            
            for (j = 1; j <= part_1->list_no; j++){
                
                k = part_1->list[j];
                
                part_2 = &part[k];
                dist = Euclid_norm (part_1->position, part_2->position);
                
                if (dist < 2 * PARTICLE_RADIUS){
                    
                    f = K_EXCLUDE * (2 * PARTICLE_RADIUS - dist) / dist;
                    
                    part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                }
                
            }
        }
        
        part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - PARTICLE_MYU * part_1->velocity[X] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / ( 2.0 * PARTICLE_MASS );
        
    }
    
    for ( i=0 ; i<NUMBER; i++) {
        
        part_1 = &part[i];
        
        part_1->position_old[X] = part_1->position[X];
        part_1->position_old[Y] = part_1->position[Y];
        part_1->position_old[Z] = part_1->position[Z];
        
        part_1->position[X] += DELTA * part_1->velocity_2[X];
        part_1->position[Y] += DELTA * part_1->velocity_2[Y];
        part_1->position[Z] += DELTA * part_1->velocity_2[Z];
    }
}

void particle_calculate( const unsigned int l, dsfmt_t *dsfmt /*, const unsigned int gene_list [CLUSTER_GENE_NUMBER] */)
//位置と速度の計算 private part_1->force dist f part_1 part_2 part_3
{
    int i, k, j, m, gene_counter = 0;
    
    double f, f_2, f_3, dist;
    double p1, p2, theta, psi;
    
    Particle *part_1, *part_2, *part_3;
    
    ////// noise /////
    for ( i=0; i<NUMBER; i++) {
        
        part_1 = &part[i];
        
        p1 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
        p2 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ));
        theta = 2.0 * PI * dsfmt_genrand_open_close(dsfmt) ;
        psi = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
        
        part_1->force[X] = p1 * sin(theta) / sqrt(DELTA);
        part_1->force[Y] = p1 * cos(theta) / sqrt(DELTA);
        part_1->force[Z] = p2 * sin(psi) / sqrt(DELTA);
        
    }
    
#pragma omp parallel for private ( j, k, m, gene_counter, p1, p2, theta, psi, dist, f, part_1, part_2, part_3, f_2, f_3) num_threads (8)
    for (i = 0; i < NUMBER; i++){
        
        part_1 = &part[i];
        
        switch (part[i].particle_type) {
            case Normal:
                
                part_2 = &part[i-1];
                part_3 = &part[i+1];
                
                
                //spring
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                membrane_exclude ( part_1 );
                
                //spb_exclude
                dist = Euclid_norm (part_1->position, spb.position);
                f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                
                if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
                    
                    part_1->force[X] += f * (part_1->position[X] - spb.position[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                }
                
                //spring2
                switch (i) {
                    case 1:
                    case 2772:
                    case 5013:
                        
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    case 2769:
                    case 5010:
                    case 6191:
                        
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-2];
                        part_3 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        dist = Euclid_norm (part_1->position, part_3->position);
                        f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        
                        break;
                }
                
                break;
                
            case Centromere:
                
                //centromere
                dist = Euclid_norm( part_1->position , spb.position);
                f = K_BOND * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                part_1->force[X] += f * (part_1->position[X] - spb.position[X]);
                part_1->force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                part_1->force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                
                //spring
                part_2 = &part[i-1];
                part_3 = &part[i+1];
                
                dist = Euclid_norm(part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //spring2
                part_2 = &part[i-2];
                part_3 = &part[i+2];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm (part_1->position, part_3->position);
                f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
 
                membrane_exclude ( part_1 );
                
                break;
                
            case Telomere:
                
                switch (i) {
                    case 0:
                    case 2771:
                    case 5012:
                        
                        part_2 = &part[i+1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 5012) { //telomere
                            
                            membrane_fix ( part_1 );
                        }
                        else { //telomere_3
                            
                            membrane_fix ( part_1 );
                        }
                        
                        //spring2
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 6192) {//telomere
                            
                            membrane_fix ( part_1 );
                        }
                        else { //telomere_3
                            
                            membrane_fix ( part_1 );
                        }
                        
                        //spring2
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                }
                break;
        }
        
        //list
        if ( l%500 == 0) {
            
            m = 0;
            
            for(j=0; j<NUMBER; j++){
                
                part_2 = &part[j];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                
                if (dist < 5.0 * PARTICLE_RADIUS && abs(i-j) > 1){
                    
                    m++;
                    part_1->list_no = m;
                    part_1->list[m] = j;
                }
            }
            if (m == 0){
                
                part_1->list_no = 0;
            }
        }
        
        //particle_exclude
        if (part_1->list_no != 0){
            
            for (j = 1; j <= part_1->list_no; j++){
                
                k = part_1->list[j];
                
                part_2 = &part[k];
                dist = Euclid_norm (part_1->position, part_2->position);
                
                if (dist < 2 * PARTICLE_RADIUS){
                    
                    f = K_EXCLUDE * (2 * PARTICLE_RADIUS - dist) / dist;
                    
                    part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                }
                
            }
        }
        
        
        part_1->velocity[X] = (2.0 * PARTICLE_MASS * part_1->velocity_2[X] + DELTA * part_1->force[X]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
        part_1->velocity[Y] = (2.0 * PARTICLE_MASS * part_1->velocity_2[Y] + DELTA * part_1->force[Y]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
        part_1->velocity[Z] = (2.0 * PARTICLE_MASS * part_1->velocity_2[Z] + DELTA * part_1->force[Z]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
        
        part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - PARTICLE_MYU * part_1->velocity[X] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / ( 2.0 * PARTICLE_MASS );
        
        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->velocity_2[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->velocity_2[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->velocity_2[Z];
        
        
    }
}

void renew () {
    
    int i;
    
    Particle *part_1;
    
    for(i = 0; i < NUMBER; i++){    //値の更新
        
        part_1 = &part[i];
        
        part_1->position_old[X] = part_1->position[X];
        part_1->position_old[Y] = part_1->position[Y];
        part_1->position_old[Z] = part_1->position[Z];
        
        part_1->position[X] = part_1->position_new[X];
        part_1->position[Y] = part_1->position_new[Y];
        part_1->position[Z] = part_1->position_new[Z];
    }
    
    spb.position_old[X] = spb.position [X];
    spb.position_old[Y] = spb.position[Y];
    spb.position_old[Z] = spb.position[Z];
    
    spb.position[X] = spb.position_new[X];
    spb.position[Y] = spb.position_new[Y];
    spb.position[Z] = spb.position_new[Z];
    
}

void membrane_to_ellipsoid (const unsigned int time) {
    
    // al_1 31.815 → 25, al_2 31.815 → 21.25, al_3 31.815 → 17.5
    const double delta = 5000;
    
    if (time <= delta) {
        
        mem.al_1 -= ( MEMBRANE_INIT_RADIUS - MEMBRANE_AXIS_1 ) / delta;
        mem.al_2 -= ( MEMBRANE_INIT_RADIUS - MEMBRANE_AXIS_2 ) / delta;
        mem.al_3 -= ( MEMBRANE_INIT_RADIUS - MEMBRANE_AXIS_3 ) / delta;
    }
}

void write_coordinate ( /*const char *number,*/ int t , int start) {
    
    int i;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "es_result_%d.txt", t + start);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++) {
        
        fprintf (fpw, "%d %d %d %lf %lf %lf %lf %lf %lf\n", i, part[i].chr_no, part[i].particle_type,
                 part[i].position_old[X], part[i].position_old[Y], part[i].position_old[Z], part[i].velocity[X], part[i].velocity[Y], part[i].velocity[Z]);
    }
    
    sprintf (str, "SPB");
    
    fprintf(fpw, "%s %s %lf %lf %lf %lf %lf %lf\n", str, str, spb.position_old[X], spb.position_old[Y], spb.position_old[Z],
            spb.velocity[X], spb.velocity[Y], spb.velocity[Z]);
    
    sprintf (str, "Axis_length");
    
    fprintf (fpw, "%s %lf %lf %lf\n", str, mem.al_1, mem.al_2, mem.al_3);
    
    fclose (fpw);
}

int main ( int argc, char **argv ) {
    
    int i, t = 0, l, membrane_flag, flag, rotate_flag;
    
    int start_number = atoi(argv[1]);
    int calculate_number = atoi(argv[2]);
    //k_expression = atof (argv[3]);
    
    unsigned int *gene_list;
    
    part = (Particle *)malloc(NUMBER * sizeof(Particle));
    
    if (part == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    for (i = 0;i < NUMBER;i++) {
        part[i].list = (int *)malloc(NUMBER * sizeof(int));
        
        if (part[i].list == NULL) {
            printf("\n error : can not secure the memory \n");
            exit(1);
        }
    }
    
    spb.list = (int *)malloc(NUMBER * sizeof(int));
    if (spb.list == NULL) {
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    /*
    gene_list = malloc (sizeof (int) * CLUSTER_GENE_NUMBER);
    if (gene_list == NULL) {
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    */

    //dSFMT
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, (unsigned)time(NULL));
    
    read_coordinate_init (start_number );
    
    //read_gene_list (gene_list);
    
    printf ("\n\t Radius of particles is %2fe-8, OK? (y:1 or n:0) : ", LENGTH * 1.0e+8);
    scanf ("%d", &flag);
    
    if (flag == 0) exit(1);
    
    printf ("\n     Membrane to ellipsoid(1) or mitigation(0) ? : ");
    scanf ("%d", &membrane_flag);
    
    printf ("\n     Rotate pos around z axis ? (y:1 or n:0) : ");
    scanf ("%d", &rotate_flag);
    
    if (rotate_flag == 1) {
        
        for (i=0; i<NUMBER; i++) {
            
            rotate_position_z (part[i].position, 107.0 * PI / 180.0);
        }
        rotate_position_z (spb.position, 107.0 * PI / 180.0);
    }
    
    
    if (membrane_flag == 1) {
        
        //核膜主成分の初期化
        mem.al_1 = MEMBRANE_INIT_RADIUS;
        mem.al_2 = MEMBRANE_INIT_RADIUS;
        mem.al_3 = MEMBRANE_INIT_RADIUS;
    }
    
    init_particle_calculate ( &dsfmt/*, gene_list*/);
    init_SPB_calculate( &dsfmt);
    
    //初期位置の出力
    write_coordinate (/* argv[3],*/ 0, start_number );
    
    for (t=1; t < calculate_number; t++) {
        
        if ( membrane_flag == 1 ) membrane_to_ellipsoid (t);
        
        for (l=1; l<=10000; l++){
            
            particle_calculate(l, &dsfmt /*, gene_list*/);
            SPB_calculate(l, &dsfmt);
            
            renew ();
        }
        
        printf("    t = %d, r = %lfe-8, al_1 = %lf, al_2 = %lf, al_3 = %lf \r", t, LENGTH*1.0e+8, mem.al_1, mem.al_2, mem.al_3);
        fflush (stdout);
        
        write_coordinate (/* argv[3],*/ t , start_number);
    }
    
    return ( 0 );
}



