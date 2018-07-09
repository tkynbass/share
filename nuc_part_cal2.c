//
//  model_cal.c
//  
//  核小体を粒子でおく　→　核内に収める
//
//  Created by tkym on 2018/05/07.
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

#define INIT_NUMBER (6193)  //rDNA領域を除いた粒子数
#define NUMBER (6743)      //rDNA領域を含めた粒子数
#define PI (M_PI)
#define MEMBRAIN_RADIUS_MINIMUM  (32.0) //膜の半径の最小値
#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )

#define K_EXCLUDE (1.0)    //排除体積効果の強さ
#define K_BOND (1.0)    //ばね定数
#define K_BOND_2 (1.0e-4)  //ひもの硬さ
#define DELTA (1.0e-11)  //刻み幅
#define MEMBRAIN_EXCLUDE (1.0)     //膜との衝突
#define MEMBRAIN_EXCLUDE_SPB (1.0) //SPBとの衝突

#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PARTICLE_MASS (M_A / N_A / 1000)      //染色体粒子の質量
#define PARTICLE_MYU (2.0 * DIMENSION * PI * PARTICLE_RADIUS * NANO * 0.000890 / 100 )    //染色体粒子粘性抵抗の強さ

#define SPB_RADIUS ( 3.0 )      //SPBの半径
#define SPB_MASS ( 27.0 * PARTICLE_MASS)      //SPBの質量
#define SPB_MYU (2.0 * DIMENSION * PI * SPB_RADIUS * NANO * 0.000890 / 100 )  //SPBの粘性

#define rDNA_RADIUS ( 2.0 )     //核小体粒子の半径
#define rDNA_MASS ( rDNA_RADIUS * rDNA_RADIUS * rDNA_RADIUS * PARTICLE_MASS)
#define rDNA_MYU (2.0 * DIMENSION * PI * rDNA_RADIUS * NANO * 0.000890 / 100 )     //核小体粒子の粘性

double nucleolus_setting_radius;

double Nucleolus_circle_center[3], k_expression;

unsigned int simulate_type;

typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere, rDNA
}TYPE;

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

Particle *part;

Particle spb;

double membrain_radius;

enum label{ X, Y, Z};

void init_particle( int start ){       //初期値設定
    
    int i, i_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    printf ("\n     reading data_name : (fission_result / 0) or (result /1) \n");
    scanf ("%d", &simulate_type);
    
    if (simulate_type == 0) {
        
        sprintf (filename, "nucleolus_particle2/fission_result_%d.dat", start);
        
        if ((fpr = fopen(filename, "r")) == NULL){
            
            printf ("\n     error : can not read \n");
            
            exit (1);
        }
        
        for (i=0; i<INIT_NUMBER; i++){
            
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
        
        fscanf (fpr, "%s %s %lf %lf %lf", dummy, dummy, &Nucleolus_circle_center[X], &Nucleolus_circle_center[Y], &Nucleolus_circle_center[Z]);
    }
    else {
        
        sprintf (filename, "nucleolus_particle2/result_%d.dat", start);
        
        if ((fpr = fopen(filename, "r")) == NULL){
            
            printf ("\n     error : can not read \n");
            
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
        
        fscanf (fpr, "%s %s %lf %lf %lf", dummy, dummy, &Nucleolus_circle_center[X], &Nucleolus_circle_center[Y], &Nucleolus_circle_center[Z]);
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

void init_SPB_calculate (dsfmt_t dsfmt) {
    
    int k, j;
    
    double origin[] = { 0.0, 0.0, 0.0};
    
    double dist = Euclid_norm ( spb.position , origin);
    double p1, p2, theta, psi;
    double f = MEMBRAIN_EXCLUDE * (membrain_radius - dist ) / dist;
    
    Particle *part_2;
    
    p1 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    
    spb.force[X] += f * (spb.position[X]);        //膜とのバネ
    spb.force[Y] += f * (spb.position[Y]);
    spb.force[Z] += f * (spb.position[Z]);
    
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
            
            if (part_2->particle_type != rDNA) {
                
                if (dist < PARTICLE_RADIUS + SPB_RADIUS){
                    
                    f = K_EXCLUDE * (PARTICLE_RADIUS + SPB_RADIUS - dist) / dist;
                    
                    spb.force[X] += f * (spb.position[X] - part_2->position[X]);
                    spb.force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                    spb.force[Z] += f * (spb.position[Z] - part_2->position[Z]);
                }
            }
            else {
                
                if (dist < rDNA_RADIUS + SPB_RADIUS){
                    
                    f = K_EXCLUDE * (rDNA_RADIUS + SPB_RADIUS - dist) / dist;
                    
                    spb.force[X] += f * (spb.position[X] - part_2->position[X]);
                    spb.force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                    spb.force[Z] += f * (spb.position[Z] - part_2->position[Z]);
                }
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

void SPB_calculate (dsfmt_t dsfmt, const unsigned int l){
    
    int k, j;
    
    double origin[] = { 0.0, 0.0, 0.0};
    
    double dist = Euclid_norm ( spb.position , origin);
    double p1, p2, theta, psi;
    double f = MEMBRAIN_EXCLUDE * (membrain_radius - dist ) / dist;
    
    Particle *part_2;
    
    p1 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    
    spb.force[X] += f * (spb.position[X]);        //膜とのバネ
    spb.force[Y] += f * (spb.position[Y]);
    spb.force[Z] += f * (spb.position[Z]);
    
    spring (&spb, &part[1880], spb.force);       //セントロメアとのバネによる力
    spring (&spb, &part[3561], spb.force);
    spring (&spb, &part[5542], spb.force);
    
    if ( l%100 == 0) spb_list (&spb);
    
    //ひも粒子との排除体積
    if (spb.list_no != 0){
        
        for (j = 1; j <= spb.list_no; j++){
            
            k = spb.list[j];
            
            part_2 = &part[k];
            dist = Euclid_norm (spb.position, part_2->position);
            
            if (part_2->particle_type != rDNA) {
                
                if (dist < PARTICLE_RADIUS + SPB_RADIUS){
                    
                    f = K_EXCLUDE * (PARTICLE_RADIUS + SPB_RADIUS - dist) / dist;
                    
                    spb.force[X] += f * (spb.position[X] - part_2->position[X]);
                    spb.force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                    spb.force[Z] += f * (spb.position[Z] - part_2->position[Z]);
                }
            }
            else {
                
                if (dist < rDNA_RADIUS + SPB_RADIUS){
                    
                    f = K_EXCLUDE * (rDNA_RADIUS + SPB_RADIUS - dist) / dist;
                    
                    spb.force[X] += f * (spb.position[X] - part_2->position[X]);
                    spb.force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                    spb.force[Z] += f * (spb.position[Z] - part_2->position[Z]);
                }
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

void init_particle_calculate( dsfmt_t dsfmt/*, const unsigned int gene_list [CLUSTER_GENE_NUMBER] */){
    
    int i, k, j, m, gene_counter=0;
    
    double f, f_2, f_3, dist, origin[] = {0.0, 0.0, 0.0};
    double p1, p2, theta, psi;
    
    Particle *part_1, *part_2, *part_3;
    
    for (i = INIT_NUMBER; i < NUMBER; i++){
        
        part_1 = &part[i];
        
        //noise
        if (simulate_type ==0 ) {
            
            part_1->force[X] = 0.0;
            part_1->force[Y] = 0.0;
            part_1->force[Z] = 0.0;
        }
        else {
            
            p1 = sqrt(2.0 * 3.0 * rDNA_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
            p2 = sqrt(2.0 * 3.0 * rDNA_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
            theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
            psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
            
            part_1->force[X] = p1 * sin(theta) / sqrt(DELTA);
            part_1->force[Y] = p1 * cos(theta) / sqrt(DELTA);
            part_1->force[Z] = p2 * sin(psi) / sqrt(DELTA);
            
            part_1->force[X] += - rDNA_MYU * part_1->velocity[X];
            part_1->force[Y] += - rDNA_MYU * part_1->velocity[Y];
            part_1->force[Z] += - rDNA_MYU * part_1->velocity[Z];
        }
        
        /*
         if (i == gene_list [gene_counter]) {    //発現量が上がる遺伝子に核中心方向の力を加える
         
         high_expression (part_1, part_1->force);
         
         gene_counter++;
         }*/
        
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
                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] /*- 0.0*/);
                    part_1->force[Y] += f * (part_1->position[Y] /*- 0.0*/);
                    part_1->force[Z] += f * (part_1->position[Z] /*- 0.0*/);
                }
                
                //spb_exclude
                dist = Euclid_norm (part_1->position, spb.position);
                f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                
                if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
                    
                    part_1->force[X] += f * (part_1->position[X] - spb.position[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                }
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                
                if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
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
                
                case Centromere:    //セントロメア粒子への相互作用
                
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
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                
                if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                }
                
                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] );
                    part_1->force[Y] += f * (part_1->position[Y] );
                    part_1->force[Z] += f * (part_1->position[Z] );
                }
                
                break;
                
                case Telomere:      //テロメア粒子との相互作用
                
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
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - origin[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                                
                                part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                                part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                                part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                            }
                            
                        }
                        else { //telomere_3 核小体との結合
                            
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
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
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - origin[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                                
                                part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                                part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                                part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                            }
                        }
                        else { //telomere_3     核小体との結合
                            
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
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
                
                case rDNA:      //核小体粒子の相互作用
                
                switch (i) {
                        
                        //
                        case 6468:
                        case 6193:
                        
                        if (i==6468) {
                            
                            part_2 = &part[5012];
                            dist = Euclid_norm (part_1->position, part_2->position);
                            
                            f_2 = K_BOND * ((PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                            part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                            part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        }
                        else{
                            
                            part_2 = &part[6192];
                            dist = Euclid_norm (part_1->position, part_2->position);
                            
                            f_2 = K_BOND * ((PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                            part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                            part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        }
                        
                        /*
                        //spring2
                        if ( i == 6468 ) part_2 = &part[5013];
                        else  part_2 = &part[i-2];
                        
                        part_3 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND_2 * ( (rDNA_RADIUS + 3.0 * PARTICLE_RADIUS) * 0.8 - dist) / dist;
                        
                        dist = Euclid_norm (part_1->position, part_3->position);
                        f_3 = K_BOND_2 * ( 2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        */
                        
                        /*
                         //核小体領域表面を滑る
                         dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                         
                         f = K_BOND * ( 4.0 * membrain_radius - dist - rDNA_RADIUS ) / dist;
                         
                         part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                         part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                         part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                         */
                        
                        //membrain_exclude
                        dist = Euclid_norm (part_1->position, origin);
                        f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - rDNA_RADIUS ) / dist;
                        
                        if ( dist + rDNA_RADIUS > membrain_radius) {
                            
                            part_1->force[X] += f * (part_1->position[X] );
                            part_1->force[Y] += f * (part_1->position[Y] );
                            part_1->force[Z] += f * (part_1->position[Z] );
                        }
                        
                        break;
                        
                        ///////核膜結合側　　テロメア/////
                        case 6742:
                        case 6467:
                        
                        part_2 = &part[i-1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * ( 2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        //Nucleolus_exclude
                        dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                        f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + rDNA_RADIUS - dist ) / dist;
                        
                        if ( dist < 4.0 * membrain_radius + rDNA_RADIUS ) {
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                        }
                        
                        /*
                        //spring2
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        */
                         
                        /*
                         //telomere
                         dist = Euclid_norm (part_1->position, origin);
                         
                         f = K_BOND * ( membrain_radius - dist - rDNA_RADIUS ) / dist;
                         
                         part_1->force[X] += f * (part_1->position[X] - origin[X]);
                         part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                         part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                         */
                        
                        //nucleolus setting
                        dist = Euclid_norm (part_1->position, origin);
                        
                        f = K_BOND * ( nucleolus_setting_radius - dist - rDNA_RADIUS ) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - origin[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-1];
                        part_3 = &part[i+1];
                        
                        //spring
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND * (2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        dist = Euclid_norm(part_1->position, part_3->position);
                        f_3 = K_BOND * (2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        
                        /*
                        //ひもの硬さポテンシャル spring2
                        switch (i) {
                                
                            case 6469:
                            case 6194:
                            
                            if ( i == 6469 ) part_2 = &part[5012];
                            else  part_2 = &part[i-2];
                            
                            part_3 = &part[i+2];
                            
                            dist = Euclid_norm (part_1->position, part_2->position);
                            f_2 = K_BOND_2 * ( (3.0 * rDNA_RADIUS + PARTICLE_RADIUS) * 0.8 - dist) / dist;
                            
                            dist = Euclid_norm (part_1->position, part_3->position);
                            f_3 = K_BOND_2 * ( 4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                            + f_3 * (part_1->position[X] - part_3->position[X]);
                            part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                            + f_3 * (part_1->position[Y] - part_3->position[Y]);
                            part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                            + f_3 * (part_1->position[Z] - part_3->position[Z]);
                            
                            break;
                            
                            case 6741:
                            case 6466:
                            
                            part_2 = &part[i-2];
                            
                            dist = Euclid_norm (part_1->position, part_2->position);
                            f = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                            
                            break;
                            
                            default:
                            
                            part_2 = &part[i-2];
                            part_3 = &part[i+3];
                            
                            dist = Euclid_norm (part_1->position, part_2->position);
                            f_2 = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                            
                            dist = Euclid_norm (part_1->position, part_3->position);
                            f_3 = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                            + f_3 * (part_1->position[X] - part_3->position[X]);
                            part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                            + f_3 * (part_1->position[Y] - part_3->position[Y]);
                            part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                            + f_3 * (part_1->position[Z] - part_3->position[Z]);
                            
                            break;
                        }
                        */
                         
                        //Nucleolus_exclude
                        dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                        f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + rDNA_RADIUS - dist ) / dist;
                        
                        if ( dist < 4.0 * membrain_radius + rDNA_RADIUS ) {
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                        }
                        
                        /*
                         //membrain_exclude
                         dist = Euclid_norm (part_1->position, origin);
                         f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - rDNA_RADIUS ) / dist;
                         
                         if ( dist + rDNA_RADIUS > membrain_radius) {
                         
                         part_1->force[X] += f * (part_1->position[X] );
                         part_1->force[Y] += f * (part_1->position[Y] );
                         part_1->force[Z] += f * (part_1->position[Z] );
                         }
                         */
                        
                        //nucleolus setting
                        dist = Euclid_norm (part_1->position, origin);
                        f = MEMBRAIN_EXCLUDE * (nucleolus_setting_radius - dist - rDNA_RADIUS ) / dist;
                        
                        if ( dist + rDNA_RADIUS > nucleolus_setting_radius) {
                            
                            part_1->force[X] += f * (part_1->position[X] );
                            part_1->force[Y] += f * (part_1->position[Y] );
                            part_1->force[Z] += f * (part_1->position[Z] );
                        }
                        
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
                
                if ( part_1->particle_type != rDNA && part_2->particle_type != rDNA) {
                    
                    if (dist < 2.0 * PARTICLE_RADIUS ){
                        
                        f = K_EXCLUDE * (2 * PARTICLE_RADIUS - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                    }
                }
                else if (part_1->particle_type == rDNA && part_2->particle_type == rDNA) {
                    
                    if (dist < 2.0 * rDNA_RADIUS){
                        
                        f = K_EXCLUDE * ( 2.0 * rDNA_RADIUS - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                    }
                }
                else {
                        
                    if (dist < (PARTICLE_RADIUS + rDNA_RADIUS)){
                        
                        f = K_EXCLUDE * ( (PARTICLE_RADIUS + rDNA_RADIUS) - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                    }
                }
            }
        }
        
        if (part_1->particle_type != rDNA) {
            
            part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - PARTICLE_MYU * part_1->velocity[X] ) / ( 2.0 * PARTICLE_MASS );
            part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / ( 2.0 * PARTICLE_MASS );
            part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / ( 2.0 * PARTICLE_MASS );
        }
        else {
            
            part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - rDNA_MYU * part_1->velocity[X] ) / ( 2.0 * rDNA_MASS );
            part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - rDNA_MYU * part_1->velocity[Y] ) / ( 2.0 * rDNA_MASS );
            part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - rDNA_MYU * part_1->velocity[Z] ) / ( 2.0 * rDNA_MASS );
        }
        
        
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

void particle_calculate( dsfmt_t dsfmt, const unsigned int l/*, const unsigned int gene_list [CLUSTER_GENE_NUMBER]*/)        //位置と速度の計算 private dist f part_1 part_2 part_3
{
    int i, k, j, m, gene_counter = 0;
    
    double force[DIMENSION], f, f_2, f_3, dist, origin[] = {0.0, 0.0, 0.0};
    double p1, p2, theta, psi;
    
    Particle *part_1, *part_2, *part_3;
    
    for (i=INIT_NUMBER; i<NUMBER; i++) {
        
        part_1 = &part[i];
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
    }
    
    for (i=INIT_NUMBER; i<NUMBER; i++) {
        
        part_1 = &part[i];
        
        //noise  核小体粒子のみ
        if (simulate_type == 1 ) {
            
            p1 = sqrt(2.0 * 3.0 * rDNA_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
            p2 = sqrt(2.0 * 3.0 * rDNA_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
            theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
            psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
            
            
            part_1->force[X] = p1 * sin(theta) / sqrt(DELTA);
            part_1->force[Y] = p1 * cos(theta) / sqrt(DELTA);
            part_1->force[Z] = p2 * sin(psi) / sqrt(DELTA);
        }
    }
    
#pragma omp parallel for private ( j, k, m, gene_counter, p1, p2, theta, psi, dist, f, part_1, part_2, part_3, f_2, f_3) num_threads (8)
    for (i = INIT_NUMBER; i < NUMBER; i++){
        
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
                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] /*- 0.0*/);
                    part_1->force[Y] += f * (part_1->position[Y] /*- 0.0*/);
                    part_1->force[Z] += f * (part_1->position[Z] /*- 0.0*/);
                }
                
                //spb_exclude
                dist = Euclid_norm (part_1->position, spb.position);
                f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                
                if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
                    
                    part_1->force[X] += f * (part_1->position[X] - spb.position[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                }
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                
                if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
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
                
            case Centromere:    //セントロメア粒子への相互作用
                
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
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                
                if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                    part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                    part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                }
                
                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    part_1->force[X] += f * (part_1->position[X] );
                    part_1->force[Y] += f * (part_1->position[Y] );
                    part_1->force[Z] += f * (part_1->position[Z] );
                }
                
                break;
                
            case Telomere:      //テロメア粒子との相互作用
                
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
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - origin[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                                
                                part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                                part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                                part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                            }
                            
                        }
                        else { //telomere_3 核小体との結合
                            
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
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
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - origin[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius - PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist + PARTICLE_RADIUS > 4.0 * membrain_radius) {
                                
                                part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                                part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                                part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                            }
                        }
                        else { //telomere_3     核小体との結合
                            
                            dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
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
                
            case rDNA:      //核小体粒子の相互作用
                
                switch (i) {
                        
                        //
                    case 6468:
                    case 6193:
                        
                        if (i==6468) {
                            
                            part_2 = &part[5012];
                            dist = Euclid_norm (part_1->position, part_2->position);
                            
                            f_2 = K_BOND * ((PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                            part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                            part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        }
                        else{
                            
                            part_2 = &part[6192];
                            dist = Euclid_norm (part_1->position, part_2->position);
                            
                            f_2 = K_BOND * ((PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 - dist) / dist;
                            
                            part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                            part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                            part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        }
                        
                        /*
                         //spring2
                         if ( i == 6468 ) part_2 = &part[5013];
                         else  part_2 = &part[i-2];
                         
                         part_3 = &part[i+2];
                         
                         dist = Euclid_norm (part_1->position, part_2->position);
                         f_2 = K_BOND_2 * ( (rDNA_RADIUS + 3.0 * PARTICLE_RADIUS) * 0.8 - dist) / dist;
                         
                         dist = Euclid_norm (part_1->position, part_3->position);
                         f_3 = K_BOND_2 * ( 2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                         
                         part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                         + f_3 * (part_1->position[X] - part_3->position[X]);
                         part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                         + f_3 * (part_1->position[Y] - part_3->position[Y]);
                         part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                         + f_3 * (part_1->position[Z] - part_3->position[Z]);
                         */
                        
                        /*
                         //核小体領域表面を滑る
                         dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                         
                         f = K_BOND * ( 4.0 * membrain_radius - dist - rDNA_RADIUS ) / dist;
                         
                         part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                         part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                         part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                         */
                        
                        //membrain_exclude
                        dist = Euclid_norm (part_1->position, origin);
                        f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - rDNA_RADIUS ) / dist;
                        
                        if ( dist + rDNA_RADIUS > membrain_radius) {
                            
                            part_1->force[X] += f * (part_1->position[X] );
                            part_1->force[Y] += f * (part_1->position[Y] );
                            part_1->force[Z] += f * (part_1->position[Z] );
                        }
                        
                        break;
                        
                        ///////核膜結合側　　テロメア/////
                    case 6742:
                    case 6467:
                        
                        part_2 = &part[i-1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * ( 2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        //Nucleolus_exclude
                        dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                        f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + rDNA_RADIUS - dist ) / dist;
                        
                        if ( dist < 4.0 * membrain_radius + rDNA_RADIUS ) {
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                        }
                        
                        /*
                         //spring2
                         part_2 = &part[i-2];
                         
                         dist = Euclid_norm (part_1->position, part_2->position);
                         f = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                         
                         part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                         part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                         part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                         */
                        
                        /*
                         //telomere
                         dist = Euclid_norm (part_1->position, origin);
                         
                         f = K_BOND * ( membrain_radius - dist - rDNA_RADIUS ) / dist;
                         
                         part_1->force[X] += f * (part_1->position[X] - origin[X]);
                         part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                         part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                         */
                        
                        //nucleolus setting
                        dist = Euclid_norm (part_1->position, origin);
                        
                        f = K_BOND * ( nucleolus_setting_radius - dist - rDNA_RADIUS ) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - origin[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - origin[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - origin[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-1];
                        part_3 = &part[i+1];
                        
                        //spring
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND * (2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        dist = Euclid_norm(part_1->position, part_3->position);
                        f_3 = K_BOND * (2.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                        
                        part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        
                        /*
                         //ひもの硬さポテンシャル spring2
                         switch (i) {
                         
                         case 6469:
                         case 6194:
                         
                         if ( i == 6469 ) part_2 = &part[5012];
                         else  part_2 = &part[i-2];
                         
                         part_3 = &part[i+2];
                         
                         dist = Euclid_norm (part_1->position, part_2->position);
                         f_2 = K_BOND_2 * ( (3.0 * rDNA_RADIUS + PARTICLE_RADIUS) * 0.8 - dist) / dist;
                         
                         dist = Euclid_norm (part_1->position, part_3->position);
                         f_3 = K_BOND_2 * ( 4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                         
                         part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                         + f_3 * (part_1->position[X] - part_3->position[X]);
                         part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                         + f_3 * (part_1->position[Y] - part_3->position[Y]);
                         part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                         + f_3 * (part_1->position[Z] - part_3->position[Z]);
                         
                         break;
                         
                         case 6741:
                         case 6466:
                         
                         part_2 = &part[i-2];
                         
                         dist = Euclid_norm (part_1->position, part_2->position);
                         f = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                         
                         part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                         part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                         part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                         
                         break;
                         
                         default:
                         
                         part_2 = &part[i-2];
                         part_3 = &part[i+3];
                         
                         dist = Euclid_norm (part_1->position, part_2->position);
                         f_2 = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                         
                         dist = Euclid_norm (part_1->position, part_3->position);
                         f_3 = K_BOND_2 * (4.0 * rDNA_RADIUS * 0.8 - dist) / dist;
                         
                         part_1->force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                         + f_3 * (part_1->position[X] - part_3->position[X]);
                         part_1->force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                         + f_3 * (part_1->position[Y] - part_3->position[Y]);
                         part_1->force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                         + f_3 * (part_1->position[Z] - part_3->position[Z]);
                         
                         break;
                         }
                         */
                        
                        //Nucleolus_exclude
                        dist = Euclid_norm (part_1->position, Nucleolus_circle_center);
                        f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + rDNA_RADIUS - dist ) / dist;
                        
                        if ( dist < 4.0 * membrain_radius + rDNA_RADIUS ) {
                            
                            part_1->force[X] += f * (part_1->position[X] - Nucleolus_circle_center[X]);
                            part_1->force[Y] += f * (part_1->position[Y] - Nucleolus_circle_center[Y]);
                            part_1->force[Z] += f * (part_1->position[Z] - Nucleolus_circle_center[Z]);
                        }
                        
                        /*
                         //membrain_exclude
                         dist = Euclid_norm (part_1->position, origin);
                         f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - rDNA_RADIUS ) / dist;
                         
                         if ( dist + rDNA_RADIUS > membrain_radius) {
                         
                         part_1->force[X] += f * (part_1->position[X] );
                         part_1->force[Y] += f * (part_1->position[Y] );
                         part_1->force[Z] += f * (part_1->position[Z] );
                         }
                         */
                        
                        //nucleolus setting
                        dist = Euclid_norm (part_1->position, origin);
                        f = MEMBRAIN_EXCLUDE * (nucleolus_setting_radius - dist - rDNA_RADIUS ) / dist;
                        
                        if ( dist + rDNA_RADIUS > nucleolus_setting_radius) {
                            
                            part_1->force[X] += f * (part_1->position[X] );
                            part_1->force[Y] += f * (part_1->position[Y] );
                            part_1->force[Z] += f * (part_1->position[Z] );
                        }
                        
                        break;
                        
                }
                break;
        }
        

        //list
        if ( l%100 == 0) {
            
            m = 0;
            
            for(j=0; j<NUMBER; j++){
                
                part_2 = &part[j];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                
                if (dist < 6.0 * PARTICLE_RADIUS && abs(i-j) > 1){
                    
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
                
                if ( part_1->particle_type != rDNA && part_2->particle_type != rDNA) {
                    
                    if (dist < 2.0 * PARTICLE_RADIUS ){
                        
                        f = K_EXCLUDE * (2 * PARTICLE_RADIUS - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                    }
                }
                else if (part_1->particle_type == rDNA && part_2->particle_type == rDNA) {
                    
                    if (dist < 2.0 * rDNA_RADIUS){
                        
                        f = K_EXCLUDE * ( 2.0 * rDNA_RADIUS - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                    }
                }
                else {
                    
                    if (dist < (PARTICLE_RADIUS + rDNA_RADIUS)){
                        
                        f = K_EXCLUDE * ( (PARTICLE_RADIUS + rDNA_RADIUS) - dist) / dist;
                        
                        part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
                        part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                    }
                }
            }
        }
        
        //速度、位置の計算
        if ( part_1->particle_type != rDNA) {
            
            part_1->velocity[X] = (2.0 * PARTICLE_MASS * part_1->velocity_2[X] + DELTA * part_1->force[X]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
            part_1->velocity[Y] = (2.0 * PARTICLE_MASS * part_1->velocity_2[Y] + DELTA * part_1->force[Y]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
            part_1->velocity[Z] = (2.0 * PARTICLE_MASS * part_1->velocity_2[Z] + DELTA * part_1->force[Z]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
            
            part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - PARTICLE_MYU * part_1->velocity[X] ) / ( 2.0 * PARTICLE_MASS );
            part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / ( 2.0 * PARTICLE_MASS );
            part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / ( 2.0 * PARTICLE_MASS );
        }
        else {
            
            part_1->velocity[X] = (2.0 * rDNA_MASS * part_1->velocity_2[X] + DELTA * part_1->force[X]) / (2.0 * rDNA_MASS + rDNA_MYU * DELTA);
            part_1->velocity[Y] = (2.0 * rDNA_MASS * part_1->velocity_2[Y] + DELTA * part_1->force[Y]) / (2.0 * rDNA_MASS + rDNA_MYU * DELTA);
            part_1->velocity[Z] = (2.0 * rDNA_MASS * part_1->velocity_2[Z] + DELTA * part_1->force[Z]) / (2.0 * rDNA_MASS + rDNA_MYU * DELTA);
            
            part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( part_1->force[X] - rDNA_MYU * part_1->velocity[X] ) / ( 2.0 * rDNA_MASS );
            part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( part_1->force[Y] - rDNA_MYU * part_1->velocity[Y] ) / ( 2.0 * rDNA_MASS );
            part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( part_1->force[Z] - rDNA_MYU * part_1->velocity[Z] ) / ( 2.0 * rDNA_MASS );
        }
        
        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->velocity_2[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->velocity_2[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->velocity_2[Z];
        
    }
}

void renew () {
    
    int i;
    
    Particle *part_1;
    
    //核小体粒子のみ値を更新

    for(i = 6193; i < NUMBER; i++){    //値の更新
        
        part_1 = &part[i];
        
        part_1->position_old[X] = part_1->position[X];
        part_1->position_old[Y] = part_1->position[Y];
        part_1->position_old[Z] = part_1->position[Z];
        
        part_1->position[X] = part_1->position_new[X];
        part_1->position[Y] = part_1->position_new[Y];
        part_1->position[Z] = part_1->position_new[Z];
    }
    
    /*
    spb.position_old[X] = spb.position[X];
    spb.position_old[Y] = spb.position[Y];
    spb.position_old[Z] = spb.position[Z];
    
    spb.position[X] = spb.position_new[X];
    spb.position[Y] = spb.position_new[Y];
    spb.position[Z] = spb.position_new[Z];
    */
}

//核小体を作るための球の初期位置決定 第３染色体テロメアの中点方向
void Nucleolus_position_init (void) {
    
    double middle_point[DIMENSION], origin[] = { 0.0, 0.0, 0.0 };
    
    middle_point[X] = ( part[5012].position[X] + part[6192].position[X] ) / 2.0;
    middle_point[Y] = ( part[5012].position[Y] + part[6192].position[Y] ) / 2.0;
    middle_point[Z] = ( part[5012].position[Z] + part[6192].position[Z] ) / 2.0;
    
    Nucleolus_circle_center[X] = 5.0 * middle_point[X] * membrain_radius / Euclid_norm( middle_point, origin );
    Nucleolus_circle_center[Y] = 5.0 * middle_point[Y] * membrain_radius / Euclid_norm( middle_point, origin );
    Nucleolus_circle_center[Z] = 5.0 * middle_point[Z] * membrain_radius / Euclid_norm( middle_point, origin );
}

//体積が核の1/4になるまで核小体領域を広げる
void Nucleolus_position (void) {
    
    double V = 0, r, t, origin[] = {0.0, 0.0, 0.0}, per, h_1, h_2, l;
    
    r = membrain_radius;
    
    t = 5.0 * r - Euclid_norm ( Nucleolus_circle_center, origin );
    
    l = ( -10.0 * r * r + 10.0 * r * t - t * t )/ ( 2.0 * ( t - 5.0 * r ) );
    
    h_1 = r - l;    //核側の球欠の長さ
    h_2 = t + l - r;    //半径4rの球に含まれる球欠の長さ
    
    /*V = PI * ( r+h ) * ( 3.0*(r*r-h*h) + ( r+h )*( r+h )) / 6.0
     + PI * ( 4.0 * r-t-h ) * ( 3.0*(r*r-h*h) + ( 4.0*r-t-h )*( 4.0*r-t-h )) / 6.0;*/
    
    V = PI * h_1 * h_1 * ( 3.0 * r - h_1) / 3.0   //核に含まれる球欠の体積
    + PI * h_2 * h_2 * ( 12.0 * r - h_2) / 3.0;
    
    if ( V < PI * r * r * r / 3.0 ) {
        
        per = ( Euclid_norm ( Nucleolus_circle_center, origin) - 0.1 )/ Euclid_norm ( Nucleolus_circle_center, origin);
        
        Nucleolus_circle_center[X] *= per;
        Nucleolus_circle_center[Y] *= per;
        Nucleolus_circle_center[Z] *= per;
    }
    
    //printf ( "\n t = %lf , h_1 = %lf. h_2 = %lf \n", t, h_1, h_2);
}

void nucleolus_particle_set (void) {
    
    int i;
    double origin [] = {0.0, 0.0, 0.0};
    
    Particle *part_1, *part_2;
    
    for (i=6193; i<6468; i++) {
        
        part_1 = &part[i];
        part_2 = &part[i-1];
        
        part_1->particle_type = rDNA;
        
        if ( i != 6193) {
            
            part_1->position[X] = part_2->position[X] * ( Euclid_norm ( part_2->position, origin ) + 2.0 * rDNA_RADIUS * 0.8 ) / Euclid_norm ( part_2->position, origin);
            part_1->position[Y] = part_2->position[Y] * ( Euclid_norm ( part_2->position, origin ) + 2.0 * rDNA_RADIUS * 0.8 ) / Euclid_norm ( part_2->position, origin);
            part_1->position[Z] = part_2->position[Z] * ( Euclid_norm ( part_2->position, origin ) + 2.0 * rDNA_RADIUS * 0.8 ) / Euclid_norm ( part_2->position, origin);
        }
        else {
            
            part_1->position[X] = part_2->position[X] * ( Euclid_norm ( part_2->position, origin ) + (PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 ) / Euclid_norm ( part_2->position, origin);
            part_1->position[Y] = part_2->position[Y] * ( Euclid_norm ( part_2->position, origin ) + (PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 ) / Euclid_norm ( part_2->position, origin);
            part_1->position[Z] = part_2->position[Z] * ( Euclid_norm ( part_2->position, origin ) + (PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 ) / Euclid_norm ( part_2->position, origin);
        }
        
        part_1->chr_no = 2;
    }
    
    for (i=6468; i<6743; i++) {
        
        part_1 = &part[i];
        part_2 = &part[i-1];
        
        part_1->particle_type = rDNA;
        
        if ( i != 6468) {
            
            part_1->position[X] = part_2->position[X] * ( Euclid_norm ( part_2->position, origin ) + 2.0 * rDNA_RADIUS * 0.8 ) / Euclid_norm ( part_2->position, origin);
            part_1->position[Y] = part_2->position[Y] * ( Euclid_norm ( part_2->position, origin ) + 2.0 * rDNA_RADIUS * 0.8 ) / Euclid_norm ( part_2->position, origin);
            part_1->position[Z] = part_2->position[Z] * ( Euclid_norm ( part_2->position, origin ) + 2.0 * rDNA_RADIUS * 0.8 ) / Euclid_norm ( part_2->position, origin);
        }
        else {
            
            part_1->position[X] = part[5012].position[X] * ( Euclid_norm ( part[5012].position, origin ) + (PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 ) / Euclid_norm ( part[5012].position, origin);
            part_1->position[Y] = part[5012].position[Y] * ( Euclid_norm ( part[5012].position, origin ) + (PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 ) / Euclid_norm ( part[5012].position, origin);
            part_1->position[Z] = part[5012].position[Z] * ( Euclid_norm ( part[5012].position, origin ) + (PARTICLE_RADIUS + rDNA_RADIUS) * 0.8 ) / Euclid_norm ( part[5012].position, origin);
        }
        
        part_1->chr_no = 2;
    }
    
    //nucleolus_setting_radius = Euclid_norm ( part[6467].position, origin) / 2.0 + Euclid_norm (part[6742].position, origin) / 2.0;
    
    if (Euclid_norm ( part[6467].position, origin) < Euclid_norm ( part[6742].position, origin)) nucleolus_setting_radius = Euclid_norm ( part[6742].position, origin) + rDNA_RADIUS;
    else nucleolus_setting_radius = Euclid_norm ( part[6467].position, origin) + rDNA_RADIUS;
    
}

void write_coordinate (int t , int start) {
    
    int i;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "nucleolus_particle2/result_%d.txt", t + start);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n    error : can not write coordinate. \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++) {
        
        fprintf (fpw, "%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", i, part[i].chr_no, part[i].particle_type,
                 part[i].position_old[X], part[i].position_old[Y], part[i].position_old[Z], part[i].velocity[X], part[i].velocity[Y], part[i].velocity[Z], membrain_radius);
    }
    
    if (simulate_type == 0) printf(" t = %d, R_n = %lf\r", t, nucleolus_setting_radius);
    else printf(" t = %d, force = {%lf, %lf, %lf} \r", t, part[6741].force[X], part[6741].force[Y], part[6741].force[Z]);
    
    sprintf (str, "SPB");
    
    fprintf(fpw, "%s %s %lf %lf %lf %lf %lf %lf %lf\n", str, str, spb.position_old[X], spb.position_old[Y], spb.position_old[Z],
            spb.velocity[X], spb.velocity[Y], spb.velocity[Z],  membrain_radius);
    
    sprintf(str, "Nucleolus");
    
    fprintf(fpw, "%s %s %lf %lf %lf\n", str, str, Nucleolus_circle_center[X], Nucleolus_circle_center[Y], Nucleolus_circle_center[Z]);
    
    fclose (fpw);
}

int main ( int argc, char **argv ) {
    
    int i, t = 0, l;
    
    int start_number = atoi(argv[1]);
    int calculate_number = atoi(argv[2]);
    //k_expression = atof (argv[3]);
    int nucleolus_flag;
    
    //unsigned int *gene_list;
    
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
    
    init_particle( start_number );
    
    printf ("nucleolus_init  on:1 off:0\n");
    scanf ("%d", &nucleolus_flag);
    
    if (nucleolus_flag == 1) Nucleolus_position_init();
    
    //read_gene_list (gene_list);
    
    if (simulate_type == 0) nucleolus_particle_set();
    else nucleolus_setting_radius = membrain_radius;
    
    init_particle_calculate (dsfmt/*, gene_list*/);
    init_SPB_calculate(dsfmt);
    
    //初期位置の出力
    write_coordinate ( 0, start_number );
    
    for (t=1; t < calculate_number; t++) {
        
        //printf ("time = %d \r", t);
        
        Nucleolus_position();
        
        for (l=1; l<=10000; l++){
            
            particle_calculate(dsfmt, l);
            //SPB_calculate(dsfmt, l);
            
            renew ();
            
            if (nucleolus_setting_radius > membrain_radius) nucleolus_setting_radius -= 1.0e-4;
            
            //write_coordinate ( l , start_number);
        }
        write_coordinate ( t , start_number);
        
    }
    
    return ( 0 );
}


