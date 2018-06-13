



//
//  nuc_cal.c
//  
//  SPBノイズあり 核小体を楕円体型に。それに伴った核小体相互作用の変更 →　中心を動かす
//  Created by 高山雄揮 on 2017/09/27.
//
//

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "dSFMT.h"
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
#define SPHERE_DIVISION (20)     //膜の分割数
#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PI (M_PI)
#define MEMBRAIN_RADIUS_MINIMUM  (32.0) //膜の半径の最小値
#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )

#define K_EXCLUDE (1.0)    //排除体積効果の強さ
#define K_BOND (1.0)    //ばね定数
#define K_BOND_2 (1.0e-4)  //ひもの硬さ
#define DELTA (1.0e-11)  //刻み幅
#define PARTICLE_MYU (2.0 * DIMENSION * PI * PARTICLE_RADIUS * NANO * 0.000890 / 100 )    //粘性抵抗の強さ
#define MEMBRAIN_EXCLUDE (1.0)     //膜との衝突
#define MEMBRAIN_EXCLUDE_SPB (1.0) //SPBとの衝突

#define SPB_RADIUS ( 3.0 )      //SPBの半径
#define SPB_MYU (2.0 * DIMENSION * PI * SPB_RADIUS * NANO * 0.000890 / 100 )  //SPBの粘性

double membrain_radius, nucleolus_mass;

//#define CLUSTER_GENE_NUMBER (146)

double Nucleolus_circle_center[3], k_expression, nucleolus_myu;
int centro_no = 5542;

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

Particle *ptr;

Particle part[NUMBER];

Particle spb;

Particle nucleolus;

enum label{ X, Y, Z};

void init_particle( int start ){       //初期値設定
    
    int i, i_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    sprintf (filename, "nuc_move_%d.txt", start);
    
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
    
    fscanf (fpr, "%s %s %lf %lf %lf %lf %lf %lf", dummy, dummy, &nucleolus.position[X], &nucleolus.position[Y], &nucleolus.position[Z], &nucleolus.velocity[X], &nucleolus.velocity[Y], &nucleolus.velocity[Z]);
    
    fclose (fpr);
    
}

/*
void read_gene_list (unsigned int gene_list[CLUSTER_GENE_NUMBER]) {     //発現量変化遺伝子の読み取り
    
    FILE *fpr;
    char filename[128];
    
    int i;
    
    sprintf (filename, "cl4_num.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf (" \n cannot read gene list \n ");
        
        exit (1);
    }
    
    for (i=0; i<CLUSTER_GENE_NUMBER; i++) {
        
        fscanf (fpr, "%d\n", &gene_list[i]);
    }
    
    fclose(fpr);
}*/

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

void nucleolus_myu_cal(void) {
    
    double nucleus_volume = 0.0, nucleolus_volume = 0.0, r, t, origin[] = {0.0, 0.0, 0.0}, per, h_1, h_2, l;
    
    r = membrain_radius;
    
    t = 5.0 * r - Euclid_norm ( nucleolus.position, origin );
    
    l = ( -10.0 * r * r + 10.0 * r * t - t * t )/ ( 2.0 * ( t - 5.0 * r ) );
    
    h_1 = r - l;    //核側の球欠の長さ
    h_2 = t + l - r;    //半径4rの球に含まれる球欠の長さ
    
    nucleus_volume = 4.0 * PI * membrain_radius * membrain_radius * membrain_radius / 3.0;
    
    nucleolus_volume = PI * h_1 * h_1 * ( 3.0 * r - h_1) / 3.0   //核に含まれる球欠の体積
    + PI * h_2 * h_2 * ( 12.0 * r - h_2) / 3.0;
    
    nucleolus_myu = 2.0 * DIMENSION * ( membrain_radius * cbrt ( nucleolus_volume / nucleus_volume )) * PI * NANO * 0.000890 / 100;
    
}

void spring (const Particle *part_1, const Particle *part_2, double force[DIMENSION]) {     //ばね
    
    double dist, dist_0;
    
    double f;
    
    Particle *part_3;
    
    part_3 = &spb;
    
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
    
    double force[3] = { 0.0, 0.0, 0.0}, origin[] = { 0.0, 0.0, 0.0};
    
    double dist = Euclid_norm ( spb.position , origin);
    double p1, p2, theta, psi;
    double f = MEMBRAIN_EXCLUDE * (membrain_radius - dist ) / dist;
    
    Particle *part_2;
    
    p1 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    
    force[X] += f * (spb.position[X]);        //膜とのバネ
    force[Y] += f * (spb.position[Y]);
    force[Z] += f * (spb.position[Z]);
    
    spring (&spb, &part[1880], force);       //セントロメアとのバネによる力
    spring (&spb, &part[3561], force);
    spring (&spb, &part[5542], force);
    
    spb_list (&spb);
    
    //ひも粒子との排除体積
    if (spb.list_no != 0){
        
        for (j = 1; j <= spb.list_no; j++){
            
            k = spb.list[j];
            
            part_2 = &part[k];
            dist = Euclid_norm (spb.position, part_2->position);
            
            if (dist < PARTICLE_RADIUS + SPB_RADIUS){
                
                f = K_EXCLUDE * (PARTICLE_RADIUS + SPB_RADIUS - dist) / dist;
                
                force[X] += f * (spb.position[X] - part_2->position[X]);
                force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                force[Z] += f * (spb.position[Z] - part_2->position[Z]);
            }
        }
    }
    
    //粘性抵抗
    force[X] += - SPB_MYU * spb.velocity[X];
    force[Y] += - SPB_MYU * spb.velocity[Y];
    force[Z] += - SPB_MYU * spb.velocity[Z];
    
    spb.velocity_2[X] = spb.velocity[X] + DELTA * ( force[X] - SPB_MYU * spb.velocity[X] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Y] = spb.velocity[Y] + DELTA * ( force[Y] - SPB_MYU * spb.velocity[Y] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Z] = spb.velocity[Z] + DELTA * ( force[Z] - SPB_MYU * spb.velocity[Z] ) / ( 2.0 * SPB_MASS );
    
    spb.position_old[X] = spb.position[X];
    spb.position_old[Y] = spb.position[Y];
    spb.position_old[Z] = spb.position[Z];
    
    spb.position[X] += DELTA * spb.velocity_2[X];
    spb.position[Y] += DELTA * spb.velocity_2[Y];
    spb.position[Z] += DELTA * spb.velocity_2[Z];

    
}

void SPB_calculate (dsfmt_t dsfmt, const unsigned int l){
    
    int k, j;
    
    double force[3] = { 0.0, 0.0, 0.0}, origin[] = { 0.0, 0.0, 0.0};
    
    double dist = Euclid_norm ( spb.position , origin);
    double p1, p2, theta, psi;
    double f = MEMBRAIN_EXCLUDE * (membrain_radius - dist ) / dist;
    
    Particle *part_2;
    
    p1 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * SPB_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    
    force[X] += f * (spb.position[X]);        //膜とのバネ
    force[Y] += f * (spb.position[Y]);
    force[Z] += f * (spb.position[Z]);
    
    spring (&spb, &part[1880], force);       //セントロメアとのバネによる力
    spring (&spb, &part[3561], force);
    spring (&spb, &part[5542], force);
    
    if ( l%1000 == 0) spb_list (&spb);
    
    //ひも粒子との排除体積
    if (spb.list_no != 0){
        
        for (j = 1; j <= spb.list_no; j++){
            
            k = spb.list[j];
            
            part_2 = &part[k];
            dist = Euclid_norm (spb.position, part_2->position);
            
            if (dist < PARTICLE_RADIUS + SPB_RADIUS){
                
                f = K_EXCLUDE * (PARTICLE_RADIUS + SPB_RADIUS - dist) / dist;
                
                force[X] += f * (spb.position[X] - part_2->position[X]);
                force[Y] += f * (spb.position[Y] - part_2->position[Y]);
                force[Z] += f * (spb.position[Z] - part_2->position[Z]);
            }
        }
    }
    
    
    spb.velocity[X] = (2.0 * SPB_MASS * spb.velocity_2[X] + DELTA * force[X]) / (2.0 * SPB_MASS + SPB_MYU * DELTA);
    spb.velocity[Y] = (2.0 * SPB_MASS * spb.velocity_2[Y] + DELTA * force[Y]) / (2.0 * SPB_MASS + SPB_MYU * DELTA);
    spb.velocity[Z] = (2.0 * SPB_MASS * spb.velocity_2[Z] + DELTA * force[Z]) / (2.0 * SPB_MASS + SPB_MYU * DELTA);
    
    spb.velocity_2[X] = spb.velocity[X] + DELTA * ( force[X] - SPB_MYU * spb.velocity[X] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Y] = spb.velocity[Y] + DELTA * ( force[Y] - SPB_MYU * spb.velocity[Y] ) / ( 2.0 * SPB_MASS );
    spb.velocity_2[Z] = spb.velocity[Z] + DELTA * ( force[Z] - SPB_MYU * spb.velocity[Z] ) / ( 2.0 * SPB_MASS );
    
    spb.position_new[X] = spb.position[X] + DELTA * spb.velocity_2[X];
    spb.position_new[Y] = spb.position[Y] + DELTA * spb.velocity_2[Y];
    spb.position_new[Z] = spb.position[Z] + DELTA * spb.velocity_2[Z];
    
    
}

void init_particle_calculate( dsfmt_t dsfmt/*, const unsigned int gene_list [CLUSTER_GENE_NUMBER] */){
    
    int i, k, j, m, gene_counter=0;
    
    double force[DIMENSION], nucleolus_force[DIMENSION], f, f_2, f_3, dist, origin[] = {0.0, 0.0, 0.0};
    double p1, p2, theta, psi;
    
    Particle *part_1, *part_2, *part_3;
    
    //noise nucleolus
    
    nucleolus_myu_cal();
    
    p1 = sqrt(2.0 * 3.0 * nucleolus_myu * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * nucleolus_myu * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    
    nucleolus_force[X] = p1 * sin(theta) / sqrt(DELTA);
    nucleolus_force[Y] = p1 * cos(theta) / sqrt(DELTA);
    nucleolus_force[Z] = p2 * sin(psi) / sqrt(DELTA);
    
    nucleolus_force[X] += - nucleolus_myu * nucleolus.velocity[X];
    nucleolus_force[Y] += - nucleolus_myu * nucleolus.velocity[Y];
    nucleolus_force[Z] += - nucleolus_myu * nucleolus.velocity[Z];
    
    for (i = 0; i < NUMBER; i++){
        
        part_1 = &part[i];
        
        //noise
        p1 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
        p2 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
        theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
        psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
        
        force[X] = p1 * sin(theta) / sqrt(DELTA);
        force[Y] = p1 * cos(theta) / sqrt(DELTA);
        force[Z] = p2 * sin(psi) / sqrt(DELTA);
        
        force[X] += - PARTICLE_MYU * part_1->velocity[X];
        force[Y] += - PARTICLE_MYU * part_1->velocity[Y];
        force[Z] += - PARTICLE_MYU * part_1->velocity[Z];

        switch (part[i].particle_type) {
            case Normal:
                
                part_2 = &part[i-1];
                part_3 = &part[i+1];
                
                
                //spring
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    force[X] += f * (part_1->position[X] /*- 0.0*/);
                    force[Y] += f * (part_1->position[Y] /*- 0.0*/);
                    force[Z] += f * (part_1->position[Z] /*- 0.0*/);
                }
                
                //spb_exclude
                dist = Euclid_norm (part_1->position, spb.position);
                f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                
                if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
                    
                    force[X] += f * (part_1->position[X] - spb.position[X]);
                    force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                    force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                }
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, nucleolus.position);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
            
                if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                    
                    force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                    force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                    force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                    
                    nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                    nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                    nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                }

                //spring2
                switch (i) {
                    case 1:
                    case 2772:
                    case 5013:
                        
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    case 2769:
                    case 5010:
                    case 6191:
                        
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-2];
                        part_3 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        dist = Euclid_norm (part_1->position, part_3->position);
                        f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        
                        break;
                }
                
                break;
                
            case Centromere:
                
                //centromere
                dist = Euclid_norm( part_1->position , spb.position);
                f = K_BOND * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                force[X] += f * (part_1->position[X] - spb.position[X]);
                force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                
                //spring
                part_2 = &part[i-1];
                part_3 = &part[i+1];
                
                dist = Euclid_norm(part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //spring2
                part_2 = &part[i-2];
                part_3 = &part[i+2];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm (part_1->position, part_3->position);
                f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, nucleolus.position);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
            
                if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                    
                    force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                    force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                    force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                    
                    nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                    nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                    nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                }

                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    force[X] += f * (part_1->position[X] );
                    force[Y] += f * (part_1->position[Y] );
                    force[Z] += f * (part_1->position[Z] );
                }
                
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
                        
                        force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 5012) { //telomere
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            force[X] += f * (part_1->position[X] - origin[X]);
                            force[Y] += f * (part_1->position[Y] - origin[Y]);
                            force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                                
                                force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                                force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                                force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                                
                                nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                                nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                                nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                            }

                        }
                        else { //telomere_3 核小体との結合
                            
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                            force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                            force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                        }
                        
                        //spring2
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 6192) {//telomere
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            force[X] += f * (part_1->position[X] - origin[X]);
                            force[Y] += f * (part_1->position[Y] - origin[Y]);
                            force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                                
                                force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                                force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                                force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                                
                                nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                                nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                                nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                            }
                        }
                        else { //telomere_3     核小体との結合
                            
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                            force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                            force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                        }
                        
                        //spring2
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
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
                    
                    force[X] += f * (part_1->position[X] - part_2->position[X]);
                    force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                    force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                }
                
            }
        }
        
        part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( force[X] - PARTICLE_MYU * part_1->velocity[X] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / ( 2.0 * PARTICLE_MASS );
        
    }
    
    nucleolus.velocity_2[X] = nucleolus.velocity[X] + DELTA * ( nucleolus_force[X] - nucleolus_myu * nucleolus.velocity[X] ) / ( 2.0 * nucleolus_mass );
    nucleolus.velocity_2[Y] = nucleolus.velocity[Y] + DELTA * ( nucleolus_force[Y] - nucleolus_myu * nucleolus.velocity[Y] ) / ( 2.0 * nucleolus_mass );
    nucleolus.velocity_2[Z] = nucleolus.velocity[Z] + DELTA * ( nucleolus_force[Z] - nucleolus_myu * nucleolus.velocity[Z] ) / ( 2.0 * nucleolus_mass );
    
    for ( i=0 ; i<NUMBER; i++) {
        
        part_1 = &part[i];
        
        part_1->position_old[X] = part_1->position[X];
        part_1->position_old[Y] = part_1->position[Y];
        part_1->position_old[Z] = part_1->position[Z];
        
        part_1->position[X] += DELTA * part_1->velocity_2[X];
        part_1->position[Y] += DELTA * part_1->velocity_2[Y];
        part_1->position[Z] += DELTA * part_1->velocity_2[Z];
    }
    
    nucleolus.position_old[X] = nucleolus.position[X];
    nucleolus.position_old[Y] = nucleolus.position[Y];
    nucleolus.position_old[Z] = nucleolus.position[Z];
    
    nucleolus.position[X] += DELTA * nucleolus.velocity_2[X];
    nucleolus.position[Y] += DELTA * nucleolus.velocity_2[Y];
    nucleolus.position[Z] += DELTA * nucleolus.velocity_2[Z];
}

void particle_calculate( dsfmt_t dsfmt, const unsigned int l/*, const unsigned int gene_list [CLUSTER_GENE_NUMBER]*/)        //位置と速度の計算 private force dist f part_1 part_2 part_3
{
    int i, k, j, m, gene_counter = 0;
    
    double force[DIMENSION], nucleolus_force[DIMENSION], f, f_2, f_3, dist, origin[] = {0.0, 0.0, 0.0};
    double p1, p2, theta, psi;
    
    Particle *part_1, *part_2, *part_3;
    
    nucleolus_myu_cal();
    
    p1 = sqrt(2.0 * 3.0 * nucleolus_myu * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    p2 = sqrt(2.0 * 3.0 * nucleolus_myu * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
    theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
    
    nucleolus_force[X] = p1 * sin(theta) / sqrt(DELTA);
    nucleolus_force[Y] = p1 * cos(theta) / sqrt(DELTA);
    nucleolus_force[Z] = p2 * sin(psi) / sqrt(DELTA);
    
    #pragma omp parallel for private ( j, k, m, gene_counter, p1, p2, theta, psi, force, dist, f, part_1, part_2, part_3, f_2, f_3) num_threads (8)
    for (i = 0; i < NUMBER; i++){
        
        part_1 = &part[i];
        
        //noise
        p1 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
        p2 = sqrt(2.0 * 3.0 * PARTICLE_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(&dsfmt) ));
        theta = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
        psi = 2.0 * PI * dsfmt_genrand_open_close(&dsfmt);
        
        force[X] = p1 * sin(theta) / sqrt(DELTA);
        force[Y] = p1 * cos(theta) / sqrt(DELTA);
        force[Z] = p2 * sin(psi) / sqrt(DELTA);
        
        
        
        /*
        if (i == gene_list [gene_counter]) {    //発現量が上がる遺伝子に核中心方向の力を加える
            
            high_expression (part_1, force);
            
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

                force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                            + f_3 * (part_1->position[X] - part_3->position[X]);
                force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                            + f_3 * (part_1->position[Y] - part_3->position[Y]);
                force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                            + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    force[X] += f * (part_1->position[X] /*- 0.0*/);
                    force[Y] += f * (part_1->position[Y] /*- 0.0*/);
                    force[Z] += f * (part_1->position[Z] /*- 0.0*/);
                }
                
                //spb_exclude
                dist = Euclid_norm (part_1->position, spb.position);
                f = K_EXCLUDE * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                
                if ( dist < PARTICLE_RADIUS + SPB_RADIUS) {
                    
                    force[X] += f * (part_1->position[X] - spb.position[X]);
                    force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                    force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                }
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, nucleolus.position);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
            
                if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                    
                    force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                    force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                    force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                    
                    nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                    nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                    nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                }
            
                //spring2
                switch (i) {
                    case 1:
                    case 2772:
                    case 5013:
                        
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                    
                    case 2769:
                    case 5010:
                    case 6191:
                        
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-2];
                        part_3 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        dist = Euclid_norm (part_1->position, part_3->position);
                        f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                        + f_3 * (part_1->position[X] - part_3->position[X]);
                        force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                        + f_3 * (part_1->position[Y] - part_3->position[Y]);
                        force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                        + f_3 * (part_1->position[Z] - part_3->position[Z]);
                        
                        break;
                }
                
                break;
            
            case Centromere:
                
                //centromere
                dist = Euclid_norm( part_1->position , spb.position);
                f = K_BOND * ( SPB_RADIUS + PARTICLE_RADIUS - dist) / dist;
                force[X] += f * (part_1->position[X] - spb.position[X]);
                force[Y] += f * (part_1->position[Y] - spb.position[Y]);
                force[Z] += f * (part_1->position[Z] - spb.position[Z]);
                
                //spring
                part_2 = &part[i-1];
                part_3 = &part[i+1];
                
                dist = Euclid_norm(part_1->position, part_2->position);
                f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm(part_1->position, part_3->position);
                f_3 = K_BOND * (INIT_DISTANCE - dist) / dist;
                
                force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //spring2
                part_2 = &part[i-2];
                part_3 = &part[i+2];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                f_2 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                dist = Euclid_norm (part_1->position, part_3->position);
                f_3 = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                
                force[X] += f_2 * (part_1->position[X] - part_2->position[X])
                + f_3 * (part_1->position[X] - part_3->position[X]);
                force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y])
                + f_3 * (part_1->position[Y] - part_3->position[Y]);
                force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z])
                + f_3 * (part_1->position[Z] - part_3->position[Z]);
                
                //Nucleolus_exclude
                dist = Euclid_norm (part_1->position, nucleolus.position);
                f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
                
                if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                    
                    force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                    force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                    force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                    
                    nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                    nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                    nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                }
                    
                //membrain_exclude
                dist = Euclid_norm (part_1->position, origin);
                f = MEMBRAIN_EXCLUDE * (membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                
                if ( dist + PARTICLE_RADIUS > membrain_radius) {
                    
                    force[X] += f * (part_1->position[X] );
                    force[Y] += f * (part_1->position[Y] );
                    force[Z] += f * (part_1->position[Z] );
                }
                
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
                        
                        force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 5012) { //telomere
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            force[X] += f * (part_1->position[X] - origin[X]);
                            force[Y] += f * (part_1->position[Y] - origin[Y]);
                            force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                                
                                force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                                force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                                force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                                
                                nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                                nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                                nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                            }
                        }
                        else { //telomere_3     核小体との結合
                            
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                            force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                            force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                        }
                        
                        //spring2
                        part_2 = &part[i+2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                        
                        break;
                        
                    default:
                        
                        part_2 = &part[i-1];
                        
                        //spring
                        dist = Euclid_norm(part_1->position, part_2->position);
                        f_2 = K_BOND * (INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f_2 * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f_2 * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f_2 * (part_1->position[Z] - part_2->position[Z]);
                        
                        if (i != 6192) {//telomere
                            
                            dist = Euclid_norm (part_1->position, origin);
                            
                            f = K_BOND * ( membrain_radius - dist - PARTICLE_RADIUS ) / dist;
                            
                            force[X] += f * (part_1->position[X] - origin[X]);
                            force[Y] += f * (part_1->position[Y] - origin[Y]);
                            force[Z] += f * (part_1->position[Z] - origin[Z]);
                            
                            //Nucleolus_exclude
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            f = MEMBRAIN_EXCLUDE * ( 4.0 * membrain_radius + PARTICLE_RADIUS - dist ) / dist;
                            
                            if ( dist < 4.0 * membrain_radius + PARTICLE_RADIUS) {
                                
                                force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                                force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                                force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                                
                                nucleolus_force[X] += f * (nucleolus.position[X] - part_1->position[X]);
                                nucleolus_force[Y] += f * (nucleolus.position[Y] - part_1->position[Y]);
                                nucleolus_force[Z] += f * (nucleolus.position[Z] - part_1->position[Z]);
                            }
                        }
                        else { //telomere_3     核小体との結合
                            
                            dist = Euclid_norm (part_1->position, nucleolus.position);
                            
                            f = K_BOND * ( 4.0 * membrain_radius + PARTICLE_RADIUS  - dist ) / dist;
                            
                            force[X] += f * (part_1->position[X] - nucleolus.position[X]);
                            force[Y] += f * (part_1->position[Y] - nucleolus.position[Y]);
                            force[Z] += f * (part_1->position[Z] - nucleolus.position[Z]);
                        }
                        
                        //spring2
                        part_2 = &part[i-2];
                        
                        dist = Euclid_norm (part_1->position, part_2->position);
                        f = K_BOND_2 * (2.0 * INIT_DISTANCE - dist) / dist;
                        
                        force[X] += f * (part_1->position[X] - part_2->position[X]);
                        force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                        force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
        
                        break;
                }
                break;
        }
        
        //list
        if ( l%2000 == 0) {
            
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
                    
                    force[X] += f * (part_1->position[X] - part_2->position[X]);
                    force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
                    force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
                }
                
            }
        }

        
        part_1->velocity[X] = (2.0 * PARTICLE_MASS * part_1->velocity_2[X] + DELTA * force[X]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
        part_1->velocity[Y] = (2.0 * PARTICLE_MASS * part_1->velocity_2[Y] + DELTA * force[Y]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
        part_1->velocity[Z] = (2.0 * PARTICLE_MASS * part_1->velocity_2[Z] + DELTA * force[Z]) / (2.0 * PARTICLE_MASS + PARTICLE_MYU * DELTA);
        
        part_1->velocity_2[X] = part_1->velocity[X] + DELTA * ( force[X] - PARTICLE_MYU * part_1->velocity[X] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * ( force[Y] - PARTICLE_MYU * part_1->velocity[Y] ) / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * ( force[Z] - PARTICLE_MYU * part_1->velocity[Z] ) / ( 2.0 * PARTICLE_MASS );
        
        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->velocity_2[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->velocity_2[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->velocity_2[Z];
    
        
    }
    
    //printf ("force[X] = %lf, force[Y] = %lf, force[Z] = %lf \n", nucleolus_force[X], nucleolus_force[Y], nucleolus_force[Z]);
    
    nucleolus.velocity[X] = (2.0 * nucleolus_mass * nucleolus.velocity_2[X] + DELTA * force[X]) / (2.0 * nucleolus_mass + nucleolus_myu * DELTA);
    nucleolus.velocity[Y] = (2.0 * nucleolus_mass * nucleolus.velocity_2[Y] + DELTA * force[Y]) / (2.0 * nucleolus_mass + nucleolus_myu * DELTA);
    nucleolus.velocity[Z] = (2.0 * nucleolus_mass * nucleolus.velocity_2[Z] + DELTA * force[Z]) / (2.0 * nucleolus_mass + nucleolus_myu * DELTA);
    
    nucleolus.velocity_2[X] = nucleolus.velocity[X] + DELTA * ( force[X] - nucleolus_myu * nucleolus.velocity[X] ) / ( 2.0 * nucleolus_mass );
    nucleolus.velocity_2[Y] = nucleolus.velocity[Y] + DELTA * ( force[Y] - nucleolus_myu * nucleolus.velocity[Y] ) / ( 2.0 * nucleolus_mass );
    nucleolus.velocity_2[Z] = nucleolus.velocity[Z] + DELTA * ( force[X] - nucleolus_myu * nucleolus.velocity[Z] ) / ( 2.0 * nucleolus_mass );
    
    nucleolus.position_new[X] = nucleolus.position[X] + DELTA * nucleolus.velocity_2[X];
    nucleolus.position_new[Y] = nucleolus.position[Y] + DELTA * nucleolus.velocity_2[Y];
    nucleolus.position_new[Z] = nucleolus.position[Z] + DELTA * nucleolus.velocity_2[Z];
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
    
    for (i = 0; i<DIMENSION; i++) {
        
        nucleolus.position_old[i] = nucleolus.position[i];
        
        nucleolus.position[i] = nucleolus.position_new[i];
    }
    
}

void Nucleolus_position_init (void) {
    
    double middle_point[DIMENSION], origin[] = { 0.0, 0.0, 0.0 };
    
    middle_point[X] = ( part[5012].position[X] + part[6192].position[X] ) / 2.0;
    middle_point[Y] = ( part[5012].position[Y] + part[6192].position[Y] ) / 2.0;
    middle_point[Z] = ( part[5012].position[Z] + part[6192].position[Z] ) / 2.0;
    
    nucleolus.position[X] = 5.0 * middle_point[X] * membrain_radius / Euclid_norm( middle_point, origin );
    nucleolus.position[Y] = 5.0 * middle_point[Y] * membrain_radius / Euclid_norm( middle_point, origin );
    nucleolus.position[Z] = 5.0 * middle_point[Z] * membrain_radius / Euclid_norm( middle_point, origin );
}

void Nucleolus_position (void) {
    
    double V = 0, r, t, origin[] = {0.0, 0.0, 0.0}, per, h_1, h_2, l;
    
    r = membrain_radius;
    
    t = 5.0 * r - Euclid_norm ( nucleolus.position, origin );
    
    l = ( -10.0 * r * r + 10.0 * r * t - t * t )/ ( 2.0 * ( t - 5.0 * r ) );
    
    h_1 = r - l;    //核側の球欠の長さ
    h_2 = t + l - r;    //半径4rの球に含まれる球欠の長さ
    
    /*V = PI * ( r+h ) * ( 3.0*(r*r-h*h) + ( r+h )*( r+h )) / 6.0
    + PI * ( 4.0 * r-t-h ) * ( 3.0*(r*r-h*h) + ( 4.0*r-t-h )*( 4.0*r-t-h )) / 6.0;*/
    
    V = PI * h_1 * h_1 * ( 3.0 * r - h_1) / 3.0   //核に含まれる球欠の体積
    + PI * h_2 * h_2 * ( 12.0 * r - h_2) / 3.0;
    
    if ( V < PI * r * r * r / 3.0 ) {
        
        per = ( Euclid_norm ( nucleolus.position, origin) - 0.1 )/ Euclid_norm ( nucleolus.position, origin);
        
        nucleolus.position[X] *= per;
        nucleolus.position[Y] *= per;
        nucleolus.position[Z] *= per;
    }
    
    //printf ( "\n t = %lf , h_1 = %lf. h_2 = %lf \n", t, h_1, h_2);
}

void write_coordinate ( int t , int start) {
    
    int i;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "nuc_move_%d.txt", t + start);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++) {
        
        fprintf (fpw, "%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", i, part[i].chr_no, part[i].particle_type,
                 part[i].position_old[X], part[i].position_old[Y], part[i].position_old[Z], part[i].velocity[X], part[i].velocity[Y], part[i].velocity[Z], membrain_radius);
        
        printf(" t = %d, i = %d \r", t, i);
        
    }
    
    sprintf (str, "SPB");
    
    fprintf(fpw, "%s %s %lf %lf %lf %lf %lf %lf %lf\n", str, str, spb.position_old[X], spb.position_old[Y], spb.position_old[Z],
            spb.velocity[X], spb.velocity[Y], spb.velocity[Z],  membrain_radius);
    
    sprintf(str, "Nucleolus");
    
    fprintf(fpw, "%s %s %lf %lf %lf %lf %lf %lf\n", str, str, nucleolus.position_old[X], nucleolus.position_old[Y], nucleolus.position_old[Z], nucleolus.velocity[X], nucleolus.velocity[Y], nucleolus.velocity[Z]);
    
    fclose (fpw);
}

int main ( int argc, char **argv ) {
    
    int i, t = 0, l;
    
    int start_number = atoi(argv[1]);
    int calculate_number = atoi(argv[2]);
    //k_expression = atof (argv[3]);
    int nucleolus_flag;
    
    //unsigned int *gene_list;
    
    ptr = (Particle *)malloc( (NUMBER +2 ) * sizeof(Particle));
    
    if (ptr == NULL) {
        
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
    
    Init_particle( start_number );
    
    nucleolus_mass = PI * membrain_radius * membrain_radius * membrain_radius / 3.0;        //核小体質量（体積）
    nucleolus_myu = 2.0 * DIMENSION * PI * PARTICLE_RADIUS * NANO * 0.000890 / 100;
    
    printf ("nucleolus_init  on:1 off:0\n");
    scanf ("%d", &nucleolus_flag);
    
    if (nucleolus_flag == 1) Nucleolus_position_init();
    
    //read_gene_list (gene_list);
    
    init_particle_calculate (dsfmt/*, gene_list*/);
    init_SPB_calculate(dsfmt);
    
    //初期位置の出力
    write_coordinate (0, start_number );
    
    for (t=1; t < calculate_number; t++) {
        
        //printf ("time = %d \r", t);
        
        //Nucleolus_position();
        
        for (l=1; l<=10000; l++){
            
            particle_calculate(dsfmt, l/*, gene_list*/);
            SPB_calculate(dsfmt, l);
            
            renew ();
        }
        
        write_coordinate ( t , start_number);
    }
    
    return ( 0 );
}
