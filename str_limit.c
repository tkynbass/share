//
//  str_limit.c
//
//  Hi-C + live imaging で構造推定
//  nucleosomeレベルでの慣性半径内にいるHMMの状態で存在確率が最高位のものをSPB側から順に選んでいく.
//  粒子はlocusデータ持ちのみ
//
//  Created by tkym on 2019/2/23.
//

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#include <omp.h>

#define DIMENSION ( 3 ) //次元
#define LENGTH ( 2.2e-6 / 75)   // 長さの単位
//#define M_A ( 1.85131596e+7 )
//#define N_A ( 6.022140857e+23 )

#define NUMBER ( 45 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )
#define RANK (7)    //HMMのランク数

#define K_BOND ( 1.0e-0 )    //ばね定数
#define K_BOND_2 ( 1.0e-1 )  //ひもの硬さ
#define K_BOND_3 ( 1.0e-1)
#define HMM_BOND (1.0e-0)

#define GYRATION_RADIUS (20 * 10e+3 / 196 * 1.4e-8 * 75 / 2.2e-6)   // GYRATION_RADIUS * 粒子数 = 慣性半径
#define PARTICLE_MYU ( 2.0 * DIMENSION * PI * PARTICLE_RADIUS * 0.000890) //粘性抵抗の強さ

#define DELTA ( 1.0e-7 )  //刻み幅

#define POTENTIAL_DELTA (1.0e-7)

unsigned int particle_number;

typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

const double nucleolus_pos[] = { 0.0, 0.0, 0.0};
const double spb_pos[] = { 1.7128e-6/LENGTH, 0.0, 0.0};
double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    //CHAIN chr_no;
    unsigned int pastis_no;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_init[DIMENSION];
    double velocity[DIMENSION];
    double nucleolus_mean[RANK];
    double spb_mean[RANK];
    double force[DIMENSION];
    
} Particle;

Particle *part;

enum label{ X, Y, Z};

void read_data ( char *cycle_dtatus, char *arm_id, unsigned int locus_list[45] ){       //初期値設定

    unsigned int loop, number = 0, locus_number = 0;
    
    char dummy[256], hmm_data[128];
    double d_dummy, enlarge_ratio;
    
    Particle *part_1;
    FILE *fpr;
    
    // Input the coordinates of Pastis //
    if ((fpr = fopen(sprintf("MDS_pos_%s.txt", arm_id), "r")) == NULL){
        
        printf ("\n\terror : cannot read coordinate.\n");
        
        exit (1);
    }
    
    fgets (dummy, 256, fpr);
    
    while (fscanf (fpr, "%d ", &part[number].pastis_no) != EOF) {
        
        part_1 = &part[number];
        fscanf (fpr, "%s %lf %lf %lf %lf %lf %lf %lf\n", &dummy,
                &part_1->position[X], &part_1->position[Y], &part_1->position[Z], &part_1->nucleolus_mean, &part_1->spb_mean, &d_dummy, &d_dummy);
        number++;
    }
    
    particle_number = number;
    
    fclose (fpr);
    
    // HMMの平均データ読み込み //
    char hmm_data[128];
    sprintf (hmm_data, "hmm_%s_%s.txt", cycle_dtatus, arm_id);
    
    if ((fpr = fopen(hmm_data, "r")) == NULL){
        
        printf ("\n\terror : cannot read coordinate.\n");
        
        exit (1);
    }

    unsigned int pastis_no;
    number = 0;
    while (fscanf (fpr, "%d ", &pastis_no) != EOF) {
        
        while (pastis_no != part[number].pastis_no) {
            
            number++;
        }
        
        part_1 = &part[number];
        
        fscanf (fpr, "%lf ", d_dummy);
        for (loop = 0; loop < RANK; loop++) {
            
            fscanf (fpr, "%lf %lf ", part_1->spb_mean[loop], part_1->nucleolus_mean[loop]);
        }
        fgets (dummy, 256, fpr);
        
        locus_list[locus_number] = number;
        locus_number ++;
    }
    
    loop = 0;
    while (locus_list != 0) {
        
        part_1 = &part[locus_list[loop]];
        printf ("\t%d status[0] %lf %lf, status[6] %lf %lf \n", part_1->pastis_no, part_1->spb_mean[0], part_1->nucleolus_mean[0],
                part_1->spb_mean[6], part_1->nucleolus_mean[6]);
        loop++;
    }
    
    /*
    for (i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        // 初期座標の保存
        part_1->position_init[X] = part_1->position[X];
        part_1->position_init[Y] = part_1->position[Y];
        part_1->position_init[Z] = part_1->position[Z];
        
        part_1->nucleolus_mean *= 1.0e-6/ LENGTH;
        part_1->spb_mean *= 1.0e-6 / LENGTH;
    
    }*/
    
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

//　ばねによる力 part_1 粒子側の力計算
void spring (Particle *part_1, const Particle *part_2, const unsigned int bond) {     //ばね
    
    double dist, dist_0;
    
    double f;
    
    //dist_0 = 自然長
    dist_0 = Euclid_norm (part_1->position_init, part_2->position_init);
    dist = Euclid_norm (part_1->position, part_2->position);
    
    f = bond * (dist_0 - dist) / dist;
    
    part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
    part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
    part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
}

void hmm_potential (Particle *part_1) {
    
    if (part_1->nucleolus_mean != 0.0) {
    
        double spb_dist = Euclid_norm (part_1->position, spb_pos);
        double nucleolus_dist = Euclid_norm (part_1->position, nucleolus_pos);
        
        double spb_f = HMM_BOND * (part_1->spb_mean - spb_dist);
        double nucleolus_f = HMM_BOND * (part_1->nucleolus_mean - nucleolus_dist);
        
        part_1->force[X] += spb_f * (part_1->position[X] - spb_pos[X]) + nucleolus_f * (part_1->position[X] - nucleolus_pos[X]);
        part_1->force[Y] += spb_f * (part_1->position[Y] - spb_pos[Y]) + nucleolus_f * (part_1->position[Y] - nucleolus_pos[Y]);
        part_1->force[Z] += spb_f * (part_1->position[Z] - spb_pos[Z]) + nucleolus_f * (part_1->position[Z] - nucleolus_pos[Z]);
        
    }
}

void calculate( unsigned int l ) {
    
    int i;
    
    Particle *part_1;
    
    for ( i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
    }
    
    
#pragma omp parallel for private (part_1) num_threads (6)
    for ( i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        // 隣同士 //
        if ( i != 0 && i != particle_number - 1 ) {
            
            spring (part_1, &part[i-1], K_BOND);
            spring (part_1, &part[i+1], K_BOND);
        }
        else if ( i == 0) spring (part_1, &part[i+1], K_BOND);
        else spring (part_1, &part[i-1], K_BOND);
        
        // 2個隣 //
        if ( 0 <= i-2 && i+2 <= particle_number-1) {
            
            spring (part_1, &part[i-2], K_BOND_2);
            spring (part_1, &part[i+2], K_BOND_2);
        }
        else if ( 0 <= i-2) spring (part_1, &part[i+2], K_BOND_2);
        else spring (part_1, &part[i-2], K_BOND_2);
        
        // 3個隣 //
        if ( 0 <= i-3 && i+3 <= particle_number-1) {
            
            spring (part_1, &part[i-3], K_BOND_3);
            spring (part_1, &part[i+3], K_BOND_3);
        }
        else if ( 0 <= i-3) spring (part_1, &part[i+3], K_BOND_3);
        else spring (part_1, &part[i-3], K_BOND_3);
        
        if (part_1->spb_mean != 0.0) {
            
            hmm_potential (part_1);
        }
        
        part_1->velocity[X] = part_1->force[X];
        part_1->velocity[Y] = part_1->force[Y];
        part_1->velocity[Z] = part_1->force[Z];
        
        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->velocity[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->velocity[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->velocity[Z];
    }
    
    // position の更新 //
    for ( i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        part_1->position[X] = part_1->position_new[X];
        part_1->position[Y] = part_1->position_new[Y];
        part_1->position[Z] = part_1->position_new[Z];
    }
    
}


void write_coordinate (int t) {
    
    int i;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "result_%d.txt", t);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    for (i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        fprintf (fpw, "%d %d %lf %lf %lf\n", i, part_1->pastis_no, part_1->position[X],
                 part_1->position[Y], part_1->position[Z]);
    }
    
    fclose (fpw);
}

int main ( int argc, char **argv ) {
    
    unsigned int loop, t = 0, l;
    unsigned int locus_list[45];
    char output_file[256];
    
    for (loop = 0; loop < 45; loop++) locus_list[loop] = 0;
     
    part = (Particle *)malloc(NUMBER * sizeof(Particle));
    
    if (part == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    read_data (argv[1], argv[2], locus_list);
    
    /*
    //初期位置の出力//
    write_coordinate (0);
    
    for (t=1; t <= calculate_number; t++) {
        
        for (l=1; l<=1.0e+5; l++){
            
            calculate(l);
        }
        
        printf ("\tt = %d \r", t);
        fflush (stdout);
        
        write_coordinate (t);
    }
    */
    
    return ( 0 );
}



