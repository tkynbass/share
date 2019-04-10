//
//  str_limit.c
//
//  Hi-C + live imaging で構造推定
//  nucleosomeレベルでの慣性半径内にいるHMMの状態で存在確率が最高位のものをSPB側から順に選んでいく.
//  粒子はlocusデータ持ちのみ

// hmm_data の　0 →　0.0
// 粒子数分のpartの確保
// malloc で part->spb_mean, nucleolus_meanのメモリ確保
// malloc を関数化　→　ポインタのポインタを使う (sub_memory は アロー演算子で渡してるからポインタ引数で大丈夫？)
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

#define NUMBER_MAX ( 184 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )
#define RANK (8)    //HMMのランク数

#define K_BOND ( 1.0e-0 )    //ばね定数
#define K_BOND_2 ( 1.0e-1 )  //ひもの硬さ
#define K_BOND_3 ( 1.0e-1)
#define HMM_BOND (1.0e-0)

#define GYRATION_RADIUS (20 * 10e+3 / 196 * 1.4e-8 * 75 / 2.2e-6)   // GYRATION_RADIUS * 粒子数 = 慣性半径
#define PARTICLE_MYU ( 2.0 * DIMENSION * PI * PARTICLE_RADIUS * 0.000890) /100 //粘性抵抗の強さ

#define DELTA ( 1.0e-5 )  //刻み幅
#define MITIGATION (1.0e+5)

#define POTENTIAL_DELTA (1.0e-7)

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
    double *nucleolus_mean;
    double *spb_mean;
    double force[DIMENSION];
    unsigned int *gr_list;
} Particle;

enum label{ X, Y, Z};

void secure_main_memory (Particle **part, unsigned int **locus_list) {   // メモリ確保 //
    
    if ( (*part = (Particle *)malloc(NUMBER_MAX * sizeof(Particle))) == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    if ( (*locus_list = (unsigned int *) calloc ( 45, sizeof (unsigned int)) ) == NULL) {
        
        printf ("\t error : can not secure the memory of locus_list\n");
        exit (1);
    }
}

void secure_sub_memory (Particle *locus) {  // locus粒子限定のメモリ確保  //
    
    if ( (locus->spb_mean = (double *) malloc (RANK * sizeof (double)))== NULL) {
        
        printf ("\t error : can not secure the memory of spb_mean \n");
        exit (1);
    }
    
    if ( (locus->nucleolus_mean = (double *) malloc (RANK * sizeof (double)))== NULL) {
        
        printf ("\t error : can not secure the memory of nucleolus_mean \n");
        exit (1);
    }
    
    if ( (locus->gr_list = (unsigned int *)malloc (RANK * sizeof (unsigned int))) == NULL) {
        
        printf ("\t error : can not secure the memory of gr_list\n");
        exit(1);
    }
}

void free_useless_memory (Particle **part, unsigned int **locus_list, const double particle_number) {
    
    
    unsigned int locus_number = 0;
    while ( (*locus_list)[locus_number] != 0) locus_number++;
    
    if ( (*locus_list = (unsigned int *) realloc ( *locus_list, locus_number * sizeof (unsigned int))) == NULL) {
        
        printf ("\t error : can not shrink the memory of locus_list\n");
        exit (1);
    }
    
    if ( (*part = (Particle *) realloc ( *part, particle_number * sizeof (Particle))) == NULL) {
        
        printf ("\t error : can not shrink the memory of part \n");
        exit (1);
    }
}

void read_data (Particle *part, char *cycle_dtatus, char *arm_id, unsigned int *locus_list, unsigned int *particle_number, unsigned int *locus_number){       //初期値設定

    unsigned int loop, number = 0, locus_count = 0, i_dummy;
    
    char dummy[256], pastis_data[128], hmm_data[128];
    double d_dummy, enlarge_ratio;
    
    Particle *part_1;
    FILE *fpr;
    
    // Input the coordinates of Pastis //
    sprintf(pastis_data, "MDS_pos_%s.txt", arm_id);
    if ((fpr = fopen(pastis_data, "r")) == NULL){
        
        printf ("\n\terror : cannot read coordinate.\n");
        
        exit (1);
    }
    
    fgets (dummy, 256, fpr);
    
    while (fscanf (fpr, "%d ", &part[number].pastis_no) != EOF) {
        
        part_1 = &part[number];
        fscanf (fpr, "%s %lf %lf %lf\n", &dummy,
                &part_1->position[X], &part_1->position[Y], &part_1->position[Z]);
        number++;
    }
    
    particle_number = &number;
    
    fclose (fpr);
    
    // HMMの平均データ読み込み //
    sprintf (hmm_data, "hmm_%s_%s.txt", cycle_dtatus, arm_id);
    
    if ((fpr = fopen(hmm_data, "r")) == NULL){
        
        printf ("\n\terror : cannot read hmm data.\n");
        
        exit (1);
    }

    unsigned int pastis_no;
    number = 0;
    while (fscanf (fpr, "%d %d", &pastis_no, &i_dummy) != EOF) {
        
        while (pastis_no != part[number].pastis_no) {
            
            number++;
        }
        
        part_1 = &part[number];
        secure_sub_memory (part_1);
        
        for (loop = 0; loop < RANK; loop++) {
            
            fscanf (fpr, " %lf %lf", &part_1->nucleolus_mean[loop], &part_1->spb_mean[loop]);
        }
        fgets (dummy, 256, fpr);
        
        locus_list[locus_count] = number;
        locus_count ++;
    }
    
    locus_number = &locus_count;
    
    fclose (fpr);
    
    loop = 0;
    while (locus_list[loop] != 0) {
        
        part_1 = &part[locus_list[loop]];
        printf ("\t%d status[0] %lf %lf, status[RANK-1] %lf %lf \n", part_1->pastis_no, part_1->spb_mean[0], part_1->nucleolus_mean[0],
                part_1->spb_mean[RANK-1], part_1->nucleolus_mean[RANK-1]);
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
    
    }
    */
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}


//　ばねによる力 part_1 粒子側の力計算//
void spring (Particle *part_1, const Particle *part_2, const unsigned int bond) {
    
    double dist, dist_0;
    
    double f;
    
    //dist_0 = 自然長 //
    dist_0 = Euclid_norm (part_1->position_init, part_2->position_init);
    dist = Euclid_norm (part_1->position, part_2->position);
    
    f = bond * (dist_0 - dist) / dist;
    
    part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
    part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
    part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
}

/*
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

void calculate() {
    
    int i;
    
    Particle *part_1;
    
    for ( i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
    }
    
    
//#pragma omp parallel for private (part_1) num_threads (6) //gdb
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

        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->force[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->force[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->force[Z];
    }
    
    // position の更新 //
    for ( i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        part_1->position[X] = part_1->position_new[X];
        part_1->position[Y] = part_1->position_new[Y];
        part_1->position[Z] = part_1->position_new[Z];
    }
    
}

void rank_optimization ( unsigned int locus, unsigned int locus_list[45]) {
    
    unsigned int time, rank_flag = 0;
    if (locus != 0) double gyration_radius = GYRATION_RADIUS * (locus_list[locus] - locus_list[locus - 1]);
    Particle *part = &part_now[locus_list[locus]], *part_old = &part[locus_list[locus - 1]];
    
    for (unsigned int rank = 0; rank < RANK; rank++) {
        
        for (time = 0; time < MITIGATION; time++) {
            
            calculate (locus_list[locus]);
        }
        
        if ( locus != 0) {
            
            if (Euclid_norm (part_now->position, part_old->position) < gyration_radius ) {
                
                part_now->gr_list[rank_flag] = rank;
                rannk_flag++;
            }
        }
    }
    
    if
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
*/


int main ( int argc, char **argv ) {
    
    unsigned int loop, t = 0, l, particle_number, locus_number;
    unsigned int *locus_list;
    char output_file[256];
    
    Particle *part, *part_1;
    
    secure_main_memory (&part, &locus_list);
    
    read_data (part, argv[1], argv[2], locus_list, &particle_number, &locus_number);
    
    //free_useless_memory (&part, &locus_list, particle_number);
    
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
/*
    for (loop=0; loop < locus_number; loop++) {
        
        rank_optimization (loop, locus_list);
    }
*/
    for (loop = 0; loop < locus_number; loop++) {
        
        part_1 = &part[locus_list[loop]];
        free (part_1->spb_mean);
        free (part_1->nucleolus_mean);
        free (part_1->gr_list);
    }
    free (part);
    
    return ( 0 );
}



