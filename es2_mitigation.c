//
//  es2_mitigation.c
//
//  Cut11-Gar2データの楕円体モデル

//  Pastis構造 5kbp のデータ補完,緩和
//
// 2019/04/26
//

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#include <omp.h>

#define DIMENSION ( 3 ) //次元
#define LENGTH ( 7.0e-8 )   //長さの単位 (粒子径)

#define M_A ( 1.85131596e+6 )   // g/mol/particle
#define N_A ( 6.022140857e+23 )
#define PASTIS_SCALING ( 1.8e-6 / 75 / LENGTH)

#define NUMBER_MAX ( 2516 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )

#define K_BOND ( 1.0e-0 )       //1つ隣　ばね定数
#define K_BOND_2 ( 1.0e-4 )     //2つ隣
#define K_BOND_3 ( 1.0e-0 )     //3つ隣

#define DELTA ( 1.0e-11 )  //刻み幅
#define MITIGATION (1.0e+6)
#define WRITE_INTERVAL (1000)

#define MEMBRANE_EXCLUDE ( 1.0 )     //膜との衝突
#define MEMBRANE_EXCLUDE_SPB ( 1.0 ) //SPBとの衝突

#define SPB_RADIUS (  3.0  )      //SPBの半径
#define SPB_MYU ( 2.0 * DIMENSION * PI * SPB_RADIUS * LENGTH * 0.000890 / 100)  //SPBの粘性

//#define POTENTIAL_DELTA (1.0e-7)

typedef enum chain {
    A, B, C
}CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

const double nucleolus_pos[] = { 0.0, 0.0, 0.0};
const double spb_pos[] = { 1.7128 * LENGTH, 0.0, 0.0};
double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    CHAIN chr_no;
    int pastis_no;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_old[DIMENSION];
    //double *nucleolus_mean;
    //double *spb_mean;
    double force[DIMENSION];
    unsigned int list_no;
    unsigned int *list;
    
} Particle;

enum label{ X, Y, Z};

void secure_main_memory (Particle **part, Particle spb) {   // メモリ確保 //
    
    if ( (*part = (Particle *)malloc(NUMBER_MAX * sizeof(Particle))) == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    for (unsigned int loop = 0; loop < NUMBER_MAX; loop++) {
        
        ( *part )[loop].pastis_no = -1;
        ( *part )[loop].list = (unsigned int *) malloc (NUMBER_MAX * sizeof (unsigned int));
        
        if ( ( *part )[loop].list = NULL) {
            
            printf ("\n error : can not secure the memory part.list \n");
            exit (1);
        }
    }
    
    if ( (spb.list = (unsigned int *) malloc (NUMBER_MAX * sizeof (unsigned int))) == NULL) {
        
        printf("\n error : can not secure the memory of spb.list \n");
        exit(1);
    }
    
}

void read_data (Particle *part){       //初期値設定

    unsigned int loop, number = 0, i_dummy;
    char chain, dummy[256], pastis_data[128], hmm_data[128], *arm_list[] = {"1long", "1short", "2short", "2long", "3short", "3long"};
    double d_dummy, enlarge_ratio;
    
    Particle *part_1;
    FILE *fpr;
    
    // Input the coordinates of Pastis //
    for ( unsigned int arm = 0; arm < 6; arm++ ){
        
        sprintf(pastis_data, "MDS_pos_%s.txt", arm_list[arm]);
        if ((fpr = fopen(pastis_data, "r")) == NULL){
            
            printf ("\n\terror : cannot read coordinate.\n");
            
            exit (1);
        }
        
        while (fscanf (fpr, "%d ", &number) != EOF) {
            
            part_1 = &part[number];
            part_1->pastis_no = number;
            fscanf (fpr, "%d %lf %lf %lf\n", &part_1->chr_no,
                    &part_1->position[X], &part_1->position[Y], &part_1->position[Z]);
            
            part_1->position[X] *= PASTIS_SCALING;
            part_1->position[Y] *= PASTIS_SCALING;
            part_1->position[Z] *= PASTIS_SCALING;
        }
    }
    
    fclose (fpr);
    
}

void completion_coordinate (Particle *part) {
    
    unsigned int loop, loop_2, division;
    double distance;
    Particle *part_1, *part_2, *part_3;
    
    // データ補完 //
    int start_list[] = { 0, 1112, 1116, 2024, 2508 };
    int end_list[] = { 2, 1115, 1120, 2036, 2515 };
    
    //　端のデータ補完 //
    for ( loop = 0; loop < sizeof (start_list) / sizeof start_list[0] ; loop++) {
        
        if ( start_list [loop] == 0 || start_list [loop] == 1116 || start_list [loop] == 2024) {
            
            for ( loop_2 = 0; loop_2 <= end_list [loop] - start_list [loop]; loop_2 ++ ) {
                
                part_1 = &part [end_list [loop] - loop_2];
                part_2 = &part [end_list [loop] - loop_2 + 1];
                part_3 = &part [end_list [loop] - loop_2 + 2];
                
                distance = Euclid_norm ( part_2->position, part_3->position );
                
                part_1->position[X] = part_2->position[X] + ( part_2->position[X] - part_3->position[X]) / distance;
                part_1->position[Y] = part_2->position[Y] + ( part_2->position[Y] - part_3->position[Y]) / distance;
                part_1->position[Z] = part_2->position[Z] + ( part_2->position[Z] - part_3->position[Z]) / distance;
            }
        }
        else if ( start_list [loop] == 1112 || start_list [loop] == 2508) {
            
            for ( loop_2 = 0; loop_2 <= end_list [loop] - start_list [loop]; loop_2 ++ ) {
                
                part_1 = &part [start_list [loop] + loop_2];
                part_2 = &part [start_list [loop] + loop_2 - 1];
                part_3 = &part [start_list [loop] + loop_2 - 2];
                
                distance = Euclid_norm ( part_2->position, part_3->position );
                
                part_1->position[X] = part_2->position[X] + ( part_2->position[X] - part_3->position[X]) / distance;
                part_1->position[Y] = part_2->position[Y] + ( part_2->position[Y] - part_3->position[Y]) / distance;
                part_1->position[Z] = part_2->position[Z] + ( part_2->position[Z] - part_3->position[Z]) / distance;
            }
        }
        /*
        else {
            
            part_2 = &part [start_list [loop] - 1];
            part_3 = &part [end_list [loop] + 1];
            
            division = end_list [loop] - start_list [loop] + 2;
            
            for ( loop_2 = 0; start_list [loop] + loop_2 <= end_list [loop]; loop_2++ ){
                
                part_1 = &part [start_list [loop] + loop_2];
    
                part_1->position [X] = part_2->position [X] + (part_2->position[X] - part_3->position[X]) / division * (loop_2 + 1);
                part_1->position [Y] = part_2->position [Y] + (part_2->position[Y] - part_3->position[Y]) / division * (loop_2 + 1);
                part_1->position [Z] = part_2->position [Z] + (part_2->position[Z] - part_3->position[Z]) / division * (loop_2 + 1);
            }
        }*/
    }
    
    unsigned int data_flag = 0, start, end;
    // 穴埋め（セントロメア等含む）のデータ補完 //
    for ( loop = 0; loop < NUMBER_MAX; loop++) {
        
        if ( data_flag == 0 && part[loop].pastis_no < 0 ) {
            
            start = loop;
            flag == 1;
        }
        else if ( data_flag == 1 && part[loop].pastis_no >= 0) {
            
            end = loop;
            flag == 0;

            part_2 = &part [start - 1];
            part_3 = &part [end];
            
            for ( loop_2 = 0; loop_2 < end - start; loop_2 ++) {
                
                part_1 = &part [start + loop_2];
                
                part_1->position [X] = part_2->position [X] + (part_3->position[X] - part_2->position[X]) / (end - start + 1) * (loop_2 + 1);
                part_1->position [Y] = part_2->position [Y] + (part_3->position[Y] - part_2->position[Y]) / (end - start + 1) * (loop_2 + 1);
                part_1->position [Z] = part_2->position [Z] + (part_3->position[Z] - part_2->position[Z]) / (end - start + 1) * (loop_2 + 1);
            }
        }
    }
    
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

/*
//　ばねによる力 part_1 粒子側の力計算//
void spring (Particle *part_1, const Particle *part_2, const unsigned int bond) {
    
    double dist, dist_0;
    
    double f;
    
    //dist_0 = 自然長 //
    
    dist = Euclid_norm (part_1->position, part_2->position);
    
    f = bond * (dist_0 - dist) / dist;
    
    part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
    part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
    part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
}

void hmm_potential (Particle *part_1, const unsigned int rank) {
    
    double spb_dist = Euclid_norm (part_1->position, spb_pos);
    double nucleolus_dist = Euclid_norm (part_1->position, nucleolus_pos);
    
    double spb_f = HMM_BOND * (part_1->spb_mean[rank] - spb_dist);
    double nucleolus_f = HMM_BOND * (part_1->nucleolus_mean[rank] - nucleolus_dist);
    
    part_1->force[X] += spb_f * (part_1->position[X] - spb_pos[X]) + nucleolus_f * (part_1->position[X] - nucleolus_pos[X]);
    part_1->force[Y] += spb_f * (part_1->position[Y] - spb_pos[Y]) + nucleolus_f * (part_1->position[Y] - nucleolus_pos[Y]);
    part_1->force[Z] += spb_f * (part_1->position[Z] - spb_pos[Z]) + nucleolus_f * (part_1->position[Z] - nucleolus_pos[Z]);
}

void calculate (Particle *part, const unsigned int target_locus, const unsigned start_number, const unsigned int rank, const unsigned int particle_number) {
    
    int loop;
    
    Particle *part_1;
    
    for ( loop = start_number; loop < particle_number; loop++) {
        
        part_1 = &part[loop];
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
    }
    
#pragma omp parallel for private (part_1) num_threads (8)
    for ( loop = start_number; loop < particle_number; loop++) {
        
        part_1 = &part[loop];
        
        // 隣同士 //
        if ( loop != 0 && loop != particle_number - 1 ) {
            
            spring (part_1, &part[loop-1], K_BOND);
            spring (part_1, &part[loop+1], K_BOND);
        }
        else if ( loop == 0) spring (part_1, &part[loop+1], K_BOND);
        else spring (part_1, &part[loop-1], K_BOND);
        
        // 2個隣 //
        if ( 2 <= loop && loop+2 <= particle_number-1) {
            
            spring (part_1, &part[loop-2], K_BOND_2);
            spring (part_1, &part[loop+2], K_BOND_2);
        }
        else if ( loop <= 2 ) spring (part_1, &part[loop+2], K_BOND_2);
        else spring (part_1, &part[loop-2], K_BOND_2);
        
        // 3個隣 //
        if ( 3 <= loop && loop+3 <= particle_number-1) {
            
            spring (part_1, &part[loop-3], K_BOND_3);
            spring (part_1, &part[loop+3], K_BOND_3);
        }
        else if ( loop <= 3) spring (part_1, &part[loop+3], K_BOND_3);
        else spring (part_1, &part[loop-3], K_BOND_3);
        
        if (loop == target_locus) hmm_potential (part_1, rank);

        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->force[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->force[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->force[Z];
    }
    
    // position の更新 //
    for ( loop = start_number; loop < particle_number; loop++) {
        
        part_1 = &part[loop];
        
        part_1->position[X] = part_1->position_new[X];
        part_1->position[Y] = part_1->position_new[Y];
        part_1->position[Z] = part_1->position_new[Z];
    }
    
}
*/
void write_coordinate (Particle *part, const unsigned int time) {
    
    unsigned int loop;
    
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
        fprintf (fpw, "%d %d %d %lf %lf %lf\n", loop, part_1->pastis_no, part_1->chr_no, part_1->position[X],
                 part_1->position[Y], part_1->position[Z]);
    }
    
    fclose (fpw);
}


int main ( int argc, char **argv ) {
    
    unsigned int loop;
    char output_file[256];
    
    Particle *part, *part_1, spb;
    
    secure_main_memory (&part, spb);
    
    read_data (part);
    
    completion_coordinate (part);
    
    write_coordinate (part, 0);
    
    for ( loop = 0 ; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part [loop];
        free (part_1->list);
    }
    free (part);
    
    return ( 0 );
}



