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
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

const double nucleolus_pos[] = { 0.0, 0.0, 0.0};
const double spb_pos[] = { 1.7128 * LENGTH, 0.0, 0.0};
double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    CHAIN chr_no;
    unsigned int pastis_no;
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
        
        part[loop].list = (unsigned int *) malloc (NUMBER_MAX * sizeof (unsigned int));
        
        if (part[loop].list = NULL) {
            
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
    char dummy[256], pastis_data[128], hmm_data[128], **arm_list = {"1long", "1short", "2short", "2long", "3short", "3long"};
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
            fscanf (fpr, "%s %lf %lf %lf\n", dummy,
                    &part_1->position[X], &part_1->position[Y], &part_1->position[Z]);
            number++;
        }
    }
    
    
    
    fclose (fpr);
    
    // HMMの平均データ読み込み //
    sprintf (hmm_data, "hmm_%s_%s.txt", cycle_status, arm_id);
    
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
            
            // spb, nucleollus mean のスケール調節 //
            part_1->nucleolus_mean[loop] *= LENGTH;
            part_1->spb_mean[loop] *= LENGTH;
        }
        fgets (dummy, 256, fpr);
        
        locus_list[ (*locus_number) ] = number;
        (*locus_number)++;
    }
    
    fclose (fpr);
    
    /*
    loop = 0;
    while (locus_list[loop] != 0) {
        
        part_1 = &part[locus_list[loop]];
        printf ("\t%d status[0] %lf %lf, status[RANK-1] %lf %lf \n", part_1->pastis_no, part_1->spb_mean[0], part_1->nucleolus_mean[0],
                part_1->spb_mean[RANK-1], part_1->nucleolus_mean[RANK-1]);
        loop++;
    }*/
    
    for (loop = 0; loop < *particle_number; loop++) {
        
        part_1 = &part[loop];
        
        // 初期座標の保存
        part_1->position_init[X] = part_1->position[X];
        part_1->position_init[Y] = part_1->position[Y];
        part_1->position_init[Z] = part_1->position[Z];
     
        part_1->position_state[X] = part_1->position[X];
        part_1->position_state[Y] = part_1->position[Y];
        part_1->position_state[Z] = part_1->position[Z];
    }
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

void write_coordinate (Particle *part, const unsigned int time, const unsigned int locus, const unsigned int particle_number, const unsigned int rank) {
    
    unsigned int loop;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "l%d_r%d_t%d.txt", locus, rank, time);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \t error : cannot write coordinate. \n");
        
        exit (1);
    }
    
    for (loop = 0; loop < particle_number; loop++) {
        
        part_1 = &part[loop];
        fprintf (fpw, "%d %d %lf %lf %lf\n", loop, part_1->pastis_no, part_1->position[X],
                 part_1->position[Y], part_1->position[Z]);
    }
    
    fclose (fpw);
}

void write_optimal_coordinate (Particle *part, const unsigned int particle_number, const unsigned int start_rank, const unsigned int order_flag) {
    
    unsigned int loop, rank;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "optimal_s%d.txt", start_rank);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \t error : cannot write optimal coordinate. \n");
        exit (1);
    }
    
    for (loop = 0; loop < particle_number; loop++) {
        
        if (order_flag == 0) part_1 = &part[loop];
        else part_1 = &part[particle_number - 1 - loop];
        
        if (part_1->gr_list != NULL) rank = part_1->gr_list[0];
        else rank = 99;
        
        fprintf (fpw, "%d %d %d %lf %lf %lf\n", loop, part_1->pastis_no, rank, part_1->position[X],
                 part_1->position[Y], part_1->position[Z]);
    }
    
    fclose (fpw);
}

void write_gr_list (Particle *part, const unsigned int locus_list[45], const unsigned int start_rank, const unsigned int locus_number, const unsigned int order_flag) {
    
    FILE *fpw;
    Particle *part_1;
    char filename[128];
    unsigned int count = 0;
    
    sprintf (filename, "gr_list_s%d.txt", start_rank);
    
    if ( (fpw = fopen (filename, "w")) == NULL ) {
        
        printf ("\t error : cannot write gr_list \n");
        exit(1);
    }
    
    for (unsigned int locus = 0; locus < locus_number; locus++) {
        
        count = 0;
        
        if (order_flag == 0) part_1 = &part[locus_list[locus]];
        else part_1 = &part[locus_list[locus_number - 1 - locus]];
        
        fprintf (fpw, "%d", part_1->pastis_no);
        
        while (part_1->gr_list[count] != 99) {
            
            fprintf (fpw, " %d", part_1->gr_list[count]);
            count++;
        }
        fprintf (fpw, "\n");
    }
    
}

void rank_optimization (Particle *part, unsigned int locus_list[45], const unsigned int locus_number, const unsigned int particle_number, const unsigned int writing_flag
                        , const unsigned int order_flag) {
    
    unsigned int time, rank_flag = 0, start_number, start_rank, rank, locus, loop;
    double gyration_radius;
    
    Particle *part_now, *part_old,  *part_1;
    
    for (unsigned int start_rank = 0; start_rank < RANK; start_rank++) {
        
        // 0番目のローカスにポテンシャルをかけ、全体を移動
        for (time = 0; time < MITIGATION; time++) {
            
            calculate (part, locus_list[0], 0, start_rank, particle_number);
        }
        
        for (loop = 0; loop < particle_number; loop++) {
            
            part_1 = &part[loop];
            
            // 0番目ローカス計算後の座標を保存 //
            part_1->position_state[X] = part_1->position[X];
            part_1->position_state[Y] = part_1->position[Y];
            part_1->position_state[Z] = part_1->position[Z];
            
            // gr_listの初期化 //
            if (part_1->gr_list != NULL) {
                
                for (rank = 0; rank < RANK; rank++) {
                    
                    part_1->gr_list[rank] = 99;
                }
            }
        }
        
        part[locus_list[0]].gr_list[0] = start_rank;
        printf ("\t [locus = %d] start_rank = %d \n", part[locus_list[0]].pastis_no, start_rank);
        
        for (locus = 1; locus < locus_number; locus++) {
            
            gyration_radius = 0.5 * NUCLEOSOME_LENGTH * sqrt( GYRATION_N * (locus_list[locus] - locus_list[locus - 1]));
            part_now = &part[locus_list[locus]];
            part_old = &part[locus_list[locus - 1]];
            
            rank_flag = 0;
            start_number = locus_list[locus - 1] + 1;
            
            for (rank = 0; rank < RANK; rank++) {
                
                if (part_now->spb_mean[rank] != 0 && part_now->nucleolus_mean[rank] != 0) {
                    
                    for (time = 0; time < MITIGATION; time++) {
                        
                        calculate (part, locus_list[locus], start_number, rank, particle_number);
                    }
                    
                    if (Euclid_norm (part_now->position, part_old->position) < gyration_radius ) {
                        
                        part_now->gr_list[rank_flag] = rank;
                        rank_flag++;
                    }
                    
                    for (loop = start_number; loop < particle_number; loop++) {
                        
                        part_1 = &part[loop];
                        
                        part_1->position[X] = part_1->position_state[X];
                        part_1->position[Y] = part_1->position_state[Y];
                        part_1->position[Z] = part_1->position_state[Z];
                    }
                }
            }
            
            if (rank_flag != 0) {
            
                // 慣性半径内に存在する状態のうち、最高ランクのもので座標決定 //
                for (time = 0; time < MITIGATION; time++) {
                    
                    calculate (part, locus_list[locus], start_number, part_now->gr_list[0], particle_number);
                    
                    if (time % WRITE_INTERVAL == 0 && writing_flag == 1) write_coordinate (part, time, part_now->pastis_no, particle_number, part_now->gr_list[0]);
                }
                
                for (loop = start_number; loop < particle_number; loop++) {
                    
                    part_1 = &part[loop];
                    
                    part_1->position_state[X] = part_1->position[X];
                    part_1->position_state[Y] = part_1->position[Y];
                    part_1->position_state[Z] = part_1->position[Z];
                }
                
                printf ("\t [locus = %d] Total %d states = { ", part_now->pastis_no, rank_flag);
                
                for ( loop = 0; loop < rank_flag; loop++) printf ("%d ", part_now->gr_list[loop]);
                printf ("}\n");
            }
            else {
                
                // 慣性半径内にいる状態がなければ、break → 0番目ローカスのランク(start_rank)を変えてやり直し //
                printf ("\t There's no state of %d in gyration radius ( start_rank = %d)\n", part_now->pastis_no, start_rank);
                break;
            }
        }
        printf ("\n");
        
        if (locus < locus_number - 1) write_optimal_coordinate (part, particle_number, start_rank, order_flag);
        write_gr_list (part, locus_list, start_rank, locus_number, order_flag);
    }
}


int main ( int argc, char **argv ) {
    
    unsigned int locus, t = 0, l, particle_number, locus_number = 0, writing_flag, order_flag = 0;
    unsigned int *locus_list;
    char output_file[256];
    
    Particle *part, *part_1. spb;
    
    secure_main_memory (&part, &spb, &locus_list);
    
    read_data (part);
    
    printf ("\t particle_number = %d \n\t locus_number = %d \n", particle_number, locus_number);
    
    //printf ("\t part[locus_list[0]].spb_mean[0] = %lf\n", part[locus_list[0]].spb_mean[0]);
    
    if ( strcmp (long_1, argv[2]) == 0 || strcmp (short_2, argv[2]) == 0 || strcmp (short_3, argv[2]) == 0 ) {
        
        reverse_order ( &part, locus_list, particle_number, locus_number);
        order_flag = 1;
    }
    
    //printf ("\t part[locus_list[locus_number -1]].spb_mean[0] = %lf\n", part[locus_list[locus_number -1]].spb_mean[0]);
    
    //free_useless_memory (&part, &locus_list, particle_number);
    
    rank_optimization (part, locus_list, locus_number, particle_number, writing_flag, order_flag);

    for (locus = 0; locus < locus_number; locus++) {
        
        part_1 = &part[locus_list[locus]];
        free (part_1->spb_mean);
        free (part_1->nucleolus_mean);
        free (part_1->gr_list);
    }
    free (part);
    free (locus_list);
    
    return ( 0 );
}



