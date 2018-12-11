//
//  str_cal.c
//  
//
//  Created by tkym on 2018/11/26.
//

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
#include <omp.h>

#define DIMENSION ( 3 ) //次元
#define LENGTH ( 1.4e-6 / 75)   // 長さの単位
//#define M_A ( 1.85131596e+7 )
//#define N_A ( 6.022140857e+23 )

#define NUMBER ( 200 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )

#define K_BOND ( 1.0e-2 )    //ばね定数
#define K_BOND_2 ( 1.0e-2 )  //ひもの硬さ
#define K_BOND_3 ( 1.0e-2)
#define HMM_BOND (1.0e-1)

#define PARTICLE_MYU ( 2.0 * DIMENSION * PI * PARTICLE_RADIUS * 0.000890) //粘性抵抗の強さ

#define DELTA ( 1.0e-9 )  //刻み幅

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

typedef struct particle {           //構造体の型宣言
    //CHAIN chr_no;
    unsigned int pastis_no;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_init[DIMENSION];
    double velocity[DIMENSION];
    double velocity_2[DIMENSION];
    double nucleolus_mean;
    //double nucleolus_var;
    double spb_mean;
    //double spb_var;
    double force[DIMENSION];
    //int list_no;
    //int *list;
    
} Particle;

Particle *part;

enum label{ X, Y, Z};

void read_data ( char *filename ){       //初期値設定
    
    unsigned int i, number = 0;
    
    char dummy[256];
    double d_dummy;
    
    Particle *part_1;
    FILE *fpr;
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
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
    
    for (i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        part_1->position_init[X] = part_1->position[X];
        part_1->position_init[Y] = part_1->position[Y];
        part_1->position_init[Z] = part_1->position[Z];
        
        part->velocity_2[X] = 0.0;
        part->velocity_2[Y] = 0.0;
        part->velocity_2[Z] = 0.0;
        
        part_1->nucleolus_mean *= 1.0e-6/ LENGTH;
        part_1->spb_mean *= 1.0e-6 / LENGTH;
        
        if (i == 1) printf ("nuc_mean , spb_mean = %lf, %lf\n", part_1->nucleolus_mean, part_1->spb_mean);
    }
}

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
        
        double spb_f = HMM_BOND * (spb_dist - part_1->spb_mean);
        double nucleolus_f = HMM_BOND * (nucleolus_dist - part_1->nucleolus_mean);
        
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
        
        /*
        part_1->velocity[X] = part_1->velocity_2[X] + DELTA * part_1->force[X] / (2.0 * PARTICLE_MASS);
        part_1->velocity[Y] = part_1->velocity_2[Y] + DELTA * part_1->force[Y] / (2.0 * PARTICLE_MASS);
        part_1->velocity[Z] = part_1->velocity_2[Z] + DELTA * part_1->force[Z] / (2.0 * PARTICLE_MASS);
        
        part_1->velocity_2[X] = part_1->velocity[X] + DELTA * part_1->force[X] / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Y] = part_1->velocity[Y] + DELTA * part_1->force[Y] / ( 2.0 * PARTICLE_MASS );
        part_1->velocity_2[Z] = part_1->velocity[Z] + DELTA * part_1->force[Z] / ( 2.0 * PARTICLE_MASS );
        
        part_1->position_new[X] = part_1->position[X] + DELTA * part_1->velocity_2[X];
        part_1->position_new[Y] = part_1->position[Y] + DELTA * part_1->velocity_2[Y];
        part_1->position_new[Z] = part_1->position[Z] + DELTA * part_1->velocity_2[Z];
        */
        
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

double calculate_potential () {
    
    unsigned int i;
    Particle *part_1;
    double potential_V = 0.0, dist_1, dist_2;
    
    for (i=0; i<particle_number; i++) {
        
        part_1 = &part[i];
        
        // 隣同士 //
        if ( i != 0 && i != particle_number - 1 ) {
            
            dist_1 = Euclid_norm (part_1->position, part[i-1].position);
            dist_2 = Euclid_norm (part_1->position, part[i+1].position);
            
            potential_V += K_BOND * (dist_1 * dist_1 + dist_2 * dist_2) / 2.0 * POTENTIAL_DELTA;
        }
        else if ( i == 0) {
            
            dist_1 = Euclid_norm (part_1->position, part[i+1].position);
            potential_V += K_BOND * dist_1 * dist_1 / 2.0 * POTENTIAL_DELTA;
        }
        else {
            
            dist_1 = Euclid_norm (part_1->position, part[i-1].position);
            potential_V += K_BOND * dist_1 * dist_1 / 2.0 * POTENTIAL_DELTA;
        }
        
        // 2個隣 //
        if ( 0 <= i-2 && i+2 <= particle_number-1) {
            
            dist_1 = Euclid_norm (part_1->position, part[i-2].position) ;
            dist_2 = Euclid_norm (part_1->position, part[i+2].position);
            
            potential_V += K_BOND_2 * (dist_1 * dist_1 + dist_2 * dist_2) / 2.0 * POTENTIAL_DELTA;
        }
        else if ( 0 <= i-2) {
            
            dist_1 = Euclid_norm (part_1->position, part[i+2].position);
            potential_V += K_BOND_2 * dist_1 * dist_1 / 2.0 * POTENTIAL_DELTA;
        }
        else{
            
            dist_1 = Euclid_norm (part_1->position, part[i-2].position);
            potential_V += K_BOND_2 * dist_1 * dist_1 / 2.0 * POTENTIAL_DELTA;
        }
        
        // 3個隣 //
        if ( 0 <= i-3 && i+3 <= particle_number-1) {
            
            dist_1 = Euclid_norm (part_1->position, part[i-3].position);
            dist_2 = Euclid_norm (part_1->position, part[i+3].position);
            
            potential_V += K_BOND_3 * (dist_1 * dist_1 + dist_2 * dist_2) / 2.0 * POTENTIAL_DELTA;
        }
        else if ( 0 <= i-3) {
            
            dist_1 = Euclid_norm (part_1->position, part[i+3].position);
            potential_V += K_BOND_3 * dist_1 * dist_1 / 2.0 * POTENTIAL_DELTA;
        }
        else {
            
            dist_1 = Euclid_norm (part_1->position, part[i-3].position);
            potential_V += K_BOND_3 * dist_1 * dist_1 / 2.0 * POTENTIAL_DELTA;
        }
        
        if (part_1->spb_mean != 0.0) {
            
            dist_1 = Euclid_norm (part_1->position, spb_pos);
            dist_2 = Euclid_norm (part_1->position, nucleolus_pos);
            
            potential_V += K_BOND * (dist_1 * dist_1 + dist_2 * dist_2) / 2.0 * POTENTIAL_DELTA;
        }
        
        potential_V += PARTICLE_MASS * (part_1->velocity[X] * part_1->velocity[X] + part_1->velocity[Y] * part_1->velocity[Y]
                                        + part_1->velocity[Z] * part_1->velocity[Z]) / 2.0 * POTENTIAL_DELTA;
    }
    
    return potential_V;
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
    
    int i, t = 0, l;
    char input_file[256], hmm_data[256], output_file[256];
    
    unsigned int calculate_number = atoi (argv[1]);
    double init_V;
    
    printf ("\t Input coordinate data : ");
    scanf ("%s", input_file);
    
    printf ("\tK_BOND1 = %lf\n", K_BOND);
    printf ("\tK_BOND2 = %lf\n", K_BOND_2);
    printf ("\tK_BOND3 = %lf\n", K_BOND_3);
    
    /*
    printf ("\t Input output_file : ");
    scanf ("%s", output_file);
    */
     
    part = (Particle *)malloc(NUMBER * sizeof(Particle));
    
    if (part == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    /*
    for (i = 0;i < NUMBER;i++) {
        part[i].list = (int *)malloc(NUMBER * sizeof(int));
        
        if (part[i].list == NULL) {
            printf("\n error : can not secure the memory \n");
            exit(1);
        }
    }
    */
     
    read_data (input_file);
    
    //init_V = calculate_potential();
    //printf ("\tinit_V = %lf\n", init_V);
    
    //初期位置の出力
    write_coordinate (0);
    
    for (t=1; t <= calculate_number; t++) {
        
        for (l=1; l<=1.0e+5; l++){
            
            calculate(l);
            //write_coordinate (/* argv[3],*/ l , start_number);
        }
        
        printf ("\tt = %d \r", t);
        //printf("    t = %d\t V - init_V = %lf \r", t, calculate_potential()-init_V);
        fflush (stdout);
        
        write_coordinate (/* argv[3],*/ t);
    }
    
    return ( 0 );
}



