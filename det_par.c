//
//  determine_spb_nucleolus.c
//  
//
//  Created by tkym on 2019/01/15.
//

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<string.h>
//#include <omp.h>

#define DIMENSION ( 3 ) //次元
//#define LENGTH ( 1.4e-6 / 75)   // 長さの単位
//#define M_A ( 1.85131596e+7 )
//#define N_A ( 6.022140857e+23 )

#define NUMBER ( 200 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )

#define K_BOND ( 1.0e+2 )    //ばね定数
#define K_BOND_2 ( 1.0e+2 )  //ひもの硬さ
#define K_BOND_3 ( 1.0e+2)
#define HMM_BOND (1.0e-0)
#define NUC_OUT_FORCE (1.0e-0)

#define PARTICLE_MYU ( 2.0 * DIMENSION * PI * PARTICLE_RADIUS * 0.000890) //粘性抵抗の強さ

#define DELTA ( 1.0e-7 )  //刻み幅

#define POTENTIAL_DELTA (1.0e-7)

unsigned int particle_number = 6;

typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    //CHAIN chr_no;
    unsigned int pastis_no;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_init[DIMENSION];
    double velocity[DIMENSION];
    //double velocity_2[DIMENSION];
    double nucleolus_mean;
    //double nucleolus_var;
    double spb_mean;
    //double spb_var;
    double force[DIMENSION];
    //int list_no;
    //int *list;
    
} Particle;

Particle *part;
Particle spb;
Particle nucleolus;

enum label{ X, Y, Z};

void read_data ( const double nuclear_radius, const unsigned int number_list[6] ){       //初期値設定
    
    unsigned int i, number = 0;
    
    char dummy[256], filename[256];
    double d_dummy, average_dist = 0.0, origin[] = {0.0, 0.0, 0.0};
    
    Particle *part_1;
    FILE *fpr;
    
    for (i=0; i<3; i++) {
        
        sprintf (filename, "MDS_data_%dshort.txt", i+1);
        
        if ((fpr = fopen(filename, "r")) == NULL){
            
            printf ("\n\terror : cannot read coordinate.\n");
            
            exit (1);
        }
        
        part_1 = &part [2*i];
        
        fgets (dummy, 256, fpr);
        
        while (fscanf (fpr, "%d ", &part_1->pastis_no) != EOF) {
            
            if (part_1->pastis_no == number_list[2*i] ) {
                
                fscanf (fpr, "%s %lf %lf %lf %lf %lf %lf %lf\n", &dummy,
                        &part_1->position[X], &part_1->position[Y], &part_1->position[Z], &part_1->nucleolus_mean, &part_1->spb_mean, &d_dummy, &d_dummy);
            }
            else fgets (dummy, 256, fpr);
        }
        
        fclose (fpr);
        
        sprintf (filename, "MDS_data_%dlong.txt", i+1);
        
        if ((fpr = fopen(filename, "r")) == NULL){
            
            printf ("\n\terror : cannot read coordinate.\n");
            
            exit (1);
        }
        
        part_1 = &part [2*i + 1];
        
        fgets (dummy, 256, fpr);
        
        while (fscanf (fpr, "%d ", &part_1->pastis_no) != EOF) {
            
            if (part_1->pastis_no == number_list[2*i + 1] ) {
                
                fscanf (fpr, "%s %lf %lf %lf %lf %lf %lf %lf\n", &dummy,
                        &part_1->position[X], &part_1->position[Y], &part_1->position[Z], &part_1->nucleolus_mean, &part_1->spb_mean, &d_dummy, &d_dummy);
            }
            else fgets (dummy, 256, fpr);
        }
        
        fclose (fpr);
    }
    
    // 初期値設定 //
    nucleolus.position[X] = 0.0;
    nucleolus.position[Y] = 0.0;
    nucleolus.position[Z] = 50.0;
    
    spb.position[X] = 0.0;
    spb.position[Y] = 0.0;
    spb.position[Z] = 0.0;
    
    for (i=0; i<6; i++) {
        
        part_1 = &part[i];
        
        part_1->nucleolus_mean *= 75 / nuclear_radius;
        part_1->spb_mean *= 75 / nuclear_radius;
        
        spb.position[X] += part_1->position[X] / 6.0;
        spb.position[Y] += part_1->position[Y] / 6.0;
        spb.position[Z] += part_1->position[Z] / 6.0;
        average_dist += part_1->spb_mean / 6.0;
    }
    
    spb.position[X] *= 2.0 + average_dist / Euclid_norm (spb.position, origin);
    spb.position[Y] *= 2.0 + average_dist / Euclid_norm (spb.position, origin);
    spb.position[Z] *= 2.0 + average_dist / Euclid_norm (spb.position, origin);

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
        
        double spb_dist = Euclid_norm (part_1->position, spb.position);
        double nucleolus_dist = Euclid_norm (part_1->position, nucleolus.position);
        
        double spb_f = HMM_BOND * (part_1->spb_mean - spb_dist);
        double nucleolus_f = HMM_BOND * (part_1->nucleolus_mean - nucleolus_dist);
        
        spb.force[X] += spb_f * (spb.position[X] - part_1->position[X]);
        spb.force[Y] += spb_f * (spb.position[Y] - part_1->position[Y]);
        spb.force[Z] += spb_f * (spb.position[Z] - part_1->position[Z]);
        
        nucleolus.force[X] += nucleolus_f * (nucleolus.position[X] - part_1->position[X]);
        nucleolus.force[Y] += nucleolus_f * (nucleolus.position[Y] - part_1->position[Y]);
        nucleolus.force[Z] += nucleolus_f * (nucleolus.position[Z] - part_1->position[Z]);
    }
}
 
void calculate( unsigned int l, const double nuclear_radius ) {
    
    int i;
    double spb_nuc_dist = 75 * 1.71 / nuclear_radius;
    
    Particle *part_1;
    
    for ( i=0; i<DIMENSION; i++) {
        
        spb.force[i] = 0.0;
        nucleolus.force[i] = 0.0;
    }
    
    
//#pragma omp parallel for private (part_1) num_threads (6)
    for ( i=0; i<6; i++) {
        
        part_1 = &part[i];
        
        if (part_1->spb_mean != 0.0) {
            
            hmm_potential (part_1);
        }
    }
    
    // spb - nucleolus spring //
    /*
    double f = HMM_BOND * ( spb_nuc_dist - Euclid_norm (spb.position, nucleolus.position));
    
    spb.force[X] += f * (spb.position[X] - nucleolus.position[X]);
    spb.force[Y] += f * (spb.position[Y] - nucleolus.position[Y]);
    spb.force[Z] += f * (spb.position[Z] - nucleolus.position[Z]);
    
    nucleolus.position[X] += f * (nucleolus.position[X] - spb.position[X]);
    nucleolus.position[Y] += f * (nucleolus.position[Y] - spb.position[Y]);
    nucleolus.position[Z] += f * (nucleolus.position[Z] - spb.position[Z]);
    */
    
    spb.velocity[X] = spb.force[X];
    spb.velocity[Y] = spb.force[Y];
    spb.velocity[Z] = spb.force[Z];
    
    nucleolus.velocity[X] = nucleolus.force[X];
    nucleolus.velocity[Y] = nucleolus.force[Y];
    nucleolus.velocity[Z] = nucleolus.force[Z];
    
    spb.position_new[X] = spb.position[X] + DELTA * spb.velocity[X];
    spb.position_new[Y] = spb.position[Y] + DELTA * spb.velocity[Y];
    spb.position_new[Z] = spb.position[Z] + DELTA * spb.velocity[Z];
    
    nucleolus.position_new[X] = nucleolus.position[X] + DELTA * nucleolus.velocity[X];
    nucleolus.position_new[Y] = nucleolus.position[Y] + DELTA * nucleolus.velocity[Y];
    nucleolus.position_new[Z] = nucleolus.position[Z] + DELTA * nucleolus.velocity[Z];
    
    // position の更新 //
    for ( i=0; i<DIMENSION; i++) {
        
        spb.position[i] = spb.position_new[i];
        nucleolus.position[i] = nucleolus.position_new[i];
    }
    
}

void write_data (const int t, const double nuclear_radius, const int calculate_number) {
    
    int i;
    double spb_strain = 0.0, nucleolus_strain = 0.0, spb_nuc_dist = 75 * 1.71 / nuclear_radius;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    /*
    for ( i=0; i<6; i++) {
        
        part_1 = &part[i];
        spb_strain += fabs (Euclid_norm (part_1->position, spb.position) - part_1->spb_mean);
        nucleolus_strain += fabs ( Euclid_norm (part_1->position, nucleolus.position) - part_1->nucleolus_mean);
    }*/
    
    for ( i=0; i<6; i++) {
        
        part_1 = &part[i];
        spb_strain += fabs ( 1.0 - Euclid_norm (part_1->position, spb.position)/part_1->spb_mean);
        nucleolus_strain += fabs ( 1.0 - Euclid_norm (part_1->position, nucleolus.position) / part_1->nucleolus_mean);
    }
    
    // SPB //
    sprintf (result, "spb_%lf.txt", nuclear_radius);
    
    if ((fpw = fopen (result, "a")) == NULL) {
        
        printf ("\n\t error : cannot write data. \n");
        
        exit (1);
    }
    
    fprintf (fpw, "%d %lf %lf %lf %lf\n", t, spb.position[X], spb.position[Y], spb.position[Z], spb_strain);
    
    fclose (fpw);
    
    // Nucleolus //
    sprintf (result, "nucleolus_%lf.txt", nuclear_radius);
    
    if ((fpw = fopen (result, "a")) == NULL) {
        
        printf ("\n\t error : cannot write data. \n");
        
        exit (1);
    }
    
    fprintf (fpw, "%d %lf %lf %lf %lf\n", t, nucleolus.position[X], nucleolus.position[Y], nucleolus.position[Z], nucleolus_strain);
    
    fclose (fpw);
    
    if (t < calculate_number) printf ("\tt = %d, spb_strain : %lf, nucleolus_strain : %lf, spb-nuc : %lf \r", t, spb_strain, nucleolus_strain, Euclid_norm (spb.position, nucleolus.position) / spb_nuc_dist);
    else printf ("\tt = %d, spb_strain : %lf, nucleolus_strain : %lf, spb-nuc : %lf \n", t, spb_strain, nucleolus_strain, Euclid_norm (spb.position, nucleolus.position) / spb_nuc_dist);
    
}


int main ( int argc, char **argv ) {
    
    int i, t = 0, l;
    char input_file[256], hmm_data[256], output_file[256];
    
    unsigned int particleGroup = atoi (argv[1]);
    unsigned int calculate_number = atoi (argv[2]);
    Particle *part_1;
    double spb_strain, nucleolus_strain, nuclear_radius, spb_nuc_dist;
    
    part = (Particle *)malloc(6 * sizeof(Particle));
    
    if (part == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
    
    // 使用する粒子リストの宣言 //
    unsigned int numberList_1[] = { 196, 177, 356, 365, 554, 568};
    unsigned int numberList_3[] = { 205, 170, 347, 374, 543, 578};
    
    //初期位置の出力
    //write_data (0, nuclear_radius, calculate_number);
    
    for (double nuclear_radius = 0.5; nuclear_radius <= 4.0; nuclear_radius += 0.1) {
        
        // 歪み値の初期化 //
        spb_strain = 0.0;
        nucleolus_strain = 0.0;
        spb_nuc_dist = 75 * 1.71 / nuclear_radius;
        
        // データの読み込み //
        if (particleGroup == 1) read_data (nuclear_radius, numberList_1);
        else read_data (nuclear_radius, numberList_3);
        
        for (t=1; t <= calculate_number; t++) {
            
            for (l=1; l<=1.0e+5; l++){
                
                calculate(l, nuclear_radius);
                //write_coordinate (/* argv[3],*/ l , start_number);
            }
            
            //write_data (t, nuclear_radius, calculate_number);
        }
        
        for ( i=0; i<6; i++) {
            
            part_1 = &part[i];
            spb_strain += fabs ( 1.0 - Euclid_norm (part_1->position, spb.position)/part_1->spb_mean);
            nucleolus_strain += fabs ( 1.0 - Euclid_norm (part_1->position, nucleolus.position) / part_1->nucleolus_mean);
        }
        
        printf ("%1.1f %lf %lf %lf\n", nuclear_radius, spb_strain, nucleolus_strain, Euclid_norm (spb.position, nucleolus.position) / spb_nuc_dist);
        
    }
    
    return ( 0 );
}


