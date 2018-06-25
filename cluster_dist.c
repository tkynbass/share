
///     各中心方向への力を加えた遺伝子の、中心との距離
///     created by tkym 2018.05.01

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define DIMENSION (3)
#define NUMBER (6193)      //粒子数
#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PI (M_PI)

#define DIS_DIVISION (10)   //距離の分割数    半径を割る
#define ANG_DIVISION (100)  //角度の分割数    180度を割る

#define TOP_NUMBER (100)
#define BOTTOM_NUMBER (100)

int GENE_NUMBER;


typedef enum chain {
    A, B, C
} CHAIN;

typedef enum type {
    Normal, Centromere, Telomere
}TYPE;

enum label{ X, Y, Z};

typedef struct particle {           //構造体の型宣言
    CHAIN chr_no;
    TYPE particle_type;
    double position[DIMENSION];
    int cluster_state;
    
} Particle;

unsigned int gene_counter, file_number;

double membrain_radius, Nucleolus_circle_center[DIMENSION], k_expression;

Particle part[NUMBER];

void read_coordinate( const char *str, int time ){       //初期値設定
    
    int i, i_dummy;
    
    double d_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    FILE *fpr;
    
    sprintf (filename, "%s/%s_c%d_%d.txt", str, str, file_number, time);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++){
        
        fscanf (fpr, "%d %d %d %lf %lf %lf %lf %lf %lf %lf",
                &i_dummy, &part[i].chr_no, &part[i].particle_type,
                &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
                &d_dummy, &d_dummy, &d_dummy, &membrain_radius);
        fgets(dummy, 128, fpr);
        
    }
    
    fgets (dummy, 128, fpr);
    
    fscanf (fpr, "%s %s %lf %lf %lf\n", dummy, dummy, &Nucleolus_circle_center[X],
            &Nucleolus_circle_center[Y], &Nucleolus_circle_center[Z]);
    
    fclose (fpr);
    
}

void read_expression_data(unsigned int gene_list[200]) {
    
    FILE *fpr;
    
    int i, i_dummy, j, number;
    
    char filename[128], dummy[128];
    
    sprintf (filename, "cl%d_num.txt", file_number);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("\n error \n");
        
        exit (1);
    }
    
    switch (file_number) {
        case 1:
            GENE_NUMBER = 86;
            break;
            
        case 2:
            GENE_NUMBER = 78;
            break;
            
        case 3:
            GENE_NUMBER = 46;
            break;
            
        case 4:
            GENE_NUMBER = 146;
            break;
    }
    
    for (i=0; i<GENE_NUMBER; i++) {
        
        fscanf (fpr,"%d\n", &gene_list[i]);
        
        part[gene_list[i]].cluster_state = 1;
    }
    
    fclose (fpr);
}

void write_hist (const char *number, const unsigned int time, const unsigned int dis_hist[32]) {
    
    unsigned int i, j;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "cluster_hist/%s_hist_%d.txt", number, time);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    for (i=0; i<64; i++) {
        
        fprintf (fpw, "%d %d\n", i, dis_hist[i]);
    }
    
    fclose (fpw);
    
}

void write_graph (const char *number, const unsigned int time, const double dist) {
    
    unsigned int i, j;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "%s_graph.txt", number);
    
    if ((fpw = fopen (result, "a")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    fprintf (fpw, "%d %lf\n", time, dist);
    
    fclose (fpw);
    
}

void write_cluster (const char *number, const unsigned int time, const double cluster_value) {
    
    unsigned int i, j;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "%s_cluster_calue.txt", number);
    
    if ((fpw = fopen (result, "a")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    fprintf (fpw, "%d %lf\n", time, cluster_value);
    
    fclose (fpw);
    
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

double inner_product (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    return ( pos_1[X] * pos_2[X] + pos_1[Y] * pos_2[Y] + pos_1[Z] * pos_2[Z]);
}

int main ( int argc, char **argv) {
    
    unsigned int i, j, t, type, width;
    
    k_expression = atof (argv[1]);
    int start_number = atoi (argv[2]);
    int calculate_number = atoi (argv[3]);
    
    Particle *ex_gene[200], *part_1;
    unsigned int gene_list[200], dis_hist[64];
    
    double origin[] = {0.0, 0.0, 0.0}, gene_dist = 0.0, counter = 0.0, ex_counter = 0.0, cluster_average;
    
    printf ("group_number : ");
    scanf ("%d", &file_number);
    

    
    printf ("type (graph:0 or hist:1 or cluster:2) : ");
    scanf ("%d", &type);
    
    printf ("width : ");
    scanf ("%d", &width);
    
    for (i=0; i<NUMBER; i++) {
        
        part[i].cluster_state = 0;
    }

    read_expression_data (gene_list);
    
    for ( t=0; t < calculate_number; t+= width) {
        
        printf ("t = %d\r", t);

        //座標の読み込み
        read_coordinate ( argv[1], t + start_number);
        
        for (i=0; i<GENE_NUMBER; i++) {
            
            ex_gene[i] = &part[gene_list[i]];
            
            //printf ("%lf %lf %lf \n", ex_gene[i]->position[X], ex_gene[i]->position[Y], ex_gene[i]->position[Z]);
        }

        switch (type) {
            case 0:
                
                gene_dist = 0.0;
                
                for (i=0; i<GENE_NUMBER; i++) {
                    
                    //printf ("%lf\n", Euclid_norm ( ex_gene[i]->position, origin));
                    gene_dist += Euclid_norm ( ex_gene[i]->position, origin);
                }
                
                write_graph( argv[1], t + start_number, gene_dist/GENE_NUMBER);
                break;
                
            case 1:
                
                for (i=0; i<64; i++) {
                    
                    dis_hist [i] = 0;
                }
                
                for (i=1; i<GENE_NUMBER; i++) {
                    
                    for (j=0; j<i; j++) {
                        
                        gene_dist = Euclid_norm ( ex_gene[i]->position, ex_gene[j]->position);
                        
                        dis_hist [ (int) gene_dist ] ++;
                    }
                }
                
                write_hist ( argv[1], t + start_number, dis_hist );
                break;
                
            case 2:
                
                cluster_average = 0.0;
                
                for (i=0; i<GENE_NUMBER; i++ ) {
                    
                    counter = 0.0;
                    ex_counter = 0.0;
                    
                    for (j = 0; j<NUMBER; j++) {
                        
                        if ( Euclid_norm (ex_gene[i]->position, part[j].position) < 10.0 * PARTICLE_RADIUS) {
                            
                            counter += 1.0;
                            
                            if (part[j].cluster_state == 1) {
                                
                                ex_counter += 1.0;
                            }
                        }
                    }
                    
                    cluster_average += ex_counter / counter / GENE_NUMBER * 100;
                }
                
                write_cluster ( argv[1], t+start_number, cluster_average);
        }

    }
    
    
    
    return (0);
}

