//
//  teritory.c
//  
//
//  Created by 高山雄揮 on 2017/10/18.
//
//

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
//#include "territory_graph.h"

#define DIMENSION (3) //次元
#define NANO (1.0e-10)   // 長さの単位
#define KBT (1.38064852e-23 / NANO / NANO) //ボルツマン
#define TEMPARTURE (300)

#define NUMBER (6193)      //粒子数
#define PARTICLE_MASS (1.0e-18)      //染色体粒子の質量
#define SPB_MASS (6.0e-18)      //SPBの質量
#define SPHERE_DIVISION (20)     //膜の分割数
#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PI (M_PI)
#define MEMBRAIN_RADIUS_MINIMUM  (32.0) //膜の半径の最小値
#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )

#define SPB_RADIUS ( 3.0 )      //SPBの半径

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
    int list_no;
    int *list;
    int territory_counter;
    
} Particle;

Particle *ptr;

Particle part[NUMBER];

double membrain_radius;

enum label{ X, Y, Z};

void read_coordinate ( int s){
    
    int i, i_dummy;
    
    double d_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    sprintf (filename, "newcal/result_%d.txt", s);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++){
        
        fscanf(fpr, "%d %d %d %lf %lf %lf %lf %lf %lf %lf \n", &i_dummy, &part[i].chr_no, &part[i].particle_type,
               &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
               &d_dummy, &d_dummy, &d_dummy,
               &membrain_radius);
        //fgets(dummy, 128, fpr);
        
        //printf ("%d %lf %lf %lf \n", i, part[i].position[X], part[i].position[Y], part[i].position[Z]);
        
    }
    
    fgets(dummy, 128, fpr);
    fgets(dummy, 128, fpr);
    
    fclose (fpr);
}

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

void write_data(int s, int t, double per_sum[3]){
    
    char filename[128];
    char strs;
    
    
    FILE *fpw;
    
    sprintf( filename, "newcal/new_ter_%d.txt", s);
    
    if ((fpw = fopen(filename, "a")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    fprintf(fpw, "%d %lf %lf %lf\n", t, per_sum[0]/2771, per_sum[1]/2241, per_sum[2]/1180);
    
    fclose(fpw);
    
}
/*
void make_graph(){
    
    FILE *gp;
    
    gp = popen ("gnuplot -persist", "w");
    
    //fprintf (gp, "set)
    
    fprintf (gp, "plot %s w l )
}
*/

void graph (int filenumber) {
    
    FILE *gp;
    
    char filename[256];
    
    if ((gp = popen ("gnuplot -persist", "w")) == NULL){
        
        printf (" cannnot open gnuplot \n");
        
        exit (1);
    }
    
    printf ("pass\n");
    
    sprintf (filename, "territory_%d.txt", filenumber);
    
    fprintf ( gp, " set yrange[0:100]\n");
    fprintf ( gp, " set linetype 1 lc rgb \"red\" lw 3\n");
    fprintf ( gp, " set linetype 2 lc rgb \"blue\" lw 3\n");
    fprintf ( gp, " set linetype 3 lc rgb \"green\" lw 3\n");
    
    fprintf ( gp, " plot \"%s\" using 1:2 w l lt 1 title \"chromosome 1\"\n", filename);
    fprintf ( gp, " replot \"%s\" using 1:3 w l lt 2 title \"chromosome 2\"\n", filename);
    fprintf ( gp, " replot \"%s\" using 1:4 w l lt 3 title \"chromosome 3\"\n", filename);
    
    fprintf ( gp, " set term postscript enhanced color\n");
    fprintf ( gp, " set output \"ter%d.eps\"\n", filenumber);
    fprintf ( gp, " replot\n");
    fprintf ( gp, " set term x11\n");
    fprintf ( gp, " set output\n");
    
    pclose ( gp );
}
int main (int argc, char **argv ) {
    
    int start_number = atoi(argv[1]);
    int finish_number = atoi(argv[2]);
    int width = atoi(argv[3]);
    int file_number = atoi(argv[4]);
    
    int i, j, m, t;
    double per_sum[3];
    
    double dist;
    
    Particle *part_1, *part_2;
    
    ptr = (Particle *)malloc(NUMBER * sizeof(Particle));
    
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
    
    //#pragma omp parallel for private (per_sum, i, j, part_1, part_2, dist, m)
    for ( t=0;  t + start_number < finish_number; t += width) {
        
        read_coordinate( t + start_number);
        
        per_sum[0] = 0;
        per_sum[1] = 0;
        per_sum[2] = 0;
        
        for ( i = 0; i<NUMBER; i++){
            
            part_1 = &part[i];
            
            part_1->territory_counter = 0;
            
            m = 0;
            
            for ( j=0; j<NUMBER; j++){
                
                part_2 = &part[j];
                
                dist = Euclid_norm (part_1->position, part_2->position);
                
                if (dist < 10.0 * PARTICLE_RADIUS && part_1 != part_2){
                    
                    m++;
                    
                    if (part_1->chr_no == part_2->chr_no){
                        
                        part_1->territory_counter ++;
                    }
                    
                }
            }
            
            per_sum[part_1->chr_no] += part_1->territory_counter * 100.0 / m;
            
        }
        
        printf ("%d %lf %lf %lf\n", t, per_sum[0]/2771, per_sum[1]/2241, per_sum[2]/1181 );
        
        write_data ( file_number , t, per_sum);
    }
    
    // graph ( file_number );
    
}




















