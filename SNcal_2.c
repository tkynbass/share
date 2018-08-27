//
//  SNcal.c
//  
//
//  Created by B141829 on 2017/10/20.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <time.h>
#include <string.h>

#define DIMENSION (3) //次元
#define NANO (4.0e-8)   // 長さの単位
#define KBT (1.38064852e-23 / NANO / NANO) //ボルツマン
#define TEMPARTURE (300)

#define NUMBER (6193)      //粒子数
#define PARTICLE_MASS (1.0e-18)      //染色体粒子の質量
#define PARTICLE_RADIUS (1.0)     //粒子の半径
#define PI (M_PI)
#define MEMBRAIN_RADIUS_MINIMUM  (32.0) //膜の半径の最小値
#define INIT_DISTANCE ((PARTICLE_RADIUS + PARTICLE_RADIUS) * 0.8 )

#define ETA (2.0 * DIMENSION * PARTICLE_RADIUS * 0.000890)       //ノイズの強度
#define K_EXCLUDE (1.0)    //排除体積効果の強さ
#define K_BOND (1.0)    //ばね定数
#define K_BOND_2 (0.05)  //ひもの硬さ
#define DELTA (1.0e-11)  //刻み幅
#define PARTICLE_MYU (2.0 * DIMESION * PARTICLE_RADIUS * NANO * 0.000890)    //粘性抵抗の強さ
#define MEMBRAIN_EXCLUDE (1.0)     //膜との衝突
#define MEMBRAIN_EXCLUDE_SPB (1.0) //SPBとの衝突

#define SPB_RADIUS ( 3.0 )      //SPBの半径
#define SPB_MASS (9.0e-18)      //SPBの質量
#define SPB_MYU (2.0 * DIMESION * SPB_RADIUS * NANO * 0.000890)

double Nucleolus_circle_center[3];
double counter[90];

typedef struct particle {           //構造体の型宣言
    double position[DIMENSION];
} Particle;

Particle *ptr;

Particle spb;

double membrain_radius;

enum label{ X, Y, Z};

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

void read_coordinate( int time ){       //初期値設定
    
    int i, i_dummy;
    
    double d_dummy;
    
    Particle *part_1;
    
    char filename[128], dummy[256];
    
    char strs;
    
    FILE *fpr;
    
    sprintf (filename, "fission_result_%d.txt", time);
    
    if ((fpr = fopen(filename, "r")) == NULL){
        
        printf ("error \n");
        
        exit (1);
    }
    
    for (i=0; i<NUMBER; i++){
        
        /*fscanf(fpr, "%d %d %d %lf %lf %lf %lf %lf %lf %lf \n", &i_dummy, &part[i].chr_no, &part[i].particle_type,
               &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
               &part[i].velocity[X], &part[i].velocity[Y], &part[i].velocity[Z],
               &membrain_radius);*/
        fgets(dummy, 128, fpr);
        
        //printf ("%d %lf %lf %lf \n", i, part[i].position[X], part[i].position[Y], part[i].position[Z]);
        
    }
    
    fscanf (fpr, "%s %s %lf %lf %lf %lf %lf %lf %lf\n", dummy, dummy, &spb.position[X], &spb.position[Y], &spb.position[Z]
            , &d_dummy, &d_dummy, &d_dummy, &membrain_radius);

    fscanf (fpr, "%s %s %lf %lf %lf", dummy, dummy, &Nucleolus_circle_center[X], &Nucleolus_circle_center[Y], &Nucleolus_circle_center[Z]);
    
    fclose (fpr);
    
}
    
void write_coordinate (const int time, const double theta ) {
    
    FILE *fpw;
    
    char result[128];

    sprintf (result, "SN.txt");
    
    if ((fpw = fopen (result, "a")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    fprintf ( fpw, "%d %lf %lf %lf\n", time, theta, -cos(theta), sin(theta));
    
    fclose (fpw);
}

void write_coordinate_heat (const unsigned int num) {
    
    int i;
    
    FILE *fpw;
    
    char result[128];
    
    sprintf (result, "SPBheat.txt");
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }
    
    for ( i=0; i<90; i++) {
        
        //角度, 円周での座標, 頻度
        //fprintf ( fpw, "%d %lf %lf %d \n", i, -cos(PI * (89-i) / 90), sin(PI * (89-i) / 90), counter[(89-i)] );
        fprintf ( fpw, "%d %lf %lf %lf\n\n", 89-i, PI - PI * (89-i) / 90, 1.0, counter[(89-i)] / num);
    }
    for ( i=1; i<90; i++) {
        
        //fprintf ( fpw, "%d %lf %lf %d \n", i, -cos(PI * i / 90), -sin(PI * i / 90), counter[i] );
        fprintf ( fpw, "%d %lf %lf %lf\n\n", i, PI * i / 90 + PI , 1.0, counter[i] / num );
    }
    
    //fprintf (fpw, "\n");
    
    fclose (fpw);
}

void write_hist (const unsigned int num) {
    
    int i;
    
    FILE *fpw;
    
    char result[128];
    
    sprintf (result, "snhist_2.txt");
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \n error \n");
        
        exit (1);
    }

    for ( i=0; i<90; i++) {
        
        fprintf ( fpw, "%d %lf\n", i, counter[i] / num );
    }
    
    //fprintf (fpw, "\n");
    
    fclose (fpw);
}

int main ( int argc, char **argv ) {
    
    int start_number = atoi ( argv[1]);
    int calculate_number = atoi ( argv[2]);
    int width = atoi ( argv[3]);
    int time_width = atoi (argv[3]);
    
    unsigned int t, i_theta, m = 0, num = 0, time = 0;
    double cos, theta, origin[] = {0.0, 0.0, 0.0}, NtoS[DIMENSION], Nucleolus[DIMENSION];
    
    for ( t=0; t<90; t++){
        
        counter [t] = 0.0;
    }
    
    for ( t=0; t < calculate_number; t += width){
        
        theta = 0.0;
        
	/*
        if ( t + start_number > 270478 && m==0) {

            width /= 10;
            m = 1;
        }
        */

        time += 1;
        num += 1;
        
        read_coordinate ( t + start_number );
        
        Nucleolus[X] = Nucleolus_circle_center[X] / Euclid_norm (Nucleolus_circle_center, origin) * ( - 18.3 );
        Nucleolus[Y] = Nucleolus_circle_center[Y] / Euclid_norm (Nucleolus_circle_center, origin) * ( - 18.3 );
        Nucleolus[Z] = Nucleolus_circle_center[Z] / Euclid_norm (Nucleolus_circle_center, origin) * ( - 18.3 );
        
        //核小体 → SPB
        NtoS[X] = spb.position[X] - Nucleolus[X];
        NtoS[Y] = spb.position[Y] - Nucleolus[Y];
        NtoS[Z] = spb.position[Z] - Nucleolus[Z];

        theta = acos (inner_product ( NtoS , Nucleolus_circle_center ) / Euclid_norm (NtoS, origin)
        / Euclid_norm (Nucleolus_circle_center, origin) ) ;
        
        i_theta = (int) ( theta / PI * 180) ;
        
        counter[i_theta] += 1.0;
        
        write_coordinate ( time * time_width, theta / PI * 180);
        
        //printf ( " theta = %d \n", counter[i_theta]);
        
        //printf (" %lf %lf %lf \n", Euclid_norm (Nucleolus_circle_center, origin), Nucleolus[Y], Nucleolus[Z]);
    }
    
    write_coordinate_heat( num );
    write_hist (num);
    
    return (0);
}

