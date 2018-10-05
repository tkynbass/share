//
//  MSD_calculate.c
//  
//
//  Created by 高山雄揮 on 2017/07/18.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PARTICLE_NUMBER (6193)
/*
#define DATA_NUMBER (30000)
#define DATA_MINIMUM (0)
#define MSD_WIDTH (1000)
*/
#define DIMENSION (3)

typedef struct particle{
    
    double position[DIMENSION];
    
}Particle;

Particle part[PARTICLE_NUMBER];

enum label { X, Y, Z };

unsigned int DATA_NUMBER, DATA_MINIMUM, MSD_WIDTH;

void read_coordinate (){
    
    int i, j, i_dummy;
    
    double d_dummy;
    
    char filename[128], dummy[256];
    
    FILE *fpr;
    
    for(j=0; j<DATA_NUMBER; j++){
        
        sprintf (filename, "himo_result_%d.txt", DATA_MINIMUM);
        
        if ((fpr = fopen(filename, "r")) == NULL){
            
            printf ("error \n");
            
            exit (1);
        }
        
        for(i=0; i<PARTICLE_NUMBER; i++){
            
            if (i==15) {
                fscanf (fpr, "%s %d %lf %lf %lf %lf %lf %lf %lf", dummy, &i_dummy,
                    &part[i].position[X], &part[i].position[Y], &part[i].position[Z],
                    &d_dummy, &d_dummy, &d_dummy, &d_dummy);
            }
            
            else {
                
                fgets( dummy, 128, fpr);
            }
            
        }
        
        //fscanf(fpr, "%s %s %lf %lf %lf %lf\n", dummy, dummy, &spb[X], &spb[Y], &spb[Z], &membrain_radius);
        
        fclose (fpr);
    }
    
}

double calculate ( const unsigned int i) {
    
    int j, m = 0;
    
    double sum = 0.0;
    
    for ( j=1 ; j < DATA_NUMBER; j += i) {
        
        sum += (part[j].position[X] - part[j-i].position[X]) * (part[j].position[X] - part[j-i].position[X])
                + (part[j].position[Y] - part[j-i].position[Y]) * (part[j].position[Y] - part[j-i].position[Y])
                + (part[j].position[Z] - part[j-i].position[Z]) * (part[j].position[Z] - part[j-i].position[Z]);
        
        m++;
    }
    
    return sum/m;
}

void write_coordinate(const unsigned int t, const double msd) {
    
    unsigned int i;
    
    char filename[128], strs[128];
    
    FILE *fpw;
    
    Particle *part_1;
    
    sprintf (filename, "msd.txt");
    
    if ((fpw = fopen(filename, "a")) == NULL ) {
        
        printf ("write error \n ");
        
        exit (1);
    }

    fprintf (fpw, "%d %lf\n", t, msd);

    fclose (fpw);
}

int main ( int argc, char **argv) {
    
    unsigned int i;
    double msd;
    
    read_coordinate();
    
    DATA_MINIMUM = atoi(argv[1]);
    
    DATA_NUMBER = atoi(argv[2]);
    
    MSD_WIDTH = atof(argv[3]);
    
    for ( i=1 ;  (double) DATA_NUMBER / i * MSD_WIDTH > 100 ; i++) {
        
        msd = calculate ( DATA_MINIMUM + i * MSD_WIDTH );
        write_coordinate ( i, msd);
    }
    
    graph();
    
    return (0);
    
}
