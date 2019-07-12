// 八面体シミュレーション

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <omp.h>

#define PI (M_PI)

#define DIMENSION (3)       // 次元
#define LENGTH (7.0e-8)     // 長さ単位
#define DELTA (1.0e-4)

#define SIZE (6)    // 各構造体を構成する粒子数

#define MEMBRANE_AXIS_1 ( 1.889011e-6 / LENGTH )
#define MEMBRANE_AXIS_2 ( 0.85 * MEMBRANE_AXIS_1 )  // ~1.6
#define MEMBRANE_AXIS_3 ( 0.75 * MEMBRANE_AXIS_1 )  // ~1.4

#define NUCLEOLUS_AXIS_1 ( 1.1426593e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.9 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.8 * NUCLEOLUS_AXIS_1 )

#define WRITE_INTERVAL (1.0e+3)

const double MEM_POS [SIZE][DIMENSION] = {
    {MEMBRANE_AXIS_1, 0.0, 0.0}, { -MEMBRANE_AXIS_1, 0.0, 0.0},
    {0.0, MEMBRANE_AXIS_2, 0.0}, { 0.0, -MEMBRANE_AXIS_2, 0.0},
    {0.0, 0.0, MEMBRANE_AXIS_3}, { 0.0, 0.0, -MEMBRANE_AXIS_3}
};

const unsigned int AXIS[3][2] = {
    {0, 1}, {2, 3}, {4,5}
};

typedef struct nucleolus {           //核小体

    double position[DIMENSION];
    double position_new[DIMENSION];
    double velocity[DIMENSION];
    double force[DIMENSION];
    unsigned int keep_list [SIZE - 1];    // 八面体を保つために相互作用する粒子リスト
    double len_list [SIZE - 1]; // 自然長リスト (八面体を保つ)
    unsigned int mem_pt[SIZE];    // 核膜主成分端点(6点)との関数番号リスト 0~3 (nn, nf, fn, ff)
    unsigned int spb_pt;
} Nuc;

typedef struct SPB {        // SPB
    
    double position[DIMENSION];
    double position_new[DIMENSION];
    double velocity[DIMENSION];
    double force[DIMENSION];
    unsigned int mem_pt [SIZE];  // 核膜の各粒子とのポテンシャルの種類 0or1 (near or far)
    unsigned int nuc_pt [SIZE];  // 核小体の各粒子とのポテンシャルの種類 0or1 (near or far)
 }Spb;

enum dimension{ X, Y, Z};
enum rl{ right, left};
enum nf {near, far};
enum pt_type { nn, nf, fn, ff};
typedef enum ellipsoid_term {
    a1_r, a1_l, a2_r, a2_l, a3_r, a3_l
}ELLIP_TERM;


double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]) {
    
    double dist = 0.0;
    
    dist += (pos_1[X] - pos_2[X]) * (pos_1[X] - pos_2[X]);
    dist += (pos_1[Y] - pos_2[Y]) * (pos_1[Y] - pos_2[Y]);
    dist += (pos_1[Z] - pos_2[Z]) * (pos_1[Z] - pos_2[Z]);
    
    return (sqrt(dist));
}

void Secure_main_memory (Nuc **nuc, Spb **spb) {   // メモリ確保 //
    
    if ( (*nuc = (Nuc *)malloc(SIZE * sizeof(Nuc))) == NULL) {
        
        printf("\n error : can not secure the memory of nuc \n");
        exit(1);
    }
    
    if ( (*spb = (Spb *)malloc(1 * sizeof (Spb))) == NULL) {
        
        printf("\n error : can not secure the memory of spb \n");
        exit(1);
    }
}

void StructInitilization (Nuc *nuc, Spb *spb) {
    
    Nuc *ncl;
    unsigned int loop, loop2;
    double gravity[] = { -0.3 * MEMBRANE_AXIS_1, 1.0, 1.0};
    
    // 核小体　初期位置 //
    double init_pos [SIZE][DIMENSION] = {
        { gravity[X], gravity[Y] + NUCLEOLUS_AXIS_1, gravity[Z]}, { gravity[X], gravity[Y] - NUCLEOLUS_AXIS_1, gravity[Z]},
        { gravity[X] + NUCLEOLUS_AXIS_2, gravity[Y], gravity[Z]}, { gravity[X] - NUCLEOLUS_AXIS_2, gravity[Y], gravity[Z]},
        { gravity[X], gravity[Y], gravity[Z] + NUCLEOLUS_AXIS_3}, { gravity[X], gravity[Y], gravity[Z] - NUCLEOLUS_AXIS_3}
    };
    
    for (loop = 0; loop < SIZE; loop++) {
        
        ncl = &nuc[loop];
        for (unsigned int dim = 0; dim < DIMENSION; dim++) ncl->position [dim] = init_pos[loop][dim];
        
//        printf (" %4.2f %4.2f %4.2f\n", ncl->position[X], ncl->position[Y], ncl->position[Z]);
    }
    
    spb->position[X] = MEMBRANE_AXIS_1;
    spb->position[Y] = 0.0;
    spb->position[Z] = 0.0;
    
    // 核小体形状保存 自然長求める
    for ( loop = 0; loop < SIZE; loop++) {
        
        ncl = &nuc[loop];
        for (loop2 = 0; loop2 < SIZE-1; loop2++){
            
            if (loop2 < loop) ncl->keep_list[loop2] = loop2;
            else ncl->keep_list[loop2] = loop2 + 1;
            
            ncl->len_list[loop2] = Euclid_norm (ncl->position, nuc[ ncl->keep_list[loop2] ].position);
            
            printf ("%4.2f ", ncl->len_list[loop2]);
        }
        printf ("\n");
    }
}

void Def_max_min (const double dist_list[4], unsigned int *max, unsigned int *min) {
    
    for (unsigned int loop = 1; loop < 4; loop++) {
        
        if ( dist_list[*max] < dist_list[loop] ) *max = loop;
        if ( dist_list[*min] > dist_list[loop] ) *min = loop;
    }
}

void Calculation (const unsigned int mitigation, Nuc *nuc, Spb *spb) {
    
    unsigned int n_ax, m_ax, dim, nf_list[2][2], loop, mem_r, mem_l, max, min;
    double mid_point[3], dist_list[4];
    Nuc *nuc_r, *nuc_l;
    
//     加えるポテンシャル種類の判別
    for (n_ax = 0; n_ax <3; n_ax++) {
        
        nuc_r = &nuc[AXIS[n_ax][right]];
        nuc_l = &nuc[AXIS[n_ax][left]];

        // 重心 (軸の中点)
        for (dim = 0; dim < DIMENSION; dim++) mid_point[dim] = (nuc_r->position[dim] + nuc_l->position[dim]) / 2.0;
        
        for (m_ax = 0; m_ax < 3; m_ax++) {
            
            max = 0;
            min = 0;
            
            mem_r = AXIS[m_ax][right];
            mem_l = AXIS[m_ax][left];
            
            dist_list[0] = Euclid_norm (nuc_r->position, MEM_POS[mem_r]);
            dist_list[1] = Euclid_norm (nuc_r->position, MEM_POS[mem_l]);
            
            dist_list[2] = Euclid_norm (nuc_l->position, MEM_POS[mem_r]);
            dist_list[3] = Euclid_norm (nuc_l->position, MEM_POS[mem_l]);
            
            Def_max_min (dist_list, &max, &min);
            
            switch (min) {
                case 0:
                    nuc_r->mem_pt[mem_r] = nn;
                    nuc_l->mem_pt[mem_r] = nf;
                    break;
                case 1:
                    nuc_r->mem_pt[mem_l] = nn;
                    nuc_l->mem_pt[mem_l] = nf;
                    break;
                case 2:
                    nuc_l->mem_pt[mem_r] = nn;
                    nuc_r->mem_pt[mem_r] = nf;
                    break;
                case 3:
                    nuc_l->mem_pt[mem_l] = nn;
                    nuc_r->mem_pt[mem_l] = nf;
                    break;
                default:
                    printf ("\t mem potential def. error \n");
                    break;
            }
            
            switch (max) {
                case 0:
                    nuc_r->mem_pt[mem_r] = ff;
                    nuc_l->mem_pt[mem_r] = fn;
                    break;
                case 1:
                    nuc_r->mem_pt[mem_l] = ff;
                    nuc_l->mem_pt[mem_l] = fn;
                    break;
                case 2:
                    nuc_l->mem_pt[mem_r] = ff;
                    nuc_r->mem_pt[mem_r] = fn;
                    break;
                case 3:
                    nuc_l->mem_pt[mem_l] = ff;
                    nuc_r->mem_pt[mem_l] = fn;
                    break;
                default:
                    printf ("\t mem potential def. error \n");
                    break;
            }
        }
        
        for (loop = 0; loop < SIZE; loop++) printf ("%d ", nuc[5].mem_pt[loop]);
        printf ("\n");
        fflush (stdout);
    }
    
//    for (int loop=0; loop<SIZE; loop++) {
//
//        for (int loop2=0; loop2<SIZE; loop2++) {
//
//            printf ("%d ", nuc[loop].mem_pt[loop2]);
//        }
//        printf ("\n");
//    }
    
}

int main ( int argc, char **argv ){
    
    Nuc *nuc;
    Spb *spb;
//    unsigned int time, mitigation, calc_max = atoi (argv[1]);
    unsigned int mitigation=0;

    Secure_main_memory (&nuc, &spb);    // 構造体のメモリ確保
    
    // 構造体の初期化
    StructInitilization (nuc, spb);
    // 計算
//    for ( time = 1; time <= calc_max; time++) {
//
//        for ( mitigation = 0; mitigation < WRITE_INTERVAL; mitigation++) {
//
//            Calculation (mitigation, nuc, spb);
//        }
//    }
    Calculation (mitigation, nuc, spb);
    
    // 書き込み
    
    return (0);
}




