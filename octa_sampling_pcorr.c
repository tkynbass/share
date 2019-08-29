// 八面体シミュレーション

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "/home/gensho/tkym/dSFMT/dSFMT.h"
#include <time.h>
//#include <omp.h>

#define PI (M_PI)

#define DIMENSION (3)       // 次元
#define LENGTH (7.0e-8)     // 長さ単位
#define DELTA (1.0e-2)

#define KBT ( 1.38064852e-23 / LENGTH / LENGTH ) //ボルツマン
#define TEMPARTURE ( 300 )

#define SIZE (6)    // 各構造体を構成する粒子数

#define MEMBRANE_AXIS_1 ( 1.889011e-6 / LENGTH )
#define MEMBRANE_AXIS_2 ( 0.85 * MEMBRANE_AXIS_1 )  // ~1.6
#define MEMBRANE_AXIS_3 ( 0.75 * MEMBRANE_AXIS_1 )  // ~1.4

#define NUCLEOLUS_AXIS_1 ( 1.1426593e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.9 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.8 * NUCLEOLUS_AXIS_1 )

#define WRITE_INTERVAL (1.0e+3)

#define MEMBRANE_EXCLUDE (1.0)
#define K_KEEP (1.0e+1)
#define K_EXPERIENCE (5.0e+1)

// SPBのノイズ用
#define DIFFUSION (5.0e-3)
#define KINEMATIC_MYU (0.000890)
//#define MYU ( 2.0 * DIMENSION * PI * 1.0 * 0.000890)
#define MYU (1.0)
#define INV_MYU (1.0 / MYU)

const double ORIGIN[] = {0.0, 0.0, 0.0};

const double MEM_POS [SIZE][DIMENSION] = {
    {MEMBRANE_AXIS_1, 0.0, 0.0}, { -MEMBRANE_AXIS_1, 0.0, 0.0},
    {0.0, MEMBRANE_AXIS_2, 0.0}, { 0.0, -MEMBRANE_AXIS_2, 0.0},
    {0.0, 0.0, MEMBRANE_AXIS_3}, { 0.0, 0.0, -MEMBRANE_AXIS_3}
};

// 軸番号と左右→粒子番号
const unsigned int AXIS[3][2] = {
    {0, 1}, {2, 3}, {4,5}
};

// 粒子番号　→　軸番号
const unsigned int AX_NUM [6] = {
    0, 0, 1, 1, 2, 2
};

typedef struct nucleolus {           //核小体

    double position[DIMENSION];
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

// サンプリング　初期配置
void RandomSetting (Nuc *nuc, Spb *spb, dsfmt_t *dsfmt) {
    
    unsigned int lp, dim;
    Nuc *nuc_r, *nuc_l;
    
    const double gravity[DIMENSION] = { MEMBRANE_AXIS_1 * dsfmt_genrand_open_open(dsfmt),
                        MEMBRANE_AXIS_2 * dsfmt_genrand_open_open (dsfmt),
                        MEMBRANE_AXIS_3 * dsfmt_genrand_open_open (dsfmt)};
    
    const double eta = PI * dsfmt_genrand_close_open (dsfmt);
    const double theta = 2 * PI * dsfmt_genrand_close_open (dsfmt);
    const double phi = 0.5 * PI * dsfmt_genrand_close_open (dsfmt);
    
    spb->position[X] = 2.0 * MEMBRANE_AXIS_1 * dsfmt_genrand_open_open (dsfmt) - MEMBRANE_AXIS_1;
    spb->position[Y] = 2.0 * MEMBRANE_AXIS_2 * dsfmt_genrand_open_open (dsfmt) - MEMBRANE_AXIS_2;
    spb->position[Z] = 2.0 * MEMBRANE_AXIS_3 * dsfmt_genrand_open_open (dsfmt) - MEMBRANE_AXIS_3;

    double nav[3][DIMENSION] = {{ cos (phi) * cos (theta) * NUCLEOLUS_AXIS_1,
                            cos (phi) * sin (theta) * NUCLEOLUS_AXIS_1,
                            sin (phi) * NUCLEOLUS_AXIS_1},
                            { (-cos(theta) * sin(eta) * sin(phi) - cos(eta) * sin(theta) ) * NUCLEOLUS_AXIS_2,
                            ( cos(theta) * cos(eta) - sin(eta) * sin(phi) * sin(theta)) * NUCLEOLUS_AXIS_2,
                            sin(eta) * cos(phi) * NUCLEOLUS_AXIS_2},
                            { ( -cos(theta) * cos(eta) * sin(phi) + sin(eta) * sin(theta) ) * NUCLEOLUS_AXIS_3,
                            ( -cos(theta) * sin(eta) - cos(eta) * sin(phi) * sin(theta) ) * NUCLEOLUS_AXIS_3,
                            cos(eta) * cos(phi) * NUCLEOLUS_AXIS_3}};
    
    for (lp = 0; lp < 3; lp++) {
        
        nuc_r = &nuc[ AXIS[lp][right] ];
        nuc_l = &nuc[ AXIS[lp][left] ];
        
        nuc_r->position[X] = gravity[X] + nav[lp][X];
        nuc_r->position[Y] = gravity[Y] + nav[lp][Y];
        nuc_r->position[Z] = gravity[Z] + nav[lp][Z];
        
        nuc_l->position[X] = gravity[X] - nav[lp][X];
        nuc_l->position[Y] = gravity[Y] - nav[lp][Y];
        nuc_l->position[Z] = gravity[Z] - nav[lp][Z];
    }
}

void StructInitilization (Nuc *nuc, Spb *spb, dsfmt_t *dsfmt) {
    
    Nuc *ncl;
    unsigned int loop, loop2;

    RandomSetting (nuc, spb, dsfmt);
    
    // 核小体形状保存 自然長求める
    for ( loop = 0; loop < SIZE; loop++) {
        
        ncl = &nuc[loop];
        for (loop2 = 0; loop2 < SIZE-1; loop2++){
            
            if (loop2 < loop) ncl->keep_list[loop2] = loop2;
            else ncl->keep_list[loop2] = loop2 + 1;
            
            ncl->len_list[loop2] = Euclid_norm (ncl->position, nuc[ ncl->keep_list[loop2] ].position);
            
//            printf ("%4.2f ", ncl->len_list[loop2]);
        }
//        printf ("\n");
    }
}

void Def_NucMem_pt (Nuc *nuc) {
    
    unsigned int dim, n_ax, m_ax, mem_r, mem_l, mem_Near, mem_Far;
    double mid_pos[DIMENSION];
    Nuc *nuc_r, *nuc_l;
    
    //     加えるポテンシャル種類の判別
    for (n_ax = 0; n_ax <3; n_ax++) {
        
        nuc_r = &nuc[AXIS[n_ax][right]];
        nuc_l = &nuc[AXIS[n_ax][left]];
        
        for (dim = 0; dim < DIMENSION; dim++) mid_pos [dim] = ( nuc_r->position[dim] + nuc_l->position[dim] ) / 2.0;
        
        for (m_ax = 0; m_ax < 3; m_ax++) {
            
            mem_r = AXIS[m_ax][right];
            mem_l = AXIS[m_ax][left];
            
//            Near point　と Far point の判別
            if (Euclid_norm (MEM_POS[mem_r], mid_pos) <= Euclid_norm (MEM_POS[mem_l], mid_pos)) {
                
                mem_Near = mem_r;
                mem_Far = mem_l;
            }
            else {
                
                mem_Near = mem_l;
                mem_Far = mem_r;
            }
            
//            nn, nfの判別
            if (Euclid_norm (nuc_r->position, MEM_POS[mem_Near]) <= Euclid_norm (nuc_l->position, MEM_POS[mem_Near])){
                
                nuc_r->mem_pt [mem_Near] = nn;
                nuc_l->mem_pt [mem_Near] = nf;
            }
            else {
                
                nuc_r->mem_pt [mem_Near] = nf;
                nuc_l->mem_pt [mem_Near] = nn;
            }
            
//            fn, ffの判別
            if (Euclid_norm (nuc_r->position, MEM_POS[mem_Far]) <= Euclid_norm (nuc_l->position, MEM_POS[mem_Far])){
                
                nuc_r->mem_pt [mem_Far] = fn;
                nuc_l->mem_pt [mem_Far] = ff;
            }
            else {
                
                nuc_r->mem_pt [mem_Far] = ff;
                nuc_l->mem_pt [mem_Far] = fn;
            }

        }
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

void Def_SpbMem_pt (Spb *spb) {
    
    unsigned int m_ax, mem_r, mem_l;
    
    for ( m_ax = 0; m_ax < 3; m_ax++) {
        
        mem_r = AXIS[m_ax][right];
        mem_l = AXIS[m_ax][left];
        
        if (Euclid_norm (spb->position, MEM_POS [mem_r]) <= Euclid_norm (spb->position, MEM_POS [mem_l])) {
            
            spb->mem_pt[mem_r] = near;
            spb->mem_pt[mem_l] = far;
        }
        else {
            
            spb->mem_pt[mem_r] = far;
            spb->mem_pt[mem_l] = near;
        }
    }
    
//    for (unsigned int loop = 0; loop<SIZE; loop++)  printf ("%d ", spb->mem_pt[loop]);
//    printf ("\n");
}

void Def_NucSpb_pt (Nuc *nuc, Spb *spb) {
    
    unsigned int n_ax;
    Nuc *nuc_r, *nuc_l;
    
    for ( n_ax = 0; n_ax < 3; n_ax++) {
        
//        主成分両端粒子を参照
        nuc_r = &nuc[ AXIS[n_ax][right] ];
        nuc_l = &nuc[ AXIS[n_ax][left] ];
        
        if ( Euclid_norm (nuc_r->position, spb->position) <= Euclid_norm (nuc_l->position, spb->position)) {
            
            nuc_r->spb_pt = near;
            nuc_l->spb_pt = far;
            
            spb->nuc_pt[ AXIS[n_ax][right] ] = near;
            spb->nuc_pt[ AXIS[n_ax][left] ] = far;
        }
        else {
            
            nuc_r->spb_pt = far;
            nuc_l->spb_pt = near;
            
            spb->nuc_pt[ AXIS[n_ax][right] ] = far;
            spb->nuc_pt[ AXIS[n_ax][left] ] = near;
        }
    }
//    for (unsigned int loop = 0; loop<SIZE; loop++)  printf ("%d ", spb->nuc_pt[loop]);
//    printf ("\n");
    
}

void Keep_ellipsoid (Nuc *nuc) {
    
    unsigned int lp, lp2;
    double f, dist;
    Nuc *ncl1, *ncl2;
    
    for ( lp = 0; lp < SIZE; lp++) {
        
        ncl1 = &nuc[lp];
        for ( lp2 = 0; lp2 < SIZE-1; lp2++ ){
            
            ncl2 = &nuc [ncl1->keep_list[lp2]];
            
            dist = Euclid_norm (ncl1->position, ncl2->position);
            f = - K_KEEP * (dist - ncl1->len_list [lp2]) / dist;
            
            ncl1->force[X] += f * (ncl1->position[X] - ncl2->position[X]);
            ncl1->force[Y] += f * (ncl1->position[Y] - ncl2->position[Y]);
            ncl1->force[Z] += f * (ncl1->position[Z] - ncl2->position[Z]);
        }
    }
    
}

void TermDIst_NucMem (Nuc *nuc, const char option) {  // 核小体-核膜間の相互作用
    
    unsigned int lp, lp2, type_lp, nuc_lp, mem_lp;
    static double para_list [4][3][3][6], pcorr_ratio[4][3][3];
    double *par, dist, f, exp1, exp2, ratio;
    char file_name[128], dummy[128], *type_list[] = {"nn", "nf", "fn", "ff"};
    Nuc *ncl;
    
//    char directory[] = "/Users/tkym/Desktop/Imaris/analysis/pcc_Cut11-Gar2_wt";
    
    FILE *fpr;
    
    if ( option == 'c') {
        
        for (nuc_lp = 0; nuc_lp < SIZE; nuc_lp++) {
            
            ncl = &nuc [nuc_lp];
            for (mem_lp = 0; mem_lp < SIZE; mem_lp++) {
                
                // type = ncl->mem_pt[mem_pt]
                par = para_list [ncl->mem_pt [mem_lp]][AX_NUM[nuc_lp]][AX_NUM[mem_lp]];
                ratio = pcorr_ratio [ncl->mem_pt [mem_lp]][AX_NUM[nuc_lp]][AX_NUM[mem_lp]];
                
                dist = Euclid_norm (ncl->position, MEM_POS [mem_lp]);
                exp1 = exp ( -0.5 * par[2]*par[2] * (dist - par[1])*(dist - par[1]));
                exp2 = exp ( -0.5 * par[5]*par[5] * (dist - par[4])*(dist - par[4]));
                
                
                f = - K_EXPERIENCE * ratio * ( par[0] * (dist - par[1]) * exp1 + par[3] * (dist - par[4]) *exp2 )
                    / ( dist * (exp1 + exp2));
                
                ncl->force[X] += f * (ncl->position[X] - MEM_POS[mem_lp][X]);
                ncl->force[Y] += f * (ncl->position[Y] - MEM_POS[mem_lp][Y]);
                ncl->force[Z] += f * (ncl->position[Z] - MEM_POS[mem_lp][Z]);
            }
        }
        
//        for (lp = 0; lp < DIMENSION; lp++) printf ("%lf\n", nuc[5].force[lp]);
    }
    else if (option == 's') { //    パラメータの読み込み
        
        for (type_lp = 0; type_lp <4; type_lp++) {
            
            sprintf (file_name, "../subdata/pote_nm_%s.txt", type_list[type_lp]);
            if ( (fpr = fopen (file_name, "r")) == NULL) {
                
                printf ("\t Cannot open %s \n", file_name);
                exit(1);
            }
            
            fgets (dummy, 256, fpr);    // column読み捨て
            
            for (nuc_lp = 0; nuc_lp < 3; nuc_lp++) {
                
                for (mem_lp = 0; mem_lp < 3; mem_lp++) {
                    
                    par = para_list [type_lp][nuc_lp][mem_lp];
                    fscanf (fpr, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", dummy,
                            &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
                    
                    par[0] /= par[2] * par[2] * par[2];
                    par[2] = 1.0 / par[2];
                    
                    par[3] /= par[5] * par[5] * par[5];
                    par[5] = 1.0 / par[5];
                }
            }
            
            fclose (fpr);
        }
        
        sprintf (file_name, "../subdata/mn_vel_pcorr.txt");
        if ( (fpr = fopen (file_name, "r")) == NULL) {
            
            printf ("\t Cannot open %s \n", file_name);
            exit(1);
        }
        fgets (dummy, 256, fpr);
        
        for (type_lp = 0; type_lp < 3; type_lp++) {
            
            for (nuc_lp = 0; nuc_lp < 3; nuc_lp++) {
                
                fscanf (fpr, "%s %lf %lf %lf\n", dummy, &pcorr_ratio[type_lp][nuc_lp][0],
                        &pcorr_ratio[type_lp][nuc_lp][1], &pcorr_ratio[type_lp][nuc_lp][2]);
            }
        }
        fclose (fpr);
        
    }
//    if (option == 'c') for (lp = 0; lp < 6; lp++) printf ("%lf\n", para_list[fn][0][0][lp]);
}
void TermDIst_SpbMem ( Spb *spb, const char option) {  // SPB-核膜間の相互作用
    
    unsigned int lp, type_lp, mem_lp;
    static double para_list [2][3][6], pcorr_ratio[2][3];
    double *par, dist, f, exp1, exp2, ratio;
    char file_name[128], dummy[128], *type_list[] = {"near", "far"};
    
//    char directory[] = "/Users/tkym/Desktop/Imaris/analysis/Cut11-Sid4";
    
    FILE *fpr;
    
    // calculate
    if ( option == 'c') {
        
        for (mem_lp = 0; mem_lp < SIZE; mem_lp++) {
            
            par = para_list [spb->mem_pt [mem_lp]][AX_NUM[mem_lp]];
            ratio = pcorr_ratio [spb->mem_pt [mem_lp]][AX_NUM[mem_lp]];
            
            dist = Euclid_norm (spb->position, MEM_POS [mem_lp]);
            exp1 = exp ( -0.5 * par[2]*par[2] * (dist - par[1])*(dist - par[1]));
            exp2 = exp ( -0.5 * par[5]*par[5] * (dist - par[4])*(dist - par[4]));
            
            f = - K_EXPERIENCE * ratio * ( par[0] * (dist - par[1]) * exp1 + par[3] * (dist - par[4]) *exp2 )
            / ( dist * (exp1 + exp2));
            
            spb->force[X] += f * (spb->position[X] - MEM_POS[mem_lp][X]);
            spb->force[Y] += f * (spb->position[Y] - MEM_POS[mem_lp][Y]);
            spb->force[Z] += f * (spb->position[Z] - MEM_POS[mem_lp][Z]);
        }
        
//        for (lp = 0; lp < DIMENSION; lp++) printf ("%lf\n", spb->force[lp]);
    }// setting
    else if (option == 's') { //    パラメータの読み込み
        
        for (type_lp = 0; type_lp < 2; type_lp++) {
            
            sprintf (file_name, "../subdata/pote_sm_%s.txt", type_list[type_lp]);
            if ( (fpr = fopen (file_name, "r")) == NULL) {
                
                printf ("\t Cannot open %s \n", file_name);
                exit(1);
            }
            
            fgets (dummy, 256, fpr);    // column読み捨て
                
            for (mem_lp = 0; mem_lp < 3; mem_lp++) {
                
                par = para_list [type_lp][mem_lp];
                fscanf (fpr, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", dummy,
                        &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
                
                par[0] /= par[2] * par[2] * par[2];
                par[2] = 1.0 / par[2];
                
                par[3] /= par[5] * par[5] * par[5];
                par[5] = 1.0 / par[5];
            }
            
            fclose (fpr);
        }
//        for (lp = 0; lp < 6; lp++) printf ("%lf\n", para_list[far][0][lp]);
        
        sprintf (file_name, "../subdata/sm_vel_pcorr.txt");
        if ( (fpr = fopen (file_name, "r")) == NULL) {
            
            printf ("\t Cannot open %s \n", file_name);
            exit(1);
        }
        
        fgets (dummy, 256, fpr);
        
        for (mem_lp = 0; mem_lp < 3; mem_lp++) {
            
            for (type_lp = 0; type_lp < 2; type_lp++){
                
                fscanf (fpr, "%s %lf\n", dummy, &pcorr_ratio[type_lp][mem_lp]);
            }
        }
        fclose (fpr);
    }
    //    if (option == 'c') for (lp = 0; lp < 6; lp++) printf ("%lf\n", para_list[fn][0][0][lp]);
}

void TermDIst_NucSpb (Nuc *nuc, Spb *spb, const char option) {  // SPB-核膜間の相互作用
    
    unsigned int lp, type_lp, nuc_lp;
    static double para_list [2][3][6], pcorr_ratio[2][3];
    double *par, dist, f, exp1, exp2, ratio;
    char file_name[128], dummy[128], *type_list[] = {"near", "far"};
    Nuc *ncl;
    
//    char directory[] = "/Users/tkym/Desktop/Imaris/analysis/06-28decon";
    
    FILE *fpr;
    
    if ( option == 'c') {
        
        for (nuc_lp = 0; nuc_lp < SIZE; nuc_lp++) {
            
            ncl = &nuc [nuc_lp];
            
            par = para_list [spb->nuc_pt [nuc_lp]][AX_NUM[nuc_lp]];
            ratio = pcorr_ratio [spb->nuc_pt [nuc_lp]][AX_NUM[nuc_lp]];
            
            dist = Euclid_norm (spb->position, ncl->position);
            exp1 = exp ( -0.5 * par[2]*par[2] * (dist - par[1])*(dist - par[1]));
            exp2 = exp ( -0.5 * par[5]*par[5] * (dist - par[4])*(dist - par[4]));
            
            f = - K_EXPERIENCE * ratio * ( par[0] * (dist - par[1]) * exp1 + par[3] * (dist - par[4]) *exp2 )
            / ( dist * (exp1 + exp2));
            
//            printf ("%c f = %lf\n", option, f);
            
            spb->force[X] += f * (spb->position[X] - ncl->position[X]);
            spb->force[Y] += f * (spb->position[Y] - ncl->position[Y]);
            spb->force[Z] += f * (spb->position[Z] - ncl->position[Z]);
            
            ncl->force[X] += f * (ncl->position[X] - spb->position[X]);
            ncl->force[Y] += f * (ncl->position[Y] - spb->position[Y]);
            ncl->force[Z] += f * (ncl->position[Z] - spb->position[Z]);
        }
    }
    else if (option == 's') { //    パラメータの読み込み
        
        for (type_lp = 0; type_lp < 2; type_lp++) {
            
            sprintf (file_name, "../subdata/pote_sn_%s.txt", type_list[type_lp]);
            if ( (fpr = fopen (file_name, "r")) == NULL) {
                
                printf ("\t Cannot open %s \n", file_name);
                exit(1);
            }
            
            fgets (dummy, 256, fpr);    // column読み捨て
            
            for (nuc_lp = 0; nuc_lp < 3; nuc_lp++) {
                
                par = para_list [type_lp][nuc_lp];
                fscanf (fpr, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", dummy,
                        &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
                
                par[0] /= par[2] * par[2] * par[2];
                par[2] = 1.0 / par[2];
                
                par[3] /= par[5] * par[5] * par[5];
                par[5] = 1.0 / par[5];
            }
        }
        
        fclose(fpr);
    
//        for (lp = 0; lp < 6; lp++) printf ("%lf\n", para_list[far][0][lp]);
        
        sprintf (file_name, "../subdata/sn_vel_pcorr.txt");
        if ( (fpr = fopen (file_name, "r")) == NULL) {
            
            printf ("\t Cannot open %s \n", file_name);
            exit(1);
        }
        
        fgets (dummy, 256, fpr);
        
        for (nuc_lp = 0; nuc_lp < 3; nuc_lp++) {
            
            for (type_lp = 0; type_lp < 2; type_lp++){
                
                fscanf (fpr, "%s %lf\n", dummy, &pcorr_ratio[type_lp][nuc_lp]);
            }
        }
        fclose (fpr);
        
    }
    //    if (option == 'c') for (lp = 0; lp < 6; lp++) printf ("%lf\n", para_list[fn][0][0][lp]);
}

void Membrane_interaction ( const double pos[DIMENSION], double force[DIMENSION], char interaction_type /* F: fix, E: exclude */) {
    
    double dist = Euclid_norm (pos, ORIGIN);
    
    double ellipsoid_dist = pos[X] * pos[X] / ( MEMBRANE_AXIS_1 * MEMBRANE_AXIS_1 )
    + pos[Y] * pos[Y] / ( MEMBRANE_AXIS_2 * MEMBRANE_AXIS_2 )
    + pos[Z] * pos[Z] / ( MEMBRANE_AXIS_3 * MEMBRANE_AXIS_3 );
    
    if ( interaction_type == 'F' || ellipsoid_dist - 1 > 0 ) {
        
        // 法線ベクトル
        double normal_vector[] = { 2.0 * pos[X] / ( MEMBRANE_AXIS_1 * MEMBRANE_AXIS_1),
            2.0 * pos[Y] / ( MEMBRANE_AXIS_2 * MEMBRANE_AXIS_2),
            2.0 * pos[Z] / ( MEMBRANE_AXIS_3 * MEMBRANE_AXIS_3) };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * ( pos[X] * normal_vector[X]
                                                                  + pos[Y] * normal_vector[Y]
                                                                  + pos[Z] * normal_vector[Z]) / normal_vector_norm;
        
        force[X] += f * normal_vector[X];
        force[Y] += f * normal_vector[Y];
        force[Z] += f * normal_vector[Z];
    }
}

void Calculation (const unsigned int mitigation, Nuc *nuc, Spb *spb) {
    
    unsigned int lp, dim;
    Nuc *ncl;
    
    // 力の初期化
    for ( lp = 0; lp < SIZE; lp++) {
        ncl = &nuc[lp];
        for (dim = 0; dim < DIMENSION; dim++) ncl->force[dim] = 0.0;
        
    }
    
    for (dim = 0; dim < DIMENSION; dim++) spb->force[dim] = 0.0;
    
    Def_NucMem_pt (nuc);
    Def_SpbMem_pt (spb);
    Def_NucSpb_pt (nuc, spb);
    
    Keep_ellipsoid (nuc);
    
    TermDIst_NucMem (nuc, 'c');
    TermDIst_SpbMem (spb, 'c');
    TermDIst_NucSpb (nuc, spb, 'c');
    
    // 核膜との排除体積　および　核膜上に固定
    for (lp = 0; lp < SIZE; lp++) Membrane_interaction (nuc[lp].position, nuc[lp].force, 'E');
    Membrane_interaction (spb->position, spb->force, 'F');
    
    for (lp = 0; lp < SIZE; lp++) {
        
        ncl = &nuc[lp];
        
        ncl->position[X] += INV_MYU * DELTA * ncl->force[X];
        ncl->position[Y] += INV_MYU * DELTA * ncl->force[Y];
        ncl->position[Z] += INV_MYU * DELTA * ncl->force[Z];
    }
    
    spb->position[X] += INV_MYU * DELTA * spb->force[X];
    spb->position[Y] += INV_MYU * DELTA * spb->force[Y];
    spb->position[Z] += INV_MYU * DELTA * spb->force[Z];
}

double Sum_force (Nuc *nuc, Spb *spb) {
    
    unsigned int lp, dim;
    double sum_force = 0.0;

    for (lp = 0; lp < SIZE; lp++) {
        
        for (dim = 0; dim < DIMENSION; dim++) sum_force += abs (nuc[lp].force[dim]);
    }
    for (dim = 0; dim < DIMENSION; dim++) sum_force += abs (spb->force[dim]);
    
//    printf ("\n\t total force = %3.2e\n", sum_force);
    
    return sum_force;
}

void Write_coordinate (Nuc *nuc, Spb *spb, const unsigned int step, const unsigned int sample_no) {
    
    unsigned int lp;
    char filename[128];
    Nuc *ncl;
    
    FILE *fpw;
    
    sprintf (filename, "%d/result_%d.txt", sample_no, step);
    
    if ( (fpw = fopen (filename, "w")) == NULL) {
        
        printf ("\t Cannot open result file.\n");
        exit(1);
    }
    
    for (lp = 0; lp < SIZE; lp++) {
        
        ncl = &nuc [lp];
        fprintf (fpw, "%d %lf %lf %lf\n", lp, ncl->position[X], ncl->position[Y], ncl->position[Z]);
    }
    fprintf (fpw, "spb %lf %lf %lf\n", spb->position[X], spb->position[Y], spb->position[Z]);
    
    fclose (fpw);
}

void Write_row_data (char *filename, Nuc *nuc, Spb *spb, const unsigned int sample_no) {
    
    unsigned int lp;
//    char filename[128];
    
    Nuc *ncl;
    FILE *fpw;
    
    if ( (fpw = fopen (filename, "a")) == NULL) {
        
        printf ("\t Cannot open result file.\n");
        exit(1);
    }
    
    fprintf (fpw, "%d ", sample_no);
    for (lp = 0; lp < SIZE; lp++) {
        
        ncl = &nuc[lp];
        fprintf (fpw, "%6.2f %6.2f %6.2f ", ncl->position[X], ncl->position[Y], ncl->position[Z]);
    }
    fprintf (fpw, "%6.2f %6.2f %6.2f\n", spb->position[X], spb->position[Y], spb->position[Z]);
    
    fclose (fpw);
}

int main ( int argc, char **argv ){
    
    Nuc *nuc;
    Spb *spb;
    unsigned int sample_no = atoi (argv[1]);
    unsigned int step, mitigation, calc_max = atoi (argv[2]);
    
    Secure_main_memory (&nuc, &spb);    // 構造体のメモリ確保
    
    //dSFMT
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, sample_no);
    
    // 構造体の初期化
    StructInitilization (nuc, spb, &dsfmt);
    
    Write_row_data ("init.txt", nuc, spb, sample_no);
    Write_coordinate (nuc, spb, 0, sample_no);
    printf ("\t DELTA = %2.1e, WRITE_INTERVAL = %2.1e \r", DELTA, WRITE_INTERVAL);
    
    // 混合ガウスポテンシャルのパラメータ読み込み
    TermDIst_NucMem (NULL, 's');
    TermDIst_SpbMem (NULL, 's');
    TermDIst_NucSpb (NULL, NULL, 's');
    // 計算
    for ( step = 1; step <= calc_max; step++) {

//        printf ("\t Calculating now ...  step = %d\r", step);
        fflush (stdout);
        
        for ( mitigation = 0; mitigation < WRITE_INTERVAL; mitigation++) {

            Calculation (mitigation, nuc, spb);
        }
        Write_coordinate (nuc, spb, step, sample_no);
        
//        if (Sum_force (nuc, spb) < 1.0e-5) break;
    }
    
    Write_row_data ("result.txt", nuc, spb, sample_no);

//    Sum_force (nuc, spb);
    fflush (stdout);
        
    free (nuc);
    free (spb);
    
    return (0);
}




