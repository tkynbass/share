
//
//
//  octa_sn_fix.c で得られた核小体・核膜・SPBの位置関係に
//  ひもを入れて緩和させる
//  各染色体の端点(テロメア)をランダムに設定 →粒子径
//  5kbpで粗視化
//  遺伝研で流す用 隠れマルコフからのポテンシャルあり
//  歪み(隣接粒子との距離の自然長からのずれ）の大きい粒子の隠れマルコフ状態を遷移させる.


// 2019/12/2
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../dSFMT/dSFMT.h"
#include <omp.h>

#define DIMENSION ( 3 ) //次元
#define LENGTH ( 9.66e-8 )   //長さの単位 (粒子径)
#define TRANS_LENGTH (1.0e-6 / LENGTH)

#define KBT ( 1.38064852e-23 / LENGTH / LENGTH ) //ボルツマン
#define TEMPARTURE ( 300 )

#define M_A ( 1.85131596e+6 )   // g/mol/particle
#define N_A ( 6.022140857e+23 )

#define NUMBER_MAX ( 2516 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )

#define K_BOND ( 1.0e+0 )       //1つ隣　ばね定数
#define K_BOND_2 ( 1.0e-1 )     //2つ隣
#define K_EXCLUDE ( 0.7 )   // 粒子間, 粒子-SPB間の排除体積効果
#define K_HMM (1.0e-1) // 隠れマルコフ状態を用いたポテンシャルの強度

#define DELTA ( 1.0e-3 )  //刻み幅
#define MITIGATION_INTERVAL (1.0e+3)
//#define LIST_INTERVAL ( 200 )   // リスト化の間隔
#define LIST_RADIUS ( 7.0 * PARTICLE_RADIUS)

#define MEMBRANE_EXCLUDE ( 1.0e+1 )     //膜との衝突
#define NUCLEOLUS_EXCLUDE ( 1.0e+0 ) //核小体との衝突

#define BOND_DISTANCE ( 1.6 * PARTICLE_RADIUS )   // １個隣ばねの自然長

#define PARTICLE_MYU ( 2.0 * DIMENSION * PI * PARTICLE_RADIUS * LENGTH * 0.000890 / 100 ) //粘性抵抗の強さ

#define KINEMATIC_MYU (0.000890)
#define DIFFUSION (5.0e-2)

#define SPB_RADIUS (  3.0  )      //SPBの半径
#define CENT_INIT_RADIUS (1.0)  // セントロメア粒子の初期半径

#define LOCUS_MAX (137)

// Ellipsoid axes parameter of nucleus & nucleolus //

#define MEMBRANE_AXIS_1 ( 1.981780e-6 / LENGTH )     // ~2.0um
#define MEMBRANE_AXIS_2 ( 0.849 * MEMBRANE_AXIS_1 )  // ~1.6um
#define MEMBRANE_AXIS_3 ( 0.737 * MEMBRANE_AXIS_1 )  // ~1.46057186 um

#define NUCLEOLUS_AXIS_1 ( 1.147035e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.895 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.783 * NUCLEOLUS_AXIS_1 )

// 核膜 排除体積計算用 (粒子が核外に出ないように半径分短く)
#define MEM_ELLIP1_EXCLUDE ( 1.0 / (MEMBRANE_AXIS_1 - 1) / (MEMBRANE_AXIS_1 - 1))
#define MEM_ELLIP2_EXCLUDE ( 1.0 / (MEMBRANE_AXIS_2 - 1) / (MEMBRANE_AXIS_2 - 1))
#define MEM_ELLIP3_EXCLUDE ( 1.0 / (MEMBRANE_AXIS_3 - 1) / (MEMBRANE_AXIS_3 - 1))

// 核小体 排除体積計算用 (粒子が核小体に入らないように半径分長く)
#define NUC_ELLIP1_EXCLUDE ( 1.0 / (NUCLEOLUS_AXIS_1 + 1) / (NUCLEOLUS_AXIS_1 + 1))
#define NUC_ELLIP2_EXCLUDE ( 1.0 / (NUCLEOLUS_AXIS_2 + 1) / (NUCLEOLUS_AXIS_2 + 1))
#define NUC_ELLIP3_EXCLUDE ( 1.0 / (NUCLEOLUS_AXIS_3 + 1) / (NUCLEOLUS_AXIS_3 + 1))

// SPBを核膜に止める際の楕円体係数
#define MEM_ELLIP1_FIX ( 1.0 / MEMBRANE_AXIS_1 / MEMBRANE_AXIS_1)
#define MEM_ELLIP2_FIX ( 1.0 / MEMBRANE_AXIS_2 / MEMBRANE_AXIS_2)
#define MEM_ELLIP3_FIX ( 1.0 / MEMBRANE_AXIS_3 / MEMBRANE_AXIS_3)

// 第３染色体末端(rDNA)を核小体境界上に止める際の楕円体係数
#define NUC_ELLIP1_FIX ( 1.0 / NUCLEOLUS_AXIS_1 / NUCLEOLUS_AXIS_1)
#define NUC_ELLIP2_FIX ( 1.0 / NUCLEOLUS_AXIS_2 / NUCLEOLUS_AXIS_2)
#define NUC_ELLIP3_FIX ( 1.0 / NUCLEOLUS_AXIS_3 / NUCLEOLUS_AXIS_3)

#define RADIUS_MITI_STEP (5.0e+2)
#define NUCLEOLUS_MITI_STEP (1.0e+3)
#define NUCLEOLUS_MITI_DELTA (1.0 / 1.0e+3)
#define STATE_MAX (8)
#define HMM_SET_INTERVAL (2000)
#define MEAN_PHASE (500)

const unsigned int CENT_LIST[] = { 754, 1440, 2244 };
const unsigned int TELO_LIST[] = { 0, 1115, 1116, 2023};
const unsigned int rDNA_LIST[] = { 2024, 2515};
const double ORIGIN[] = { 0.0, 0.0, 0.0};

const unsigned int LIST_INTERVAL[] = {500, 500, 500};

//#define POTENTIAL_DELTA (1.0e-7)

typedef enum chain {
    A, B, C
}CHAIN;

typedef enum type {
    Normal, Centromere, Telomere, rDNA
}TYPE;

// 染色体末端粒子のリスト UP:上流 DOWN:下流 //
typedef enum chr_end {
    
    TELO1_UP = 0, TELO1_DOWN = 1115,
    TELO2_UP = 1116, TELO2_DOWN = 2023,
    rDNA_UP = 2024, rDNA_DOWN = 2515
    
}CHR_END;

typedef enum hmm_option {
    FIRST, RANDOM, CHANGE
}HMM_OPTION;

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    unsigned int chr_no;
    unsigned int particle_type;
    double position[DIMENSION];
    double force[DIMENSION];
    double radius;
    unsigned int list_no;
    unsigned int *list;
    
    // Locus 対応粒子 (137個)のみ
    double *nuc_mean;
    double *spb_mean;
    double *hmm_prob;
    double *nuc_var;
    double *spb_var;
    double **coefficient;
    int hmm_count;
    int hmm_state;
    int hmm_state_old;
    unsigned int **state_candidate;
    
} Particle;

typedef struct nucleolus{
    double position[DIMENSION]; // 重心座標
    double theta;
    double phi;
    double eta;
    double al1;
    double al2;
    double al3;
} Nuc;

enum label{ X, Y, Z};

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

// rotate function about z axis //
void rotate_about_z ( double pos[DIMENSION], const double theta ) {
    
    double pos_new[DIMENSION];
    
    pos_new[X] = cos (theta) * pos[X] - sin (theta) * pos[Y];
    pos_new[Y] = sin (theta) * pos[X] + cos (theta) * pos[Y];
    
    pos[X] = pos_new[X];
    pos[Y] = pos_new[Y];
}

// rotate function about y axis x-z平面での回転 //
void rotate_about_y ( double pos[DIMENSION], const double theta ) {
    
    double pos_new[DIMENSION];
    
    pos_new[X] = cos (theta) * pos[X] - sin (theta) * pos[Z];
    pos_new[Z] = sin (theta) * pos[X] + cos (theta) * pos[Z];
    
    pos[X] = pos_new[X];
    pos[Z] = pos_new[Z];
}

// rotate function about x axis y-z平面での回転 //
void rotate_about_x ( double pos[DIMENSION], const double theta ) {
    
    double pos_new[DIMENSION];
    
    pos_new[Y] = cos (theta) * pos[Y] - sin (theta) * pos[Z];
    pos_new[Z] = sin (theta) * pos[Y] + cos (theta) * pos[Z];
    
    pos[Y] = pos_new[Y];
    pos[Z] = pos_new[Z];
}

void secure_main_memory (Particle **part, Nuc **nuc, Particle **spb) {   // メモリ確保 //
    
    if ( (*part = (Particle *)malloc(NUMBER_MAX * sizeof(Particle))) == NULL) {
        
        printf("\n error : can not secure the memory part\n");
        exit(1);
    }
    
    for (unsigned int loop = 0; loop < NUMBER_MAX; loop++) {

        ( *part )[loop].list = (unsigned int *) malloc (NUMBER_MAX * sizeof (unsigned int));
        
        if ( ( *part )[loop].list == NULL) {
            
            printf ("\n error : can not secure the memory part.list \n");
            exit (1);
        }
    }
    
    if ( ( *nuc = (Nuc *)malloc (sizeof (Nuc))) == NULL) {
        
        printf("\n error : can not secure the memory nuc\n");
        exit(1);
    }
    
    if ( ( *spb = (Particle *)malloc (sizeof (Particle))) == NULL) {
        
        printf("\n error : can not secure the memory spb\n");
        exit(1);
    }
}

// 途中から計算する場合の 座標データ読みこみ //
void Read_coordinate (Particle *part, const unsigned int start, const char *dir) {
    
    unsigned int loop, number, i_dummy;
    char input_file[256];
    FILE *fpr;
    Particle *part_1;
    
    sprintf (input_file, "%s/result_%d.txt", dir, start);
    if ((fpr = fopen (input_file, "r")) == NULL){
        
        printf ("\t error : cannot read coordinate data.\n");
        exit (1);
    }
    
    for ( loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        
        fscanf (fpr, "%d %d %d %lf %lf %lf %lf %d\n", &i_dummy, &part_1->chr_no, &part_1->particle_type,
                &part_1->position[X], &part_1->position[Y], &part_1->position[Z], &part_1->radius, &part_1->hmm_state);
    }
    
    fclose (fpr);
}

void Read_structure (Nuc *nuc, Particle *spb, const unsigned int stable_no) {
    
    FILE *fpr;
    char filename[128], dummy[256];
    unsigned int loop, i_dummy;
    
    sprintf (filename, "stable_status.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open stable_status.txt \n");
        exit (1);
    }
    
    for (loop = 0; loop < stable_no + 1; loop++) fgets (dummy, 256, fpr);
    
    fscanf (fpr, "%d %lf %lf %lf %lf %lf %lf\n", &i_dummy, &nuc->position[X], &nuc->position[Y], &nuc->position[Z],
            &nuc->theta, &nuc->phi, &nuc->eta);
    
    fclose (fpr);
    
    sprintf (filename, "stable_spb.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open stable_spb.txt \n");
        exit (1);
    }
    
    fgets (dummy, 256, fpr);
    
    fscanf (fpr, "%d %lf %lf %lf\n", &i_dummy, &spb->position[X], &spb->position[Y], &spb->position[Z]);
    
    fclose (fpr);
    
}

void Read_hmm_state (Particle *part, int *hmm_list) {
    
    FILE *fpr;
    char filename[128], dummy[512];
    unsigned int loop, number;
    Particle *part_1;
    
    sprintf (filename, "../subdata/G2_status_5k_all_80per.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open %s \n", filename);
        exit (1);
    }
    
    fgets (dummy, 512, fpr);
    
    hmm_list [0] = 0; // hmm_stateを持つ粒子 (locus)の個数
    
    while (fscanf (fpr, "%d", &number) != EOF) {
        
        part_1 = &part[number];
        
        hmm_list[0]++; // locus数の更新
        hmm_list [hmm_list [0]] = number;
        
        if ( (part_1->spb_mean = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL
            || (part_1->nuc_mean = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL
            || (part_1->hmm_prob = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL) {
            
            printf ("\t Cannot secure memories related to hmm_state.\n");
        }
        part_1->hmm_count = 0;
        
        for (loop = 0; loop < STATE_MAX; loop++) {
        
            fscanf (fpr, "\t%lf\t%lf\t%lf", &part_1->spb_mean[loop], &part_1->nuc_mean[loop], &part_1->hmm_prob[loop]);
            
            if (part_1->spb_mean [loop] != 0.0) {
                
                part_1->hmm_count++;
            }
            part_1->spb_mean[loop] *= 1.0e-6 / LENGTH;   // シミュレーションの長さ単位に合わせる
            part_1->nuc_mean[loop] *= 1.0e-6 / LENGTH;
        }
        fgets (dummy, 256, fpr);
    }
}

void Read_hmm_state_all (Particle *part, int *hmm_list) {
    
    FILE *fpr;
    char filename[128], dummy[1024];
    unsigned int loop, number, state;
    double /*spb_var [STATE_MAX], nuc_var [STATE_MAX],*/ covar [STATE_MAX];
    Particle *part_1;
    
    sprintf (filename, "../subdata/G2_status_5k_all_80per.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open %s \n", filename);
        exit (1);
    }
    
    fgets (dummy, 1024, fpr);
    
    hmm_list [0] = 0; // hmm_stateを持つ粒子 (locus)の個数
    
    while (fscanf (fpr, "%d", &number) != EOF) {
        
        part_1 = &part[number];
        
        hmm_list[0]++; // locus数の更新
        hmm_list [hmm_list [0]] = number;
        
        // 平均,分散,共分散を格納するリストのメモリ確保
        if ( (part_1->spb_mean = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL
            || (part_1->nuc_mean = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL
            || (part_1->hmm_prob = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL
            || (part_1->spb_var = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL
            || (part_1->nuc_var = (double *)malloc (sizeof (double) * STATE_MAX)) == NULL) {
            
            printf ("\t Cannot secure memories related to hmm_state.\n");
        }
        
        part_1->hmm_count = 0;
        
        for (loop = 0; loop < STATE_MAX; loop++) {
            
            fscanf (fpr, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &part_1->spb_mean[loop], &part_1->nuc_mean[loop], &part_1->spb_var[loop],
                    &part_1->nuc_var[loop], &covar[loop], &part_1->hmm_prob[loop]);
            
            if (part_1->spb_mean [loop] != 0.0) {
                
                part_1->hmm_count++;
                
                part_1->spb_mean[loop] *= TRANS_LENGTH;   // シミュレーションの長さ単位に合わせる
                part_1->nuc_mean[loop] *= TRANS_LENGTH;
                part_1->spb_var[loop] *= TRANS_LENGTH * TRANS_LENGTH;
                part_1->nuc_var[loop] *= TRANS_LENGTH * TRANS_LENGTH;
                covar[loop] *= TRANS_LENGTH * TRANS_LENGTH;
            }
        }
//        fgets (dummy, 512, fpr);
        
        // ガウスポテンシャル計算時に使う係数のメモリ確保
        if ( (part_1->coefficient = (double **)malloc (sizeof (double *) * part_1->hmm_count)) == NULL ||
            (part_1->state_candidate = (unsigned int **)malloc (sizeof (unsigned *) * part_1->hmm_count)) == NULL ) {
            
            printf ("\t Cannot secure memories related to hmm coefficient.\n");
        }
        
        for (loop = 0; loop < part_1->hmm_count; loop++) {
            
            if ( (part_1->coefficient [loop] = (double *)malloc (sizeof (double) * 3)) == NULL ||
                (part_1->state_candidate [loop] = (unsigned int *)malloc (sizeof (unsigned) * part_1->hmm_count)) == NULL ) {
                
                printf ("\t Cannot secure memories related to hmm coefficient.\n");
            }
            
            part_1->coefficient [loop][0] = 1.0 / part_1->nuc_var [loop];
            part_1->coefficient [loop][1] = 1.0 / part_1->spb_var [loop];
            part_1->coefficient [loop][2] = part_1->coefficient[loop][0] * part_1->coefficient [loop][1] * covar [loop];
        }
        
    }
    
    sprintf (filename, "../subdata/trans_candidate.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open %s \n", filename);
        exit (1);
    }
    
    fgets (dummy, 256, fpr);
    
    while ( fscanf (fpr, "%d,%d", &number, &state) != EOF) {
        
        part_1 = &part [number];
        fscanf (fpr, ",%d", &part_1->state_candidate [state][0]);
        
        for ( loop = 1; loop <= part_1->state_candidate [state][0]; loop++ ) {
            
            fscanf (fpr, ",%d", &part_1->state_candidate [state][loop]);
        }
        fgets (dummy, 256, fpr); // 改行読み込み
    }
    
}

void Set_hmm_state (Particle *part_1, dsfmt_t *dsfmt, const int option) {
    
    unsigned int state = 0, candidate_idx;
    double prob_value;
    
    switch (option) {
        case RANDOM:
            
            prob_value = dsfmt_genrand_close_open (dsfmt);
            while (prob_value > part_1->hmm_prob [state]) state++;
            part_1->hmm_state = state;
            part_1->hmm_state_old = state;
            break;
            
        case CHANGE:
           // 存在比率による重み付きありの状態決定
//            do {
//                state = 0;
//                prob_value = dsfmt_genrand_close_open (dsfmt);
//                while (prob_value > part_1->hmm_prob [state]) state++;
//
//            } while (part_1->hmm_state_old == state);
            
            // 重みなし
//            do {
//                // state (0~hmm_countを乱数で振る)
//                state = dsfmt_genrand_uint32( dsfmt ) % part_1->hmm_count ;
//            } while (part_1->hmm_state_old == state);
//
//            part_1->hmm_state = state;
            
            state = part_1->hmm_state;
            // 遷移確率を元に抽出した候補の中から次の状態を決める
            if (part_1->state_candidate [state][0] > 0) {
            
                part_1->hmm_state_old = state; // oldにstateを保存して今とは違う状態に変更
                candidate_idx = dsfmt_genrand_uint32 (dsfmt) % part_1->state_candidate [state][0] + 1;
                part_1->hmm_state = part_1->state_candidate [state][candidate_idx];
            }
            else {  // 遷移する状態候補がなければ１つ前の状態に戻す
                
                part_1->hmm_state = part_1->hmm_state_old;
                part_1->hmm_state_old = state;
            }
            
            break;
        
        case FIRST:
            
            part_1->hmm_state = 0;
            part_1->hmm_state_old = 0;
            break;
            
        default:
            printf ("\t Cannot set hmm_state\n");
            break;
    }

}

// ひも粒子の初期座標決定
void Particle_initialization (Particle *part, Nuc *nuc, Particle *spb, dsfmt_t *dsfmt) {
    
    unsigned int loop, loop2, arm_num, dim, territory_flag[2] = {0, 0};
    double theta[2], init_theta, phi, chr_dist, vector[DIMENSION], telo_dist[2], telo_center[DIMENSION];
    Particle *part_1, *cent;
    
    // テロメアの初期位置を染色体ごとでほぼ同位置に配置する。
    for (loop = 0; loop < 2; loop++) {
        
        init_theta = 2 * PI * dsfmt_genrand_close_open (dsfmt);
        if (spb->position[Z] > 0 ) phi = 0.5 * PI * dsfmt_genrand_close_open (dsfmt) + 0.5 * PI;
        else phi = 0.5 * PI * dsfmt_genrand_close_open (dsfmt);
        
        theta[0] = init_theta + PI / 36;
        theta[1] = init_theta - PI / 36;
        
        for (loop2 = 0; loop2 < 2; loop2++) {
         
            part_1 = &part[ TELO_LIST [2 * loop + loop2] ];
            
            part_1->position[X] = MEMBRANE_AXIS_1 * sin (phi) * cos (theta [loop2]);
            part_1->position[Y] = MEMBRANE_AXIS_2 * sin (phi) * sin (theta [loop2]);
            part_1->position[Z] = MEMBRANE_AXIS_3 * cos (phi);
            
            part_1->particle_type = Telomere;
            part_1->chr_no = loop;
        }
    }
    
    //rDNA末端粒子の初期値を核小体表面上でランダムに設定 //
    init_theta = 2 * PI * dsfmt_genrand_close_open (dsfmt);
    phi = PI * dsfmt_genrand_close_open (dsfmt);
    
    theta[0] = init_theta + PI / 36;
    theta[1] = init_theta - PI / 36;
    
    for (loop = 0; loop < 2; loop++) {
        
        part_1 = &part[ rDNA_LIST[loop]];
        
        part_1->position[X] = nuc->al1 * sin (phi) * cos (theta [loop]);
        part_1->position[Y] = nuc->al2 * sin (phi) * sin (theta [loop]);
        part_1->position[Z] = nuc->al3 * cos (phi);
        
        rotate_about_x (part_1->position, nuc->eta);
        rotate_about_y (part_1->position, nuc->phi);
        rotate_about_z (part_1->position, nuc->theta);
        
        part_1->position[X] += nuc->position[X];
        part_1->position[Y] += nuc->position[Y];
        part_1->position[Z] += nuc->position[Z];
        
        part_1->chr_no = 2;
        part_1->particle_type = rDNA;
    }
    

    double init_radius, arm_dist, unit_vector[DIMENSION];
    const unsigned int chr_term[3][2] = {{TELO1_UP, TELO1_DOWN}, {TELO2_UP, TELO2_DOWN}, {rDNA_UP, rDNA_DOWN}};
    for (loop = 0; loop < 3; loop++) {
        
        part_1 = &part[ chr_term [loop][0]];
        cent = &part [ CENT_LIST[loop]];
        
        cent->chr_no = loop;
        cent->particle_type = Centromere;
        cent->radius = CENT_INIT_RADIUS;
        
        chr_dist = Euclid_norm (spb->position, part_1->position );
        
        for (dim = 0; dim < DIMENSION; dim++) {
            
            vector [dim] = part_1->position[dim] - spb->position[dim];
            // セントロメアの座標決定
            cent->position [dim] = spb->position[dim] + vector [dim] * ( CENT_INIT_RADIUS + SPB_RADIUS ) / chr_dist;
        }
        
        // 上流側　//
        // 対象腕に含まれる粒子数
        arm_num = abs (chr_term[loop][0] - CENT_LIST[loop]) - 1;
        // セントロメア-テロメアの距離
        arm_dist = Euclid_norm (cent->position, part_1->position);
        // 粒子の初期径
        init_radius = (arm_dist - BOND_DISTANCE * 0.5 * CENT_INIT_RADIUS) / (arm_num * BOND_DISTANCE + BOND_DISTANCE * 0.5) ;

        // セントロメア→テロメア方向への単位ベクトル
        for (dim = 0; dim < DIMENSION; dim++){
            
            unit_vector [dim] = (part_1->position[dim] - cent->position[dim]) / arm_dist;
        }
        // テロメア粒子の半径
        part_1->radius = init_radius;
        
        // 上流側粒子の初期座標決定
        for (loop2 = 1; loop2 <= arm_num; loop2++) {
            
            part_1 = &part [CENT_LIST [loop] - loop2];
            part_1->radius = init_radius;
            
            for (dim = 0; dim < DIMENSION; dim++) {
                
                part_1->position [dim] = cent->position [dim] + unit_vector[dim] * ( (CENT_INIT_RADIUS + init_radius ) * BOND_DISTANCE * 0.5 + (loop2 - 1) * init_radius * BOND_DISTANCE );
            }
            
            part_1->chr_no = loop;
            part_1->particle_type = Normal;
        }
        
        // 同様に下流側も決める //
        part_1 = &part[ chr_term [loop][1] ];
        arm_num = abs (chr_term[loop][1] - CENT_LIST[loop]) - 1;
        arm_dist = Euclid_norm (cent->position, part_1->position);
        init_radius = (arm_dist - BOND_DISTANCE * 0.5 * CENT_INIT_RADIUS) / (arm_num * BOND_DISTANCE + BOND_DISTANCE * 0.5);

        // セントロメア→テロメア方向への単位ベクトル
        for (dim = 0; dim < DIMENSION; dim++){
            
            unit_vector [dim] = (part_1->position[dim] - cent->position[dim]) / arm_dist;
        }
        // テロメア粒子の半径
        part_1->radius = init_radius;
        
        // 下流側粒子の初期座標決定
        for (loop2 = 1; loop2 <= arm_num; loop2++) {
            
            part_1 = &part [CENT_LIST [loop] + loop2];
            part_1->radius = init_radius;
            
            for (dim = 0; dim < DIMENSION; dim++) {
                
                part_1->position [dim] = cent->position [dim] + unit_vector [dim] * ( (CENT_INIT_RADIUS + init_radius ) * BOND_DISTANCE * 0.8 + (loop2 - 1) * init_radius * BOND_DISTANCE);
            }
            part_1->chr_no = loop;
            part_1->particle_type = Normal;
        }
        
    }
    
}

//　ばねによる力 part_1側の力計算//
void spring (Particle *part_1, const Particle *part_2, unsigned int interval) {
    
    // 線形バネの強さ　0:spb-centromere, 1,2,3: n個隣 //
    const double bonding_power[] = { K_BOND, K_BOND, K_BOND_2};
    double dist, f, dist_0;
    
    switch (interval) {
        case 0:
            
            dist_0 = part_1->radius + SPB_RADIUS;
            
            dist = Euclid_norm (part_1->position, part_2->position);
            
            f = bonding_power[interval] * (dist_0 - dist) / dist;
            
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
            break;
            
        case 1:
        case 2:
            
            dist_0 = part_1->radius * BOND_DISTANCE * interval;
            dist = Euclid_norm (part_1->position, part_2->position);
            
            f = bonding_power[interval] * (dist_0 - dist) / dist;
            
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
            break;
            
        default:
            
            break;
    }
}

void Membrane_interaction ( Particle *part_1, const char interaction_type) {
    
    double ellipsoid_dist = part_1->position[X] * part_1->position[X] * MEM_ELLIP1_EXCLUDE
    + part_1->position[Y] * part_1->position[Y] * MEM_ELLIP2_EXCLUDE
    + part_1->position[Z] * part_1->position[Z] * MEM_ELLIP3_EXCLUDE;
    
    if ( interaction_type == 'F' || ellipsoid_dist - 1 > 0 ) {
        
        // 法線ベクトル
        double normal_vector[] = { 2.0 * part_1->position[X] * MEM_ELLIP1_EXCLUDE,
            2.0 * part_1->position[Y] * MEM_ELLIP2_EXCLUDE,
            2.0 * part_1->position[Z] * MEM_ELLIP3_EXCLUDE };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * ( part_1->position[X] * normal_vector[X]
                                                                  + part_1->position[Y] * normal_vector[Y]
                                                                  + part_1->position[Z] * normal_vector[Z]) / normal_vector_norm;
        
        part_1->force[X] += f * normal_vector[X];
        part_1->force[Y] += f * normal_vector[Y];
        part_1->force[Z] += f * normal_vector[Z];
    }
}

// 核小体との排除体積効果　//
void Nucleolus_interaction ( Particle *part_1, Nuc *nuc, const char interaction_type) {
    
    //核小体中心から粒子へのベクトル
    double nuc_to_pos[DIMENSION] = { part_1->position[X] - nuc->position[X],
        part_1->position[Y] - nuc->position[Y],
        part_1->position[Z] - nuc->position[Z]};
    
    // 核小体座標系に変換 //
    //位置座標をx-y平面で-theta回転
    rotate_about_z (nuc_to_pos, -nuc->theta );
    // x-z平面
    rotate_about_y (nuc_to_pos, -nuc->phi);
    // y-z平面
    rotate_about_x (nuc_to_pos, -nuc->eta);
    
    double ellipsoid_dist =  nuc_to_pos[X] * nuc_to_pos[X] / ( (nuc->al1 + 1) * (nuc->al1 + 1) )
    + nuc_to_pos[Y] * nuc_to_pos[Y] / ( (nuc->al2 + 1) * (nuc->al2 + 1) )
    + nuc_to_pos[Z] * nuc_to_pos[Z] / ( (nuc->al3 + 1) * (nuc->al3 + 1) );
    
    if ( interaction_type == 'F' || ellipsoid_dist < 1.0 ) {
        
        // 法線ベクトル @核小体座標系
        double normal_vector[] = { 2.0 * nuc_to_pos[X] / ( (nuc->al1 + 1) * (nuc->al1 + 1) ),
            2.0 * nuc_to_pos[Y] / ( (nuc->al2 + 1) * (nuc->al2 + 1) ),
            2.0 * nuc_to_pos[Z] / ( (nuc->al3 + 1) * (nuc->al3 + 1) ) };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * NUCLEOLUS_EXCLUDE * Inner_product (nuc_to_pos, normal_vector)
        / normal_vector_norm;
        
        rotate_about_x (normal_vector, nuc->eta);
        rotate_about_y (normal_vector, nuc->phi);
        rotate_about_z (normal_vector, nuc->theta);
        
        part_1->force[X] += f * normal_vector[X];
        part_1->force[Y] += f * normal_vector[Y];
        part_1->force[Z] += f * normal_vector[Z];
    }
}

void spb_exclusion (Particle *part_1, Particle *spb) {
    
    // spb_exclusion //
    double dist = Euclid_norm (part_1->position, spb->position);
    double f = K_EXCLUDE * ( SPB_RADIUS + part_1->radius - dist) / dist;
    
    if ( dist < part_1->radius + SPB_RADIUS) {
        
        part_1->force[X] += f * (part_1->position[X] -  spb->position[X]);
        part_1->force[Y] += f * (part_1->position[Y] -  spb->position[Y]);
        part_1->force[Z] += f * (part_1->position[Z] -  spb->position[Z]);
    }
}

//　ひも同士のリスト化 //
void make_ve_list (Particle *part, Particle *part_1, const unsigned int target) {
    
    unsigned int loop, list_count = 0;
    double dist;
    Particle *part_2;
    
    part_1->list_no = 0;
    
    for ( loop = 0; loop < NUMBER_MAX; loop++ ){
        
        part_2 = &part[loop];
        dist = Euclid_norm ( part_1->position, part_2->position);
        
        // 1個隣 そうでなくてもテロメア同士 (番号は隣だが染色体No.が異なる) //
        if ( dist < LIST_RADIUS && (abs (target - loop) > 1 || part_1->chr_no != part_2->chr_no ) ) {
            
            list_count++;
            part_1->list_no = list_count;
            part_1->list [list_count] = loop;
        }
    }
}

// Volume exclusion between particles //
void particle_exclusion (Particle *part, Particle *part_1) {
    
    unsigned int loop;
    double dist, f;
    Particle *part_2;
    
    for ( loop = 1; loop <= part_1->list_no; loop++) {
        
        part_2 = &part [part_1->list [loop]];
        dist = Euclid_norm (part_1->position, part_2->position);
        
        if ( dist < part_1->radius + part_2->radius ){
            
            f = K_EXCLUDE * (part_1->radius + part_2->radius - dist) / dist;
            
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
        }
    }
}

void Hmm_spring_potential (Particle *part_1, Nuc *nuc, Particle *spb) {
    
    unsigned int loop;
    double spb_f, nuc_f, dist;
    
    dist = Euclid_norm (part_1->position, spb->position);
    spb_f = K_HMM * (part_1->spb_mean[ part_1->hmm_state ] - dist) / dist;
    
    dist = Euclid_norm (part_1->position, nuc->position);
    nuc_f = K_HMM * (part_1->nuc_mean[ part_1->hmm_state ] - dist) / dist;
    
    for (loop = 0; loop < DIMENSION; loop++) {
     
        part_1->force[loop] += spb_f * (part_1->position[loop] - spb->position[loop]) + nuc_f * (part_1->position[loop] - nuc->position[loop]);
    }
}

void Hmm_gauss_potential (Particle *part_1, Nuc *nuc, Particle *spb) {
    
    unsigned int loop, state;
    double spb_dist, nuc_dist, dist, *coef;
    
    nuc_dist = Euclid_norm (part_1->position, nuc->position);
    spb_dist = Euclid_norm (part_1->position, spb->position);
    
    state = part_1->hmm_state;
    coef = part_1->coefficient [state];
    
    for (loop = 0; loop < DIMENSION; loop++) {
        
        part_1->force[loop] += - K_HMM * (coef [0] * (nuc_dist - part_1->nuc_mean[state]) - coef[2] * (spb_dist - part_1->spb_mean[state]))
                            * ( part_1->position[loop] - nuc->position[loop]) / nuc_dist
                            - K_HMM * (coef[1] * (spb_dist - part_1->spb_mean[state]) - coef[2] * (nuc_dist - part_1->nuc_mean[state]))
                               * ( part_1->position[loop] - spb->position[loop]) / spb_dist;

    }
}

void Noise (double *force, dsfmt_t *dsfmt) {
    
    //noise dsfmt
//    double p1 = sqrt(2.0 * 3.0 * KINEMATIC_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ))
//    / sqrt (DELTA);
//    double p2 = sqrt(2.0 * 3.0 * KINEMATIC_MYU * KBT * TEMPARTURE) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ))
//    / sqrt (DELTA);
    
    double p1 = sqrt(2.0 * DIFFUSION) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ))
    / sqrt (DELTA);
    double p2 = sqrt(2.0 * DIFFUSION) * sqrt(-2.0 * log( dsfmt_genrand_open_close(dsfmt) ))
    / sqrt (DELTA);
    double theta = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
    double psi = 2.0 * PI * dsfmt_genrand_open_close(dsfmt);
    
    force[X] += p1 * sin(theta);
    force[Y] += p1 * cos(theta);
    force[Z] += p2 * sin(psi);
}

// 各stepごとの座標計算 //
void calculation (Particle *part, Nuc *nuc, Particle *spb, const unsigned int mitigation, dsfmt_t *dsfmt, const int *hmm_list, unsigned int calc_phase) {
    
    unsigned int loop;
    Particle *part_1, *part_2, *part_3;
    
    // 位置の計算 & 力の初期化 //
    for ( loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
    }
    
    if (calc_phase > 0) {
        
        for (loop = 0; loop < NUMBER_MAX; loop++) Noise (part[loop].force, dsfmt);
        for (loop = 1; loop <= hmm_list[0]; loop++) Hmm_spring_potential (&part [hmm_list [loop]], nuc, spb);
    }

    #pragma omp parallel for private (part_1) num_threads (8)
    for ( loop = 0; loop < NUMBER_MAX; loop++ ){

        part_1 = &part [loop];

        switch (part_1->particle_type) {
            case Normal:
                
                // spring 1 //
                spring (part_1, &part [loop + 1], 1);
                spring (part_1, &part [loop - 1], 1);

                // spring 2 & 3 //
                switch (loop)  {

                    case TELO1_UP + 1:
                    case TELO2_UP + 1:
                    case rDNA_UP + 1:

                        spring (part_1, &part[loop + 2], 2);
                        
                        break;

                    case TELO1_UP + 2:
                    case TELO2_UP + 2:
                    case rDNA_UP + 2:

                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);

                        
                        break;

                    case TELO1_DOWN - 1:
                    case TELO2_DOWN - 1:
                    case rDNA_DOWN - 1:

                        spring (part_1, &part[loop - 2], 2);

                        
                        break;

                    case TELO1_DOWN - 2:
                    case TELO2_DOWN - 2:
                    case rDNA_DOWN - 2:

                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);

                        
                        break;

                    default:

                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);

                        
                        
                        break;
                }

                spb_exclusion (part_1, spb);
                Nucleolus_interaction (part_1, nuc, 'E');
                Membrane_interaction (part_1, 'E');

                if ( mitigation % LIST_INTERVAL[calc_phase] == 0 ) make_ve_list (part, part_1, loop);
                particle_exclusion (part, part_1);

                break;

            case Centromere:    // calc_phase >= 1のときのみ

                if (calc_phase > 0) {
                    
                    spring (part_1, &part[loop + 1], 1);
                    spring (part_1, &part[loop - 1], 1);

                    spring (part_1, &part[loop + 2], 2);
                    spring (part_1, &part[loop - 2], 2);

                    
                    

                    Nucleolus_interaction (part_1, nuc, 'E');
                    Membrane_interaction (part_1, 'E');
                    spring (part_1, spb, 0);

                    if ( mitigation % LIST_INTERVAL[calc_phase] == 0 ) make_ve_list (part, part_1, loop);
                    particle_exclusion (part, part_1);
                    
                } else {
                    
                    // 粒子に一括で与えたノイズをリセット
                    part_1->force[X] = 0.0;
                    part_1->force[Y] = 0.0;
                    part_1->force[Z] = 0.0;
                }
                

                break;

            case Telomere:
                
                switch (loop) {
                    case TELO1_UP:
                    case TELO2_UP:

                        spring (part_1, &part[loop + 1], 1);

                        spring (part_1, &part[loop + 2], 2);

                        break;

                    case TELO1_DOWN:
                    case TELO2_DOWN:

                        spring (part_1, &part[loop - 1], 1);

                        spring (part_1, &part[loop - 2], 2);
                        
                        break;
                }
                
                spb_exclusion (part_1, spb);
                Membrane_interaction (part_1, 'F');
                Nucleolus_interaction (part_1, nuc, 'E');

                if ( mitigation % LIST_INTERVAL[calc_phase] == 0 ) make_ve_list (part, part_1, loop);
                particle_exclusion (part, part_1);
                
                break;

            case rDNA:

                switch (loop) {
                    case rDNA_UP:

                        spring (part_1, &part[loop + 1], 1);

                        spring (part_1, &part[loop + 2], 2);

                        

                        break;

                    case rDNA_DOWN:

                        spring (part_1, &part[loop - 1], 1);

                        spring (part_1, &part[loop - 2], 2);

                        

                        break;
                }

                spb_exclusion (part_1, spb);
                Membrane_interaction (part_1, 'E');
                Nucleolus_interaction (part_1, nuc, 'F');

                if ( mitigation % LIST_INTERVAL[calc_phase] == 0 ) make_ve_list (part, part_1, loop);
                particle_exclusion (part, part_1);

                break;

            default:
                printf ("\t Labeling error occured. \n");
                exit(1);
        }

    }

    for ( loop = 0; loop < NUMBER_MAX; loop++ ) {
        
        part_1 = &part[loop];
        
        part_1->position[X] += DELTA * part_1->force[X];
        part_1->position[Y] += DELTA * part_1->force[Y];
        part_1->position[Z] += DELTA * part_1->force[Z];
    }
}

void update_radius (Particle *part, const char operation) {
    
    unsigned int loop;
    static double diff_radius[6];
    
    if (operation == 'c' && part[0].radius < 1.0) {
        
        for (loop = TELO_LIST[0]; loop < CENT_LIST[0]; loop++) part [loop].radius += diff_radius [0];
        for (loop = CENT_LIST[0] + 1; loop <= TELO_LIST[1]; loop++) part [loop].radius += diff_radius [1];
    
        for (loop = TELO_LIST[2]; loop < CENT_LIST[1]; loop++) part [loop].radius += diff_radius [2];
        for (loop = CENT_LIST[1] + 1; loop <= TELO_LIST[3]; loop++) part [loop].radius += diff_radius [3];
    
        for (loop = rDNA_LIST[0]; loop < CENT_LIST[2]; loop++) part [loop].radius += diff_radius [4];
        for (loop = CENT_LIST[2] + 1; loop <= rDNA_LIST[1]; loop++) part [loop].radius += diff_radius [5];
        
        if (part[0].radius + diff_radius[0] > 1.0 ) {
            
            for (loop = 0; loop < NUMBER_MAX; loop++) {
                
                part [loop].radius = 1.0;
            }
        }
    }
    else if (operation == 's') {
        
        for (loop = 0; loop < 4; loop++) {
            
            diff_radius [loop] = (PARTICLE_RADIUS - part[ TELO_LIST [loop]].radius) / RADIUS_MITI_STEP;
        }
        for (loop = 0; loop < 2; loop++) {
            
            diff_radius [loop + 4] = (PARTICLE_RADIUS - part[ rDNA_LIST [loop]].radius) / RADIUS_MITI_STEP;
        }
    }
}



void write_coordinate (Particle *part, const unsigned int sample_no, const char *dir) {
    
    unsigned int loop;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "%s/result_%d.txt", dir, sample_no);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \t error : cannot write coordinate. \n");
        
        exit (1);
    }
    
    for (loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        fprintf (fpw, "%d %d %d %lf %lf %lf %lf %d\n", loop, part_1->chr_no, part_1->particle_type, part_1->position[X],
                 part_1->position[Y], part_1->position[Z], part_1->radius, part_1->hmm_state);
    }
    
//    fprintf (fpw, "Radius %1.1e\n", LENGTH);
    fclose (fpw);
}

void Save_settings (const char *dir, const int start, const int phase) {
    
    FILE *fpw;
    char filename[256];
    
    sprintf (filename, "%s/readme_%s.txt", dir, dir);
    
    if ((fpw = fopen (filename, "a")) == NULL) {
        
        printf (" \t error : cannot write coordinate. \n");
        
        exit (1);
    }
    
    fprintf (fpw, "start = %d\t phase %d\n\n", start, phase);
    
    fprintf (fpw, "K_BOND = %2.1e\nK_BOND_2 = %2.1e\nK_BOND_HMM = %2.1e\nLIST_INTERVAL = {%d, %d, %d}\nK_EXCLUDE = %2.1e\n",
             K_BOND, K_BOND_2, K_HMM, LIST_INTERVAL[0], LIST_INTERVAL[1], LIST_INTERVAL[2], K_EXCLUDE);
    fprintf (fpw, "DELTA = %2.1e\nMITIGATION_INTERVAL = %2.1e\n", DELTA, MITIGATION_INTERVAL);
    fprintf (fpw, "DIFFUSION = %2.1e\n", DIFFUSION);
    
    fclose (fpw);
}

double Max (double a, double b) {
    return a > b ? a:b;
}

void Calculate_strain (Particle *part, int *hmm_list, double *total_strain_mean, double *strain_max_list) {
    
    unsigned int loop;
    double strain [2], strain_max, strain_mean = 0.0;
    Particle *part_1;
    
    for (loop = 1; loop <= hmm_list[0]; loop++) {
        
        part_1 = &part [ hmm_list [loop]];
        
        // locus対応粒子の隣接粒子とのばねのずれ　0:上流側 1:下流側
        strain [0] = fabs (Euclid_norm (part_1->position, part [ hmm_list [loop] - 1].position) - part_1->radius * BOND_DISTANCE);
        strain [1] = fabs (Euclid_norm (part_1->position, part [ hmm_list [loop] + 1].position) - part_1->radius * BOND_DISTANCE);
        
        // 自然長とのずれの総和を求める
        strain_mean += ( strain [0] + strain [1] ) * 0.5;
        strain_max = Max (strain [0], strain [1]);  // 上流側下流側の歪みの最大値
        strain_max_list [loop] = Max (strain_max_list [loop], strain_max);  // 最大値の更新
    }
    
    *total_strain_mean += strain_mean / hmm_list[0]; // 歪みの粒子平均
}

int Evaluate_gauss (Particle *part, Nuc *nuc, Particle *spb, int *hmm_list, unsigned int *eval_list) {
    
    unsigned int loop, state, eval_count = 0;
    double nuc_diff, spb_diff, ellipsoid_dist;
    Particle *part_1;
    
    
    for (loop = 1; loop <=hmm_list [0]; loop++) {
        
        part_1 = &part [ hmm_list [loop] ];
        state = part_1->hmm_state;
        
        nuc_diff = Euclid_norm (part_1->position, nuc->position) - part_1->nuc_mean [state];
        spb_diff = Euclid_norm (part_1->position, spb->position) - part_1->spb_mean [state];
        
        if ( nuc_diff * nuc_diff / part_1->nuc_var [state] + spb_diff * spb_diff / part_1->spb_var [state] < 1.0 ) {
            
            eval_list [loop] = 1;
            
        }
        else {
            
            eval_count++;
            eval_list [loop] = 0;
        }
    }
    
    return eval_count;
}

void Hmm_set_mitigation (Particle *part, Nuc *nuc, Particle *spb, unsigned int *hmm_list, dsfmt_t *dsfmt,
                         unsigned int calc_phase, unsigned int *total_time, const char *directory) {
    
    unsigned int loop, time, mitigation;
    unsigned int try_count = 0, eval_list [138], change_count, eval_count;
    double total_strain_mean, strain_max_list [138];
    strain_max_list [0] = 137;
    eval_list [0] = 137;
    
    do {
        
        try_count++;
        change_count = 0;
        total_strain_mean = 0.0;
        for (loop = 1; loop <= 137; loop++) strain_max_list [loop] = 0.0;
        
        if (calc_phase != 2) {
            
            // 緩和
            for ( unsigned int time = 1; time <= HMM_SET_INTERVAL - MEAN_PHASE; time++) {
                
                ++ *total_time;
                printf ("\t Now calculating... phase %d / try_count = %d, time = %d \r", calc_phase, try_count, time);
                fflush (stdout);
                
                for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){
                    
                    calculation (part, nuc, spb, mitigation, dsfmt, hmm_list, calc_phase);
                }
                write_coordinate (part, *total_time, directory);
                
            }
            
            // 緩和 + 歪みの平均・最大値計算
            for ( unsigned int time = HMM_SET_INTERVAL - MEAN_PHASE + 1; time <= HMM_SET_INTERVAL; time++) {
                
                ++ *total_time;
                printf ("\t Now calculating... phase %d / try_count = %d, time = %d \r", calc_phase, try_count, time);
                fflush (stdout);
                
                for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){
                    
                    calculation (part, nuc, spb, mitigation, dsfmt, hmm_list, calc_phase);
                }
                write_coordinate (part, *total_time, directory);
                
                Calculate_strain (part, hmm_list, &total_strain_mean, strain_max_list);
            }
        }else {     // 核小体領域　増大時
            
            // 緩和
            for ( unsigned int time = 1; time <= HMM_SET_INTERVAL - MEAN_PHASE; time++) {
                
                ++ *total_time;
                printf ("\t Now calculating... phase %d / try_count = %d, time = %d \r", calc_phase, try_count, time);
                fflush (stdout);
                
                for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){
                    
                    calculation (part, nuc, spb, mitigation, dsfmt, hmm_list, calc_phase);
                }
                write_coordinate (part, *total_time, directory);
                
                if ( time % 10 == 0 || nuc->al1 < NUCLEOLUS_AXIS_1) {        // 核小体領域増大
                    
                    nuc->al1 += NUCLEOLUS_AXIS_1 * NUCLEOLUS_MITI_DELTA * 10;
                    nuc->al2 += NUCLEOLUS_AXIS_2 * NUCLEOLUS_MITI_DELTA * 10;
                    nuc->al3 += NUCLEOLUS_AXIS_3 * NUCLEOLUS_MITI_DELTA * 10;
                }
            }
            
            // 緩和 + 歪みの平均・最大値計算
            for ( unsigned int time = HMM_SET_INTERVAL - MEAN_PHASE + 1; time <= HMM_SET_INTERVAL; time++) {
                
                ++ *total_time;
                printf ("\t Now calculating... phase %d / try_count = %d, time = %d \r", calc_phase, try_count, time);
                fflush (stdout);
                
                for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){
                    
                    calculation (part, nuc, spb, mitigation, dsfmt, hmm_list, calc_phase);
                }
                write_coordinate (part, *total_time, directory);
                
                if ( time % 10 == 0 || nuc->al1 < NUCLEOLUS_AXIS_1) {      // 核小体領域増大
                    
                    nuc->al1 += NUCLEOLUS_AXIS_1 * NUCLEOLUS_MITI_DELTA * 10;
                    nuc->al2 += NUCLEOLUS_AXIS_2 * NUCLEOLUS_MITI_DELTA * 10;
                    nuc->al3 += NUCLEOLUS_AXIS_3 * NUCLEOLUS_MITI_DELTA * 10;
                }
                
                Calculate_strain (part, hmm_list, &total_strain_mean, strain_max_list);
            }
        }
    
        total_strain_mean /= MEAN_PHASE; // 時間方向の平均
        
        eval_count = Evaluate_gauss (part, nuc, spb, hmm_list, eval_list);
        
        for (loop = 1; loop <= hmm_list [0]; loop++) {
            
            if ( eval_list [loop] == 0 || strain_max_list [loop] > 0.5 ) {
                
                change_count++;
                Set_hmm_state (&part [hmm_list[loop]], dsfmt, CHANGE);
            }
        }
        
        printf ("\t try_count = %d, strain_mean = %lf change_count = %d  eval_count = %d \n", try_count, total_strain_mean, change_count, eval_count);
        
        if (calc_phase == 2 && nuc->al1 >= NUCLEOLUS_AXIS_1 ) {
            
            nuc->al1 = NUCLEOLUS_AXIS_1;
            nuc->al2 = NUCLEOLUS_AXIS_2;
            nuc->al3 = NUCLEOLUS_AXIS_3;
            
            break;
        }
        
    } while ( total_strain_mean > 0.11 || change_count > 30);
    
}

int main ( int argc, char **argv ) {
    
    unsigned int loop, mitigation, start, total_time, stable_no, sample_no, calc_phase, hmm_set_option, init_time;
    char output_file[256], directory[256], hmm_select;
    double mem_al[3];
    
    Particle *part, *part_1, *spb, *part_cp;
    Nuc *nuc;
    
    int hmm_list [LOCUS_MAX + 1];
    for (loop = 0; loop < LOCUS_MAX + 1; loop++) hmm_list [loop] = -1;
    
    secure_main_memory (&part, &nuc, &spb);
    
    dsfmt_t dsfmt;
    
    if (argc == 5 ) {
        
        stable_no = atoi (argv[1]);
        sample_no = atoi (argv[2]);
        total_time = atoi (argv[3]);
        calc_phase = atoi (argv[4]);
    }
    else {
        
        printf ("\t error : Number of arguments error \n");
        exit (1);
    }
    
    init_time = total_time;
    dsfmt_init_gen_rand (&dsfmt, sample_no);
    sprintf (directory, "%d_%d", stable_no, sample_no);
    
    Read_structure (nuc, spb, stable_no);

    nuc->al1 = 0.1 * NUCLEOLUS_AXIS_1;
    nuc->al2 = 0.1 * NUCLEOLUS_AXIS_2;
    nuc->al3 = 0.1 * NUCLEOLUS_AXIS_3;
    
    
    hmm_set_option = FIRST;
    Read_hmm_state_all (part, hmm_list);   // 隠れマルコフ状態のデータを読み込み
    
    if (total_time == 0) {
        
        Particle_initialization (part, nuc, spb, &dsfmt);
        write_coordinate (part, 0, directory);
        update_radius (part, 's');
        
        // 隠れマルコフ状態の設定
        for (loop = 0; loop < NUMBER_MAX; loop++) part [loop].hmm_state = -1;
        for (loop = 1; loop <= hmm_list[0]; loop++) {
            Set_hmm_state (&part[ hmm_list[loop]], &dsfmt, hmm_set_option);
        }
    }
    else if (total_time > RADIUS_MITI_STEP) {
        
        Read_coordinate (part, total_time, directory);
    }
    else {
        
        Read_coordinate (part, total_time, directory);
        // 隠れマルコフ状態の設定
        for (loop = 0; loop < NUMBER_MAX; loop++) part [loop].hmm_state = -1;
        for (loop = 1; loop <= hmm_list[0]; loop++) {
            Set_hmm_state (&part[ hmm_list[loop]], &dsfmt, hmm_set_option);
        }
    }
    
    Save_settings (directory, total_time, calc_phase);
    
    printf ("\t K_BOND = %2.1e, K_BOND2 = %2.1e, K_HMM = %2.1e, DIFFUSE = %2.1e\n ", K_BOND, K_BOND_2, K_HMM, DIFFUSION);
    
    if (calc_phase == 0) {  // 粒子径 増加
        
        for ( unsigned int time = 1; time <= RADIUS_MITI_STEP; time++) {
            
            total_time++;
            printf ("\t Now calculating... phase 0 / time = %d \r", time);
            fflush (stdout);
            
            for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){
                
                calculation (part, nuc, spb, mitigation, &dsfmt, hmm_list, calc_phase);
            }
            
            write_coordinate (part, total_time, directory);
            
            update_radius (part, 'c');
        }
        calc_phase++;
    }
    
    
    
    // 隠れマルコフ状態セットの最適化
    // 隣接間のバネのずれの最大値 < 0.1 && ずれ平均が隠れマルコフポテンシャル無のときとの同じくらい
    // 1回の緩和時間5000step 最後の1000stepでずれ平均・最大値を計算　→評価

    if (calc_phase == 1) {
        
        Hmm_set_mitigation (part, nuc, spb, hmm_list, &dsfmt, 1, &total_time, directory);
        calc_phase++;
    }
    if (calc_phase == 2) {  // 核小体 増大
        
        Hmm_set_mitigation (part, nuc, spb, hmm_list, &dsfmt, 2, &total_time, directory);
        calc_phase++;
    }
    if (calc_phase == 3) {
        
        Hmm_set_mitigation (part, nuc, spb, hmm_list, &dsfmt, 3, &total_time, directory);
    }
    
    
//    write_coordinate (part, sample_no, directory);
    
    // メモリ解放 //
    for ( loop = 0 ; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part [loop];
        free (part_1->list);
    }
    free (part);
    free (nuc);
    free (spb);
    
    return ( 0 );
}



