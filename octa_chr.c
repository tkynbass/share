//
//
//  octa_sn_fix.c で得られた核小体・核膜・SPBの位置関係に
//  ひもを入れて緩和させる
//  各染色体の端点(テロメア)をランダムに設定 →粒子径
//  5kbpで粗視化
//


// 2019/09/05
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "/home/gensho/tkym/dSFMT/dSFMT.h"
#include <omp.h>

#define DIMENSION ( 3 ) //次元
#define LENGTH ( 7.0e-8 )   //長さの単位 (粒子径)

#define M_A ( 1.85131596e+6 )   // g/mol/particle
#define N_A ( 6.022140857e+23 )
#define PASTIS_SCALING ( 1.8e-6 / 75 / LENGTH)

#define NUMBER_MAX ( 2516 )    //粒子数
#define PARTICLE_MASS ( 1.0)    //染色体粒子の質量 Kg
#define PARTICLE_RADIUS ( 1.0 )     //粒子の半径
#define PI ( M_PI )

#define K_BOND ( 1.0e+0 )       //1つ隣　ばね定数
#define K_BOND_2 ( 1.0e+0 )     //2つ隣
#define K_BOND_3 ( 1.0e+0 )     //3つ隣
#define K_EXCLUDE ( 1.0e+1 )

#define DELTA ( 1.0e-3 )  //刻み幅
#define MITIGATION_INTERVAL (1.0e+3)
#define LIST_INTERVAL ( 50 )   // リスト化の間隔
#define LIST_RADIUS ( 10.0 * PARTICLE_RADIUS)

#define MEMBRANE_EXCLUDE ( 1.0 )     //膜との衝突
#define MEMBRANE_EXCLUDE_SPB ( 1.0 ) //SPBとの衝突

#define BOND_DISTANCE ( 2.0 * PARTICLE_RADIUS * 0.8 )   // １個隣ばねの自然長

#define PARTICLE_MYU ( 2.0 * DIMENSION * PI * PARTICLE_RADIUS * LENGTH * 0.000890 / 100 ) //粘性抵抗の強さ

#define SPB_RADIUS (  3.0  )      //SPBの半径
#define CENT_INIT_RADIUS (0.5)  // セントロメア粒子の初期半径

// Ellipsoid axes parameter of nucleus & nucleolus //

#define MEMBRANE_AXIS_1 ( 1.889011e-6 / LENGTH )
#define MEMBRANE_AXIS_2 ( 0.85 * MEMBRANE_AXIS_1 )
#define MEMBRANE_AXIS_3 ( 0.75 * MEMBRANE_AXIS_1 )

#define NUCLEOLUS_AXIS_1 ( 1.1426593e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.9 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.8 * NUCLEOLUS_AXIS_1 )

#define RADIUS_MITI_STEP (1.0e+3)

const unsigned int CENT_LIST[] = { 754, 1440, 2244 };
const unsigned int TELO_LIST[] = { 0, 1115, 1116, 2023};
const unsigned int rDNA_LIST[] = { 2024, 2515};
const double ORIGIN[] = { 0.0, 0.0, 0.0};

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

double Euclid_norm (const double pos_1[DIMENSION], const double pos_2[DIMENSION]);

typedef struct particle {           //構造体の型宣言
    unsigned int chr_no;
    unsigned int particle_type;
    double position[DIMENSION];
    double velocity[DIMENSION];
    double force[DIMENSION];
    double radius;
    unsigned int list_no;
    unsigned int *list;
    
} Particle;

typedef struct nucleolus{
    
    double position[DIMENSION]; // 重心座標
    double theta;
    double phi;
    double eta;
} Nuc;

enum label{ X, Y, Z};

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
void read_coordinate (Particle *part, const unsigned int start) {
    
    unsigned int loop, number, i_dummy;
    char pastis_data[256], input_file[256], *arm_list[] = {"1long", "1short", "2short", "2long", "3short", "3long"};
    FILE *fpr;
    Particle *part_1;
    
    sprintf (input_file, "result_%d.txt", start);
    if ((fpr = fopen (input_file, "r")) == NULL){
        
        printf ("\t erroe : cannot read coordinate data.\n");
        exit (1);
    }
    
    for ( loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        
        fscanf (fpr, "%d %d %d %lf %lf %lf\n", &i_dummy, &part_1->chr_no, &part_1->particle_type,
                &part_1->position[X], &part_1->position[Y], &part_1->position[Z]);
    }
    
    fclose (fpr);
}

void Read_structure (Nuc *nuc, Particle *spb, const unsigned int stable_no) {
    
    FILE *fpr;
    char filename[128], dummy[256];
    unsigned int loop, i_dummy;
    
    sprintf (filename, "../subdata/stable_status.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open stable_status.txt \n");
        exit (1);
    }
    
    for (loop = 0; loop < stable_no + 1; loop++) fgets (dummy, 256, fpr);
    
    fscanf (fpr, "%d %lf %lf %lf %lf %lf %lf\n", &i_dummy, &nuc->position[X], &nuc->position[Y], &nuc->position[Z],
            &nuc->theta, &nuc->phi, &nuc->eta);
    
    fclose (fpr);
    
    sprintf (filename, "../subdata/stable_spb.txt");
    
    if ( (fpr = fopen (filename, "r")) == NULL) {
        
        printf ("\t Cannot open stable_spb.txt \n");
        exit (1);
    }
    
    fgets (dummy, 256, fpr);
    
    fscanf (fpr, "%d %lf %lf %lf\n", &i_dummy, &spb->position[X], &spb->position[Y], &spb->position[Z]);
    
    fclose (fpr);
}

// ひも粒子の初期座標決定
void Particle_initialization (Particle *part, Nuc *nuc, Particle *spb, dsfmt_t *dsfmt) {
    
    unsigned int loop, loop2, arm_num, dim;
    double theta, phi, chr_dist, vector[DIMENSION];
    Particle *part_1, *cent;
    
    //テロメア rDNA末端粒子の初期値を核膜上または核小体表面上でランダムに設定 //
    for (loop = 0; loop < 4; loop++) {
        
        part_1 = &part[ TELO_LIST[loop]];
        
        theta = 2 * PI * dsfmt_genrand_close_open (dsfmt);
        phi = 0.5 * PI * dsfmt_genrand_close_open (dsfmt);
        
        part_1->position[X] = MEMBRANE_AXIS_1 * sin (phi) * cos (theta);
        part_1->position[Y] = MEMBRANE_AXIS_2 * sin (phi) * sin (theta);
        part_1->position[Z] = MEMBRANE_AXIS_3 * cos (phi);
        
        part_1->particle_type = Telomere;
        part_1->chr_no = loop / 2;
    }
    for (loop = 0; loop < 2; loop++) {
        
        part_1 = &part[ rDNA_LIST[loop]];
        
        theta = 2 * PI * dsfmt_genrand_close_open (dsfmt);
        phi = 0.5 * PI * dsfmt_genrand_close_open (dsfmt);
        
        part_1->position[X] = NUCLEOLUS_AXIS_1 * sin (phi) * cos (theta) + nuc->position[X];
        part_1->position[Y] = NUCLEOLUS_AXIS_2 * sin (phi) * sin (theta) + nuc->position[Y];
        part_1->position[Z] = NUCLEOLUS_AXIS_3 * cos (phi) + nuc->position [Z];
        
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
        init_radius = (arm_dist - 0.9 * CENT_INIT_RADIUS) / (arm_num * 1.8 + 0.9) ;

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
                
                part_1->position [dim] = cent->position [dim] + unit_vector[dim] * ( (CENT_INIT_RADIUS + init_radius ) * 0.9 + (loop2 - 1) * init_radius * 1.8 );
            }
            
            part_1->chr_no = loop;
            part_1->particle_type = Normal;
        }
        
        // 同様に下流側も決める //
        part_1 = &part[ chr_term [loop][1] ];
        arm_num = abs (chr_term[loop][1] - CENT_LIST[loop]) - 1;
        arm_dist = Euclid_norm (cent->position, part_1->position);
        init_radius = (arm_dist - 0.9 * CENT_INIT_RADIUS) / (arm_num + 0.9);

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
                
                part_1->position [dim] = cent->position [dim] + unit_vector [dim] * ( (CENT_INIT_RADIUS + init_radius ) * 0.9 + (loop2 - 1) * init_radius * 1.8);
            }
            part_1->chr_no = loop;
            part_1->particle_type = Normal;
        }
        
    }
    
}

// particle_type labeling セントロメアorテロメアor rDNA末端 //
//void type_labeling (Particle *part) {
//
//    unsigned int loop;
//
//    for ( loop = 0; loop < NUMBER_MAX; loop++ ) part [loop].particle_type = Normal;
//
//    for ( loop = 0; loop < sizeof (CENT_LIST) / sizeof (CENT_LIST[0]); loop++ ) part [CENT_LIST [loop]].particle_type = Centromere;
//    for ( loop = 0; loop < sizeof (TELO_LIST) / sizeof (TELO_LIST[0]); loop++ ) part [TELO_LIST [loop]].particle_type = Telomere;
//    for ( loop = 0; loop < sizeof (rDNA_LIST) / sizeof (rDNA_LIST[0]); loop++ ) part [rDNA_LIST [loop]].particle_type = rDNA;
//}

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

//　ばねによる力 part_1側の力計算//
void spring (Particle *part_1, const Particle *part_2, unsigned int interval) {
    
    // 線形バネの強さ　0:spb-centromere, 1,2,3: n個隣 //
    const double bonding_power[] = { K_BOND, K_BOND, K_BOND_2, K_BOND_3 };
    double dist, f, dist_0;
    
    switch (interval) {
        case 0:
            
            dist_0 =  part_1->radius + SPB_RADIUS;
            
            dist = Euclid_norm (part_1->position, part_2->position);
            
            f = bonding_power[interval] * (dist_0 - dist) / dist;
            
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
            break;
            
        case 1:
            
            dist_0 = BOND_DISTANCE * interval;
            dist = Euclid_norm (part_1->position, part_2->position);
            
            f = bonding_power[interval] * (dist_0 - dist) / dist;
            
            part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
            part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
            part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
            break;
            
        default:
            
//            if ( part_1->pastis_no >= 0 && part_2->pastis_no >= 0 ) {
//
//                dist_0 = Euclid_norm (part_1->position_init, part_2->position_init);
//                dist = Euclid_norm (part_1->position, part_2->position);
//
//                f = bonding_power[interval] * (dist_0 - dist) / dist;
//
//                part_1->force[X] += f * (part_1->position[X] - part_2->position[X]);
//                part_1->force[Y] += f * (part_1->position[Y] - part_2->position[Y]);
//                part_1->force[Z] += f * (part_1->position[Z] - part_2->position[Z]);
//            }
            break;
    }
}


void membrane_interaction ( Particle *part_1, char interaction_type /* F: fix, E: exclude */) {
    
    double dist = Euclid_norm (part_1->position, ORIGIN);
    
    double ellipsoid_dist = part_1->position[X] * part_1->position[X] / ( MEMBRANE_AXIS_1 * MEMBRANE_AXIS_1 )
    + part_1->position[Y] * part_1->position[Y] / ( MEMBRANE_AXIS_2 * MEMBRANE_AXIS_2 )
    + part_1->position[Z] * part_1->position[Z] / ( MEMBRANE_AXIS_3 * MEMBRANE_AXIS_3 );
    
    if ( interaction_type == 'F' || ellipsoid_dist - 1 > 0 ) {
        
        // 法線ベクトル
        double normal_vector[] = { 2.0 * part_1->position[X] / ( MEMBRANE_AXIS_1 * MEMBRANE_AXIS_1),
            2.0 * part_1->position[Y] / ( MEMBRANE_AXIS_2 * MEMBRANE_AXIS_2),
            2.0 * part_1->position[Z] / ( MEMBRANE_AXIS_3 * MEMBRANE_AXIS_3) };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * ( part_1->position[X] * normal_vector[X]
                                                                  + part_1->position[Y] * normal_vector[Y]
                                                                  + part_1->position[Z] * normal_vector[Z]);
        
        part_1->force[X] += f * normal_vector[X] / normal_vector_norm;
        part_1->force[Y] += f * normal_vector[Y] / normal_vector_norm;
        part_1->force[Z] += f * normal_vector[Z] / normal_vector_norm;
    }
}

// nucleolus interaction //
void nucleolus_interaction ( Particle *part_1, Nuc *nuc, const char interaction_type ) {
    
    //核小体中心から粒子へのベクトル
    double nuc_to_pos[DIMENSION] = { part_1->position[X] - nuc->position[X],
        part_1->position[Y] - nuc->position[Y],
        part_1->position[Z] - nuc->position[Z]};
    
    // 核小体座標系に変換 //
    //位置座標をx-y平面で-theta回転
    rotate_about_z (nuc_to_pos, nuc->theta );
    // x-z平面
    rotate_about_y (nuc_to_pos, nuc->phi);
    // y-z平面
    rotate_about_x (nuc_to_pos, nuc->eta);
    
    double ellipsoid_dist =  nuc_to_pos[X] * nuc_to_pos[X] / ( NUCLEOLUS_AXIS_1 * NUCLEOLUS_AXIS_1 )
    + nuc_to_pos[Y] * nuc_to_pos[Y] / ( NUCLEOLUS_AXIS_3 * NUCLEOLUS_AXIS_3 )
    + nuc_to_pos[Z] * nuc_to_pos[Z] / ( NUCLEOLUS_AXIS_2 * NUCLEOLUS_AXIS_2 );
    
    if ( interaction_type == 'F' || ellipsoid_dist < 1.0 ) {
        
        // 法線ベクトル @核小体座標系
        double normal_vector[] = { 2.0 * nuc_to_pos[X] / ( NUCLEOLUS_AXIS_1 * NUCLEOLUS_AXIS_1),
            2.0 * nuc_to_pos[Y] / ( NUCLEOLUS_AXIS_3 * NUCLEOLUS_AXIS_3),
            2.0 * nuc_to_pos[Z] / ( NUCLEOLUS_AXIS_2 * NUCLEOLUS_AXIS_2) };
        
        double normal_vector_norm = Euclid_norm (normal_vector, ORIGIN);
        
        double f = - ( ellipsoid_dist - 1 ) * MEMBRANE_EXCLUDE * Inner_product (nuc_to_pos, normal_vector)
        / normal_vector_norm;
        
        rotate_about_x (normal_vector, -nuc->eta);
        rotate_about_y (normal_vector, -nuc->phi);
        rotate_about_z (normal_vector, -nuc->theta);
        
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

// 各stepごとの座標計算 //
void calculation (Particle *part, Nuc *nuc, Particle *spb, const unsigned int mitigation ) {
    
    unsigned int loop;
    Particle *part_1, *part_2, *part_3;
    
    // 位置の計算 & 力の初期化 //
    for ( loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        
        part_1->force[X] = 0.0;
        part_1->force[Y] = 0.0;
        part_1->force[Z] = 0.0;
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
                        spring (part_1, &part[loop + 3], 3);
                        break;
                    
                    case TELO1_UP + 2:
                    case TELO2_UP + 2:
                    case rDNA_UP + 2:
                        
                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);
                        
                        spring (part_1, &part[loop + 3], 3);
                        break;
                    
                    case TELO1_DOWN - 1:
                    case TELO2_DOWN - 1:
                    case rDNA_DOWN - 1:
                        
                        spring (part_1, &part[loop - 2], 2);
                        
                        spring (part_1, &part[loop - 3], 3);
                        break;
                        
                    case TELO1_DOWN - 2:
                    case TELO2_DOWN - 2:
                    case rDNA_DOWN - 2:
                        
                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);
                        
                        spring (part_1, &part[loop - 3], 3);
                        break;
                    
                    default:
                        
                        spring (part_1, &part[loop + 2], 2);
                        spring (part_1, &part[loop - 2], 2);
                        
                        spring (part_1, &part[loop + 3], 3);
                        spring (part_1, &part[loop - 3], 3);
                        break;
                }
                
                spb_exclusion (part_1, spb);
                nucleolus_interaction (part_1, nuc, 'E');
                membrane_interaction (part_1, 'E');
                
                break;
            
            case Centromere:
                
                spring (part_1, &part[loop + 1], 1);
                spring (part_1, &part[loop - 1], 1);
                
                spring (part_1, &part[loop + 2], 2);
                spring (part_1, &part[loop - 2], 2);
                
                spring (part_1, &part[loop + 3], 3);
                spring (part_1, &part[loop - 3], 3);
                
                nucleolus_interaction (part_1, nuc, 'E');
                membrane_interaction (part_1, 'E');
                spring (part_1, spb, 0);
                
                break;
                
            case Telomere:
                
                switch (loop) {
                    case TELO1_UP:
                    case TELO2_UP:
                        
                        spring (part_1, &part[loop + 1], 1);
                        
                        spring (part_1, &part[loop + 2], 2);
                        
                        spring (part_1, &part[loop + 3], 3);
                        
                        break;
                        
                    case TELO1_DOWN:
                    case TELO2_DOWN:
                        
                        spring (part_1, &part[loop - 1], 1);
                        
                        spring (part_1, &part[loop - 2], 2);
                        
                        spring (part_1, &part[loop - 3], 3);
                        
                        break;
                }
                
                spb_exclusion (part_1, spb);
                membrane_interaction (part_1, 'F');
                nucleolus_interaction (part_1, nuc, 'E');
                
                break;
                
            case rDNA:
                
                switch (loop) {
                    case rDNA_UP:
                        
                        spring (part_1, &part[loop + 1], 1);
                        
                        spring (part_1, &part[loop + 2], 2);
                        
                        spring (part_1, &part[loop + 3], 3);
                        
                        break;
                        
                    case rDNA_DOWN:
                        
                        spring (part_1, &part[loop - 1], 1);
                        
                        spring (part_1, &part[loop - 2], 2);
                        
                        spring (part_1, &part[loop - 3], 3);
                
                        break;
                }
                
                spb_exclusion (part_1, spb);
                membrane_interaction (part_1, 'E');
                nucleolus_interaction (part_1, nuc, 'F');
                
                break;
                
            default:
                printf ("\t Labeling error occured. \n");
                exit(1);
        }
        
        if ( mitigation % LIST_INTERVAL == 0 ) make_ve_list (part, part_1, loop);
        particle_exclusion (part, part_1);
        
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
    
    if (operation == 'c') {
        
        for (loop = TELO_LIST[0]; loop < CENT_LIST[0]; loop++) part [loop].radius += diff_radius [0];
        for (loop = CENT_LIST[0] + 1; loop <= TELO_LIST[1]; loop++) part [loop].radius += diff_radius [1];
        
        for (loop = TELO_LIST[2]; loop < CENT_LIST[1]; loop++) part [loop].radius += diff_radius [2];
        for (loop = CENT_LIST[1] + 1; loop <= TELO_LIST[3]; loop++) part [loop].radius += diff_radius [3];
        
        for (loop = rDNA_LIST[0]; loop < CENT_LIST[2]; loop++) part [loop].radius += diff_radius [4];
        for (loop = CENT_LIST[2] + 1; loop <= rDNA_LIST[1]; loop++) part [loop].radius += diff_radius [5];
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

// 初期データの出力 //
void write_init_coordinate (Particle *part) {
    
    unsigned int loop;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "result_init.txt");
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \t error : cannot write coordinate. \n");
        
        exit (1);
    }
    
    for (loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        fprintf (fpw, "%d %d %d %lf %lf %lf %lf\n", loop, part_1->chr_no, part_1->particle_type, part_1->position[X],
                 part_1->position[Y], part_1->position[Z], part_1->radius);
    }
    
    fclose (fpw);
}

void write_coordinate (Particle *part, const unsigned int time) {
    
    unsigned int loop;
    
    Particle *part_1;
    
    FILE *fpw;
    
    char result[128], str[128];
    
    sprintf (result, "result_%d.txt", time);
    
    if ((fpw = fopen (result, "w")) == NULL) {
        
        printf (" \t error : cannot write coordinate. \n");
        
        exit (1);
    }
    
    for (loop = 0; loop < NUMBER_MAX; loop++) {
        
        part_1 = &part[loop];
        fprintf (fpw, "%d %d %d %lf %lf %lf %lf\n", loop, part_1->chr_no, part_1->particle_type, part_1->position[X],
                 part_1->position[Y], part_1->position[Z], part_1->radius);
    }
    
//    fprintf (fpw, "Radius %1.1e\n", LENGTH);
    fclose (fpw);
}

int main ( int argc, char **argv ) {
    
    unsigned int loop, mitigation, start, calculation_max, stable_no, sample_no = 0;
    char output_file[256];
    double mem_al[3];
    
    Particle *part, *part_1, *spb;
    Nuc *nuc;
    
    secure_main_memory (&part, &nuc, &spb);
    
    dsfmt_t dsfmt;
    
    
    if ( argc == 3 ) {
        
        start = 0;
        calculation_max = atoi (argv[1]);
        stable_no = atoi (argv[2]);
        
        dsfmt_init_gen_rand (&dsfmt, sample_no);
        
//        type_labeling (part);
        Read_structure (nuc, spb, stable_no);
        Particle_initialization (part, nuc, spb, &dsfmt);
        
        write_coordinate (part, 0);
    }
    else if (argc == 4 ) {
        
        start = atoi (argv[1]);
        calculation_max = atoi (argv[2]);
        stable_no = atoi (argv[3]);
    }
    else {
        
        printf ("\t error : Number of arguments error \n");
        exit (1);
    }
    
    update_radius (part, 's');
    
    for ( unsigned int time = 1; time <= calculation_max; time++) {

        printf ("\t Now calculating...  time = %d \r", time);
        fflush (stdout);

        for ( mitigation = 0; mitigation < MITIGATION_INTERVAL; mitigation++ ){

            calculation (part, nuc, spb, mitigation);
        }

        write_coordinate (part, start + time);

        update_radius (part, 'c');
    }
    
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



