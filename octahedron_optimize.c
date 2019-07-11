// 八面体シミュレーション

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define PI (M_PI)

#define DIMENSION (3)       // 次元
#define LENGTH (7.0e-8)     // 長さ単位
#define DELTA (1.0e-4)

#define SIZE (6)    // 各構造体を構成する粒子数

#define MEMBRANE_AXIS_1 ( 1.889011e-6 / LENGTH )
#define MEMBRANE_AXIS_2 ( 0.85 * MEMBRANE_AXIS_1 )
#define MEMBRANE_AXIS_3 ( 0.75 * MEMBRANE_AXIS_1 )

#define NUCLEOLUS_AXIS_1 ( 1.1426593e-6 / LENGTH )
#define NUCLEOLUS_AXIS_2 ( 0.9 * NUCLEOLUS_AXIS_1 )
#define NUCLEOLUS_AXIS_3 ( 0.8 * NUCLEOLUS_AXIS_1 )

typedef struct object {           //構造体の型宣言
//    int pastis_no;
//    unsigned int chr_no;
//    unsigned int particle_type;
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_init[DIMENSION];
    double velocity[DIMENSION];
    //double velocity_h[DIMENSION];
    double force[DIMENSION];
//    unsigned int list_no;
//    unsigned int *list;
    double equi_length[6]   // 自然長リスト (八面体を保つ)
    
} Object;

typedef struct spot {
    
    double position[DIMENSION];
    double position_new[DIMENSION];
    double position_init[DIMENSION];
    double velocity[DIMENSION];
    double force[DIMENSION];
};


enum label{ X, Y, Z};
typedef enum ellipsoid_term {
    
    a1_r, a1_l, a2_r, a2_l, a3_r, a3_l
}ELLIP_TERM;

void secure_main_memory (Object **obj) {   // メモリ確保 //
    
    if ( (*obj = (Object *)malloc(SIZE * sizeof(Object))) == NULL) {
        
        printf("\n error : can not secure the memory \n");
        exit(1);
    }
}

void def_membrane_pos (Object *membrane) {
    
    
}

int main ( int argc, char **argv ){
    
    Object *membrane, *nucleolus;
    
    
    secure_main_memory (&membrane);
    secure_main_memory (&nucleolus);
    
    // 核小体初期配置を与える or 座標読み込み (核膜粒子の座標を与える)
    def_membrane_pos ( membrane);
    // 計算
    
    // 書き込み
    
}
