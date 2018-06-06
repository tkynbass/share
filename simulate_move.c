//
//  simulate_move.c
//  
//
//  Created by tkym on 2018/06/06.
//

#include "common.h"
#include "func_move.h"

int main (int argc, char **argv) {
    
    printf ("\n     Simulation      \n");
    fflush (stdout);
    
    Init_particle (start_no);
    
    //dSFMT seed
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, (unsigned)time(NULL));
    
    //first calculation
    init_particle_calculate (dsfmt);
    init_SPB_calculate (dsfmt);
    
    // output result of first calculation
    write_coordinate (0, start_no);
    
    unsigned int time, loop_count;
    
    for ( time = 1; time < save_max; time++ ) {
        
        printf(" now time = %d\r", time);
        fflush (stdout);
        
        for (loop_count = 1; loop_count <= loop; loop_count++) {
            
            particle_calculate (dsfmt, loop_count);
            SPB_calculate (dsfmt, loop_count);
            
            renew ();
        }
        
        write_coordinate (time, start_no);
    }
    
    return (0);
}
