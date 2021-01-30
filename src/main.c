/* main.c */

#include <time.h>
#include "init_algo.h"
#include "step_1.h"
#include "step_2.h"
#include "fact.h"
#include "test.h"

int main(int argc, char **argv) {
    char *file_name = "70_s_fb.txt";
    test_s_fb(file_name, 0, 50); 
    
    /*
    int lp_var  = 0; 
    int eas     = 0; 
    fact_F7(lp_var, eas);
    */

    return 0; 
}
