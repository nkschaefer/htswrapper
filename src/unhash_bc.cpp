#include <stdio.h>   
#include <stdlib.h> 
#include <ctype.h>
#include <string.h>
#include <map>
#include <string>
#include "bc.h"

using namespace std;

int main(int argc, char *argv[]) {   
    
    if (argc < 2){
        fprintf(stderr, "USAGE: unhash_bc [unsigned long int]\n");
        exit(1);
    }    
    unsigned long ul = atol(argv[1]);
    string bc_str = bc2str(ul);
    fprintf(stdout, "%s\n", bc_str.c_str());
}

