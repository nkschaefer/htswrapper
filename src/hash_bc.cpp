#include <stdio.h>   
#include <stdlib.h> 
#include <ctype.h>
#include <string.h>
#include <map>
#include <string>
#include <iostream>
#include "bc.h"

using namespace std;

int main(int argc, char *argv[]) {   
    /*    
    if (argc < 2){
        fprintf(stderr, "USAGE: hash_bc [barcode sequence]\n");
        exit(1);
    }    
    string bc = argv[1];
    unsigned long ul = bc_ul(bc);
    fprintf(stdout, "%ld\n", ul);
    */
    string str;
    while (cin >> str){
        unsigned long ul = bc_ul(str);
        fprintf(stdout, "%ld\n", ul);
    }
}

