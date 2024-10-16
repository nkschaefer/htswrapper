#include <stdio.h>   
#include <stdlib.h> 
#include <ctype.h>
#include <string.h>
#include <map>
#include <string>
#include <algorithm>
#include "bc.h"

using namespace std;

int main(int argc, char *argv[]) {   
    bc_whitelist wlx;
    wlx.init("whitelist.txt");

    bool exact;
    unsigned long ulx;
    
    string str1 = "AAACCCAAGAAACCAT";
    string str2 = "AAACCCAAGAAAGCAT";
    string str3 = "AAGATAGAGTTACACT";
    string str4 = "TAGATAGAGTTACACT";
    string str5 = "CTCGAGGTCCTAATGA";
    string str6 = "CTCGAGGTCCGAATGA";
    string str7 = "CTNGAGGTCCTAATGA";
    
    vector<string> strs{ str1, str2, str3, str4, str5, str6, str7 };
    for (int i = 0; i < strs.size(); ++i){
        bool success = wlx.lookup(strs[i].c_str(), ulx, exact);
        fprintf(stderr, "%s: %d %d\n", strs[i].c_str(), success, exact);
    }
}
