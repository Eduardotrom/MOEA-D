#include <iostream>
#include "moea_d.hpp"
#include "aux_structs.hpp"

int main(){
    srand(time(NULL));
    int nvars=24,pobs=100;
    mat_db lims;
    lims.reserve(nvars);
    for(int i=0;i<nvars;i++)
        lims.emplace_back(vec_db{0,2.0*(1+i)});
    MOEA_D A(nvars,2,std::ceil(pobs*0.2),pobs,lims);
    A.populate();
    A.Optimiza(100);
    std::vector<solution> PF=A.get_front();
    fitprint("pareto.txt",PF);
    #ifdef DEBUG
    mat_int doms(10,vec_int(10,0));
    PF.resize(10);
    
    for(int i=0;i<10;i++){
        for(int j=0;j<10;j++){
            doms[i][j]=dominancia(PF[i],PF[j]);
        }
    }
    for(int i=0;i<10;i++){
        for(int j=0;j<10;j++){
            std::string bla;
            if(doms[i][j]==LEFT_DOMINATES)bla="First_Dominates";
            if(doms[i][j]==RIGHT_DOMINATES)bla="Second_Dominates";
            if(doms[i][j]==EQUAL)bla="Equal";
            if(doms[i][j]==NON_DOMINATION)bla="Non Domination";
            printf("%17s ",bla.c_str());
        }
        printf("\n");
    }
    printf("\n\n\n\n");
    printf("Last one\n");
    std::string bla;
    int domss=dominancia(PF[0],PF[1]);
    if(domss==LEFT_DOMINATES)bla="First_Dominates";
    if(domss==RIGHT_DOMINATES)bla="Second_Dominates";
    if(domss==EQUAL)bla="Equal";
    if(domss==NON_DOMINATION)bla="Non Domination";
    printf("%17s ",bla.c_str());   
    #endif
    return 0;
}