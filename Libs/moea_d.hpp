#ifndef MOEA_D_HPP
#define MOEA_D_HPP
#include<iostream>
#include "aux_structs.hpp"
#include <vector>
#include <algorithm>
#include "Toolkit/ExampleProblems.h"
#include "Toolkit/TransFunctions.h"
#include <fstream>
#define RANDU (rand()/(double)RAND_MAX)
#define LEFT_DOMINATES 0
#define RIGHT_DOMINATES 1 
#define NON_DOMINATION 2
#define EQUAL 3
#define EPS 1e-6
class solution{
    private:
        int modSize,functionTotal;
        vec_db fitness;
    public:
        vec_db mod;
        solution(){};
        solution(int dims,int sub_p);
        void set_mod(vec_db &mod);
        void set_fit(vec_db &nfit);
        vec_db get_mod();
        vec_db get_fit();
        int get_size();
        bool operator == (solution &s2);
        bool operator <  (solution &s2);
};

class MOEA_D{
    private:
        int T,dims,sub_p,pobs;
        mat_db pesos,distancias,limits;
        mat_db B,FV;
        mat_int Bi;
        std::vector<solution> poblacion,EP;
    public:
        MOEA_D(){};
        MOEA_D(int dims,int sbp,int T,int pobs,mat_db limits);
        void populate();
        void B_fill(int sub_p);
        solution cruza(solution &s1, solution&s2);
        void Optimiza(int max_gen);
        void Dist_calc();
        void Show_B();
        void Show_Pop();
        void range_check(solution &sol);
        std::vector<solution> get_pob();
        double gte(vec_db &z,vec_db &lamb, vec_db &y);
        std::vector<solution> get_front();
};


int dominancia(solution & s1, solution &s2);
void fitprint(std::string nom,std::vector<solution>& vec);
void pareto_front(std::vector<solution> &s1,std::vector<solution> &pareto);
void divide(int init,int end,std::vector<solution> &solutions,std::vector<solution> &pareto);
void check_last(std::vector<solution> &pareto,int last_size);
#endif