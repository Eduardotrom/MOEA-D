#ifndef PARETO_F_HPP
#define PARETO_F_HPP

#include"aux_structs.hpp"
#define FIRST_DOMINATES 0
#define SECOND_DOMINATES 1
#define EQUAL 2
#define NON_DOMINATION 3

int dominancia(const tii &a,const tii &b);
int dominancia(const tiii &a,const tiii &b);
std::vector<tii> dominancia_ineficiente( std::vector<tii> &a);
std::vector<tiii> dominancia_ineficiente( std::vector<tiii> &a);
void dominancia_eficiente(std::vector<tii> solutions,std::vector<tii> &pareto);
void dominancia_eficiente(std::vector<tiii> &solutions, std::vector<tiii> &pareto);
void check_last(std::vector<tii> &pareto,int last_size);
void divide(int init,int end,const std::vector<tii> &solutions,std::vector<tii> &pareto);
bool comp_2(const tiii &a, const tiii &b);
bool comp_3(const tiii &a, const tiii &b);
#endif