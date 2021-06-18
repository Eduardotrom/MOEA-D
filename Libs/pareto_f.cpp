#include "pareto_f.hpp"
#include<algorithm>
#include<cmath>
#include<iostream>

bool comp_2(const tiii &a, const tiii &b){
    return (std::get<1>(a)<std::get<1>(b));
}

bool comp_3(const tiii &a, const tiii &b){
    return (std::get<2>(a)<std::get<2>(b));
}

int dominancia(const tii &a,const tii &b){
    if(a==b)return EQUAL;
    bool p1= std::get<0>(a)<std::get<0>(b);
    bool p2= std::get<1>(a)<std::get<1>(b);
    bool p4= std::get<0>(a)==std::get<0>(b);
    bool p5= std::get<1>(a)==std::get<1>(b);
    if((p1||p4)&&(p2||p5))return FIRST_DOMINATES;
    if(((!p1||p4)&&(!p2||p5)))return SECOND_DOMINATES;
    return NON_DOMINATION;    
}
int dominancia(const tiii &a,const tiii &b){
    if(a==b)return EQUAL;
    bool p1= std::get<0>(a)<std::get<0>(b);
    bool p2= std::get<1>(a)<std::get<1>(b);
    bool p3= std::get<2>(a)<std::get<2>(b);
    bool p4= std::get<0>(a)==std::get<0>(b);
    bool p5= std::get<1>(a)==std::get<1>(b);
    bool p6= std::get<2>(a)==std::get<2>(b);
    if((p1||p4)&&(p2||p5)&&(p3||p6))return FIRST_DOMINATES;
    if((!(p1||p4)&&!(p2||p5)&&!(p3||p6)))return SECOND_DOMINATES;
    return NON_DOMINATION;
}


std::vector<tii> dominancia_ineficiente(std::vector<tii> &a){
    std::vector<tii> res,aux=a;
    std::sort(aux.begin(),aux.end());
    tii temp;
    int n=aux.size();
    res.reserve(n);
    res.push_back(aux[0]);
    for(int i=1; i<n ; i++){ 
        int indicador=FIRST_DOMINATES;
        for(int j=0; j<n ;j++){ 
            if(i==j)continue;
            indicador=dominancia(aux[i],aux[j]);
            if(indicador==EQUAL)break;
            if(indicador==SECOND_DOMINATES)break;
        }
        if(SECOND_DOMINATES==indicador||indicador==EQUAL)continue;
        res.push_back(aux[i]);
    }
    res.shrink_to_fit();
    return res;
}

void dominancia_eficiente(std::vector<tii> solutions,std::vector<tii> &pareto){
    pareto.clear();
    pareto.reserve(solutions.size());
    std::sort(solutions.begin(),solutions.end());
    divide(0,solutions.size()-1,solutions,pareto);
    return;
}

void check_last(std::vector<tii> &pareto,int last_size){
    if(last_size!=0){
        auto starp=pareto.begin()+last_size;
        auto test=pareto[last_size-1];
        for(int i=pareto.size()-1;i>last_size-1;i--){
            int a=dominancia(test,pareto[i]);
            if(a==FIRST_DOMINATES||a==EQUAL){
                std::swap(pareto[i],*(pareto.end()-1));
                pareto.pop_back();
            }
        }
        std::sort(starp,pareto.end());
    }
}

void divide(int init,int end,const std::vector<tii> &solutions,std::vector<tii> &pareto){
    if(end-init==1){
        int dm=dominancia(solutions[init],solutions[end]);
        switch (dm)
        {
        case FIRST_DOMINATES:
            pareto.push_back(solutions[init]);
            break;
        case SECOND_DOMINATES:
            pareto.push_back(solutions[end]);
        case EQUAL:
            pareto.push_back(solutions[init]);
            break;  
        case NON_DOMINATION:
            pareto.push_back(solutions[init]);
            pareto.push_back(solutions[end]);
            break;
        default:
            break;
        }
    }else if(end-init==0){
        pareto.push_back(solutions[init]);
    }else{
        int n=init+std::ceil((end-init)/2.0);
        divide(init,n,solutions,pareto);
        int last=pareto.size();
        divide(n,end,solutions,pareto);
        check_last(pareto,last);
    }
}

void Coord_compression(std::vector<tiii> &a){
    std::sort(a.begin(),a.end(),comp_2);
    int cont=0;
    int temp=std::get<1>(a[0]);
    for(int i=0;i<a.size();i++){ 
        if(temp!=std::get<1>(a[i])){
            cont++;
            temp=std::get<1>(a[i]);
        }
        std::get<1>(a[i])=cont;
    }
        std::sort(a.begin(),a.end());
}
std::vector<tiii> dominancia_ineficiente(std::vector<tiii> &a){
    std::vector<tiii> res,aux=a;
    std::sort(aux.begin(),aux.end());
    tii temp;
    int n=aux.size();
    res.reserve(n);
    res.push_back(aux[0]);
    for(int i=1; i<n ; i++){ 
        int indicador=FIRST_DOMINATES;
        for(int j=0; j<n ;j++){ 
            if(i==j)continue;
            indicador=dominancia(aux[i],aux[j]);
            if(indicador==EQUAL)break;
            if(indicador==SECOND_DOMINATES)break;
        }
        if(SECOND_DOMINATES==indicador||indicador==EQUAL)continue;
        res.push_back(aux[i]);
    }
    res.shrink_to_fit();
    return res;
    
}

void dominancia_eficiente(std::vector<tiii> &solutions, std::vector<tiii> &pareto){
    pareto.clear();
    int n=solutions.size();
    pareto.reserve(n);
    Coord_compression(solutions);
    vec_int aux(n,INT8_MAX);
    SegmentTree<int> tree(aux);
    for(auto x:solutions){
        int temp=std::get<1>(x);
        if(std::get<2>(x)<tree.rmq_v(0,temp)){;
            tree.update(temp,std::get<2>(x));
            pareto.push_back(x);
        }
    }
}
