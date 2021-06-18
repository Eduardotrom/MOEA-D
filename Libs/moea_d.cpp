#include "moea_d.hpp"
using namespace std;
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
void solution::set_mod(vec_db &mod){
    this->mod=mod;
}

solution::solution(int dims,int sub_p){
    this->mod=vec_db(dims,0);
    this->fitness=vec_db(sub_p,0);
    this->modSize=dims;
    this->functionTotal=sub_p;
}
void solution::set_fit(vec_db &nfit){
    this->fitness=nfit;
}
int solution::get_size(){
    return this->mod.size();
}
vec_db solution::get_mod(){
    return this->mod;
}
vec_db solution::get_fit(){
    return this->fitness;
}
bool solution::operator<(solution &s2){
    return (this->fitness<s2.fitness);
}
bool solution::operator == (solution &s2){
    vec_db fit=s2.get_fit();
    if(this->fitness.size()!=fit.size()){
        printf("La igualdad no se puede calcular\n");
        exit(0);
    }  
    bool test=true;
    for(int i=0;i<fit.size();i++){ 
        if(abs(this->fitness[i]-fit[i])>EPS)return false;
    }
    return true;
}

MOEA_D::MOEA_D(int dims,int sbp,int T,int pobs,mat_db limits){
    this->T=T;
    this->dims=dims;
    this->sub_p=sbp;
    this->limits=limits;
    this->pobs=pobs;
    //std::vector<solution>(pobs,solution(dims,this->sub_p));
    this->B_fill(sub_p);
    this->Dist_calc();
}


void MOEA_D::populate(){
    this->poblacion.resize(this->pobs);
    vec_db temporal(this->dims);
    for(int i=0;i<this->pobs;i++){
        for(int j=0;j<this->dims;j++){
            double dist=this->limits[j][1]-this->limits[j][0];
            temporal[j]=this->limits[j][0]+RANDU*dist;
        }
        poblacion[i].set_mod(temporal);
    }
}

void MOEA_D::B_fill(int sub_p){
    this->B=mat_db(this->pobs,vec_db(sub_p,0));
    for(int i=0;i<B.size();i++){
        for(int j=0;j<B[i].size();j++){ 
            B[i][j]=rand();
        }
        normalize(B[i]);
    }
}

void MOEA_D::Dist_calc(){
    int n=this->pobs;
    this->distancias=mat_db(n,vec_db(n,0));
    for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
            distancias[i][j]=dist(B[i],B[j]);
            distancias[j][i]=distancias[i][j];
        }
    }
    #ifdef DEBUG
        std::cout<<"Las distancias son: "<<std::endl;
        show(distancias);
    #endif
    this->Bi=mat_int(this->pobs,vec_int(this->T,0));
    vec_db temp=distancias[0];
    for(int i=0;i<this->pobs;i++){
        temp=distancias[i];
        std::sort(temp.begin(),temp.end());
        Bi[i].clear();
        double dmax=temp[this->T];
        #ifdef DEBUG
            printf("El vector ordenado es:\n");
            show(temp);
        #endif
        for(int j=0;j<distancias[i].size();j++){
            if(i==j)continue;
            if(dmax>=distancias[i][j])Bi[i].push_back(j);
        }
        #ifdef DEBUG
            printf("El dmax es %lf\n",dmax);
            printf("los indices para lamb%d son:\n",i);
            show(Bi[i]);
        #endif
        Bi[i].shrink_to_fit();
    }
}

//Cruza utilizando operador de Evolucion diferencial
solution MOEA_D::cruza(solution &s1, solution&s2){
    int n=s1.get_size();
    double a=-0.5+2*RANDU;
    vec_db x(n,0),y(n,0);
    for(int i=0;i<n;i++){
        x[i]=a*s1.mod[i]+(1.0-a)*s2.mod[i];
    } 
    solution res(n,s1.get_fit().size());
    res.set_mod(x);
    range_check(res);
    #ifdef DEBUG
        printf("n=%d\n",n);
        printf("El tamanio del hijo es de %d\n",res.get_mod().size());
    #endif
    return res;
}
double MOEA_D::gte(vec_db &z,vec_db &lamb, vec_db &y){
    double res=MAXFLOAT;
    for(int i=0;i<y.size();i++)
        if(res<lamb[i]*y[i]-z[i])res=lamb[i]*y[i]-z[i];
    return res;
}
std::vector<solution> MOEA_D::get_pob(){
    return this->poblacion;
}
void MOEA_D::range_check(solution &sol){
    int n=sol.get_size();
    for(int i=0;i<n;i++){ 
        bool rl=(this->limits[i][0]>=sol.mod[i]);
        bool ll=(this->limits[i][1]<=sol.mod[i]);
        if(!(rl&&ll)){
            sol.mod[i]=limits[i][0]+RANDU*(limits[i][1]-limits[i][0]);
        }
    }
}

void MOEA_D::Optimiza(int max_gen){
    this->EP.reserve(this->pobs);
    this->FV.resize(this->pobs);
    vec_db z(this->sub_p);
    int k=4;
    auto funcion=&Problems::WFG2;
    #ifdef WFG8F
        funcion=&Problems::WFG8;
    #endif
    vec_db fft,fft2;
    for(int i=0;i<this->pobs;i++){
        fft=poblacion[i].get_mod();
        #ifdef DEBUG 
            printf("Funcion 1\n"); 
        #endif
        FV[i]=funcion(fft,k,this->sub_p);
        this->poblacion[i].set_fit(FV[i]);
        for(int j=0;j<this->sub_p;j++){
            if(z[i]>FV[i][j])z[i]=FV[i][j];
        }
    }
    solution newsol;
    for(int gen=0;gen<max_gen;gen++){
        for(int i=0;i<this->pobs;i++){
            int i1=rand()%(this->T-1),i2=rand()%(this->T-1);
            i1=Bi[i][i1];
            i2=Bi[i][i2];
            newsol=cruza(this->poblacion[i1],this->poblacion[i2]);
            fft2=newsol.get_mod();
            #ifdef DEBUG 
                printf("Funcion 2\n"); 
                printf("El tamanio de fft es de %d\n",fft2.size());
            #endif
            fft=funcion(fft2,k,this->sub_p);
            newsol.set_fit(fft);
            for(int j=0;j<this->sub_p;j++)
                if(z[j]>fft[j])z[j]=fft[j];
            for(int j:Bi[i]){
                auto fftt=this->poblacion[j].get_fit();
                if(gte(z,B[j],fft)<=gte(z,B[j],fftt)){
                    this->poblacion[j]=newsol;
                    FV[j]=fft;
                }
            }
            int dom=NON_DOMINATION,domm;
            for(int j=0;j<EP.size();j++){
                domm=dominancia(newsol,EP[j]);
                if(domm==LEFT_DOMINATES){
                    std::swap(EP[j],EP[EP.size()-1]);
                    EP.pop_back();
                }
                if(domm==RIGHT_DOMINATES)dom=RIGHT_DOMINATES;
                if(domm==EQUAL)dom=EQUAL;
            }
            if(dom==NON_DOMINATION)EP.push_back(newsol);        }
    }
    auto EP_T=EP;
    pareto_front(EP_T,EP);
}

vector<solution> MOEA_D::get_front(){
    return this->EP;
}
void MOEA_D::Show_B(){
    show(this->B);
}

void MOEA_D::Show_Pop(){
    for(auto x:this->poblacion)
        show(x.mod);
}

int dominancia(solution & s1, solution &s2){
    if(s1.get_size()!=s2.get_size()){
        printf("Los tamaños de los elementos a revisar son incompatibles}n");
        exit(0);
    }
    if (s1==s2)return EQUAL;
    bool p1=true,p2=true;
    vec_db fit1=s1.get_fit(),fit2=s2.get_fit();
    int n=fit1.size();
    for(int i=0;(i<n)&&(p1||p2);i++){
        bool t1=(fit1[i]<fit2[i]);
        bool t2=(fit1[i]==fit2[i]);
        p1=p1&&(t1||t2);
        p2=p2&&(!t1);
        #ifdef DEBUG
            printf("Para el elemento %d los bool son:\n",i);
            printf("Los valores son %lf %lf\n",fit1[i],fit2[i]);
            printf("t1=%d t2=%d p1=%d p2=%d\n",t1,t2,p1,p2);
        #endif
    }
    if(p1)return LEFT_DOMINATES;
    if(p2)return RIGHT_DOMINATES;
    return NON_DOMINATION;
}
void fitprint(std::string nom,std::vector<solution>& vec){
    ofstream f(nom);
    std::sort(vec.begin(),vec.end());
    for(solution x:vec){
        vec_db fit=x.get_fit();
        for(int i=0;i<fit.size();i++){
            f<<fit[i]<<" ";
        }f<<std::endl;
    }
    f.close();
}

void pareto_front(std::vector<solution> &s1,std::vector<solution> &pareto){
    pareto.clear();
    pareto.reserve(s1.size());
    std::sort(s1.begin(),s1.end());
    divide(0,s1.size()-1,s1,pareto);
}

void divide(int init,int end,std::vector<solution> &solutions,std::vector<solution> &pareto){
    if(end-init==1){
        int dm=dominancia(solutions[init],solutions[end]);
        switch(dm){
            case LEFT_DOMINATES:
                pareto.push_back(solutions[init]);
                break;
            case RIGHT_DOMINATES:
                pareto.push_back(solutions[end]);
                break;
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


void check_last(std::vector<solution> &pareto,int last_size){
      if(last_size!=0){
        auto starp=pareto.begin()+last_size;
        auto test=pareto[last_size-1];
        for(int i=pareto.size()-1;i>last_size-1;i--){
            int a=dominancia(test,pareto[i]);
            if(a==LEFT_DOMINATES||a==EQUAL){
                std::swap(pareto[i],*(pareto.end()-1));
                pareto.pop_back();
            }
        }
        std::sort(starp,pareto.end());
    }
}