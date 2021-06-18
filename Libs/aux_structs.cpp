#include "aux_structs.hpp"
template<typename T>
int SegmentTree<T>:: left(int p) {return p<<1;}

template<typename T>
int SegmentTree<T>::right(int p) {return (p<<1)+1;}
        
template<typename T>
void SegmentTree<T>::build(int p, int L,int R){
    if(L==R)
        st[p]=L;
    else{
        build(left(p),L,(L+R)/2);
        build(right(p),(L+R)/2+1,R);
        int p1=st[left(p)],p2=st[right(p)];
        st[p]=(A[p1]<=A[p2])?p1:p2;
    }
}

template<typename T>
int SegmentTree<T>:: rmq(int p,int L,int R, int i,int j){
    if(i>R|| j<L )return -1;
    if(L>=i && R<=j)return st[p];

    int p1=rmq(left(p), L, (L+R)/2,i,j);
    int p2=rmq(right(p),(L+R)/2+1,R,i,j);
    if(p1==-1)return p2;
    if(p2==-1)return p1;
    return (A[p1]<=A[p2])?p1:p2;
}
template<typename T>
SegmentTree<T>::SegmentTree(std::vector<T> & A_i){
    this->A=A_i;this->n=A_i.size();
    st.assign(4*this->n,0);
    build(1,0,n-1);
}
template<typename T>
int SegmentTree<T>::rmq(int i,int j){return rmq(1,0,this->n-1,i,j);}

template<typename T>
void SegmentTree<T>::update(int p,T val){
    this->A[p]=val;
    this->build(1,0,p);
}
template<typename T>
T SegmentTree<T>::rmq_v(int i,int j){
    return A[this->rmq(i,j)];
}

void printf_file(std::string nom,const std::vector<tii> &v){
    std::ofstream f(nom);
    if(!f.is_open()){
        printf("El archivo %s no se pudo abrir\n",nom.c_str());
        exit(0);
    }
    for(auto x:v)
        f<<std::get<0>(x)<<" "<<std::get<1>(x)<<std::endl;
    f.close();
}
void printf_file(std::string nom,const std::vector<tff> &v){
    std::ofstream f(nom);
    if(!f.is_open()){
        printf("El archivo %s no se pudo abrir\n",nom.c_str());
        exit(0);
    }
    for(auto x:v)
        f<<std::get<0>(x)<<" "<<std::get<1>(x)<<std::endl;
    f.close();
}
void printf_file(std::string nom,const std::vector<tiii> &v){
    std::ofstream f(nom);
    if(!f.is_open()){
        printf("El archivo %s no se pudo abrir\n",nom.c_str());
        exit(0);
    }
    for(auto x:v)
        f<<std::get<0>(x)<<" "<<std::get<1>(x)<<" "<<std::get<2>(x)<<std::endl;
    f.close();
}


//Se crea un conjunto de soluciones aleatorias limitadas por un plano
std::vector<tiii> curv_dat(int n){
    std::vector<tiii> a(n);
    for(int i=0;i<n;i++){
        int x=rand()%n,y=rand()%n;
        int kk=abs(n-curv(x,y,n));
        kk=(kk==0)?0:rand()%kk;
        kk=(abs((int)curv(x,n)+kk)>n)?n:abs((int)curv(x,n)+kk);
        a[i]=std::make_tuple(x,y,kk);
    }
    return a;
}
//Se crea un conjunto de soluciones completamente aleatorias
std::vector<tiii> curv_rand(int n){
    std::vector<tiii> a(n);
    for(int i=0;i<n;i++){
        int x=rand()%n,y=rand()%n,z=rand()%n;
        a[i]=std::make_tuple(x,y,z);
    }
    return a;
}
//Se crea un conjunto de soluciones de no dominados
std::vector<tiii> curv_nd(int n){
    std::vector<tiii> a(n);
    for(int i=0;i<n;i++){
        int x=n-i,y=rand()%n,z=-4*x-3*y-n;
        a[i]=std::make_tuple(x,y,z);
    }
    return a;
}
//Se crea un conjunto de soluciones aleatorias limitadas por un plano
std::vector<tii> curv_dat2(int n){
    std::vector<tii> a(n);
    for(int i=0;i<n;i++){
            int x=rand()%n;
            int kk=abs(n-curv(x,n));
            kk=(kk==0)?0:rand()%kk;
            a[i]=std::make_tuple(x,abs((int)curv(x,n)+kk));
        }
    return a;
}
//Se crea un conjunto de soluciones completamente aleatorias
std::vector<tii> curv_rand2(int n){
    std::vector<tii> a(n);
    for(int i=0;i<n;i++){
        int x=rand()%n,y=rand()%n;
        a[i]=std::make_tuple(x,y);
    }
    return a;
}
//Se crea un conjunto de soluciones de no dominados
std::vector<tii> curv_nd2(int n){
    std::vector<tii> a(n);
    for(int i=0;i<n;i++){
        int x=n-i,y=x*-2+n;
        a[i]=std::make_tuple(x,y);
    }
    return a;
}

int curv(int x,int n){
    return -4*x+n;
}
int curv(int x, int y, int n){
    return -4*x-3*y+n;
}

template class SegmentTree<int>;
template class SegmentTree<float>;
template class SegmentTree<double>;
template class SegmentTree<tdd>;
template class SegmentTree<tff>;
template class SegmentTree<tii>;
template class SegmentTree<tiii>;
//Modifica el vector para que la suma de sus elementos sea 1
void normalize(vec_db &v){
    double sum=0;
    #ifndef DEBUG
    for(int i=0;i<v.size();i++)
        sum+=v[i];
    #else
    printf("El vector a normalizar es: \n");
    for(int i=0;i<v.size();i++){
        sum+=v[i];
        printf("%.3lf ",v[i]);
    }
    printf("\n La suma es %.5lf\n",sum);
    #endif
    for(int i=0;i<v.size();i++)
        v[i]/=sum;
}

void show(mat_db const &mat){
    for(int i=0;i<mat.size();i++){
        for(double x:mat[i])
            printf("%.3lf ",x);
        printf("\n");
    }
}

void show(vec_db const &vec){
    for(double x:vec)
        printf("%.3lf ",x);
    printf("\n");
}
void show(vec_int const &vec){
    for(int x:vec)
        printf("%d ",x);
    printf("\n");
}

double dist(vec_db const &v1,vec_db const &v2){
    int n=std::min(v1.size(),v2.size());
    double dist=0;
    for(int i=0;i<n;i++)
        dist+=(v1[i]-v2[i])*(v1[i]-v2[i]);
    dist=sqrt(dist);
    return dist;
}
