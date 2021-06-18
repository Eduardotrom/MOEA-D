#ifndef AUX_STRUCTS_HPP
#define AUX_STRUCTS_HPP
#include<vector>
#include<tuple>
#include<fstream>
#include<cstring>
#include<cmath>
typedef std::tuple<int,int> tii;
typedef std::tuple<int,int,int> tiii;
typedef std::tuple<float,float> tff;
typedef std::tuple<double,double> tdd;
typedef std::vector<int> vec_int;
typedef std::vector<double> vec_db;
typedef std::vector<vec_int> mat_int;
typedef std::vector<std::vector<double>> mat_db;

template <class T>
class SegmentTree{
    private:
        int n;
        std::vector<T>  vi,A;
        std::vector<int> st;
        int left(int p);
        int right(int p);
        void build(int p, int L,int R);
        int rmq(int p,int L,int R, int i,int j);
    public:
        SegmentTree(std::vector<T> & A_i);
        int rmq(int i,int j);
        void update(int p,T val);
        void update(int p, T val, int i,int j);
        T rmq_v(int i,int j);
};

void printf_file(std::string nom,const std::vector<tii> &v);
void printf_file(std::string nom,const std::vector<tiii> &v);
void printf_file(std::string nom,const std::vector<tff> &v);
int curv(int x,int n);
int curv(int x, int y, int n);
std::vector<tiii> curv_nd(int n);
std::vector<tiii> curv_rand(int n);
std::vector<tiii> curv_dat(int n);
std::vector<tii> curv_nd2(int n);
std::vector<tii> curv_rand2(int n);
std::vector<tii> curv_dat2(int n);
void normalize(vec_db &v);
double dist(vec_db const &v1,vec_db const &v2);
void show(mat_db const &mat);
void show(vec_db const &vec);
void show(vec_int const &vec);
vec_db algo(vec_db &x);
#endif