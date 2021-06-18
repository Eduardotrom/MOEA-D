#include <bits/stdc++.h>
#include "Toolkit/ExampleProblems.h"
#include "Toolkit/TransFunctions.h"
using namespace std;
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;

int main()
{
   srand(time(NULL));
   int nobj=2, nvar =24;
   vector<double> y_obj, x_var(nvar);
   for(int i =0; i < 100; i++) //random points...
   {
      for(int i = 0; i < nvar; i++) x_var[i] = (rand()/(double)RAND_MAX)*2.0*(i+1.0);
      int k=4; //position variables
       y_obj = Problems::WFG2( x_var, k, nobj);
       //y_obj = Problems::WFG8( x_var, k, nobj);
      for(auto obj:y_obj) cout << obj <<" ";
      cout <<endl;
   }
  return 0;
}

