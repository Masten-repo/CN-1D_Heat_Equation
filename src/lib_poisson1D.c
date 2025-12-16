/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){

int k = 0;
int j = 0;
int m = *la ; // colonnes
int n = *lab ; // lignes
int v = *kv ; // index diagonale 

// On parcours les colonne de la matrice 
for(int i = 0; i<m;++i){
  k=i*n;
  if(v >= 0 ){
    for(int j = 0 ; j< v ; ++j){
      AB[k+j] =0.0;
    }
  }

  AB[k+v]=-1.0; // On met les élements de la surdiagonale à -1.0
  AB[k+v+1]= 2.0; // On met les éléments de la diagonale principale a 2.0
  AB[k+v+2]=-1.0; // On met les éléments de la sous diagonale a -1.0

}

// On remplit le reste de la matrice avec des 0
AB[0] = 0.0;
  if(v == 1){
    AB[1]=0 ;
  }
AB[n*m-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){

int k = 0;
int m = *la ;
int n = *lab ;
int v = *kv ;
  // On parcours les colonnes de la matrice 
  for(int i = 0 ; i<m ;++i){
    k=i*n ;

    if(v >=0){
      for(int j = 0; j < v ; j++){
        AB[k+j] = 0.0 ;
      }
    }
    AB[k+v]=0.0; // On met les élements de la surdiagonale a 0.0
    AB[k+v+1]=1.0; // On met les élement de la diagonale a 1.0
    AB[k+v+2]=0.0; //On met les élement de la sousdiagonale a 0.0
  }
  AB[1] = 0.0; // deuxieme élément de la matrice a 0
  AB[n*m-1] = 0.0 ; //dernier élement a 0
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int m = *la ;
  // On met les élement (entre premier et dernier) du vecteur RHS a 0
  for (int i = 0 ; i < m ;++i){
    RHS[i] = 0.0 ;
  }
  RHS[0] = *BC0 ; //on met la premiere valeur du vecteur RHS a la valeur BC0
  RHS[m-1] = *BC1 ; // on met la derniere valeur du vecteur RHS a la valeur BC1
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int m = *la ;
  double d = (*BC1)-(*BC0); // calcul du delta 
  //calcul de la solution analytique
  for(int i = 0 ; i < m ;++i){
    EX_SOL[i] = *BC0 + X[i] * d ;
  }
}  

void set_grid_points_1D(double* x, int* la){
  // cas d'erreur
  if (x == NULL || la == NULL || *la <= 0) {
    return;
  }
  int m = *la ;
  double a = 1.0 /(1.0 * m+1);
  for(int i = 0 ; i < m ; ++i){
    x[i] = (i+1)* a ;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  // norme de x
  double norm = cblas_dnrm2(*la,x,1);
  cblas_daxpy(*la,-1.0,y,1,x,1);

  // calcul de l'erreur avant 
  return (cblas_dnrm2(*la,x,1)/norm);
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  double pivot = 0;

  for (int i = 1; i < (*la); i++)
  {
    pivot = AB[i*(*lab)+1]*AB[(i-1)*(*lab)+3]; 
    pivot /= AB[(i-1)*(*lab)+2]; 
    AB[i*(*lab)+2] -= pivot; 
    AB[(i-1)*(*lab)+3] /= AB[(i-1)*(*lab)+2];
  }
  return *info;  return *info;
}
