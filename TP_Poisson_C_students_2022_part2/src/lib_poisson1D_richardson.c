/**********************************************/
/* lib_poisson1D_richardson.c                 */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <stdlib.h>
#include <math.h>
#include <string.h> 

void eig_poisson1D(double* eigval, int *la){
}

double eigmax_poisson1D(int *la){
  double max = sin(*la * M_PI_2 * (1.0 / (*la + 1)));
  return 4 * max * max;
}

double eigmin_poisson1D(int *la){
  double min = sin(M_PI_2 * (1.0 / (*la + 1)));
  return 4 * min * min;
}

double richardson_alpha_opt(int *la){
  double opt = eigmax_poisson1D(la) + eigmin_poisson1D(la);
  return (2.0/opt);
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double norm_b;
  double norm_r;
  double *r = (double *)malloc(sizeof(double)*(*la));
  
  norm_b = cblas_dnrm2(*la, RHS, 1);
  if (norm_b == 0.0) norm_b = 1.0;

  // r = b - Ax
  cblas_dcopy(*la, RHS, 1, r, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1);
  
  norm_r = cblas_dnrm2(*la, r, 1);
  *nbite = 0;
  resvec[*nbite] = norm_r / norm_b;

  while((resvec[*nbite] > *tol) && (*nbite < *maxit))
  {
    cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
    
    cblas_dcopy(*la, RHS, 1, r, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1);
    
    norm_r = cblas_dnrm2(*la, r, 1);
    (*nbite)++;
    resvec[*nbite] = norm_r / norm_b;
  }
  free(r);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // Jacobi : M = D
  // Pour dgbtrf avec KL=1, KU=0, la diagonale doit être à la ligne 1.
  // Ligne 0 est pour le fill-in.
  int n = *la;
  int ld = *lab; // ld doit être au moins 3 (2*kl + ku + 1 = 2*1+0+1 = 3)

  for(int j=0; j<n; j++){
      int k_ab = j*ld; // AB a KU=1 -> Diag en ligne 1
      int k_mb = j*ld; // MB a KU=0, KL=1 -> Diag en ligne 1 (après ligne 0 de garde)
      
      MB[k_mb + 0] = 0.0;          // Ligne 0 : Garde pour dgbtrf
      MB[k_mb + 1] = AB[k_ab + 1]; // Ligne 1 : Diagonale
      MB[k_mb + 2] = 0.0;          // Ligne 2 : Sous-diagonale (0 pour Jacobi)
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // Gauss-Seidel : M = D - E (Diag + Sous-Diag)
  // Pour dgbtrf avec KL=1, KU=0 :
  // Ligne 0 : Fill-in
  // Ligne 1 : Diagonale
  // Ligne 2 : Sous-diagonale
  int n = *la;
  int ld = *lab;
  
  for(int j=0; j<n; j++){
      int k_ab = j*ld;
      int k_mb = j*ld;
      
      MB[k_mb + 0] = 0.0;          // Ligne 0 : Garde
      MB[k_mb + 1] = AB[k_ab + 1]; // Ligne 1 : Diagonale
      MB[k_mb + 2] = AB[k_ab + 2]; // Ligne 2 : Sous-diagonale
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double norm_b, norm_r;
  int info;
  int nrhs = 1;
  int ku_precond = 0; 

  double *r = (double *)malloc(sizeof(double)*(*la));
  int *ipiv = (int *)malloc(sizeof(int)*(*la));

  dgbtrf_(la, la, kl, &ku_precond, MB, lab, ipiv, &info);
  
  if(info != 0){
      *nbite = 0;
      resvec[0] = 1.0; 
      free(r); free(ipiv);
      return;
  }
  
  norm_b = cblas_dnrm2(*la, RHS, 1);
  if (norm_b == 0.0) norm_b = 1.0;

  // Calcul résidu initial
  cblas_dcopy(*la, RHS, 1, r, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1);
  
  norm_r = cblas_dnrm2(*la, r, 1);
  *nbite = 0;
  resvec[*nbite] = norm_r / norm_b;

  while((resvec[*nbite] > *tol) && (*nbite < *maxit))
  {
    // Résolution M*z = r
    dgbtrs_("N", la, kl, &ku_precond, &nrhs, MB, lab, ipiv, r, la, &info);
    
    // x = x + z
    cblas_daxpy(*la, 1.0, r, 1, X, 1);
    
    // r = b - Ax
    cblas_dcopy(*la, RHS, 1, r, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, r, 1);
    
    norm_r = cblas_dnrm2(*la, r, 1);
    (*nbite)++;
    resvec[*nbite] = norm_r / norm_b;
  }
  
  free(r);
  free(ipiv);
}