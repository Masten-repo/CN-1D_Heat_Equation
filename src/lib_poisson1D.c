/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <math.h> // Nécessaire pour sqrt dans l'erreur relative manuelle

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int k = 0;
  int m = *la; // colonnes (taille du problème)
  int n = *lab; // lignes (leading dimension)
  int v = *kv; // index diagonale (ku)

  // On parcours les colonnes de la matrice 
  for(int i = 0; i < m; ++i){
    k = i * n; // Début de la colonne i

    // Mise à zéro des éléments au dessus de la première sur-diagonale (fill-in)
    if(v >= 0){
      for(int j = 0; j < v; ++j){
        AB[k+j] = 0.0;
      }
    }

    // Remplissage Poisson 1D (-1, 2, -1)
    // AB[k+v] correspond à la ligne 'ku' (sur-diagonale)
    AB[k+v] = -1.0; 
    // AB[k+v+1] correspond à la ligne 'ku+1' (diagonale)
    AB[k+v+1] = 2.0; 
    // AB[k+v+2] correspond à la ligne 'ku+2' (sous-diagonale)
    AB[k+v+2] = -1.0; 
  }

  // Correction des bords (Conditions aux limites)
  AB[0] = 0.0; // Pas de voisin gauche pour la première colonne
  if(v == 1){
    AB[1] = 0.0; // Sécurité pour la première sur-diagonale
  }
  AB[n*m - 1] = 0.0; // Pas de voisin droit pour la dernière colonne (sous-diagonale)
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int k = 0;
  int m = *la;
  int n = *lab;
  int v = *kv;

  for(int i = 0; i < m; ++i){
    k = i * n;
    if(v >= 0){
      for(int j = 0; j < v; j++){
        AB[k+j] = 0.0;
      }
    }
    AB[k+v] = 0.0;     // Sur-diagonale
    AB[k+v+1] = 1.0;   // Diagonale (Identité)
    AB[k+v+2] = 0.0;   // Sous-diagonale
  }
  // Nettoyage bords
  AB[1] = 0.0; 
  AB[n*m - 1] = 0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int m = *la;
  // Initialisation à 0 (Source g=0)
  for (int i = 0; i < m; ++i){
    RHS[i] = 0.0;
  }
  // Application des conditions aux limites Dirichlet
  RHS[0] = *BC0; 
  RHS[m-1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int m = *la;
  double d = (*BC1) - (*BC0); // Delta T
  
  for(int i = 0; i < m; ++i){
    EX_SOL[i] = *BC0 + X[i] * d; // T(x) = T0 + x(T1-T0)
  }
}  

void set_grid_points_1D(double* x, int* la){
  if (x == NULL || la == NULL || *la <= 0) return;
  
  int m = *la;
  double h = 1.0 / (double)(m + 1); 
  
  for(int i = 0; i < m; ++i){
    x[i] = (i + 1) * h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  
  double norm_diff = 0.0;
  double norm_x = 0.0;
  double diff;
  int m = *la;

  for(int i=0; i<m; i++){
      diff = x[i] - y[i];
      norm_diff += diff * diff;
      norm_x += x[i] * x[i];
  }

  norm_diff = sqrt(norm_diff);
  norm_x = sqrt(norm_x);

  if (norm_x == 0.0) return 0.0;
  return norm_diff / norm_x;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab) + i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  *info = 0;
  
  int m = *la;
  int ld = *lab;
  
  // 1. Initialisation de ipiv (Obligatoire pour dgbtrs)
  // dgbtrs suppose que des pivots ont été choisis. Si on ne pivote pas, 
  // on remplit ipiv avec l'identité (1, 2, 3...).
  for (int i = 0; i < m; i++) {
      ipiv[i] = i + 1; // Indices 1-based pour LAPACK
  }

  // 2. Factorisation LU sans pivotage (Algorithme de Thomas simplifié)
  // On parcourt les colonnes de 0 à m-2 (l'avant dernière)
  for (int i = 0; i < m - 1; i++) {
      // Indices dans le stockage bande (kv=1) :
      // Diagonale actuelle : colonne i, ligne 2
      double pivot = AB[i * ld + 2]; 

      if (pivot == 0.0) {
          *info = i + 1; // Erreur pivot nul
          return *info;
      }

      // Sous-diagonale actuelle (à éliminer) : colonne i, ligne 3
      double sub_diag = AB[i * ld + 3];

      // Calcul du facteur L
      double factor = sub_diag / pivot;
      
      // Stockage du facteur L à la place de la sous-diagonale
      AB[i * ld + 3] = factor;

      // Mise à jour de la diagonale suivante (colonne i+1)
      // Diagonale suivante = Diagonale suivante - (facteur * sur-diagonale suivante)
      // Sur-diagonale suivante : colonne i+1, ligne 1
      double super_diag_next = AB[(i + 1) * ld + 1];
      
      AB[(i + 1) * ld + 2] -= factor * super_diag_next;
  }
  
  return *info;
}