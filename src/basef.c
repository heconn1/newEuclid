/*
Copyright (C) 2012  Pierre Lezowski.

This file is part of the euclid package.

euclid is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/* Base functions, which are repeatedly used */

#include "header.h"




void init_numberfield(numberfield *k, int n, int r, int s, double **a, double **b, double **u, int o, double *g_1, int no, int h_K){
  k->dim = n;
  k->r1 = r;
  k->r2 = s;
  k->sigma=a;
  k->sigma_inv = b;
  k->unit = u;
  k->ordre_1 = o;
  k->gene_1 = g_1;
  k->norme_ideal = no;
  k->h = h_K;
}


void free_numberfield(numberfield *k){
  freeMat(k->sigma,k->dim);
  freeMat(k->sigma_inv,k->dim);
  freeMat(k->unit,k->r1+k->r2-1);
  free (k->gene_1);
  return;
}


inline int partie_entiere(double x, double eps){
	int a = ceil(x);
	if ( (a -x) < eps){
		return a;
	}
	else{
		return a-1;
	}
}




double *
allocVec(const long dim) { return (double *)malloc(dim * sizeof(double)); }

void freeVec(double *v){ free (v); v= NULL;}

int *
allocVecEnt(const long dim) { return (int *)malloc(dim * sizeof(int)); }

void freeVecEnt(int *v){ free (v) ; v = NULL ;}

double **
allocMat(const long J, const long I)
{
  double **sigma = (double **) malloc(J * sizeof(double *));
  long j; for (j = 0; j < J; j++) sigma[j] = allocVec(I);
  return sigma;
}



void freeMat(double **m, int dim){
  int i;
  for (i=0 ; i< dim ; i++) freeVec(m[i]);
  free (m);
  m = NULL;
}

void freeMatEnt(int **m, int dim){
  int i;
  for (i=0 ; i< dim ; i++) freeVecEnt(m[i]);
  free (m);
  m = NULL;
}



void
reallocMat(double **sigma, const long J_nouv, const long J_anc, const long I){
   long j;
   for (j=0 ; j < J_anc ; j++) freeVec(sigma[j]);
   sigma = allocMat(J_nouv, I);
}



int **
allocMatEnt(const long J, const long I)
{
  int **sigma = (int **) malloc(J * sizeof(int *));
  long j; for (j = 0; j < J; j++) sigma[j] = allocVecEnt(I);
  return sigma;
}


double fmin(double a, double b){
  if (a<b) {return a;} else {return b;}
}

double fmax(double a, double b){
  if (a<b) {return b;} else {return a;}
}

inline double norme(double *v,int r1, int r2)
{
  int i,j ;
  double norm = 1. ;
  for (i=0 ; i < r1; i++) norm *= v[i] ;
  for (j=r1 ; j < r1+r2; j++) norm *= (v[j] * v[j] + v[j+r2] * v[j+r2]) ;
  return norm;
}

inline double module(double *v, int i, int r1, int r2){
  if (i < r1){
     return fabs(v[i]);
  }
  else{
    if (i < r1+r2){
       return sqrt(v[i]*v[i] + v[i+r2] *v[i+r2]);
    }
    else{
       return sqrt(v[i-r2]*v[i-r2] + v[i] * v[i]);
    }
  }
}

inline double module_diff(double *v1, double *v2, int i, int r1, int r2){ 
  if (i < r1){
     return fabs(v1[i] - v2[i]);
  }
  else{
    if(i<r1+r2){
       return sqrt(   (v1[i] - v2[i]) * (v1[i] - v2[i]) + (v1[i+r2] - v2[i+r2]) * (v1[i+r2] - v2[i+r2]));
    }
    else{
      return sqrt(   (v1[i] - v2[i]) * (v1[i] - v2[i]) + (v1[i-r2] - v2[i-r2]) * (v1[i-r2] - v2[i-r2]));
    }
  }
}



double norml2(double *v, int dim)
{
  double t = 0;
  int i;
  for (i = 0; i < dim; i++) t += v[i] * v[i] ;
  return t;
}

/* affichage d'une matrice*/
void affiche_matrix(double **M, int dim1, int dim2){  
  int i, j; 
  if ((dim1 ==0)||(dim2==0)) {
    if (LANGUAGE){
      fprintf(stdout,"de rang 0\n"); 
    }
    else{
      fprintf(stdout,"of rank 0\n"); 
    }
    return;
  }
  fprintf(stdout,"["); 
  for (i = 0; i < dim1 - 1; i++) 
    {
      for (j = 0; j < dim2 - 1; j++) 
	{ 
	  fprintf(stdout,"%f, ", M[i][j]); 
	}
      fprintf(stdout,"%f;\n ", M[i][dim2 - 1]); 
    }

  for (j = 0; j < dim2 - 1; j++) 
    fprintf(stdout,"%f, ", M[dim1-1][j]); 
  fprintf(stdout,"%f]\n", M[dim1-1][dim2-1]);
  fflush(stdout);
}

void affiche_matrixEnt(int **M, int dim1, int dim2){  
  int i, j; 
  if ((dim1 ==0)||(dim2==0)) {printf("de rang 0\n"); return;}
  fprintf(stdout,"["); 
  for (i = 0; i < dim1 - 1; i++) {
    for (j = 0; j < dim2 - 1; j++) { 
      fprintf(stdout,"%d, ", M[i][j]); 
    }
    fprintf(stdout,"%d;\n ", M[i][dim2 - 1]); 
  }
  for (j = 0; j < dim2 - 1; j++) 
    fprintf(stdout,"%d, ", M[dim1-1][j]); 
  fprintf(stdout,"%d]\n", M[dim1-1][dim2-1]);
  fflush(stdout);
}

/* affichage d'un vecteur*/
void affiche_vecteur(double *v, int dim){
  int i;
  fprintf(stdout,"[");
  for (i=0; i< dim-1 ; i++) 
    fprintf(stdout,"%f; ",v[i]);
  fprintf(stdout,"%f]\n",v[dim-1]);
  fflush(stdout);
}

void affiche_vecteurEnt(int *v, int dim){
  int i;
  fprintf(stdout,"[");
  for (i=0; i< dim-1 ; i++) 
    fprintf(stdout,"%d; ",v[i]);
  fprintf(stdout,"%d]\n",v[dim-1]);
  fflush(stdout);
}

/* rounds a double to an integer, provided it is close enough to an integer */ 
int arrondi(double x){
  int num; 
  double diff;
  num = floor(x+0.5);
  diff = fabs(x - num);
  /*f (diff > 0.5) { num++ ; diff = 1 - diff;}*/
  if (diff> EPS4) { 
    if (LANGUAGE){
      fprintf(stderr, "Rounding problem - round : x= %f \n", x);
    }
    else{
      fprintf(stderr, "Problème d'arrondi - arrondi : x= %f \n", x);
    }
  }
  return num;
}


/* when x = y/n^DIM, writes x = z/t */

double partie_fractionnaire(double x){
  return x-floor(x);
}

inline void produit_tordu(double *v1, double *v2,double *prod, int R1, int R2)
{
  int i,j ;
  for (i=0 ; i < R1; i++) prod[i]= v1[i] * v2[i] ;
  for (j=R1 ; j < R1+R2; j++) 
    {
     prod[j] = (v1[j] * v2[j] - v1[j+R2] * v2[j+R2]) ;
     prod[j+R2] = ( v1[j] * v2[j+R2] + v1[j+R2] * v2[j]) ;
    }
}


/* multiplication of an integral vector by a matrix */
void MatVecEnt(double *v, double **ma, int *w, int m, int p){
  int i, k;
  double x; 
  for (i = 0; i < m; i++) 
    {
      x = 0; 
      for (k = 0; k < p; k++) 
	  x += ma[i][k] * w[k]; 
      v[i] = x; 
    }
}


/* multiplication of a float vector by a matrix */
void MatVec(double *v, double **ma, double *w, int m, int p) {
  int i, k;
  double x; 
  for (i = 0; i < m; i++){
    x = 0; 
    for (k = 0; k < p; k++) 
      x += ma[i][k] * w[k]; 
    v[i] = x; 
  }
}


/* copy of a matrix */
double **copie_matrice(double **po,int dim,int pb){
  double **res;  
  int i,j;
  res = allocMat(dim,pb);
  for (i=0; i< dim ; i++){
    for (j=0 ; j< pb ;j++){
     res[i][j] = po[i][j];
    }
  }
  return res;
}


/* copy of a matrix without allocating memory */
void copie_matrice_sans_alloc(double **po, double **copie, int dim,int pb){
  int i,j;
  for (i=0; i< dim ; i++){
    for (j=0 ; j< pb ;j++){
     copie[i][j] = po[i][j];
    }
  }
  return;
}

/* copy of a vector */
double *copie_vecteur(double *v,int dim){
  double *res; 
  int i;
  res = allocVec(dim);
  for (i=0; i< dim ; i++){
     res[i]= v[i];
  }
  return res;
}


/* copy of a vector without allocation */
void copie_vecteur_sans_alloc(double *v, double *res, int dim){ 
  int i;
  for (i=0; i< dim ; i++){
     res[i]= v[i];
  }
}



/* list of integers */

void initlist(struct listvec *L){
  L->ind = 0;
  L->max = 512;
  L->list = (double **) malloc(L->max * sizeof(double *));
}

void copie_liste(struct listvec *L, struct listvec *M, int dim){
  int a;
  L->ind = M->ind;
  L->max = M->ind;
  L->list=malloc(M->ind *sizeof(double*));
  for(a=0 ; a < L->ind ; a++){
    L->list[a] = malloc(dim*sizeof(double));
    memcpy(L->list[a],M->list[a],dim*sizeof(double));
  }
}

void append(struct listvec *L, double *v){
  if (L->ind >= L->max-1)
  {
    L->max *= 2;
    L->list = (double**)realloc(L->list, L->max * sizeof(double*));
  }
  L->list[ L->ind++ ] = v;
}


void append_2(struct listvec *L, double *v,double *w){
  {
    append(L,v);
    append(L,w);
  }
  return;
}



void optimiser(struct listvec *L){
  if(L->max > L->ind){
    L->max = L->ind;
    L->list = (double**)realloc(L->list, L->max*sizeof(double*)); 
  }
}

#ifdef OPENMP
#else
void freelist(struct listvec *L)
{
  long i;
  for (i = 0; i < L->ind; i++) free((void*)L->list[i]);
  free((void*)L->list);
}
#endif

/* colle M au bout de L */
#ifdef OPENMP
void merge(struct listvec *L, struct listvec *M)
{ 
#pragma omp critical
  {
    long i;
    L->max = L->ind + M->ind +2;
    L->list = (double**)realloc(L->list,L->max*sizeof(double*));
    for(i=0; i < M->ind; i++){
      L->list[L->ind+i] = M->list[i];
    }
    L->ind += M->ind;
  }
}
#endif /*OPENMP*/




void freeListOfIntegers(ListOfIntegers e, int nb_ent,int n){
#ifdef OPENMP
  int i;
  for(i=0 ; i < n ; i++){
    freeMat(e[i],nb_ent);
  }
  free (e);
#else
  freeMat(e,nb_ent);
#endif /*OPENMP*/
}



/* end list of integers */

inline double carre(double x){
  return x*x;
}


/* a custom norm for the absorption test */
inline double norme_bricolee(double *v, double *probleme, double *pas,int r1, int r2)
{ 
  int i ;
  double norm = 1. ;
  for (i=0 ; i < r1; i++) norm *= fabs(v[i] - probleme[i])+pas[i] ;
  for (i=r1 ; i < r1+r2; i++) norm *= carre(fabs(v[i] - probleme[i])+pas[i]) + carre(fabs(v[i+r2] - probleme[i+r2])+pas[i+r2]) ;
  return norm;
}



/* stack for Tarjan's algorithm */

void initpile(struct pile_ent *L)
{
  L->ind = 0;
  L->max = 2;
  L->list = (int *) malloc(L->max * sizeof(int));
}

void push(struct pile_ent *L, int e)
{
  if (L->ind >= L->max-1)
  {
    L->max *= 2;
    L->list = (int *)realloc(L->list, L->max * sizeof(int));
  }
  L->list[ L->ind++ ] = e;
}

void freepile(struct pile_ent *L)
{
  free((void*)L->list);
}

void pop(struct pile_ent *L, int *e){
  if(L->ind >0){
    *e = L->list[-- L->ind];
  }
  else{
    fprintf(stderr,"erreur : pop appliqué à pile vide");
  }
}

/* end stack for Tarjan's algorithm */




inline void echange(double **matrice, int i, int j, int dim){
  double t;
  int a;
  for (a=0 ; a< dim ; a++) {
    t = matrice[i][a];
    matrice[i][a] = matrice[j][a];
    matrice[j][a] = t;
  }
  return;
}


inline void ramene_en_tete3(double **matrice, int i, int j, int dim){
  int a;
  double *t;
#ifdef OPENMP
#pragma omp critical
#endif /*OPENMP*/
  {
	if (i<j){
	  t = matrice[j];  /* the coefficient we will bring back to head */
		for (a = j ; a > i ; a--){
			matrice[a] = matrice[a-1];
		}
		matrice[i] = t;  
	}
  }
  return;
}


int min(int a, int b){
  if(a<b) return a;
    else return b;
}
