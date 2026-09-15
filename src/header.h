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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <pari/pari.h>
#ifdef OPENMP
#include <omp.h>
#endif
#ifdef SQLITE
#include <sqlite3.h>
#endif

struct listvec {
  double **list; /* list of vectors */
  long max; /* L[0..max-1] already allocated */
  long ind; /* < max; L[0..ind-1] filled */
};


typedef struct {
  int dim;
  int r1;
  int r2;
  double **sigma;
  double **sigma_inv;
  double **unit;
  int ordre_1;
  double *gene_1;
  int norme_ideal;
  int h;
} numberfield;


/* List of Integers (of the number fields) */
#ifdef OPENMP
typedef double *** ListOfIntegers;
#else
typedef double ** ListOfIntegers;
#endif /*OPENMP*/



/* Stack of integers */
struct pile_ent {
  int *list; /* éléments de la pile */
  long max; /* longueur allouée */
  long ind; /* < max, hauteur de pile  */
};

struct listvec liste_ent;



int LANGUAGE,NIV_AFF,PRINCIPAL;
long PARI_STACK_SIZE;
double MINIMAL_VALUE_K,EPS4;
GEN calcul_min_liste_cycles(GEN cyc, GEN unite, GEN k, GEN mini, long petite_n,long prec);
GEN calculs_corps_de_nombres(GEN P, long petite_n, long prec);
GEN minimum_quad_im(GEN Q);
#ifdef SQLITE
void creation_table_numberfield(sqlite3 *db);
void ajout_numberfield(sqlite3 *db, int n, int r1, int r2, char *disc, char *disc2, int clnb, char *pol, char *M1, int T1, char * C1, double approx, long euclidean);
void ajout_numberfield2(sqlite3 *db, int n, int r1, int r2, char *disc, char *disc2, int clnb, char *pol, long euclidean);
#endif /*SQLITE*/
int liminf(int dim);
int limsup(int dim);
int nb_de_dec(numberfield *c);
void init_numberfield(numberfield *k, int n, int r, int s, double **a, double **b, double **u, int o, double *g_1, int no, int h_K);
void free_numberfield(numberfield *k);
int partie_entiere(double x, double eps);
double *allocVec(const long dim);
void freeVec(double *v);
int *allocVecEnt(const long dim);
void freeVecEnt(int *v);
double **allocMat(const long J, const long I);
void freeMat(double **m, int dim);
void freeMatEnt(int **m, int dim);
void reallocMat(double **sigma, const long J_nouv, const long J_anc, const long I);
int **allocMatEnt(const long J, const long I);
double fmin(double a, double b);
double fmax(double a, double b);
double norme(double *v,int r1, int r2);
double module(double *v, int i, int r1, int r2);
double module_diff(double *v1, double *v2, int i, int r1, int r2);
double norml2(double *v, int dim);
void affiche_matrix(double **M, int dim1, int dim2);
void affiche_matrixEnt(int **M, int dim1, int dim2);
void affiche_vecteur(double *v, int dim);
void affiche_vecteurEnt(int *v, int dim);
int arrondi(double x);
double partie_fractionnaire(double x);
void produit_tordu(double *v1, double *v2,double *prod, int R1, int R2);
void MatVecEnt(double *v, double **ma, int *w, int m, int p);
void MatVec(double *v, double **ma, double *w, int m, int p);
double **copie_matrice(double **po,int dim,int pb);
void copie_matrice_sans_alloc(double **po, double **copie, int dim,int pb);
double *copie_vecteur(double *v,int dim);
void copie_vecteur_sans_alloc(double *v, double *res, int dim);
void initlist(struct listvec *L);
void copie_liste(struct listvec *L, struct listvec *M, int dim);
void append(struct listvec *L, double *v);
void append_2(struct listvec *L, double *v,double *w);
void optimiser(struct listvec *L);
void freelist(struct listvec *L);
void merge(struct listvec *L, struct listvec *M);
void freeListOfIntegers(ListOfIntegers e, int nb_ent,int n);
double carre(double x);
double norme_bricolee(double *v, double *probleme, double *pas,int r1, int r2);
void initpile(struct pile_ent *L);
void push(struct pile_ent *L, int e);
void freepile(struct pile_ent *L);
void pop(struct pile_ent *L, int *e);
void ramene_en_tete3(double **matrice, int i, int j, int dim);
int min(int a, int b);
#ifdef SQLITE
GEN calcul_pari2(GEN P, long prec, numberfield* k , int c_i, double *K2, double *K4,sqlite3 *db);
GEN lancement_pari(int argc, char **argv, numberfield *k, int c_i, double *K2, double *K4,sqlite3 *db);
void recherche_pol(sqlite3 *db, char *pol);
#else
GEN calcul_pari2(GEN P, long prec, numberfield* k , int c_i, double *K2, double *K4);
GEN lancement_pari(int argc, char **argv, numberfield *k, int c_i, double *K2, double *K4);
#endif
int *cfc(int nb_pb, int *nb_imag, int **pb_imag, int *k);
void affichage_graphe(int nb_sommets, int *nb_aretes, int **liste_but_aretes, int ***liste_etiqu_aretes, int dim);
int premier_indice_faux(bool *tab , int lg);
int unique_image_cfc(int i,int *pb_imag, int nb_imag, int *comp,int k, int *ind);

