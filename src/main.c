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


#include "header.h"

#define DISPLAY_INF_MIN "E" /* what to display if we find less than the minimal value */
#define DISPLAY_NON_PRINCIPAL "H" /* what to display if the number field is not principal */
#define DATABASE_NAME "db/table_euc.db"
#define CONFIG_FILE "euclid.cfg"


/* Global variables read in config file, LANGUAGE, NIV_AFF,
   MINIMAL_VALUE_K,PRINCIPAL, EPS4 are in header*/
double INITIAL_VALUE_K;
int ITERATIONS_WITHOUT_UNITS;
double EPS0,EPS2,EPS3;
int MAXITERATIONS,MAX_NUMBER_PB,MAX_NUMBER_PB2,MAX_PB_CONV,NB_UNITS;
int BOUNDSM[15],BOUNDSP[15],CUTTING[15];


void read_config(FILE *file){
  char line[256],val[256];
  char **v2 = malloc(15*sizeof(char *));
  int i;
  for(i=0 ; i< 15 ; i++)
    v2[i] = malloc(256*sizeof(char));
  fgets(line, 256,file);
  sscanf(line, "%*s %s",  val);
  INITIAL_VALUE_K = atof(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  MINIMAL_VALUE_K = atof(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  PRINCIPAL = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  ITERATIONS_WITHOUT_UNITS = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  NIV_AFF = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  EPS0 = atof(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  EPS2 = atof(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  EPS3 = atof(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  EPS4 = atof(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  MAXITERATIONS = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  MAX_NUMBER_PB = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  MAX_NUMBER_PB2 = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  MAX_PB_CONV = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  NB_UNITS = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  LANGUAGE = atoi(val);
  fgets(line, 256,file);
  sscanf(line, "%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6],v2[7],v2[8],v2[9],v2[10],v2[11],v2[12],v2[13],v2[14]);
  for(i=0 ; i < 14 ; i++)
    BOUNDSM[i] = atoi(v2[i]);
  fgets(line, 256,file);
  sscanf(line, "%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6],v2[7],v2[8],v2[9],v2[10],v2[11],v2[12],v2[13],v2[14]);
  for(i=0 ; i < 14 ; i++)
    BOUNDSP[i] = atoi(v2[i]);
  fgets(line, 256,file);
  sscanf(line, "%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", v2[0],v2[1],v2[2],v2[3],v2[4],v2[5],v2[6],v2[7],v2[8],v2[9],v2[10],v2[11],v2[12],v2[13],v2[14]);
  for(i=0 ; i < 14 ; i++)
    CUTTING[i] = atoi(v2[i]);
  fgets(line, 256,file);
  sscanf(line, "%*s %s", val);
  PARI_STACK_SIZE = atoi(val);
}

int liminf(int dim){
  if(dim < 14){
    return(BOUNDSM[dim-1]);
  }
  else{
    return(BOUNDSM[14]);
  }   
}

int limsup(int dim){
  if(dim < 14){
    return(BOUNDSP[dim-1]);
  }
  else{
    return(BOUNDSP[14]);
  }   
}

int nb_de_dec(numberfield *c){
  if(c->dim < 14){
    return(CUTTING[c->dim-1]);
  }
  else{
    return(CUTTING[14]);
  }   
}

/* Integers (of the number field) used*/

ListOfIntegers small_elts(numberfield *corps, double K, int liminf,int limsup, int n,long int *nb_po){
  long i, k, ind;
  double *max,*min,x,c;
#ifdef OPENMP
  int aaa;
  struct listvec *r = malloc(n*sizeof(struct listvec));
  struct listvec ** q = malloc(n*sizeof(struct listvec *));
  for(aaa=0 ; aaa < n ; aaa++){
    q[aaa] = & (r[aaa]);
    initlist(q[aaa]);
  }
#else
  int *vars;
  double *v;
  struct listvec e;
  struct listvec *po;
  po= &e;
  initlist(po);
  v=allocVec(corps->dim);
  vars =allocVecEnt(corps->dim);
  for (i = 0; i < corps->dim; i++) vars[i] = liminf;
#endif /*OPENMP*/
  max = malloc((corps->dim)*sizeof(double)); 
  min = malloc((corps->dim)*sizeof(double));
  x = exp(log(K)/corps->dim) + EPS0;
  for (i=0;i<corps->dim;i++)
    {
      max[i] = x; min[i] = -x;
      for (k = 0; k < corps->dim; k++)
	{
	  c = corps->sigma[i][k];
	  if (c > 0) max[i] += c; else min[i] += c;
	}
    }
#ifdef OPENMP
  { /* selection part */
    int bbb,th_id, **vars2;
    double **v2;
    vars2=allocMatEnt(n,corps->dim);
    for(th_id=0 ; th_id < n ; th_id ++){
      for(bbb=0 ; bbb < corps->dim ; bbb++){
	vars2[th_id][bbb] = liminf;
      }
    }
    v2 = allocMat(n,corps->dim);
#pragma omp parallel for private(i,aaa,th_id,ind)
    for(bbb= liminf ; bbb <= limsup ; bbb++){
      th_id = omp_get_thread_num();
      /*fprintf(stdout,"Entier %d (thread %d).\n",bbb,th_id);*/
      vars2[th_id][0] = bbb;
      /*    fprintf(stdout,"Entier %d\n",vars2[0]);*/
      for (;;) {
	MatVecEnt(v2[th_id],corps->sigma,vars2[th_id],corps->dim,corps->dim);
	for (i = 0; i < corps->dim; i++)
	  if (v2[th_id][i] < max[i] && v2[th_id][i] > min[i]) break;
	if(i < corps->dim){
	  append(q[th_id],v2[th_id]);
	  v2[th_id]=malloc(corps->dim*sizeof(double));
	}
	ind = 1;
	while(ind < corps->dim && vars2[th_id][ind] == limsup) vars2[th_id][ind++] = liminf;
	if (ind >= corps->dim) break;
	vars2[th_id][ind]++;
      }
    }
    freeMatEnt(vars2,n);
    free (v2);
    optimiser(q[0]);
    for(aaa=1 ; aaa <n ; aaa++){
      optimiser(q[aaa]);
      merge(q[0],q[aaa]);
      free (q[aaa]->list);
    }
    optimiser(q[0]);
    if (NIV_AFF >0){
      if (LANGUAGE){
	fprintf(stdout,"We will use %ld integers, copy in progress.\n",q[0]->ind);
      }
      else{
	fprintf(stdout,"On va utiliser %ld entiers, copie en cours.\n",q[0]->ind);
      }    
    }
    /*#pragma omp parallel for*/
    for(aaa=1 ; aaa <n ; aaa++){
      copie_liste(q[aaa],q[0],corps->dim);
      /*fprintf(stdout,"Copie %d finie.\n",aaa);*/
    }
  }
#else /*OPENMP*/
  for (;;) {
    /* FIXME: most of the time, only 1 coord of var changes and adding a column
     * would be enough */
    MatVecEnt(v,corps->sigma,vars,corps->dim,corps->dim);
    for (i = 0; i < corps->dim; i++)
      if (v[i] < max[i] && v[i] > min[i]) break;
    if(i < corps->dim){
      append(po,v); v=allocVec(corps->dim);
    }
    ind = 0;
    while(ind < corps->dim && vars[ind] == limsup) vars[ind++] = liminf;
    if (ind >= corps->dim) break;
    vars[ind]++;
  }
#endif
#ifdef OPENMP
  { /* optimization and merge */
    double *** p;
    p = malloc(n*sizeof(double**));
    for (aaa = 0 ; aaa < n ; aaa++){
      optimiser(q[aaa]);
      /*fprintf(stdout,"%ld entiers.\n",q[aaa]->ind);*/
      p[aaa] = q[aaa]->list;
    }
    *nb_po = q[0]->ind;
    return (p);
  }
#else /*OPENMP*/
  optimiser(po);
  *nb_po = po->ind;
  freeVec(max);
  freeVec(min);
  return (po->list);
#endif
}




/* Cutting of the fundamental domain and covering by parallelotopes */

void bornes_pour_parallelotopes_utiles(int J, double **oo,double **y, const int *v, int *l1, int *l2, int dim, int *bornes, const double pas_J){
  long p, j;
  double m,M;
  double ma=bornes[J];
  double mi=0;
  double o;
  double ssmin; 
  double ssmax;
  *l1 = 0;
  *l2 = bornes[J]; 
  for (p = 0; p < dim; p++){
    double Min, Max;
    ssmin = 0; 
    ssmax = 0;
    for (j = 0; j < J ; j++)
      {
	o = oo[p][j];
	if (o > 0)
	  {
	    Min=  o * y[j][v[j]];
	    Max=  o * y[j][v[j]+1];
	  }
	else
	  {
	    Max = o * y[j][v[j]];
	    Min =  o * y[j][v[j]+1];
	  }
	ssmin += Min;
	ssmax += Max;
      }
    for (j = J+1; j < dim ; j++)
      {
	o = oo[p][j];
	if (o > 0)
	  {
	    Min = o * y[j][0];            
	    Max = o * y[j][bornes[j]+1];
	  }
	else
	  {
	    Max = o * y[j][0];
	    Min = o * y[j][bornes[j]+1];
	  }
	ssmin += Min;
	ssmax += Max;
      }
    o = oo[p][J];
    if (o != 0){  /* dividing by 0 should probably be avoided */
      M = -1/(2*pas_J) * (y[J][0]+ssmax/o);
      m = -1/(2*pas_J) * (y[J][0]+(ssmin-1)/o);
      if (o > 0) { M -= 1; if (M >mi) mi = M; if (m < ma) ma = m; }
      else       { m -= 1; if (M < ma) ma = M; if (m > mi) mi = m; }
      /*fprintf(stdout,"o = %f, p= %ld, m=%f, M=%f, mi = %f, ma = %f, pas = %f   ", o,p,m,M,mi, ma, pas_J);*/
    }
  }
  if (mi - floor(mi) < EPS3) mi = floor(mi); else mi = ceil(mi);
  if (ceil(ma) - ma < EPS3) ma = ceil(ma); else ma = floor(ma);
  *l1 = (mi >0  ? mi : 0) ;
  *l2 = (ma < bornes[J] ? ma : bornes[J]);
  /*fprintf(stdout,"l1 = %d, l2 = %d\n", *l1, *l2);*/
  return;
}



/* includes the fundamental domain into a parallelogram of R^dim and computes subdivision */

double **calcul_y(double **sigma, const int dim, int* nob, double **pas){
  int i,k;
  double fmax1,fmin1, **y;
  *pas = allocVec(dim);
  y = malloc(dim * sizeof(double *));
  for (i=0 ; i < dim; i++) y[i] = malloc((nob[i]+2) *sizeof(double));
  for (i = 0; i < dim; i++){
    fmax1 = 0.; fmin1 = 0.; 
    for (k = 0; k < dim; k++){
      if (sigma[i][k] > 0) fmax1 += sigma[i][k]; 
      else fmin1 += sigma[i][k]; 
    }
    (*pas)[i] = (fmax1 - fmin1) / 2 / (nob[i]+1);
    if (NIV_AFF >0){
      if (LANGUAGE){
	fprintf(stdout, "[%f, %f], step = %f\n", fmin1, fmax1, (*pas)[i]);
      }
      else{
	fprintf(stdout, "[%f, %f], pas = %f\n", fmin1, fmax1, (*pas)[i]);
      } 
    }
    for (k = 0; k <= nob[i]+1; k++){
      y[i][k] = fmin1 + 2 * k * (*pas)[i]; 
    }
  }
  return y;
}


/* Exploration bounds in any direction and first test*/



bool test_absorption(double *probleme, double ** entiers, int nb_entiers, int dim, int r1, int r2, int *entier, double k, double *pas){
  int j=0;
  bool b = false; /* dit si on a réussi à absorber */
  while ((j < nb_entiers)&& (!b)){
    b = (norme_bricolee(entiers[j],probleme,pas,r1,r2) < k-EPS2);
    j++;
  }
  *entier = j-1;
  return b;
}


/* treatment of the problem "vec"*/
void traitement_vec(int *vec, numberfield *corps,double **y, double *pas, double **entiers, const int nb_entiers, const double K,struct listvec *problemes){  
  double *pb ; 
  int i,entier;
  pb= allocVec(corps->dim);
  for (i=0; i < corps->dim; i++) 
    pb[i] = y[i][vec[i]] + pas[i];
  entier=0;
  if (test_absorption(pb, entiers,nb_entiers, corps->dim, corps->r1, corps->r2, &entier, K, pas)){
    ramene_en_tete3(entiers,0,entier,corps->dim); /*place l'entier qui a absorbé en tête.*/
    freeVec(pb);
  }
  else{
    double *pb2 = allocVec(corps->dim);
    int l1,l2;
    double t;
    for (l1=0 ; l1< corps->dim ; l1++){
      t = 0;
      for (l2 = 0 ; l2 < corps->dim ; l2++){
	t += corps->sigma[l1][l2];
      }
      pb2[l1] = t - pb[l1];
    }
    append_2(problemes,pb,pb2);
  }
}

/* Recursive treatment */
void parcours_coordonnees(int j, int *w,numberfield *corps, double **y, double *pas, double **entiers, const int nb_entiers, const double K,struct listvec *problemes, int *nob, FILE* fich){  
  /* fonction de parcours des coordonnées, celles avant j sont fixées*/
  if (j>= corps->dim) {
    traitement_vec(w,corps, y, pas, entiers,nb_entiers, K,problemes); 
  }
  else{
    int b1,b2;
    bornes_pour_parallelotopes_utiles(j, corps->sigma_inv,y, w, &b1, &b2, corps->dim, nob, pas[j]);
    for (w[j] = b1 ; w[j] <= b2 ; w[j] = w[j]+1){
      /*		
      if (j==1){
	fprintf(stdout,"%d,%d: %ld problèmes.\n", w[0],w[1],problemes->ind); fflush(stdout);
	}
      */
      /*
	if (j==2){
	fprintf(stdout,"%d,%d,%d: %ld problèmes.\n", w[0],w[1],w[2],problemes->ind); fflush(stdout);
	}
      */
      parcours_coordonnees(j+1,w,corps,y,pas,entiers,nb_entiers,K,problemes,nob,fich);
    }      
  }
  return;
}


void decoupage_initial5(double **y, numberfield* corps, int* nob, ListOfIntegers entiers, int nb_entiers, struct listvec * problemes, int *nombre, double *pas, double K, FILE* fich){
  /*fprintf(stdout,"K=%lf\n",K);*/
  /*fprintf(fich,"Test pour K=%lf\n",K);*/
  int aa;
#ifdef OPENMP
  int th_id;
#pragma omp parallel for private(th_id)
  for (aa=0 ; aa< (nob[0]+1)/2 ; aa ++){
    struct listvec problemes2;
    struct listvec *pb2;
    int *z;
    pb2 = &problemes2;
    th_id = omp_get_thread_num();
    initlist(pb2);
    /*fprintf(stdout,"début de %d : %ld problèmes.\n", aa,problemes->ind); fflush(stdout);*/
    z = allocVecEnt(corps->dim);
    z[0] = aa;
    parcours_coordonnees(1,z,corps,y,pas,entiers[th_id],nb_entiers,K,pb2,nob,fich);
    optimiser(pb2);
    merge(problemes,pb2);
    /*fprintf(stdout,"fin de %d : %ld problèmes (dont %ld).\n", aa,problemes->ind,pb2->ind); fflush(fich);*/
  }
#else /*OPENMP*/
  for (aa=0 ; aa< (nob[0]+1)/2 ; aa ++){
    /*fprintf(fich,"début de %d : %ld problèmes.\n", aa,problemes->ind); fflush(fich);*/
    int *z;
    z = allocVecEnt(corps->dim);
    z[0] = aa;
    parcours_coordonnees(1,z,corps, y,pas,entiers,nb_entiers,K,problemes,nob,fich);
    /*fprintf(fich,"fin de %d : %ld problèmes.\n", aa,problemes->ind); fflush(fich);*/
  }
#endif /*OPENMP*/
  optimiser(problemes);
  return;
}





/* display the problems */
/* meant not to be used, debugging function */

void affichage_problemes(double **problemes, double *pas, double **sigma_inv, int nb_pb , int dim, int b){
  if (b==0){
    return;
  }
  else{
    if (nb_pb >0){
      int affich = b;
      if (b != 1){
        fprintf(stdout, "Voulez-vous afficher le(s) problème(s) ? ");
        fscanf(stdin,"%d", &affich);
      }
      if (affich==1){
        double * probleme=allocVec(dim);
        int i,j;
        for (j=0 ; j < nb_pb ; j++){
          MatVec(probleme,sigma_inv,problemes[j],dim,dim);
          fprintf(stdout,"problème n°%d : [ ", j); 
          for (i=0 ; i < dim ; i++){
            fprintf(stdout,"%f", probleme[i]);
            if (i < dim -1) fprintf(stdout," , ");
          }
          fprintf(stdout,"]\n");
        }
        freeVec(probleme);
	if (NIV_AFF >0){
	  fprintf(stdout," ---------------- \n");
	  for (j=0 ; j < nb_pb ; j++){
	    fprintf(stdout,"problème n°%d : [ ", j); 
	    for (i=0 ; i < dim ; i++){
	      fprintf(stdout,"%f", problemes[j][i]);
	      if (i < dim -1) fprintf(stdout," , ");
	    }
	    fprintf(stdout,"]\n");
	  }
	}
	if (NIV_AFF >0){
	  fprintf(stdout,"---------------- \n");
	  for (i=0 ; i< dim ; i++) fprintf(stdout, "pas[%d] = %f\n", i, pas[i]);
	}
      }
    }
  }
}


/* LOOP : action of the units and further cutting */

/* cutting + absorption */

void decoupage_problemes_abs(double ***liste_pb, int *nb_pb, double *pas,int dim, int r1, int r2, ListOfIntegers entiers, int nb_entiers, double K){ 
  /* un découpage bête avec absorption directe */
  int i,j,k,ij;
  int nb_pb2;
  nb_pb2= pow(2,dim) * (*nb_pb);
  for (i=0 ; i< dim ; i++) pas[i] = pas[i]/2;
#ifdef OPENMP
  {
    int n;
    struct listvec *n_pb ;
    double **liste_pb_bis =allocMat(nb_pb2,dim);
    bool *vecteur = malloc(dim * sizeof(bool));
    k=0;
    for (j=0 ; j< *nb_pb ; j=j+2){
      for (i=0 ; i< dim ; i++) vecteur[i]=false;
      ij=0;
      do{ 
	for (i=0 ; i< dim ; i++){
	  if (vecteur[i]){
	    liste_pb_bis[k][i] = (*liste_pb)[j][i] - pas[i] ;
	    liste_pb_bis[k+1][i] = (*liste_pb)[j+1][i] + pas[i] ;		   
	  } 
	  else{
	    liste_pb_bis[k][i] = (*liste_pb)[j][i] + pas[i] ;
	    liste_pb_bis[k+1][i] = (*liste_pb)[j+1][i] - pas[i] ;
	  }
	}
	k=k+2;
	if (ij < dim){
	  while((ij < dim) && vecteur[ij]){
	    vecteur[ij] = false;
	    ij++;
	  }
	  if(ij<dim){
	    vecteur[ij] =true;
	    ij = 0;
	  }
	}
      }
      while (ij< dim );
    }
    if (NIV_AFF >0){
      if (LANGUAGE){
	fprintf(stdout,"%d problems are being processed.\n",nb_pb2);
      }
      else{
	fprintf(stdout,"On traite %d problèmes.\n",nb_pb2);
      }
    }
    n= omp_get_max_threads();
    n_pb = malloc(n*sizeof(struct listvec));
    for(ij=0 ; ij < n ; ij++){
      initlist(&(n_pb[ij]));
    }
#pragma omp parallel for 
    for(k=0 ; k <nb_pb2 ; k=k+2){
      int entier;
      int i_t = omp_get_thread_num();
      if (test_absorption(liste_pb_bis[k], entiers[i_t], nb_entiers, dim, r1, r2, &entier, K ,pas)){
	ramene_en_tete3(entiers[i_t],0,entier,dim);
	freeVec(liste_pb_bis[k]);
	freeVec(liste_pb_bis[k+1]);
      }
      else{
	append(&(n_pb[i_t]),liste_pb_bis[k]);
	append(&(n_pb[i_t]),liste_pb_bis[k+1]);
      }
    }
    for(ij=1 ; ij < n ;ij++){
      optimiser(&(n_pb[ij]));
      /*fprintf(stdout,"Processeur %d, %ld problèmes.\n",ij,(&(n_pb[ij]))->ind);*/
      merge(&(n_pb[0]),&(n_pb[ij]));
    }
    optimiser(&(n_pb[0]));
    *nb_pb=(&(n_pb[0]))->ind;
    *liste_pb=(&(n_pb[0]))->list;
  }
#else /*OPENMP*/
  {
    int nb_pb_a_traiter = 0; /* compteur de vérification*/
    double **liste_pb_bis = malloc(2*sizeof(double*));
    bool *vecteur = malloc(dim * sizeof(bool));
    int entier ; /* l'indice d'un entier */
    k=0; /* l'indice du nouveau problème traité.*/
    liste_pb_bis[0] = malloc(dim *sizeof(double)); /* on alloue les deux premières cases */
    liste_pb_bis[1] = malloc(dim *sizeof(double)); /* pour les deux premiers problèmes */
    for (j=0 ; j< *nb_pb ; j=j+2){
      for (i=0 ; i< dim ; i++) vecteur[i]=false;
      ij=0;
      do{ 
	for (i=0 ; i< dim ; i++){
	  if (vecteur[i]){
	    liste_pb_bis[k][i] = (*liste_pb)[j][i] - pas[i] ; 
	  } 
	  else{
	    liste_pb_bis[k][i] = (*liste_pb)[j][i] + pas[i] ; 
	  }
	}
	if (test_absorption(liste_pb_bis[k], entiers, nb_entiers, dim, r1, r2, &entier, K ,pas)){
	  ramene_en_tete3(entiers,0,entier,dim);
	}
	else{
	  for (i=0 ; i< dim ; i++){
	    if (vecteur[i]){ 
	      liste_pb_bis[k+1][i] = (*liste_pb)[j+1][i] + pas[i] ;
	    }
	    else{
	      liste_pb_bis[k+1][i] = (*liste_pb)[j+1][i] - pas[i] ;
	    }
	  }
	  k=k+2;
	  liste_pb_bis = realloc(liste_pb_bis,(k+2)*sizeof(double*));
	  liste_pb_bis[k] = allocVec(dim);    /* on prend de la place pour les deux problèmes suivants*/
	  liste_pb_bis[k+1] = allocVec(dim);
	}
	nb_pb_a_traiter = nb_pb_a_traiter+2;
	if (ij < dim){
	  while((ij < dim) && vecteur[ij]){
	    vecteur[ij] = false;
	    ij++;
	  }
	  if(ij<dim){
	    vecteur[ij] =true;
	    ij = 0;
	  }
	}
      }
      while (ij< dim );
    }  
    if (nb_pb_a_traiter != nb_pb2) { 
      if (LANGUAGE){
	fprintf(stderr,"error in cutting, k= %d, nb_pb2 = %d.\n", nb_pb_a_traiter, nb_pb2);
      }
      else{
	fprintf(stderr,"erreur dans découpage, k= %d, nb_pb2 = %d.\n", nb_pb_a_traiter, nb_pb2);
      }
    }
    free (vecteur);
    vecteur = NULL;  
    if (NIV_AFF >0){
      if (LANGUAGE){
	fprintf(stdout,"%d possible problems (on %d a priori)\n",k, nb_pb_a_traiter);
      }
      else{
	fprintf(stdout,"%d problèmes possibles (sur %d a priori)\n",k, nb_pb_a_traiter);  
      }
    }
    freeMat(*liste_pb,*nb_pb);
    *nb_pb = k;
    (*liste_pb) = copie_matrice(liste_pb_bis,k,dim);
    freeMat(liste_pb_bis,k+2); 
  }
#endif /*OPENMP*/     
}


/* action of the units */

void vecteurs_possibles2(double *probleme,double *unite, double **sigma, double **sigma_inv, double *pas, int r1, int r2, int dim, double **c2, int **minorant, int **majorant){
  /* version plus compacte*/
  int i,j;
  double temp_j_p,temp_j_m ;
  double *c1 = allocVec(dim);
  double *T=allocVec(dim);
  int *t = allocVecEnt(dim);
  double *X0 = allocVec(dim);
  double *minorantR=allocVec(dim);
  double *majorantR=allocVec(dim);
  /*description des centres des problèmes*/
  *c2 = allocVec(dim);
  produit_tordu(unite,probleme,c1,r1,r2); 
  MatVec(T, sigma_inv, c1, dim, dim);
  for (i=0 ; i<dim ; i++) { if (ceil(T[i])-T[i] < EPS3) { t[i] = ceil(T[i]) ;} else  {t[i] = floor(T[i]);} }
  MatVecEnt(X0,sigma,t,dim,dim);
  for (i=0 ; i<dim ; i++) (*c2)[i] = c1[i] - X0[i];
  /* calcul des bornes à partir du centre translaté*/
  for (i=0 ; i<dim; i++){
    minorantR[i] =0;
    majorantR[i] =0;
    for (j=0 ; j< dim; j++){
      temp_j_p= (*c2)[j];
      temp_j_m= (*c2)[j];
      if (sigma_inv[i][j] >0){
	temp_j_m -= module(unite,j,r1,r2) * module(pas,j,r1,r2);
	temp_j_p += module(unite,j,r1,r2) * module(pas,j,r1,r2);
      }
      else{
	temp_j_m += module(unite,j,r1,r2) * module(pas,j,r1,r2);
	temp_j_p -= module(unite,j,r1,r2) * module(pas,j,r1,r2);
      }
      minorantR[i] += sigma_inv[i][j] * temp_j_m;
      majorantR[i] += sigma_inv[i][j] * temp_j_p;
    }
    (*minorant)[i] = partie_entiere(minorantR[i],EPS3);
    (*majorant)[i] = partie_entiere(majorantR[i],EPS3);
  }
  freeVec(X0);
  freeVecEnt(t);
  freeVec(T);
  freeVec(minorantR);
  freeVec(majorantR);
  freeVec(c1);
}



bool intersection_unite( double *c2, int *vec, double *unite, double **sigma, double **liste_pb, int nb_pb, double *pas, int dim, int r1, int r2, bool *deja_vu) {
  int i,j;
  double *vect = allocVec(dim);
  /*fprintf(stdout,"vec[0] = %d, vec[1] = %d, vec[2] = %d\n",vec[0],vec[1],vec[2]);*/
  double *translate = allocVec(dim);
  bool trouve;
  bool intersection = false;
  MatVecEnt(vect,sigma,vec,dim,dim);
  /*affiche_vecteur(vect,dim);*/
  for (i=0; i< dim ; i++) translate[i] = c2[i] -vect[i];
  /*fprintf(stdout,"t[0] = %f, t[1] = %f, t[2] = %f\n",translate[0],translate[1],translate[2]);*/
  j=0 ;
  while ( (j< nb_pb) && (!intersection)){
    if(!deja_vu[j]){
      i=0;
      trouve = true; /* dit si on intersecte avec un problème */
      while ( (i < r1+r2) && trouve){ /* les composantes complexes sont inutiles, elles donnent la même chose */
	trouve = (module_diff(translate,liste_pb[j],i,r1,r2) <= ( (1 + module(unite,i,r1,r2)) * module(pas,i,r1,r2) + EPS2) );
	i++; 
      }
      j++;
      if (trouve) intersection=true;
    }
    else{
      j++;
    }
  }
  freeVec(vect);
  freeVec(translate);
  return intersection; /* 'true' si intersection, 'false' si aucune intersection*/
}


void test_des_unites(double ***liste_pb, int *nb_pb, double **sigma, double **sigma_inv, double *unite, double *pas, int dim, int r1, int r2){
  if (*nb_pb ==0){
    return ; /*pas de problème, ce n'est pas la peine de faire le test*/
  }
  else{
#ifdef OPENMP
    int i,n,**minorant,**majorant,**vec,ik;
    struct listvec *n_prob;
    bool *deja_vu=malloc(*nb_pb * sizeof(bool));
    for(i=0 ; i < *nb_pb ; i++){
      deja_vu[i] = false;
    }
    n = omp_get_max_threads();
    minorant = allocMatEnt(n,dim);
    majorant = allocMatEnt(n,dim);
    vec = allocMatEnt(n,dim);
    n_prob = malloc(n*sizeof(struct listvec));
    for(ik=0 ; ik < n ; ik++){
      initlist(&(n_prob[ik]));
    }
#else
    int i,j;
    bool intersection;
    int *minorant = allocVecEnt(dim);
    int *majorant = allocVecEnt(dim);
    int *vec = allocVecEnt(dim);
    bool *deja_vu=malloc(*nb_pb * sizeof(bool));
    for(i=0 ; i < *nb_pb ; i++){
      deja_vu[i] = false;
    }  
#endif /*OPENMP*/
#ifdef OPENMP
#pragma omp parallel for
    for(i=0 ; i < *nb_pb ; i=i+2){
      double *c2;
      bool intersection;
      int j,ij;
      intersection=0;
      ij=omp_get_thread_num();
      /*fprintf(stdout,"Thread %d\n",ij);*/
      vecteurs_possibles2((*liste_pb)[i],unite, sigma, sigma_inv, pas, r1, r2,dim, &c2, &(minorant[ij]),
			  &(majorant[ij]));
      for (j=0 ; j< dim ; j++) vec[ij][j] = minorant[ij][j];
      vec[ij][0] =vec[ij][0]-1;
      j=0;
      while (!intersection && j<dim){
	if(vec[ij][j]>= majorant[ij][j]){
	  vec[ij][j] = minorant[ij][j];
	  j++;
	}
	else{
	  vec[ij][j] = vec[ij][j]+1;
	  j=0;
	  intersection =  intersection_unite(c2, vec[ij], unite, sigma, *liste_pb, *nb_pb, pas, dim, r1, r2,deja_vu);
	}
      }
      freeVec(c2);
      if (intersection){
	append(&(n_prob[ij]),(*liste_pb)[i]);
	append(&(n_prob[ij]),(*liste_pb)[i+1]);
      }
      else{
	deja_vu[i] = true;
	deja_vu[i+1] = true;
	/*freeVec(liste_pb[i]);*/
      /*freeVec(liste_pb[i+1]);*/
      }
    }
    for(ik=1 ; ik < n ; ik++){
      optimiser(&(n_prob[ik]));
      /*fprintf(stdout,"Thread %d, %ld problèmes.\n",ik,(&(n_prob[ik]))->ind);*/
      merge(&(n_prob[0]),&(n_prob[ik]));
    }
    for(ik=0 ; ik < *nb_pb ; ik++){
      if(deja_vu[ik]){
	freeVec((*liste_pb)[ik]);
      }
    }
    optimiser(&(n_prob[0]));
    *nb_pb = (&(n_prob[0]))->ind;
    *liste_pb = (&(n_prob[0]))->list;
#else /*OPENMP*/
    i=0;
    while ( (i < *nb_pb) ){
      double *c2;
      intersection = 0;
      vecteurs_possibles2((*liste_pb)[i],unite, sigma, sigma_inv, pas, r1, r2,dim, &c2, &minorant,
			  &majorant);
      /*fprintf(stdout,"vec_poss = [");*/
      /*for(j=0 ; j < dim ; j++) fprintf(stdout,"(%d,%d) ", minorant[j],majorant[j]);*/
      /*fprintf(stdout,"]\n");*/
      for (j=0 ; j< dim ; j++) vec[j] = minorant[j];
      vec[0] =vec[0]-1;
      j=0;
      while (!intersection && j<dim){
	if(vec[j]>= majorant[j]){
	  vec[j] = minorant[j];
	  j++;
	}
	else{
	  vec[j] = vec[j]+1;
	  j=0;
	  intersection =  intersection_unite(c2, vec, unite, sigma, *liste_pb, *nb_pb, pas, dim, r1, r2,deja_vu);
	}
      }
      freeVec(c2);
      if (intersection){
	i = i+2;
      }
      else{
	/*{ affichage de l'élimination
	  fprintf(stdout,"Elimination du problème :\n");
	  affiche_vecteur(pb_aux[i],dim);
	  fprintf(stdout,"par l'unité :\n");
	  affiche_vecteur(unites[r],dim);
	  }*/
	for (j=0 ; j < dim ; j++){
	  (*liste_pb)[i][j] = (*liste_pb)[*nb_pb-2][j];
	  (*liste_pb)[i+1][j] = (*liste_pb)[*nb_pb-1][j];
	}
	freeVec((*liste_pb)[*nb_pb-1]);
	freeVec((*liste_pb)[*nb_pb-2]);
	(*nb_pb) = (*nb_pb) -2;
      }
    }
#endif /*OPENMP*/    
#ifdef OPENMP
    for(ik=0 ; ik < n ;ik++){
      freeVecEnt(minorant[ik]);
      freeVecEnt(majorant[ik]);
      freeVecEnt(vec[ik]);
    }
#else /*OPENMP*/
    freeVecEnt(minorant);
    freeVecEnt(majorant);
    freeVecEnt(vec);
#endif /*OPENMP*/
    free (deja_vu);
  }
  return;
}


/* end action of units */  


/* returns the max(r1+r2-1,ordre) first fundamental units*/
int generateur_unites(int nombre, double **unit, double *gen_rac_unit, int ordre, int r1, int r2, int dim, double ***unites){
  int i,j;
  int nb_calc=r1+r2;
  if(nb_calc>nombre){
    nb_calc = nombre;
  }
  *unites = allocMat(nb_calc,dim);
  i=0;
  while (i < nb_calc){
    if( i < r1+ r2-1){
      for (j=0 ; j< dim ; j++) (*unites)[i][j] = unit[i][j];
    }
    else{
      for (j=0 ; j< dim ; j++) (*unites)[i][j] = gen_rac_unit[j];
    }
    i++;
  }
  if (NIV_AFF >0){
    if (LANGUAGE){
      fprintf(stdout,"Units used :\n");
    }
    else{
      fprintf(stdout,"Unités utilisées :\n");
    }
    affiche_matrix(*unites,nb_calc,dim);
  }
  return nb_calc;
}



/* loop: 1/ action of units -> repeated until the number of problems is unchanged
   2/ exit condition: no more problem or too many iterations without improvement
   3/ cutting
*/


/* boucle : 1/ action des unités -> répétition tant que le nombre de problèmes diminue
   2/ Condition de sortie : plus de problèmes ou trop d'itérations sans rien améliorer          
   3/ découpage


   on compte le nombre de découpages qui augmentent le nombre de problèmes
   on se fixe une borne maximale d'augmentation du nombre de problèmes
*/

/* version itérative*/
void boucle_unites_et_decoupage_v4(double **liste_pb, int *nb_pb, double **sigma, double **sigma_inv, double *pas, double **meilleure_liste_pb, int *meilleur_nb_pb, double *meilleur_pas, double **unites, int nb_unites , int dim, int r1, int r2, int nb_max_iterations, int *nb_iterations_en_cours, ListOfIntegers petits_elts, int nb_p_e, double K, int *nb_passages ){
  bool continue_iteration = true;
  int r;
  while ((*nb_pb >0)&&(continue_iteration)){
    int nb_pb_depart ; 
    if (*nb_pb > MAX_NUMBER_PB2){ continue_iteration = false ; } /* on s'arrête s'il y a beaucoup trop de problèmes */
    /* 1ère étape : action des unités */ 
    if(*nb_passages == 0){
      do{
	nb_pb_depart = *nb_pb;
	if (NIV_AFF >0){
	  if (LANGUAGE){
	    fprintf(stdout, "Action of units: ");	    
	  }
	  else{
	    fprintf(stdout, "Action des unités : ");
	  }
	}
        for(r=0 ; r < nb_unites ; r++){
	  test_des_unites(&liste_pb, nb_pb, sigma, sigma_inv, unites[r],pas,dim,r1,r2);
	}
        if (NIV_AFF >0){
	  if (LANGUAGE){
	    fprintf(stdout, "there remains %d problems.\n", *nb_pb);
	  }
	  else{
	    fprintf(stdout, "il reste %d problèmes.\n", *nb_pb);
	  }
	}
	/* si on a éliminé au moins un problème, on relance la boucle "tant que"*/
      }
      while(*nb_pb && *nb_pb < nb_pb_depart);
    }
    else{
      if (NIV_AFF>0){
	if (LANGUAGE){
	  fprintf(stdout,"Pas d'action des unités.\n");
	}
	else{
	  fprintf(stdout,"Pas d'action des unités.\n");
	}
      }
      (*nb_passages)--;
    }
    if ((*nb_pb > *meilleur_nb_pb) &&(*nb_iterations_en_cours >= nb_max_iterations)){
      continue_iteration = false;
    }
    if ((*nb_pb >0)&&(continue_iteration)){ /* on ne fait le découpage que s'il reste des problèmes*/
      if (*nb_pb <= *meilleur_nb_pb){ /* version moins économique, on copie même quand on retrouve autant de problèmes, on suppose
					 qu'ainsi le découpage est plus fin et meilleur*/
	copie_matrice_sans_alloc(liste_pb, meilleure_liste_pb, *nb_pb, dim);
	copie_vecteur_sans_alloc(pas, meilleur_pas, dim);
	if (*nb_pb < *meilleur_nb_pb) *nb_iterations_en_cours = 0; 
	*meilleur_nb_pb = *nb_pb;
      }
      /* on teste si on doit s'arrêter parce qu'on a fait trop d'itérations*/
     
      /* il reste des problèmes, on découpe les problèmes*/
      if (NIV_AFF >0){
	if (LANGUAGE){
	  fprintf(stdout, "Cutting of problems: ");
	}
	else{
	  fprintf(stdout, "Découpage des problèmes : ");
	}
      }
      decoupage_problemes_abs(&liste_pb, nb_pb, pas, dim, r1, r2, petits_elts, nb_p_e, K);
      if (NIV_AFF >0){
	if (LANGUAGE){
	  fprintf(stdout,"%d problems.\n",*nb_pb);
	}
	else{
	  fprintf(stdout,"%d problèmes.\n",*nb_pb);
	}
      }
      if ((*nb_pb >0)){
	if (*nb_pb <= *meilleur_nb_pb){ /* version moins économique, on copie même quand on retrouve autant de problèmes, on suppose
					   qu'ainsi le découpage est plus fin et meilleur */
	  copie_matrice_sans_alloc(liste_pb,meilleure_liste_pb, *nb_pb, dim);
	  copie_vecteur_sans_alloc(pas,meilleur_pas, dim);
	  if (*nb_pb < *meilleur_nb_pb){
            (*nb_iterations_en_cours) = 0;
            *meilleur_nb_pb = *nb_pb;
	  }
	  else{ /* on n'a pas amélioré le nombre de problèmes */
            (*nb_iterations_en_cours) ++; 
	  }
	}
	else{
	  (*nb_iterations_en_cours) ++; 
	}
      }
    }
  }
  if (*nb_pb == 0){
    *meilleur_nb_pb = 0;
  }
  freeMat(liste_pb, *nb_pb);
  freeVec(pas);
  return ;
}
/* fin boucle */


/* FIN BOUCLE ACTION UNITES + DECOUPAGE */

/* AFFICHAGE DU "GRAPHE" */


void vecteurs_possibles3(double *probleme,double *unite, double **sigma, double **sigma_inv, double *pas, int r1, int r2, int dim, double **c2, int **t,int **minorant, int **majorant){
  int i,j;
  double temp_j_p ;
  double temp_j_m ;
  double *c1 = allocVec(dim);
  double *T=allocVec(dim);
  double *X0 = allocVec(dim);
  double *minorantR=allocVec(dim);
  double *majorantR=allocVec(dim);
  /*description des centres des problèmes */
  produit_tordu(unite,probleme,c1,r1,r2); 
  MatVec(T, sigma_inv, c1, dim, dim);
  for (i=0 ; i<dim ; i++) { if (ceil(T[i])-T[i] < EPS3) { (*t)[i] = ceil(T[i]) ;} else  {(*t)[i] = floor(T[i]);} }
  freeVec(T);
  MatVecEnt(X0,sigma,*t,dim,dim);
  for (i=0 ; i<dim ; i++) (*c2)[i] = c1[i] - X0[i];
  /* calcul des bornes à partir du centre translaté */
  for (i=0 ; i<dim; i++){
    minorantR[i] =0;
    majorantR[i] =0;
    for (j=0 ; j< dim; j++){
      temp_j_p= (*c2)[j];
      temp_j_m= (*c2)[j];
      if (sigma_inv[i][j] >0){
	temp_j_m -= module(unite,j,r1,r2) * module(pas,j,r1,r2);
	temp_j_p += module(unite,j,r1,r2) * module(pas,j,r1,r2);
      }
      else{
	temp_j_m += module(unite,j,r1,r2) * module(pas,j,r1,r2);
	temp_j_p -= module(unite,j,r1,r2) * module(pas,j,r1,r2);
      }
      minorantR[i] += sigma_inv[i][j] * temp_j_m;
      majorantR[i] += sigma_inv[i][j] * temp_j_p;
    }
    (*minorant)[i] = partie_entiere(minorantR[i],EPS3);
    (*majorant)[i] = partie_entiere(majorantR[i],EPS3);
  }
  freeVec(X0);
  freeVec(minorantR);
  freeVec(majorantR);
  freeVec(c1);
}

void intersection_unite3( double *c2, int *t, int *vec, double *unite, double **sigma, double **liste_pb, int nb_pb, double *pas, int dim, int r1, int r2, int *pb_imag, int *nb_pb_imag, int **unite_imag, bool *est_testable) {
  int i,j;
  bool trouve;
  double *vect = allocVec(dim);
  double *translate = allocVec(dim);
  MatVecEnt(vect,sigma,vec,dim,dim);
  for (i=0; i< dim ; i++) translate[i] = c2[i] -vect[i];
  /*fprintf(stdout,"t[0] = %f, t[1] = %f, t[2] = %f\n",translate[0],translate[1],translate[2]);*/
  j=0 ;
  for (j= 0 ; j < nb_pb ; j++){
    i=0;
    trouve =true;
    while ( (i < r1+r2)  && trouve ){ /* les autres composantes complexes sont inutiles*/
      trouve = (module_diff(translate,liste_pb[j],i,r1,r2) <= ( (1 + module(unite,i,r1,r2)) * module(pas,i,r1,r2) + EPS2) );
      i++; 
    }
    if (trouve){
      (*nb_pb_imag) ++;
      
      if(*nb_pb_imag > nb_pb){
	*est_testable = false;
      }
      else{
	/*fprintf(stdout,"%d -",*nb_pb_imag -1);*/
	/*fflush(stdout);*/
	(pb_imag)[(*nb_pb_imag) -1]  = j;
	/*fprintf(stdout,"%d via l'entier ",j); */ /*affiche qu'on intersecte le problème j*/
	for(i=0 ; i< dim ; i++){
	  (unite_imag)[(*nb_pb_imag) -1][i] = vec[i]+t[i];
	  /*fprintf(stdout,"%d ", vec[i]+t[i]);*/
	}
	/*fprintf(stdout,", ");*/
      }
    }
  }
  /*
   *pb_imag = realloc(*pb_imag, (*nb_pb_imag) * sizeof(int));
   *unite_imag = realloc(*unite_imag, (*nb_pb_imag)*sizeof(int*));
   for (i=*nb_pb_imag ; i < nb_pb ;i++){
   freeVecEnt( (*unite_imag)[i]);
   }
  */
  freeVec(translate);
  freeVec(vect);
}

bool affichage_action_unite(double **liste_pb, int nb_pb, double **sigma, double **sigma_inv, double *unite, double *pas, int dim, int r1, int r2, int ***pb_imag, int **nb_pb_imag, int ****unite_imag){
  if (nb_pb ==0){
    return true; /*pas de problème, ce n'est pas la peine de faire le test*/
  }
  else{
    /*fprintf(stdout,"Affichage de l'action de l'unité 0 :\n");*/
    int i,j;
    int *minorant = allocVecEnt(dim);
    int *majorant = allocVecEnt(dim);
    int *vec = allocVecEnt(dim);
    double *c2 = allocVec(dim);
    int *t= allocVecEnt(dim);
    bool est_testable = true;
    *pb_imag = malloc(nb_pb*sizeof(int*)); /* les problèmes atteints*/
    *nb_pb_imag = allocVecEnt(nb_pb); 
    *unite_imag = malloc(nb_pb * sizeof(int**));
    for( i=0 ; i < nb_pb ; i++ ){
      (*nb_pb_imag)[i] = 0;
      (*pb_imag)[i] = malloc(nb_pb *sizeof(int));
      (*unite_imag)[i] = malloc(nb_pb *sizeof(int*));
      for (j=0 ; j< nb_pb ; j++){
	(*unite_imag)[i][j] = malloc(dim *sizeof(int));		
      }
      /*fprintf(stdout,"problème %d envoyé sur problème(s) ",i);*/
      vecteurs_possibles3(liste_pb[i],unite, sigma, sigma_inv, pas, r1, r2,dim, &c2,&t, &minorant,
			  &majorant);
      for (j=0 ; j< dim ; j++) vec[j] = minorant[j];
      j=0;	
      vec[0] =vec[0]-1;
      while ((j < dim) && est_testable){
	if (vec[j] >= majorant[j]){
	  vec[j] = minorant[j];
	  j++;
	}
	else{
	  vec[j] = vec[j]+1;
	  j= 0;
	  intersection_unite3(c2, t, vec, unite, sigma, liste_pb, nb_pb, pas, dim, r1, r2, ((*pb_imag)[i]), & ((*nb_pb_imag)[i]),  ((*unite_imag)[i]),&est_testable);
	}
      }
      (*pb_imag)[i] = realloc( (*pb_imag)[i] , (*nb_pb_imag)[i] * sizeof(int));
      for (j = (*nb_pb_imag)[i] ; j < nb_pb ; j++){
	freeVecEnt((*unite_imag)[i][j]);
      }
      (*unite_imag)[i] = realloc((*unite_imag)[i], (*nb_pb_imag)[i] *sizeof(int*));
    }    
    freeVecEnt(minorant);
    freeVecEnt(majorant);
    freeVecEnt(vec);
    freeVec(c2);
    freeVecEnt(t);
    return est_testable;
  }
}


/* FIN AFFICHAGE DU "GRAPHE" */



/* CALCUL DU MINIMUM */

/* entrée à la main des cycles */
int entree_manuelle_cycles(int **lung, int ****beta, double ****betta, double **sigma, int dim){
  int nb_cycles;
  int i,j,k;
  fprintf(stdout,"Entrez le nombre de cycles :");
  fscanf(stdin,"%d",&nb_cycles);
  *lung = allocVecEnt(nb_cycles);
  *beta = malloc(nb_cycles *sizeof(int**));
  *betta = malloc(nb_cycles * sizeof(double **));
  for(i=0 ; i < nb_cycles ; i++){
    fprintf(stdout,"Cycle %d\nEntrez le nombre d'éléments :",i); 
    fscanf(stdin,"%d",&((*lung)[i]));
    (*beta)[i] = allocMatEnt( (*lung)[i],dim );
    (*betta)[i] = allocMat((*lung)[i],dim );
    for(j=0 ; j< (*lung)[i];j++){
      for(k=0 ; k< dim ; k++){
	fprintf(stdout,"beta[%d][%d][%d] =",i,j,k);
	fscanf(stdin,"%d",&((*beta)[i][j][k]));
      }
      MatVecEnt((*betta)[i][j],sigma,(*beta)[i][j],dim,dim);
    }
  }
  return nb_cycles;
}

/* fin entrée à la main des cycles */


/*  ELIMINATION DES PROBLEMES EN 1            */
/* Il se peut que des problèmes se trouvent "au bord"
   Pour simplifier, on préfère tous les voir près de 0
*/

void traitement_pb_1(double **liste_pb, double *pas, double **sigma_inv, double **sigma, int nb_pb, int dim){
  double ** probleme=allocMat(nb_pb,dim);
  double total;
  double total_pas;
  int i,j,k;
  for (k=0 ; k < nb_pb ; k++){
    for (i=0 ; i< dim ; i++){
      total = 0;
      total_pas = 0;
      for (j=0 ; j < dim ; j++){
	if (sigma_inv[i][j] > 0){
	  total_pas += sigma_inv[i][j] * (liste_pb[k][j] + pas[j] ); 
	}
	else{
	  total_pas += sigma_inv[i][j] * (liste_pb[k][j] - pas[j] ); 
	}
	total += sigma_inv[i][j] * liste_pb[k][j];
      }
      if (total_pas >1){
	probleme[k][i] = total -1;
      }
      else{
	probleme[k][i] = total; 
      }
    }
  }
  for (k=0 ; k< nb_pb ; k++){
    MatVec(liste_pb[k],sigma,probleme[k],dim,dim);
  }
}

void action_unite_sur_entier(int *resultat, double* unite, int* entier, int dim, double **sigma_inv, double **sigma, int r1,int r2){
  double *w,*v;
  int i;
  w = allocVec(dim);
  v = allocVec(dim);
  MatVecEnt(w,sigma,entier,dim,dim);
  produit_tordu(unite,w,v,r1,r2);
  MatVec(w,sigma_inv,v,dim,dim);
  for (i=0 ; i< dim ;i++) resultat[i]  = arrondi(w[i]);
  freeVec(v);
  freeVec(w);
}

void traitement_pb_1_bis(double **liste_pb, double *pas, double **sigma_inv, double **sigma, int nb_pb, int dim, int *nb_imag , int **pb_imag, int ***entier_imag, double *unite, int r1,int r2){
  double **probleme,total,total_pas;
  int i,j,k,l,*vecteur,*v_aux;
  vecteur = allocVecEnt(dim); /* vecteur de translation des problemes */
  v_aux = allocVecEnt(dim);
  probleme = allocMat(nb_pb,dim);
  for (k=0 ; k < nb_pb ; k++){
    for (i=0 ; i< dim ; i++){
      total = 0;
      total_pas = 0;
      for (j=0 ; j < dim ; j++){
	if (sigma_inv[i][j] > 0){
	  total_pas += sigma_inv[i][j] * (liste_pb[k][j] + pas[j] ); 
	}
	else{
	  total_pas += sigma_inv[i][j] * (liste_pb[k][j] - pas[j] ); 
	}
	total += sigma_inv[i][j] * liste_pb[k][j];
      }
      if (total_pas >1){
	probleme[k][i] = total -1;
	vecteur[i] = 1;
      }
      else{
	probleme[k][i] = total; 
	vecteur[i] = 0;
      }
    }
    action_unite_sur_entier(v_aux, unite, vecteur, dim, sigma_inv, sigma, r1,r2);
    for (l= 0 ; l< nb_imag[k] ; l++){
      for (i=0 ; i< dim ;i++){
	entier_imag[k][l][i] =  entier_imag[k][l][i] - v_aux[i];
      }			
    }
    for (l=0 ; l < nb_pb ; l++){
      for (j=0 ; j < nb_imag[l]; j++){
	if (pb_imag[l][j] == k){
	  for (i=0 ; i< dim ; i++){
	    entier_imag[l][j][i] = entier_imag[l][j][i] + vecteur[i];
	  }			
	}			
      }
    }
  }
  freeVecEnt(vecteur);
  freeVecEnt(v_aux);
  for (k=0 ; k< nb_pb ; k++){
    MatVec(liste_pb[k],sigma,probleme[k],dim,dim);
  }
  freeMat(probleme,nb_pb);
}



/* FIN  ELIMINATION DES PROBLEMES EN 1   */


/* SIMPLIFICATION DU GRAPHE DES PROBLEMES */


/* REGROUPEMENT PAR IMAGE */

bool vecteurs_egaux(int *v1,int *v2, int dim){
  int i;
  bool b;
  i=0;
  b=true;
  while (b&& (i<dim)){
    b = (v1[i]==v2[i]);
    i++;
  }
  return b;
}

bool est_regroupable(int i, int nb_prob, int *nb_pb_imag, int ***entier_imag, int dim){
  /*fprintf(stdout,"i=%d\n",i);*/
  if ((nb_pb_imag[i] <2)||(i >=nb_prob)){
    return false;
  }
  else{
    int k=1;
    bool b=true;
    while ((k<nb_pb_imag[i])&&b){
      b= (vecteurs_egaux(entier_imag[i][0], entier_imag[i][k],dim));
      k++;
    }
    return b;	
  }
}

bool appartient_liste_pb(int i, int j, int *nb_imag, int **pb_imag){
  /* dit si 'i' appartient à pb_imag[j] */
  int k;
  bool b;
  k=0;
  b=false;
  while ((k< nb_imag[j])&& !b){
    b = (i== pb_imag[j][k]);
    k++;
  }
  return b;
}

int plus_petit_indice(int i, int j, int *nb_imag, int **pb_imag){
  /* plus petit indice >= i n'appartenant pas à pb_imag[j]*/
  int k=i;
  while ((k< nb_imag[j])&& appartient_liste_pb(i,j,nb_imag,pb_imag)){
    k++;
  }
  return k;	
}

bool test_appartenance(int a, int* l, int lg, int *ind){
  int k=0;
  bool b = true;
  while ((k< lg ) && b){
    b = (l[k] != a);
    k++;
  }
  if (!b){
    *ind = k-1;
  }
  else{
    *ind = lg;
  }
  return !b;	
}

bool regroupement_prob(int nb_prob, int *nb_imag, int **pb_imag, int ***entier_imag, int dim, int *n_nb_prob, int **n_nb_imag, int ***n_pb_imag, int **** n_entier_imag){
  /* essaie de regrouper les problèmes, en commençant par l'image de 0*/
  int i=0;
  while ((i < nb_prob) && !(est_regroupable(i, nb_prob, nb_imag, entier_imag, dim))  ){
    i++;
  }
  if (i== nb_prob){ 
    return false;
  }
  else{
    /* à ce stade, on tente de regrouper les images de 'i'*/
    /* 
       'i' est envoyé sur 0
       les images regroupées forment le problème '1'
       de sorte que 'i' a une seule image, 1
       il faut faire attention si "0 = 1"
    */
    int l;
    int *tab_images ; /* un tableau qui décrit les images des problèmes */
    /* ie i -> 0 , ... */
    int indice_possible = 2;
    int un1 = 1;
    bool zero_un ; 
    *n_nb_prob = nb_prob - nb_imag[i] + 1;
    zero_un = appartient_liste_pb(i,i,nb_imag,pb_imag);
    tab_images = malloc(nb_prob * sizeof(int));
    if (zero_un){
      indice_possible = 1;
      un1 = 0;
    }
    /*fprintf(stdout,"i = %d, un1 = %d\n",i,un1); */
    /*fprintf(stdout,"tab_images = [ ");*/
    for (l= 0 ; l < nb_prob ; l++){
      if (appartient_liste_pb(l, i, nb_imag, pb_imag)){
	tab_images[l] = un1;
      }
      else{
	if (l==i){
	  tab_images[l] = 0;
	}
	else{
	  tab_images[l] = indice_possible;
	  indice_possible++;
	}
      }
      /*fprintf(stdout,"%d ",tab_images[l]);*/
    }
    /*fprintf(stdout,"]\n");*/
    
    *n_nb_imag = malloc( (*n_nb_prob) *sizeof(int));
    *n_pb_imag = malloc( (*n_nb_prob) *sizeof(int*));
    *n_entier_imag = malloc( (*n_nb_prob) *sizeof(int**));
    if (!zero_un){
      (*n_nb_imag)[0] = 1;
      (*n_pb_imag)[0] = malloc(sizeof(int));
      (*n_pb_imag)[0][0] = 1;
      (*n_entier_imag)[0] = malloc(sizeof(int*));
      (*n_entier_imag)[0][0] = malloc(dim*sizeof(int));
      for (l=0 ; l< dim ; l++){
	(*n_entier_imag)[0][0][l] = entier_imag[i][0][l];
      }
    }
    {
      int m,n,nb_pb_reel;
      int ind;
      int pb;
      bool b;
      for (l=0 ; l< nb_prob ; l++){
	nb_pb_reel = 0;
	/*if (!((appartient_liste_pb(l, i, nb_imag, pb_imag))|| l==i)){*/
	if ((l != i)&& (tab_images[l] != un1)){
	  (*n_nb_imag)[tab_images[l]] = nb_imag[l];
	  (*n_pb_imag)[tab_images[l]] = malloc( nb_imag[l]*sizeof(int));
	  (*n_entier_imag)[tab_images[l]] = malloc(nb_imag[l] *sizeof(int*));
	  for (m = 0 ; m < nb_imag[l] ; m++){
	    /*if (!appartient_liste_pb(pb_imag[l][m],i,nb_imag,pb_imag)){*/
	    if (tab_images[pb_imag[l][m]] != un1){
	      (*n_pb_imag)[tab_images[l]][nb_pb_reel] = tab_images[pb_imag[l][m]];
	      (*n_entier_imag)[tab_images[l]][nb_pb_reel] = malloc(dim*sizeof(int));
	      for (n = 0; n < dim ; n++){
		(*n_entier_imag)[tab_images[l]][nb_pb_reel][n] = entier_imag[l][m][n];	
	      }
	      nb_pb_reel++;
	    }
	    else{
	      b =true;
	      n=0;
	      while (b && n < nb_pb_reel){
		b= ((*n_pb_imag)[tab_images[l]][n] != un1);
		n++;
	    }
	      if (b){
		(*n_pb_imag)[tab_images[l]][nb_pb_reel] = un1;
		(*n_entier_imag)[tab_images[l]][nb_pb_reel] = malloc(dim*sizeof(int));
		for (n = 0; n < dim ; n++){
		  (*n_entier_imag)[tab_images[l]][nb_pb_reel][n] = entier_imag[l][m][n];	
		}
		nb_pb_reel++;
	      }
	    }
	  }
	  (*n_nb_imag)[tab_images[l]] = nb_pb_reel;
	  (*n_pb_imag)[tab_images[l]] = realloc((*n_pb_imag)[tab_images[l]], nb_pb_reel*sizeof(int));
	  (*n_entier_imag)[tab_images[l]] = realloc((*n_entier_imag)[tab_images[l]],nb_pb_reel *sizeof(int*));
	}
      }
      /*reste à traiter le cas de 1*/
      (*n_nb_imag)[un1] = (*n_nb_prob); /* a priori*/
      (*n_pb_imag)[un1] = malloc( (*n_nb_prob)*sizeof(int));
      (*n_entier_imag)[un1] = malloc((*n_nb_prob) *sizeof(int*));
      nb_pb_reel = 0;
      for (l=0 ; l < nb_imag[i] ; l++){
	for (n= 0 ; n < nb_imag[pb_imag[i][l]] ; n++){
	  pb = pb_imag[pb_imag[i][l]][n];
	  if (!(test_appartenance(tab_images[pb], (*n_pb_imag)[un1] , nb_pb_reel, &ind))){
	    (*n_pb_imag)[un1][nb_pb_reel] = tab_images[pb];
	    (*n_entier_imag)[un1][nb_pb_reel] = malloc (dim * sizeof(int));
	    for (m=0 ; m < dim ;m++){
	      (*n_entier_imag)[un1][nb_pb_reel][m] = entier_imag[pb_imag[i][l]][n][m];
	      /*fprintf(stdout, "e_i[%d][%d][%d] = %d. ", i,l,m,entier_imag[i][l][m]);*/
	    }
	    nb_pb_reel++;
	  }
	  else{
	    if (!(vecteurs_egaux((*n_entier_imag)[un1][ind], entier_imag[pb_imag[i][l]][n], dim))){
	      if (LANGUAGE){
		fprintf(stdout,"Simplification problem: different vectors.\n");
	      }
	      else{
		fprintf(stdout,"Problème de simplification : vecteurs différents.\n");
	      }
	      return false;		
	    }
	  }
	}
      }
      (*n_nb_imag)[un1] = nb_pb_reel;
      (*n_pb_imag)[un1] = realloc((*n_pb_imag)[un1],nb_pb_reel * sizeof(int));
      (*n_entier_imag)[un1] = realloc((*n_entier_imag)[un1],nb_pb_reel *sizeof(int*));
      return true;
    }
  }
}


/* SUPPRESSION DES SOMMETS SANS IMAGE */

void suppression_sans_image(int i, int *nb_prob, int **nb_imag, int ***pb_imag, int ****entier_imag, int dim){
  int a,b,c;
  for (a = 0 ; a < *nb_prob; a++){
    for (b= 0 ; b < (*nb_imag)[a] ; b++){
      if ((*pb_imag)[a][b] == i){
	((*nb_imag)[a]) --;
	(*pb_imag)[a][b] = (*pb_imag)[a][(*nb_imag)[a]];
	(*pb_imag)[a] = realloc((*pb_imag)[a], ((*nb_imag)[a]) *sizeof(int*));
	for (c=0 ; c < dim ; c++){
	  (*entier_imag)[a][b][c] = (*entier_imag)[a][(*nb_imag)[a]][c];
	}
	(*entier_imag)[a] = realloc( (*entier_imag)[a], ((*nb_imag)[a]) *sizeof(int*));
      }
      else{
	if ((*pb_imag)[a][b] == *nb_prob -1){
	  (*pb_imag)[a][b] = i;
					
	}				
      }
    }		
  }
  (*nb_prob)--;
  for (b=0 ; b < (*nb_imag)[i] ; b++){
    free ((*entier_imag)[i][b]);
  }
  (*nb_imag)[i] = (*nb_imag)[*nb_prob];
  freeVecEnt((*pb_imag)[i]);
  (*pb_imag)[i] = malloc(((*nb_imag)[*nb_prob]) *sizeof(int));
  free ((*entier_imag)[i]);
  (*entier_imag)[i] = malloc(((*nb_imag)[*nb_prob]) *sizeof(int*));
  for (b=0 ; b < (*nb_imag)[i] ; b++){
    (*pb_imag)[i][b] = (*pb_imag)[*nb_prob][b];
    (*entier_imag)[i][b] = malloc( dim *sizeof(int));
    for (c=0; c < dim ; c++){
      (*entier_imag)[i][b][c] = (*entier_imag)[*nb_prob][b][c];	
    }
  }
  for (b=0 ; b < (*nb_imag)[*nb_prob] ; b++){
    freeVecEnt((*entier_imag)[*nb_prob][b]);
  }
  freeVecEnt((*pb_imag)[*nb_prob]);
  free ((*entier_imag)[*nb_prob]);
  (*nb_imag) = realloc( (*nb_imag), (*nb_prob) *sizeof(int));
  (*pb_imag) = realloc( (*pb_imag), (*nb_prob) *sizeof(int*));
  (*entier_imag) = realloc( (*entier_imag), (*nb_prob) *sizeof(int**));
  return;
}

void suppression_sommets_sans_image(int *nb_prob, int **nb_imag, int ***pb_imag, int ****entier_imag, int dim){
  /*fprintf(stdout,"Nettoyage des sommets sans image...");*/
  int i= 0;
  while (i < *nb_prob){
    if ((*nb_imag)[i] == 0){
      suppression_sans_image(i,nb_prob, nb_imag, pb_imag, entier_imag, dim);
      i=0;
    }
    else{
      i++;
    }
  }
  /*fprintf(stdout,"-> %d problèmes.\n",*nb_prob);*/
  return;
}


/* FIN SUPPRESSION DES SOMMETS SANS IMAGE*/



/* SUPPRESSION DES SOMMETS NON ATTEINTS*/

void sommets_atteints(int nb_prob, int *nb_imag, int **pb_imag, bool *t){
  int i,j;
  for (i=0 ; i < nb_prob ; i++) t[i] = false;
  for (i=0 ; i < nb_prob ; i++){
    for (j= 0 ; j < nb_imag[i] ; j++){
      t[pb_imag[i][j]]= true;			
    }		
  }
  return;	
}

void suppression_sommet_non_atteint(int i,int *nb_prob, int **nb_imag, int ***pb_imag, int ****entier_imag, int dim){
  /* fprintf(stdout,"Là-");
     fprintf(stdout,"Suppression du sommet %d.\n", i);
  */	
  int a,b,j,k;
  /* suppression du sommet i
     1ère étape : i <-> nb_prob -1
  */	
  (*nb_prob)--;
	
  for (a= 0 ; a < *nb_prob ; a++){
    for (b= 0 ; b < (*nb_imag)[a] ; b++){
      if ((*pb_imag)[a][b] == *nb_prob ){
	(*pb_imag)[a][b] = i;	
      }
    }
  }
  if (i != *nb_prob){
    for (a= 0 ; a < (*nb_imag)[i] ; a++){
      freeVecEnt((*entier_imag)[i][a]);
    } 
    free ((*entier_imag)[i]);
    (*nb_imag)[i] = (*nb_imag)[*nb_prob];
    freeVecEnt((*pb_imag)[i]),
      (*pb_imag)[i] = malloc( (*nb_imag)[i] *sizeof(int));
    (*entier_imag)[i] = malloc( (*nb_imag)[i] *sizeof(int*));
    for (j=0 ; j < 	(*nb_imag)[i] ;j++){
      /*fprintf(stdout, "p[%d][%d] = %d\n", *nb_prob, j, (*pb_imag)[*nb_prob][j]);*/
      if ( (*pb_imag)[*nb_prob][j] != *nb_prob){ 
	(*pb_imag)[i][j] = (*pb_imag)[*nb_prob][j];
      }
      else{
	(*pb_imag)[i][j] = i;
      }
      (*entier_imag)[i][j] = malloc(dim *sizeof(int));
      for (k= 0 ; k < dim ;k++){
	(*entier_imag)[i][j][k] = (*entier_imag)[*nb_prob][j][k];
      }
    }
  }	
  /* 2nde étape : suppression*/
  for (j=0 ; j < (*nb_imag)[*nb_prob] ; j++){
    freeVecEnt((*entier_imag)[*nb_prob][j]);
  }
  freeVecEnt((*pb_imag)[*nb_prob]);
  free ((*entier_imag)[*nb_prob]);
  /*fprintf(stdout,"Ici--\n");*/
  (*nb_imag) = realloc((*nb_imag), (*nb_prob) *sizeof(int));
  (*pb_imag) = realloc((*pb_imag), (*nb_prob) *sizeof(int*));
  (*entier_imag) = realloc((*entier_imag), (*nb_prob) *sizeof(int**));
  return;
}

void nettoyage_sommets_non_atteints(int *nb_prob, int **nb_imag, int ***pb_imag, int ****entier_imag, int dim){
  /*return; 
    fprintf(stdout,"Nettoyage des sommets non atteints...");
  */
  bool *t = malloc( (*nb_prob) * sizeof(bool));
  bool b =false;
  int i;
  while (!b){
    /*fprintf(stdout,"nb_prob = %d\n",*nb_prob);*/
    sommets_atteints(*nb_prob, *nb_imag, *pb_imag, t);
    i=0;
    b=true;
    while ((i < *nb_prob)&& b){
      b = t[i];
      if (b){ 
	i++;
      }
    }
    if (!b){
      suppression_sommet_non_atteint(i,nb_prob, nb_imag, pb_imag, entier_imag, dim);
      /*affichage_graphe(*nb_prob, *nb_imag, *pb_imag, *entier_imag, dim);*/
    }
  }
  free (t);
  /*fprintf(stdout," -> %d problèmes.\n", *nb_prob);*/
  /*affichage_graphe(*nb_prob, *nb_imag, *pb_imag, *entier_imag, dim);*/
  return;
	
}

/* FIN SOMMETS ATTEINTS*/


void regroupement_des_problemes(int nb_pb, int *nb_imag, int **pb_imag, int ***entier_imag, int dim, int *c1_nb_prob, int **c1_nb_imag, int ***c1_pb_imag, int ****c1_entier_imag  ){ 
  bool reg=true;
  int nombre;
  int aa, ab, ac;
  /* copie 1*/
  *c1_nb_prob = nb_pb;
  *c1_nb_imag = malloc(*c1_nb_prob * sizeof(int));
  *c1_pb_imag = malloc(*c1_nb_prob * sizeof(int*));
  *c1_entier_imag = malloc(*c1_nb_prob * sizeof(int**));
  for (aa = 0 ; aa < (*c1_nb_prob) ; aa++){
    (*c1_nb_imag)[aa] = nb_imag[aa];
    (*c1_pb_imag)[aa] = malloc(nb_imag[aa] * sizeof(int));
    (*c1_entier_imag)[aa] = malloc(nb_imag[aa] * sizeof(int*));
    for (ab = 0 ; ab < nb_imag[aa]; ab++){
      (*c1_pb_imag)[aa][ab] = pb_imag[aa][ab];
      (*c1_entier_imag)[aa][ab] = malloc(dim * sizeof(int));
      for (ac = 0 ; ac < dim ; ac++){
	(*c1_entier_imag)[aa][ab][ac] = entier_imag[aa][ab][ac];
      }
    }
  }
	
  /* premier nettoyage*/
  suppression_sommets_sans_image(c1_nb_prob, c1_nb_imag, c1_pb_imag, c1_entier_imag, dim);
  nettoyage_sommets_non_atteints(c1_nb_prob, c1_nb_imag, c1_pb_imag, c1_entier_imag, dim);
  /*affichage_graphe(c1_nb_prob, c1_nb_imag, c1_pb_imag, c1_entier_imag, dim);*/
	
  /* copie 2 -> déclaration*/
  {
    int c2_nb_prob;
    int *c2_nb_imag;
    int **c2_pb_imag;
    int ***c2_entier_imag;
    int l,m;
    
    while (reg){
      reg =  regroupement_prob(*c1_nb_prob, *c1_nb_imag, *c1_pb_imag, *c1_entier_imag, dim, &c2_nb_prob, &c2_nb_imag, &c2_pb_imag, &c2_entier_imag);
      /*fprintf(stdout,"Simplification.\n");*/
      if (reg){
	/* libération*/
	for (l=0 ; l < (*c1_nb_prob) ;l++){
	  for (m=0 ; m < (*c1_nb_imag)[l] ; m++){
	    freeVecEnt((*c1_entier_imag)[l][m]);
	  }
	  freeVecEnt((*c1_pb_imag)[l]);
	  free ((*c1_entier_imag)[l]);
	}
	freeVecEnt((*c1_nb_imag));
	free ((*c1_pb_imag));
	free ((*c1_entier_imag));
	/* copie*/
	*c1_nb_prob = c2_nb_prob;
	*c1_nb_imag = malloc( *c1_nb_prob * sizeof(int));
	*c1_pb_imag = malloc(*c1_nb_prob * sizeof(int*));
	*c1_entier_imag = malloc(*c1_nb_prob * sizeof(int**));
	for (aa = 0 ; aa < *c1_nb_prob ; aa++){
	  (*c1_nb_imag)[aa] = c2_nb_imag[aa];
	  (*c1_pb_imag)[aa] = malloc(c2_nb_imag[aa] * sizeof(int));
	  (*c1_entier_imag)[aa] = malloc(c2_nb_imag[aa] * sizeof(int*));
	  for (ab = 0 ; ab < c2_nb_imag[aa]; ab++){
	    (*c1_pb_imag)[aa][ab] = c2_pb_imag[aa][ab];
	    (*c1_entier_imag)[aa][ab] = malloc(dim * sizeof(int));
	    for (ac = 0 ; ac < dim; ac++){
	      (*c1_entier_imag)[aa][ab][ac] = c2_entier_imag[aa][ab][ac];
	    }
	    freeVecEnt(c2_entier_imag[aa][ab]);
	  }
	  freeVecEnt(c2_pb_imag[aa]);
	  free (c2_entier_imag[aa]);
	}
	freeVecEnt(c2_nb_imag);
	free (c2_pb_imag);
	free (c2_entier_imag);			
      }
      else{
	nombre = *c1_nb_prob;
	/*fprintf(stdout,"Simplification--- -> %d problèmes\n", c1_nb_prob);*/
	suppression_sommets_sans_image(c1_nb_prob, c1_nb_imag, c1_pb_imag, c1_entier_imag, dim);
	/*affichage_graphe(c1_nb_prob, c1_nb_imag, c1_pb_imag, c1_entier_imag, dim);*/
	nettoyage_sommets_non_atteints(c1_nb_prob, c1_nb_imag, c1_pb_imag, c1_entier_imag, dim);
	if (*c1_nb_prob == nombre){
	  if (NIV_AFF >0){
	    affichage_graphe(*c1_nb_prob, *c1_nb_imag, *c1_pb_imag, *c1_entier_imag, dim);
	  }
	}
	else{
	  reg = true;
	}
			
			
      }
    }
  }
  /*
    for(l=0 ; l < c1_nb_prob ; l++){
    for (m=0 ; m < c1_nb_imag[l] ; m++){
    freeVecEnt(c1_entier_imag[l][m]);
    }
    freeVecEnt(c1_pb_imag[l]);
    free (c1_entier_imag[l]);				
    }
    free (c1_pb_imag);
    freeVecEnt(c1_nb_imag);	
    free (c1_entier_imag);
  */
  return;
}





/* FIN SIMPLIFICATION DU GRAPHE DES PROBLEMES */

/* GRAPHE CONVENABLE */
/* On cherche à déterminer si on a un graphe convenable.
   En ce cas, on donne les cycles disjoints du graphe
*/


/* test de graphe convenable */



bool est_convenable_cfc(int nb_pb, int *nb_imag, int **pb_imag, int ***entier_imag, int *nb_cycles, int **lg_cycles, int ****cycles, int dim, double ****betta, double **sigma){
  if (nb_pb > MAX_PB_CONV){
    return false;
  }
  else{
    int i=0;
    int j,k,l;
    int dep = 0;
    int *comp;
    bool nouveau_cycle = true;
    int cycle_en_cours = -1;
    int ind=0;
    int nb;
    bool b= true;
    bool* deja_vu = malloc(nb_pb *sizeof(bool));
    *nb_cycles = 0;
    comp = cfc(nb_pb, nb_imag, pb_imag,nb_cycles);
    *cycles = malloc(*nb_cycles *sizeof(int**));
    *lg_cycles = malloc(nb_pb *sizeof(int));
    for (j=0 ; j< nb_pb ; j++){ 
      if(comp[j]==0){
	deja_vu[j] = true;
      }
      else{
	deja_vu[j] = false;
      }
    }
    while (b && (dep <nb_pb)){
      if (nouveau_cycle){
	dep = premier_indice_faux(deja_vu,nb_pb);
	if (dep < nb_pb){
	  nouveau_cycle = false;
	  cycle_en_cours++ ;
	  (*cycles)[cycle_en_cours] =allocMatEnt(nb_pb, dim);
	  i= dep;
	  j=0;
	  deja_vu[dep] = true;
	}
      }
      else{
	nb = unique_image_cfc(i,pb_imag[i],nb_imag[i],comp,comp[dep],&ind);
	if (nb ==0){
	  nouveau_cycle = true;
	}
	else{
	  b=  b && (nb==1);
	  for (k=0 ; k < dim ; k++){
	    (*cycles)[cycle_en_cours][j][k] = entier_imag[i][ind][k];
	  }
	      
	  if ((pb_imag)[i][ind] == dep){
	    nouveau_cycle = true;
	    (*lg_cycles)[cycle_en_cours] = j+1;
	    for (l = j+1 ; l < nb_pb ; l++){
	      freeVecEnt((*cycles)[cycle_en_cours][l]);			
	    }
	    (*cycles)[cycle_en_cours] = realloc( (*cycles)[cycle_en_cours], (j+1) * sizeof(int *));
	  }
	  else{
	    i = pb_imag[i][ind];
	    b = b && (!(deja_vu[i]));
	    deja_vu[i] = true;
	    j++;	
	  }
	}

      }
    }
    if (NIV_AFF > 0){
      if(b){
	if (LANGUAGE){
	  fprintf(stdout,"The graph is convenient.\n");
	}
	else{
	  fprintf(stdout,"Le graphe est convenable.\n");
	}
      }
      else{
	if (LANGUAGE){
	  fprintf(stdout,"The graph is not convenient.\n");
	}
	else{
	  fprintf(stdout,"Le graphe n'est pas convenable.\n");
	}
      }
    }
    fflush(stdout);
    if (!b){/* la décomposition a échoué, on libère la mémoire*/
      /*TODO*/
    }
    else{/* calcul de betta */
      *betta = malloc(*nb_cycles *sizeof(double**));
      for(i=0 ; i < *nb_cycles ; i++){
	(*betta)[i] = allocMat((*lg_cycles)[i],dim );
	for(j=0 ; j< (*lg_cycles)[i];j++){
	  MatVecEnt((*betta)[i][j],sigma,(*cycles)[i][j],dim,dim);
	}
      }
    }
    free (deja_vu);
    free (comp);
    return b;
  }
}

GEN
traitement_resultat(GEN v, double *res)
{
  long i, j = 0, l = lg(v); 
  GEN M = gen_0, w = cgetg(l, t_VEC);
  for (i = 1; i < l; i++) {
    GEN m, z = gel(v, i);
    if (lg(z) < 4) continue;
    m = gel(z,3);
    switch(gcmp(m, M))
      {
      case 0: j++; gel(w,j) = gel(v,i); break;
      case 1: M = m; j = 1; gel(w,j) = gel(v,i); break;
      }
  }
  setlg(w, j+1); *res = gtodouble(M); return w;
}
#ifdef SQLITE
void calcul_minimum_liste_cycles(sqlite3* db, int ***cycl, int nb_cyc, int* lg_cyc, int ind_unite, char *pol, int dim, double min,double *res, char* fichier,GEN Q,char *fichier_resultat, double K,bool *success){
#else
void calcul_minimum_liste_cycles(int ***cycl, int nb_cyc, int* lg_cyc, int ind_unite, char *pol, int dim, double min,double *res, char* fichier,GEN Q,char *fichier_resultat, double K,bool *success){
#endif
  GEN k;
  GEN cyc;
  GEN mini;
  GEN resultat;
  int i,j,ij;
  int petite_norme = 0;
  char* mini1;
  int nb_orb;
  int T1;
  char *C1 ;
#ifdef SQLITE
  int r1,r2; /* just used to add the number field to the database */
  char *poly;
#endif
  pari_sp ltop=avma;
  /*fprintf(stdout,"ici");
    fprintf(stdout,"%s",pol);
    P = gp_read_str(pol);
    if (typ(P) == t_VEC) P = gtopoly(P, 0);
    pariprintf("%Ps\n",Q);*/
  k = gel(Q,13);
  /*pariprintf("%Ps\n",k);*/
  mini = dbltor(min);
  cyc = cgetg(nb_cyc+1, t_VEC);
  for(i=1; i <= nb_cyc; i++){
    gel(cyc,i) = cgetg(lg_cyc[i-1]+1,t_VEC);
    for(j=1; j <= lg_cyc[i-1] ; j++){
      gel(gel(cyc,i),j) = cgetg(dim+1, t_VEC);
      for(ij=1; ij <= dim ; ij++){
	gel(gel(gel(cyc,i),j),ij) = stoi(cycl[i-1][j-1][ij-1]);
      }
    }
  }
  if(PRINCIPAL==2){
    petite_norme=1;
  }
  resultat=calcul_min_liste_cycles(cyc,gel(member_fu(k),ind_unite), k, mini,petite_norme,DEFAULTPREC);
  resultat = traitement_resultat(resultat, res);
  /*pari_printf("%Ps\n",resultat);*/
  /*
  switchout(fichier_resultat);
  pari_printf("%Ps",resultat);
  */
  /* insertion des résultats obtenus dans la base de données en cas de minimum >= K.*/
  if(gcmp(gel(gel(resultat,1),3),dbltor(K))>=0){
    mini1 = GENtostr(gel(gel(resultat,1),3));
    nb_orb = glength(resultat);
    T1 = itos(gel(gel(resultat,1),2));
    C1 = malloc( 40000 );
    strcpy(C1,GENtostr(geval(gel(gel(resultat,1),1))));
    for(i=2 ; i <= nb_orb ; i++){
      T1+= itos(gel(gel(resultat,i),2));
      strcat(C1, " ; ");
      strcat(C1,GENtostr(geval(gel(gel(resultat,i),1))));
    }
#ifdef SQLITE
    poly = GENtostr(member_pol(k));
    r1 = itos(member_r1(k));
    r2 = itos(member_r2(k));
#endif
    if (LANGUAGE){
      fprintf(stdout,"The minimum is %s, reached at %d ",mini1,T1);
    }
    else{
      fprintf(stdout,"Le minimum est %s, atteint en %d ",mini1,T1);
    }
    if(T1>1){
      if(LANGUAGE){
	fprintf(stdout,"critical points");
      }
      else{
	fprintf(stdout,"points critiques ");
      }
    }
    else{
      if(LANGUAGE){
	fprintf(stdout,"critical point");
      }
      else{
	fprintf(stdout,"point critique ");
      }
    }
    fprintf(stdout,": %s.\n",C1);
    fflush(stdout);
#ifdef SQLITE
    ajout_numberfield(db,r1+2*r2,r1,r2,GENtostr(member_disc(k)),GENtostr(gabs(member_disc(k),DEFAULTPREC)),itos(member_no(member_clgp(k))),poly,mini1,T1,C1,*res,gcmpgs(gel(gel(resultat,1),3),1)); 
    sqlite3_close(db);
#endif
    avma=ltop;
    pari_close();
    *success=true;
  } 
  else{
    if (NIV_AFF >0){
      if (LANGUAGE){
	pari_printf("The minimum found (%Ps) is smaller than %f.\n",gel(gel(resultat,1),3),K);
      }
      else{
	pari_printf("Le minimum trouvé %Ps est plus petit que %f.\n",gel(gel(resultat,1),3),K);
      }
    }
    fflush(stdout);
    avma=ltop;
    *success = false;
  }
}

/* FIN GRAPHE CONVENABLE */


#ifdef SQLITE
int main_loop(sqlite3* db, numberfield *corps,double *minim, int classes_ideaux, double K, int nb_de_decoupage, int nombre_unites_utilisees, double *pas, int *nob, double **y, ListOfIntegers uu, int nombre, char *fichier, char* pol , double K4,GEN Q, FILE *fich, pari_timer T, long temps, char *fichier_resultat, int nb_p){
#else /*SQLITE*/
int main_loop(numberfield *corps,double *minim, int classes_ideaux, double K, int nb_de_decoupage, int nombre_unites_utilisees, double *pas, int *nob, double **y, ListOfIntegers uu, int nombre, char *fichier, char* pol , double K4,GEN Q, FILE *fich, pari_timer T, long temps, char *fichier_resultat, int nb_p){
#endif /*SQLITE*/
  int i,j ; /* des indices dont on se servira*/
  
  /* 0 -> calcul réussi
     1 -> pas de problème
     2 -> trop de problèmes (sûr)
     3 -> trop de problèmes (mais peut-être à regarder plus tard)
  */
  int* nb_passages = &nb_p;
  int ii;
  int nb_pb, nb_pb_initial;
  long temps2;
  double **problemes;  
  double **meilleure_liste_pb;
  double *meilleur_pas;
  int meilleur_nb_pb;
  int nb_it_en_cours;
  int nb_max_iterations;

  /*int nombre_unites_utilisees;*/
  int nb_unites;
  double **unites;

  struct listvec liste_pb;
  struct listvec *pblemes;
  bool alamain = false; /* for hand treatment, should NOT be used */
  
  int *lg;
  int ***beta;
  double ***betta;
  int nb_cycles ;

  int **pb_imag; 
  int *nb_imag; 
  int ***entier_imag;


  int nb_pb2;
  int **pb_imag2; 
  int *nb_imag2; 
  int ***entier_imag2;

  bool est_testable ;
  bool success ; 

  K= K*corps->norme_ideal;
  if (alamain){
    goto ALAMAIN;
  }
  /* premier découpage */
  if (NIV_AFF >0){
    if (LANGUAGE){
      fprintf(stdout, "First cutting and absorption test.\n");
    }
    else{
      fprintf(stdout, "Premiers découpage et test d'absorption...\n");
    }
    fflush(stdout);
  }

  /*decoupage_initial4(y, sigma, sigma_inv, DIM, R1,R2,nob, uu, nombre, &problemes, &nb_pb, pas, K,fich);*/
  temps2  = TIMER(&T);
  temps+=temps2;
  if (NIV_AFF >0){
    if (LANGUAGE){
      fprintf(stdout,"Time : %ld ms.\n",temps2);
    }
    else{
      fprintf(stdout,"Temps : %ld ms.\n",temps2);
    }
  }
  pblemes = & liste_pb;
  initlist(pblemes);
  decoupage_initial5(y,corps,nob,uu,nombre,pblemes,&nb_pb,pas,K,fich);
  problemes = pblemes ->list;
  nb_pb = pblemes ->ind;
  temps2 = TIMER(&T);
  temps+=temps2;
  if (NIV_AFF >0){
    if(LANGUAGE){
      fprintf(stdout,"Time (initial cutting): %ld ms.\n",temps2);
    }
    else{
      fprintf(stdout,"Temps (découpage initial) : %ld ms.\n",temps2);
    }
  }
  /*fprintf(stdout,"(%d,%ld)\n", nb_pb,pblemes->ind);
    freeMat(y,DIM);
    freeVecEnt(nob);
    affichage_problemes(problemes,pas,sigma_inv,nb_pb,DIM);
  */
  /* fin premier découpage*/
  /* boucle unités et découpage*/
  nb_pb_initial = nb_pb;
  if (NIV_AFF >0){
    if (LANGUAGE){
      fprintf(stdout,"Initial number of problems: %d\n", nb_pb);
    }
    else{
      fprintf(stdout,"Nombre initial de problèmes : %d\n", nb_pb);
    }
  }
  if (nb_pb >= MAX_NUMBER_PB){
    freeMat(problemes,nb_pb);
    /*freeMat(uu,nombre+1);*/
    return 3;
  }
  meilleure_liste_pb=copie_matrice(problemes,nb_pb,corps->dim);
  meilleur_pas=copie_vecteur(pas,corps->dim);
  meilleur_nb_pb = nb_pb;
  nb_it_en_cours = 0;
  nb_max_iterations =MAXITERATIONS;
  
  /*int nombre_unites_utilisees;*/
  /*
    if (argc <=4){
    fprintf(stdout, "entrer le nombre d'unités utilisées : ");
    fscanf(stdin, "%d", &nombre_unites_utilisees);
    }
    else{
    nombre_unites_utilisees = atoi(argv[4]);
    }
  */
  nb_unites = generateur_unites(nombre_unites_utilisees, corps->unit, corps->gene_1, corps->ordre_1,corps->r1,corps->r2, corps->dim, &unites);
  if (NIV_AFF >0){
    if (LANGUAGE){
      fprintf(stdout,"%d unit(s) will be used.\n",nb_unites);
    }
    else{
      fprintf(stdout,"On utilise %d unité(s).\n",nb_unites);
    }
  }
  boucle_unites_et_decoupage_v4(problemes, &nb_pb, corps->sigma, corps->sigma_inv, pas, meilleure_liste_pb, &meilleur_nb_pb, meilleur_pas, unites, nb_unites , corps->dim, corps->r1,corps->r2, nb_max_iterations, &nb_it_en_cours, uu, nombre, K, nb_passages );
  
  /* translation retardée
       translation en 0
       traitement_pb_1(meilleure_liste_pb,meilleur_pas,sigma_inv,sigma,meilleur_nb_pb,DIM);
       fin translation en 0
  */
  
  
  /* fin boucle unités et découpages*/
  if (meilleur_nb_pb==0){
    if (NIV_AFF >0){
      if (LANGUAGE){
	fprintf(stdout, "No problem, we will try with k < %f.\n",K/corps->norme_ideal);
      }
      else{
	fprintf(stdout, "Pas de problème, on va essayer k < %f.\n",K/corps->norme_ideal);
      }
    }
    freeMat(meilleure_liste_pb,nb_pb_initial);
    freeVec(meilleur_pas);
    /*freeMat(uu,nombre+1);*/
    freeMat(unites,nb_unites);
    return 1;
  }
  else{
    if (NIV_AFF > 0){
      if (LANGUAGE){
	fprintf(stdout, "%d problem(s) were obtained in the best case.\n", meilleur_nb_pb);
      }
      else{
	fprintf(stdout, "On a obtenu %d problème(s) dans le meilleur cas.\n", meilleur_nb_pb);
      }
    }
  }

  /* translation de 1 en 0*/
  est_testable = affichage_action_unite(meilleure_liste_pb, meilleur_nb_pb, corps->sigma, corps->sigma_inv, unites[0], meilleur_pas, corps->dim,  corps->r1, corps->r2,&pb_imag,&nb_imag, &entier_imag);
  
  if (est_testable){
    traitement_pb_1_bis(meilleure_liste_pb, meilleur_pas, corps->sigma_inv, corps->sigma, meilleur_nb_pb, corps->dim, nb_imag , pb_imag, entier_imag, unites[0],corps->r1,corps->r2);
    /* fin translation de 1 en 0*/
    
    if (NIV_AFF >0){
      if (LANGUAGE){
	fprintf(stdout,"Translation from 1 to 0. ");
	fprintf(stdout,"Simplification of the graph.\n");
      }
      else{
	fprintf(stdout,"Translation de 1 en 0. ");
	fprintf(stdout,"Simplification du graphe.\n");
      }
    }
    /*affichage_problemes(meilleure_liste_pb,meilleur_pas,sigma_inv,meilleur_nb_pb,DIM,aff);*/
    regroupement_des_problemes(meilleur_nb_pb, nb_imag, pb_imag, entier_imag, corps->dim,  &nb_pb2, &nb_imag2, &pb_imag2, &entier_imag2);
  }
  /* suppression du graphe*/
  if (est_testable){
    for (i=0 ; i < meilleur_nb_pb ; i++){
      for (j=0 ; j < nb_imag[i] ; j++){
	freeVecEnt(entier_imag[i][j]);
      }
      freeVecEnt(pb_imag[i]);
      free (entier_imag[i]);
    }
  }
  freeVecEnt(nb_imag);
  free (pb_imag);
  free (entier_imag);
  
  /* fin suppression du graphe*/
  
  if (!est_testable){
    if (NIV_AFF > 0 ){
      if (LANGUAGE){
	fprintf(stdout,"The decomposition of the graph into disjoint circuits failed (too many arrows).\n");
      }
      else{
	fprintf(stdout,"La décomposition du graphe en cycles disjoints a échoué (trop de flèches).\n");
      }
    }    
    return(2);
  }

 ALAMAIN:;
 

  if(alamain){
    fprintf(stdout,"Décomposition à la main.\n");
    nb_cycles = entree_manuelle_cycles(&lg, &beta, &betta, corps->sigma, corps->dim);
  }
  else{
    if (! (est_convenable_cfc(nb_pb2, nb_imag2, pb_imag2, entier_imag2, &nb_cycles, &lg, &beta, corps->dim, &betta,corps->sigma))){
      if( NIV_AFF >0){
	if (LANGUAGE){
	  fprintf(stdout,"The decomposition of the graph into disjoint circuits failed.\n");
	}
	else{
	  fprintf(stdout,"La décomposition du graphe en cycles disjoints a échoué.\n");
	}
      }
      return(2);  	
      /*nb_cycles = entree_manuelle_cycles(&lg, &beta, &betta, sigma, DIM);*/
    }
    else{
      if( NIV_AFF >0){
	if (LANGUAGE){
	  fprintf(stdout,"Le graphe se décompose en %d cycle", nb_cycles);
	}
	else{
	  fprintf(stdout,"The graph has %d circuit", nb_cycles);
	}
	if (nb_cycles >1){
	  fprintf(stdout,"s");
	}
	fprintf(stdout,".\n");
	for (i = 0 ; i< nb_cycles ; i++){
	  if (LANGUAGE){
	    fprintf(stdout,"cycle %d (%d élément", i , lg[i]);
	  }
	  else{
	    fprintf(stdout,"circuit %d (%d element", i , lg[i]);
	  }
	  if (lg[i] > 1){
	    fprintf(stdout,"s");
	  }
	  fprintf(stdout,") :\n");
	  for (j=0 ; j < lg[i] ; j++){
	    if (LANGUAGE){
	      fprintf(stdout,"integer[%d] = [", j);
	    }
	    else{
	      fprintf(stdout,"entier[%d] = [", j);
	    }
	    for (ii=0 ; ii < corps->dim ; ii++){
	      fprintf(stdout,"%d ", beta[i][j][ii]);
	      if (ii < corps->dim-1){
		fprintf(stdout,"; ");
	      }
	      else{
		fprintf(stdout,"]\n");
	      }
	    }
	  }
	}
      }
    }
  }

#ifdef SQLITE
  calcul_minimum_liste_cycles(db,beta,nb_cycles, lg, 1, pol , corps->dim, K4, minim, fichier,Q,fichier_resultat, K, &success);
#else
  calcul_minimum_liste_cycles(beta,nb_cycles, lg, 1, pol , corps->dim, K4, minim, fichier,Q,fichier_resultat, K, &success);
#endif  
  if(!success){
    /* On a trouvé un minimum plus petit que K*/
    return(1);
  }

  /* suppression du graphe*/
  for (i=0 ; i < nb_pb2 ; i++){
    for (j=0 ; j < nb_imag2[i] ; j++){
      freeVecEnt(entier_imag2[i][j]);
    }
    freeVecEnt(pb_imag2[i]);
    free (entier_imag2[i]);
  }
  freeVecEnt(nb_imag2);
  free (pb_imag2);
  free (entier_imag2);
  
  /* fin suppression du graphe */
  

  /*libérer proprement beta, betta et lg*/
  
  freeMat(meilleure_liste_pb,nb_pb_initial);
  freeVec(meilleur_pas);
  for (i=0 ; i< nb_cycles ; i++){
    for (j= 0 ; j < lg[i] ; j++){
      freeVecEnt(beta[i][j]);
      freeVec(betta[i][j]);
    }
    free (beta[i]);
    free (betta[i]);
  }
  free (beta);
  free (betta);
  freeVecEnt(lg);
  /*freeMat(uu,nombre+1);*/
  freeMat(unites,nb_unites);
  return 0; 
}



void affichage_donnees_corps(GEN Q, numberfield* corps){
  GEN P = member_pol(gel(Q,13));
  if (LANGUAGE){
    pariprintf("polynomial =%Ps, n = %d , r1 = %d , r2 = %d.\n",P,corps->dim,corps->r1,corps->r2);
    fprintf(stdout,"norm of the ideal considered = %d.\n", corps->norme_ideal);
  }
  else{
    pariprintf("polynôme =%Ps, n = %d , r1 = %d , r2 = %d.\n",P,corps->dim,corps->r1,corps->r2);
    fprintf(stdout,"norme de l'idéal considéré = %d.\n", corps->norme_ideal);
  }
  /*if (NIV_AFF >0)FIXME*/{
    fprintf(stdout,"sigma = \n");
    affiche_matrix(corps->sigma, corps->dim, corps->dim);  
    fprintf(stdout,"sigma^(-1) = \n");
    affiche_matrix(corps->sigma_inv, corps->dim,corps->dim);
    if (LANGUAGE){
      fprintf(stdout,"units (transposed)= \n");
    }
    else{
      fprintf(stdout,"unités (transposées)= \n");
    }
    affiche_matrix(corps->unit, corps->r1+corps->r2-1,corps->dim);
    if (LANGUAGE){
      fprintf(stdout,"There are %d roots of unity in the number field, a generator is \n",corps->ordre_1);
    }
    else{
      fprintf(stdout,"Les racines de l'unité sont d'ordre %d, un générateur est \n",corps->ordre_1);
    }
    affiche_vecteur(corps->gene_1,corps->dim); 
  }
  return;
}



int main(int argc, char ** argv){
  int nb_threads=1;
  pari_timer T;
  long temps = 0;
  double min;
  int i;
  int classes_ideaux = 0;
  /*fprintf(stdout,"nombre de découpage =%d\n",nb_de_decoupage);*/
  int nb_unites;
  double K;
  int nb_passages;
  double K2;
  double K3;
  double K4;
#ifdef SQLITE
  sqlite3 *db;
  int rc;
#endif /*SQLITE*/

  char **ecriture = malloc(sizeof(char*));
  char ecr[] = "result";

  char** fichier_resultat = malloc(sizeof(char*));
  
  FILE* fichier;
  double temp;
	
  numberfield corps_de_nombres;
  numberfield *corps;
  GEN Q;

  int nb_de_decoupage;

  int *nob;
  
  double *pas;
  double *sauv_pas ;
  double **y ;

  long int nb_po;
  long temps2 ;

  ListOfIntegers po;

  int resultat;
  bool reussi;

  GEN P;
#ifdef SQLITE
  GEN k; /* just used to add to the database */
#endif
  FILE* config_file;

  /* reading the config file */
  config_file = fopen(CONFIG_FILE,"r");
  read_config(config_file);
  fclose(config_file);
  /* done */
  nb_unites = NB_UNITS ;
  K=INITIAL_VALUE_K;
  nb_passages = ITERATIONS_WITHOUT_UNITS ; 
  TIMERstart(&T);
#ifdef OPENMP
  {
    nb_threads= omp_get_max_threads();
    if (NIV_AFF > 1){
      if (LANGUAGE){
	fprintf(stdout,"We use openmp (%d threads will be used at most).\n",nb_threads);
      }
      else{
	fprintf(stdout,"On utilise open mp (%d threads au plus).\n",nb_threads);
      }
    }
  }
#endif /*OPENMP*/
  if(PRINCIPAL ==2){
    classes_ideaux=1;
  }
#ifdef SQLITE
  rc = sqlite3_open(DATABASE_NAME, &db);
  if( rc ){
    if (LANGUAGE){
      fprintf(stderr, "Failed to open the database: %s\n", sqlite3_errmsg(db));
    }
    else{
      fprintf(stderr, "Impossible d'ouvrir la base de données : %s\n", sqlite3_errmsg(db));
    }
    sqlite3_close(db);
    exit(1);
  }
  creation_table_numberfield(db);
#endif /*SQLITE*/
  /*fprintf(stdout,"%s\n",argv[1]);*/
  if(argc>4){
    K= atof(argv[4]);
  }
  if(argc>3){
    *ecriture = argv[3];
  }
  else{
    *ecriture = ecr ;
  }

  if(argc>2){
    *fichier_resultat = argv[2];
  }
  else{
    *fichier_resultat = "stdout";
  }
  fichier =fopen(*ecriture,"a");

  corps = &corps_de_nombres;
#ifdef SQLITE
  Q=lancement_pari(argc,argv,corps,classes_ideaux,&K2,&K4,db);
#else /*SQLITE*/
  Q=lancement_pari(argc,argv,corps,classes_ideaux,&K2,&K4);
#endif /*SQLITE*/	
  if((PRINCIPAL == 0) &&(corps->h) >1){
    if (LANGUAGE){
      fprintf(stdout,"The field is not principal.\n");
    }
    else{
      fprintf(stdout,"Le corps n'est pas principal.\n");
    }
    fprintf(fichier,DISPLAY_NON_PRINCIPAL);
    return 0;
  }
  if((PRINCIPAL == 2) && (corps->h == 1)){
    if (LANGUAGE){
      fprintf(stdout,"The field is principal principal.\n");
    }
    else{
      fprintf(stdout,"Le corps est principal.\n");
    }
    fprintf(fichier,DISPLAY_NON_PRINCIPAL);
    return 0;
  }

  nb_de_decoupage = nb_de_dec(corps);
  if (argc>5){
    nb_de_decoupage = atoi(argv[5]);
  }

  if ( corps->r1+corps->r2 <2){
    if(LANGUAGE){
      fprintf(stdout,"The algorithm does not work in this case, but you can do without it.\n");
    }
    else{
      fprintf(stdout,"L'algorithme ne fonctionne pas en ce cas, mais on peut s'en passer.\n");
    }
    if(corps->r1==1){
      if(LANGUAGE){
	fprintf(stdout,"The Euclidean minimum is 1/2.\n");
      }
      else{
	fprintf(stdout,"Le minimum euclidien vaut 1/2.\n");
      }
    }
    else{
      GEN P=gel(Q,4);
      GEN mini = minimum_quad_im(P);
      if(LANGUAGE){
	pariprintf("The Euclidean minimum is %Ps.\n",mini);
      }
      else{
	pariprintf("Le minimum euclidien est %Ps.\n",mini);
      }
      pari_close();
    }
    fprintf(stdout,"Temps : %ld ms.\n",TIMER(&T));
    return 0;
  }
	
  /*fprintf(stdout,"Affichage des information sur le corps : \n");*/
  if(NIV_AFF >0){
    affichage_donnees_corps(Q, corps); 
  }
  nob =  malloc(corps->dim*sizeof(int));
	  
  for (i=0 ; i< corps->dim ; i++){
    /*fprintf(stdout,"ATTENTION, nb de découpage personnalisé.\n");*/   
    nob[i] = nb_de_decoupage-1;
  }
  /*
    fprintf(stdout,"ATTENTION, nb de découpage personnalisé.\n");
    nob[0]=49;
    nob[1]=699;
    nob[2]=699;
  */

  sauv_pas=malloc(corps->dim*sizeof(double));
  y = calcul_y(corps->sigma, corps->dim, nob, &pas);
  for(i=0 ; i< corps->dim ; i++){
    sauv_pas[i] = pas[i];
  }

  /* calcul de nb_elements et small_elements*/
  /*if (NIV_AFF >0) fprintf(stdout,"Calcul du nombre d'entiers utilisés...");*/
  temps2 = TIMER(&T);
  temps+=temps2;
  po= small_elts(corps, K, liminf(corps->dim),limsup(corps->dim),nb_threads,&nb_po);
  temps2=TIMER(&T);
  /*fprintf(stdout,"On utilisera %ld entiers.\n",nb_po);*/
  if (NIV_AFF >0){
    if(LANGUAGE){
      fprintf(stdout,"Time of computation of the integers: %ld ms.\n",temps2);
    }
    else{
      fprintf(stdout,"Temps de calcul des entiers : %ld ms.\n",temps2);
    }
  }
  /* fin calcul de nb_elements et small_elements*/

  /* K2 est fixé à 1/norme minimale, mais cette valeur est possible, on diminue donc K2*/
  K2 = 0.9*K2;
  K3=K2;
  K4=K; /*FIXME donner une vraie borne*/
#ifdef SQLITE	  
  resultat= main_loop(db,corps,&min, classes_ideaux,K,nb_de_decoupage,nb_unites,pas,nob,y,po,nb_po,*ecriture, argv[1],K4,Q,fichier,T,temps,*fichier_resultat,nb_passages);
#else
  resultat= main_loop(corps,&min, classes_ideaux,K,nb_de_decoupage,nb_unites,pas,nob,y,po,nb_po,*ecriture, argv[1],K4,Q,fichier,T,temps,*fichier_resultat,nb_passages);
#endif
  reussi = true;
  while ((resultat !=0) && (K4>MINIMAL_VALUE_K)){	
    if (K4-K2 < 0.000000001){
      reussi = false;
      resultat = 0;
    }
    if (resultat != 0 ){
      if (resultat == 1){
	if(K<K4) {K4=K;}
	K = (5*K+K3)/6 ;
	K3 = (K2+8*K3)/9;
      }
      if (resultat == 2){
	if(K==K4){
	  K4+=5;
	  K=K4;
	}
	else{
	  temp = K;
	  K = (2*K+5*K4)/7;
	  if(temp > K2) { K2 = temp; }
	  if(K2>K3) { K3 = K2 ;}
	} 
      }
      else{
	if (resultat == 3){
	  if(K<K4){
	    K3 = K;
	    K = (K+K4)/2;
	  }
	  else{
	    K4 = K4 + 2;
	    K = K +2;
	  }
	}
      }
      pas=malloc(corps->dim*sizeof(double));
      for(i=0 ; i< corps->dim ; i++){
	pas[i] = sauv_pas[i];
      }
      if(K4>MINIMAL_VALUE_K){
	if(K<MINIMAL_VALUE_K){
	  K = fmax( ((float) MINIMAL_VALUE_K) -0.01,K);
	}
	if( NIV_AFF >0)
	  fprintf(stdout,"K2=%f, K3 = %f, K=%f, K4 = %f\n",K2,K3,K,K4);
#ifdef SQLITE
	resultat= main_loop(db,corps, &min, classes_ideaux,K,nb_de_decoupage,nb_unites,pas,nob,y,po,nb_po,*ecriture,argv[1],K4,Q,fichier,T,temps,*fichier_resultat, nb_passages);
#else
	resultat= main_loop(corps, &min, classes_ideaux,K,nb_de_decoupage,nb_unites,pas,nob,y,po,nb_po,*ecriture,argv[1],K4,Q,fichier,T,temps,*fichier_resultat, nb_passages);
#endif
      }		
    }
  }
  if (reussi){
    if(K4<=MINIMAL_VALUE_K){
      if (NIV_AFF >0){
	fprintf(stdout,"Le corps a un minimum plus petit que la valeur max %f.\n",MINIMAL_VALUE_K);
      }
      if(MINIMAL_VALUE_K == 1){
	/* Le corps est euclidien, on va le stocker. */
	P = gp_read_str(argv[1]);
	if (typ(P) == t_VEC) P = gtopoly(P, 0);
#ifdef SQLITE
	k = bnfinit0(P,0,NULL,DEFAULTPREC);
	ajout_numberfield2(db, corps->dim, corps->r1, corps->r2, GENtostr(member_disc(k)), GENtostr(gabs(member_disc(k),DEFAULTPREC)), itos(member_no(member_clgp(k))), GENtostr(member_pol(k)), -1);
#endif 
      }
      fprintf(fichier,DISPLAY_INF_MIN);
    }
  }
  else{
    if( LANGUAGE){
      fprintf(stdout,"The algorithm failed, the minimum is roughly %f, you should try again with other parameters.\n",K2);
    }
    else{
      fprintf(stdout,"L'algorithme a échoué, mais le minimum doit valoir environ %f.\n",K2);
    }
  }

  freeVec(sauv_pas);
  freeVecEnt(nob);
  freeListOfIntegers(po,nb_po,nb_threads);
  freeMat(y,corps->dim);
  /*fprintf(fichier,"\n");*/
  free_numberfield(corps);
  fclose(fichier);
  free (ecriture);
  if(LANGUAGE){
    fprintf(stdout,"Time: %ld ms.\n",TIMER(&T));
  }
  else{
    fprintf(stdout,"Temps : %ld ms.\n",TIMER(&T));
  }
  return 0;
}
