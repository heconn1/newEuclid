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

/* functions for graph handling */

#include "header.h"

/* display graph of problems */
void affichage_graphe(int nb_sommets, int *nb_aretes, int **liste_but_aretes, int ***liste_etiqu_aretes, int dim){
  int i,j,k;
  if (LANGUAGE){
    fprintf(stdout,"We found %d problem",nb_sommets);
  }
  else{
    fprintf(stdout,"On trouve %d problème",nb_sommets);
  }
  if (nb_sommets >1){
    fprintf(stdout,"s");
  }
  fprintf(stdout,".\n");
  for (i=0 ; i < nb_sommets ; i++){
    if (LANGUAGE){
      fprintf(stdout,"\n PROBLEM N°%d sent on\n", i);
    }
    else{
      fprintf(stdout,"\n PROBLEME N°%d envoyé sur \n", i);
    }
    for (j=0 ; j < nb_aretes[i] ; j++){
      if (LANGUAGE){
	fprintf(stdout,"%d through the integer [", liste_but_aretes[i][j]);
      }
      else{
	fprintf(stdout,"%d via l'entier [", liste_but_aretes[i][j]);
      }
      for (k= 0 ; k< dim-1 ; k++){
	fprintf(stdout,"%d ,", liste_etiqu_aretes[i][j][k]);	
      }
      fprintf(stdout,"%d]\n", liste_etiqu_aretes[i][j][dim-1]);
    }
  }
  fflush(stdout);
}

/* FIN AFFICHAGE GRAPHE PROBLEMES */


int tarjan(int *i, struct pile_ent *p,int nb_pb, int *nb_imag, int **pb_imag, int v, int *marque, bool *succ, int *k, int *comp){
  int j,w,m,lg;
  *i = *i +1;  
  marque[v] = *i;
  m = marque[v];
  push(p,v);
  for(j=0 ; j < nb_imag[v] ; j++){
    w= pb_imag[v][j];
    succ[v] = (v==w) || succ[v];
    if(marque[w] == 0){
      m = min(m,tarjan(i,p,nb_pb,nb_imag,pb_imag,w,marque,succ,k,comp));
    }
    else{
      if(comp[w]==0){
	m = min(m,marque[w]);
      }
    }
  }
  /*fprintf(stdout,"(m,v,marque[v]) = (%d,%d,%d)\n",m,v,marque[v]);*/
  if(marque[v] == m){
    *k= *k+1;
    /*    fprintf(stdout,"CFC n°%d:\n",*k);*/
    lg=0;
    do{
      pop(p,&j);
      /*if(indice[j] == indice[v]){*/
      comp[j] = *k;
      lg++;
      /*}*/
      /*fprintf(stdout,"%d,",j);*/
    }
    while (j != v);
    /*fprintf(stdout,"\n");*/
    if (lg==1){
      if(!(succ[v])){
	*k = *k - 1;
	comp[v] = -1;
      }
    }
  }
  return m; 
}


int *cfc(int nb_pb, int *nb_imag, int **pb_imag, int *k){
  struct pile_ent pile_pb;
  struct pile_ent *p = & pile_pb;
  int i, j,*marque, *comp;
  bool *emp;
  initpile(p);
  marque = allocVecEnt(nb_pb);
  emp = malloc(nb_pb * sizeof(bool));
  comp = allocVecEnt(nb_pb);
  i=0;
  for(j=0 ; j < nb_pb ; j++){
    marque[j] = 0;
    comp[j] = 0;
    emp[j] = false;
  }
  i=0;
  for(j= 0 ; j < nb_pb ; j++){
    if ((marque[j] ==0) && (comp[j]==0)){
      tarjan(&i,p,nb_pb,nb_imag,pb_imag,j,marque,emp,k,comp);
    }
  }
  for(j=0 ; j < nb_pb ; j++){
    if(comp[j] < 0){
      comp[j] = 0;
    }
    if (NIV_AFF >1)
      fprintf(stdout,"comp[%d] = %d\n", j, comp[j]);
  }
  /*
    for(j=0 ; j < nb_pb ; j++){
    fprintf(stdout,"indice[%d] = %d\n", j, indice[j]);
    }
    for(j=0 ; j < nb_pb ; j++){
    fprintf(stdout,"marque[%d] = %d\n", j, marque[j]);
    }
  */
  free (marque);
  free (emp);
  freepile (p);
  if (NIV_AFF >0){
    if (LANGUAGE){
      fprintf(stdout,"We obtained %d strongly connected components.\n",*k);
    }
    else{
      fprintf(stdout,"On obtient %d composantes fortement connexes.\n",*k);
    }
  }
  return comp;
}


/* pour le test de graphe convenable*/

int premier_indice_faux(bool *tab , int lg){
  int i =0;
  bool b =true;
  while ( (i< lg) && b){
    b = tab[i];
    i++;
  }
  if (b){
    return i;
  }
  else{
    return i-1;
  }
}

int unique_image_cfc(int i,int *pb_imag, int nb_imag, int *comp,int k, int *ind){
  /* dit si i a au plus une image dans la cfc.*/
  int nb = 0;
  int j;
  for(j=0 ; j < nb_imag ; j++){
    if(comp[pb_imag[j]] == k){
      nb++;
      *ind = j;
    }
  }
  return(nb);
}
