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



/* This file contains basic functions to get data on number fields  */
/* as well as functions to compute the local Euclidean minimum      */
/* given circuits of reunions of parellelotopes.                    */
/* At the beginning, some functions to go from pari to standard     */
/* objects.                                                         */


#include <pari/pari.h>
#include "header.h"


/* Protoype */
void GENtoMat(double*** matrice, GEN M);
double* affGENtoVec(double *v, GEN w);
int* affGENtoVecEnt(int *v, GEN w);
#ifdef SQLITE
GEN calcul_pari2(GEN P, long prec, numberfield* k , int c_i, double *K2, double *K4,sqlite3 *db);
GEN lancement_pari(int argc, char **argv, numberfield *k, int c_i, double *K2, double *K4,sqlite3 *db);
#else
GEN calcul_pari2(GEN P, long prec, numberfield* k , int c_i, double *K2, double *K4);
GEN lancement_pari(int argc, char **argv, numberfield *k, int c_i, double *K2, double *K4);
#endif
GEN matrice_plongement(GEN k);
GEN calcul_plongement(GEN xi, long r1, long r2, long n, GEN m, GEN l, long prec);
GEN calcul_gamma(GEN k, GEN mini, GEN m, GEN l, long prec);
GEN calcul_bornes_grossieres(GEN gamma_k, GEN m, long prec);
long rebnfprinc(GEN k, GEN a);
long est_principal(GEN k, GEN a, GEN b);
GEN calcul_pt_pb(GEN k, GEN ent, GEN unite, GEN base_entiers);
GEN calcul_norme(GEN point, GEN ent, GEN base_entiers_t, long prec);
GEN matrice_des_entiers(GEN l, GEN ll);
GEN polmod_vecteur(GEN p, long m);
GEN vecteur_pol(GEN v, long m);
GEN reduction_mod_entiers(GEN p, GEN po, GEN a, GEN b, long n);
GEN orbite3(GEN p, GEN k, GEN ll);
void bornes_fines(long j, GEN borne, GEN xi2, GEN v, GEN gamma_k, GEN k, GEN m, GEN ll, long n, long prec,GEN **t);
long est_entier_non_nul_ni_unite(GEN q, long prec);
GEN calcul_minimum_un_point_orbite(GEN point, GEN borne_g, GEN gam, long dim, GEN base_entiers, GEN k, GEN m, GEN ll, long prec);
void parcours(long i, GEN point, GEN borne_g, GEN gam, GEN base_entiers, GEN k, GEN dim, GEN m, GEN ll, GEN v, GEN xi2, GEN nor, long prec);
GEN calcul_minimum(GEN point, GEN k, GEN borne_g, GEN gam, long dim, GEN base_entiers, GEN m, GEN ll, long prec);
GEN calcul_min_cycle(GEN ent, GEN unite, GEN borne_g, GEN gam, long dim, GEN base_entiers, GEN k, GEN m, GEN ll, GEN pts_deja_vus, long prec);
GEN calcul_min_point(GEN p, GEN k, GEN mini, long prec);
GEN ideal_plus_petit(GEN k, long petite_n);
GEN sqfp(GEN n);
GEN calcul_minimum_un_point_orbite_dim_qcq(GEN point, GEN borne_g, GEN gam, GEN base_entiers, GEN k, GEN m, GEN ll, int n,long prec);
GEN liste_inv_entiers_pas_unites(GEN k, long n);
GEN calcul_rapide(GEN k, long n, GEN s);
/* End of prototype */



void GENtoMat(double*** matrice, GEN M)
{
  //  long i, j, l = lg(M), h = lg(M[1]); 
  long i, j, l = lg(M), h = lg(gel(M,1)); 
  for (i = 1; i < h; i++)
    for (j = 1; j < l; j++) {(*matrice)[i-1][j-1] = gtodouble(gcoeff(M,i,j)); }
}

double* affGENtoVec(double *v, GEN w){
  long i, l = lg(w);
  for (i = 1; i < l; i++) v[i-1] = gtodouble(gel(w,i));
  return v;
}

int* affGENtoVecEnt(int *v, GEN w){
  long i, l = lg(w);
  for (i = 1; i < l; i++) v[i-1] = gtolong(gel(w,i));
  return v;
}






/* computation of the matrix of the embedding      */
/* 1 .. r_1: real ones                             */
/* r_1+1 .. r_1+r_2: real part of the complex ones */
/* r_1+r_1+1 .. n: imag part of the complex ones   */
GEN matrice_plongement(GEN k){
  int n,r1,r2,ii,jj;
  GEN m, l, ret;
  pari_sp ltop,lbot;
  ltop=avma;
  k = checkbnf(k);
  r1 = nf_get_r1(bnf_get_nf(k));
  r2 = nf_get_r2(bnf_get_nf(k));
  n = r1+2*r2;
  m = cgetg(n+1, t_MAT);
  for (ii = 1; ii <= n; ii++)
    gel(m, ii) = cgetg(n+1, t_COL);
  for(ii=1 ; ii <= r1 ; ii++){
    for(jj=1 ; jj <= n ;jj++)
      gcoeff(m,ii,jj)= gcoeff(gel(gel(bnf_get_nf(k), 5), 1),ii,jj);
  }
  for(ii=r1+1 ; ii <= r1+r2 ; ii++){
    for(jj=1 ; jj <= n ; jj++){
      gcoeff(m,ii,jj)= greal(gcoeff(gel(gel(bnf_get_nf(k), 5), 1),ii,jj));
      gcoeff(m,ii+r2,jj)= gimag(gcoeff(gel(gel(bnf_get_nf(k), 5), 1),ii,jj));
    }
  }
  l = lll(m);
  /*pari_printf("determinant=%Ps\n", det(m));*/
  lbot=avma;
  ret = cgetg(3, t_VEC);
  gel(ret, 1) = gmul(m, l);
  gel(ret, 2) = gcopy(l);
  return gerepile(ltop,lbot,ret);
}

/* computation of the embedding of xi in a nf                 */
/* the first case deals with the case where xi is rational    */
/* i.e. xi is a constant polynomial                           */
GEN calcul_plongement(GEN xi, long r1, long r2, long n, GEN m, GEN l, long prec){
  GEN a, b;
  long i,j;
  pari_sp ltop,lbot;
  ltop=avma;
  if (degree(lift(xi)) <= 0)
  {
    b = cgetg(n+1, t_VEC);
    for (i = 1; i <= n; ++i){
        gel(b, i) = lift(xi);
    }
  }
  else{
    b = conjvec(xi, prec);
  }
  a = cgetg(2, t_MAT);
  gel(a, 1) = cgetg(n+1, t_COL); 
  for (j = 1; j <= r1; j++){
    gcoeff(a, j, 1) = gcopy(gel(b,j));
  }
  for(j=r1+1; j<n ;j=j+2){
    gcoeff(a,r1+(j-r1-1)/2+1,1)= greal(gel(b,j));
    gcoeff(a,r1+(j-r1-1)/2+1+r2,1)= gimag(gel(b,j+1));
  }
  /*pari_printf("plongement : %Ps\n", gmul(gmul(gmul(m, l), ginv(m)), a));*/
  lbot=avma;
  return gerepile(ltop,lbot,gmul(gmul(gmul(m, l), ginv(m)), a));
}


/* computation of gamma */
GEN calcul_gamma(GEN k, GEN mini, GEN m, GEN l, long prec){
  long i,j,r1,r2,r,n;
  GEN gamma_i,unit,unites,temp,res;
  pari_sp ltop,lbot;
  ltop=avma;
  r1 = nf_get_r1(bnf_get_nf(k));
  r2 = nf_get_r2(bnf_get_nf(k));
  r = (r1 + r2) - 1;
  n = r1 + (2*r2);
  gamma_i = cgetg(maxss(r,r1)+1,t_VEC);
  for(i=1 ; i<= maxss(r,r1) ; i++)
    gel(gamma_i,i) = gen_1;
  unit = conjvec(gmodulo(lift(bnf_get_fu(k)), nf_get_pol(bnf_get_nf(k))), prec);
  unites = cgetg(r+1,t_MAT);
  for(j=1 ; j<= r ; j++){
    gel(unites,j) = cgetg(n+1,t_COL);
    for(i=1 ; i<= r1 ; i++)
      gcoeff(unites,i,j) = gcopy(gcoeff(unit,i,j));
    for(i=r1+1 ; i < n ; i=i+2){
      gcoeff(unites,r1+(i-r1-1)/2+1,j) = greal(gcoeff(unit,i,j));
      gcoeff(unites,r1+(i-r1-1)/2+1+r2,j)=gimag(gcoeff(unit,i+1,j));
    }
  }
  for(i=1 ; i<= r1 ; i++){
    for(j=1 ; j<= r ; j++){
      if (gcmpgs(gabs(gcoeff(unites, i, j), prec), 1) > 0)
	gel(gamma_i, i) = gmul(gel(gamma_i, i), gabs(gcoeff(unites, i, j), prec));
      else
	gel(gamma_i, i) = gdiv(gel(gamma_i, i), gabs(gcoeff(unites, i, j), prec));
    }
  }
  for(i=r1+1 ; i<= r ; i++){
    for(j=1 ; j<= r ; j++){
      temp = gadd(gsqr(gcoeff(unites, i, j)), gsqr(gcoeff(unites, i+r2, j)));
      if (gcmpgs(temp, 1) > 0)
          gel(gamma_i, i) = gmul(gel(gamma_i, i), temp);
        else
          gel(gamma_i, i) = gdiv(gel(gamma_i, i), temp);
    }
  }
  res = mini;
  if(r2 == 0){
    for(j=1 ; j < r1 ; j++)
      res = gmul(res,gel(gamma_i,j));
  }
  else{
    for(j=1 ; j<= r1 ; j++){
      res = gmul(res,gel(gamma_i,j));
    }
    for(j=r1+1 ; j<= r ; j++){
      res = gmul(res,gsqr(gel(gamma_i,j)));
    }
  }
  lbot=avma;
  return gerepile(ltop,lbot,gpow(res, ginv(stoi(n)), prec));
}


/* rough bounds on each coordinate using gamma */
GEN calcul_bornes_grossieres(GEN gamma_k, GEN m, long prec){
  GEN mp, borne,temp;
  long i,j,n;
  mp = ginv(m);
  n = lg(m)-1;
  borne = cgetg(n+1,t_VEC);
  for(i=1 ; i <= n ; i++){
    temp = gen_0;
    for(j=1 ; j <= n ; j++){
      temp=gadd(temp,gabs(gcoeff(mp,i,j),prec));
    }
    gel(borne,i)=gceil(gmul(gamma_k,temp));
  }
  return borne;
}


/* tells if k is principal */
long rebnfprinc(GEN k, GEN a)
{
  GEN b;	  /* vec */
  k = checkbnf(k);
  b = gcopy(gel(bnfisprincipal0(k, a, 1), 1));
  /*pari_printf("principal : %Ps\n", b);*/
  if ((lg(b)-1) == 0)
    return 1;
  else
    return gequal0(gel(b, 1));
  return 0;
}


/* tells if the ideal (a,b) is principal */
long est_principal(GEN k, GEN a, GEN b)
{
  k = checkbnf(k);
  return rebnfprinc(k, idealadd(k, idealhnf0(k, a, NULL), idealhnf0(k, b, NULL)));
}


/* given a unit (unite) and a circuit of integers (ent), */
/* computes a critical point associated in the bnf k     */
GEN calcul_pt_pb(GEN k, GEN ent, GEN unite, GEN base_entiers)
{
  long l;
  GEN Omega;
  pari_sp ltop,lbot;
  ltop=avma;
  k = checkbnf(k);
  l = lg(ent)-1;
  Omega = gen_0;
  /*pari_printf("La base entiers est %Ps\n",base_entiers);*/
  {
    long i;
    for (i = 1; i <= l; i++){
      Omega = gadd(gmul(gel(ent, i), base_entiers), gmul(Omega, unite));
    }
  }
  /* to check Generalized Euclideanity */
  /* 
    if (est_principal(k, Omega, gsubgs(gpowgs(unite, l), 1)))
    pari_printf(" est principal\n");
    else
    pari_printf(" n'est pas principal.\n");
  */
  lbot=avma;
  return gerepile(ltop,lbot,gdiv(Omega, gsubgs(gpowgs(unite, l), 1)));
}


/* norm of (point minus ent) */
GEN calcul_norme(GEN point, GEN ent, GEN base_entiers_t, long prec){
  return gabs(gnorm(gsub(point, gmul(ent, base_entiers_t))), prec);
}


/* matrix of a Z-basis of the integers of l */
/* multiplied by ll (obtained through LLL)  */
GEN matrice_des_entiers(GEN l, GEN ll){
  long i,j,n;
  GEN a, z;
  pari_sp ltop,lbot;
  ltop=avma;
  l = checkbnf(l);
  n = glength(nf_get_pol(bnf_get_nf(l))) - 1;
  z = gmul(nf_get_zk(bnf_get_nf(l)), ll);
  lbot=avma;
  a= cgetg(n+1,t_MAT);
  for(i=1 ; i <=n ; i++){
    gel(a,i ) = cgetg(n+1,t_COL);
    for(j=1 ; j <= n ; j++){
      gcoeff(a,j,i) = polcoeff0(gel(z,i),j-1,-1);
    }
  }
  return gerepile(ltop,lbot,a);
}


/* converts a polynomial into a vector of length m */
/* (in fact a matrix with one line)                */
GEN polmod_vecteur(GEN p, long m){
  GEN v;
  long j;
  v=cgetg(2,t_MAT);
  gel(v,1) = cgetg(m+1,t_COL);
  for(j=1 ; j<= m ;j++){
    gcoeff(v,j,1) = polcoeff0(lift(p), j-1, -1);
  }
  return v;
}


/* converts a vector into a polynomial of degree < m */
/* in fact, the vector is a matrix with one line     */
GEN vecteur_pol(GEN v, long m)
{
  GEN q;
  long i;
  q = gen_0;
  for (i = 1; i <= m; ++i)
    q = gadd(q, gmul(gpowgs(pol_x(fetch_user_var("x")), i - 1), gcoeff(v, i, 1)));
  return q;
}


/* reduction of a point modulo integers */
GEN reduction_mod_entiers(GEN p, GEN po, GEN a, GEN b, long n){
  GEN v, v2;
  pari_sp ltop,lbot;
  long i;
  ltop=avma;
  v = polmod_vecteur(gmodulo(lift(p), po), n);
  v2 = gmul(b, v);
  for (i = 1; i <= n; ++i)
    gcoeff(v2, i, 1) = gfrac(gcoeff(v2, i, 1));
  v = gmul(a, v2);
  lbot=avma;
  return gerepile(ltop,lbot,vecteur_pol(v, n));
}


/* orbit of p (in the nf k) under the action of units */
/* modulo the integers of k                           */
GEN orbite3(GEN p, GEN k, GEN ll){
  GEN q, qdep,s,s1,lng,ae,be,u,po;
  long b,i,j,r,n,sdep;
  pari_sp ltop,lbot;
  ltop=avma;
  k = checkbnf(k);
  ae = matrice_des_entiers(k, ll);
  be = ginv(ae);
  n = glength(nf_get_pol(bnf_get_nf(k))) - 1;
  po = gcopy(nf_get_pol(bnf_get_nf(k)));
  q = reduction_mod_entiers(p, po, ae, be, n);
  qdep = q;
  s = gtoset(q);
  r = (nf_get_r1(bnf_get_nf(k)) + nf_get_r2(bnf_get_nf(k))) - 1;
  u = cgetg(r+2,t_VEC);
  for (i = 1; i <= r; i++)
    gel(u, i) = gcopy(gel(bnf_get_fu(k), i));
  gel(u, r + 1) = gcopy(gel(member_tu(k), 2));
  lng = cgetg(r+2,t_VEC);
  gel(lng, r + 1) = gsubgs(gel(member_tu(k), 1), 1);
  s1=gtoset(q); /* just in case r=0*/
  for(i=r ; i>0 ; i--){
    gel(lng,i)=gen_0;
    b=0;
    q = qdep;
    s1 = gtoset(qdep);
    while(b==0){
      q = reduction_mod_entiers(gmul(q, gel(u, i)), po, ae, be, n);
      b = setsearch(s1, q, 0);
      if (b == 0) {
	s1 = setunion(s1, gtoset(q));
	gel(lng, i) = gaddgs(gel(lng, i), 1);
      }
    }
    if (NIV_AFF >1){
      pari_printf("Lg[%ld]=%Ps\n", i, gel(lng, i));
    }
  }
  i = r + 1;
  s = setunion(s1, s);
  for(i=2 ; i <= r+1 ; i++){
    /*pari_printf("%Ps\n", s);*/
      sdep = (lg(s) - 1)-1;
      while ((lg(s)-1) > sdep){
	sdep = lg(s)-1;
	for(j=1 ; j < lg(s) ; j++){
	  q = reduction_mod_entiers(gmul(geval(gel(s, j)), gel(u, i)), po, ae, be, n);
	  /*pari_printf("%Ps\n", q);*/
	  s = setunion(s, gtoset(q));
	}
      }
  }
  lbot=avma;
  return gerepile(ltop,lbot,gcopy(s));
}

/* sharp bounds on the j-th coordinate assuming the previous */
/* ones are computed                                         */
void bornes_fines(long j, GEN borne, GEN xi2, GEN v, GEN gamma_k, GEN k, GEN m, GEN ll, long n, long prec,GEN **t){
  GEN borneg, borned, temp, b1, b2, xi3;
  long i,ii;
  pari_sp ltop;
  ltop=avma;
  k = checkbnf(k);
  b1 = gneg(gel(borne, j));
  b2 = gcopy(gel(borne, j));
  xi3 = cgetg(n+1,t_VEC);
  for(i=1 ; i <= n ; i++)
    gel(xi3,i) = gcopy(gcoeff(xi2,i,1));
  for(i=1 ; i <= n ; i++){
    borneg = gmul(gcoeff(m,i,j),gel(xi3,j));
    for(ii=1; ii<j; ii++)
	borneg = gadd(borneg,gmul(gcoeff(m,i,ii),gsub(gel(xi3,ii),gel(v,ii))));
    borned=borneg;
    borneg = gsub(borneg,gamma_k);
    borned = gadd(borned,gamma_k);
    for(ii=j+1 ; ii<=n ; ii++){
      if (gcmpgs(gcoeff(m, i, ii), 0) > 0){
	borneg = gadd(borneg,gmul( gcoeff(m,i,ii), gsub(gel(xi3,ii),gel(borne,ii))));
	borned = gadd(borned,gmul( gcoeff(m,i,ii), gadd(gel(xi3,ii),gel(borne,ii))));
      }
      else{
	borneg = gadd(borneg,gmul( gcoeff(m,i,ii), gadd(gel(xi3,ii),gel(borne,ii))));
	borned = gadd(borned,gmul( gcoeff(m,i,ii), gsub(gel(xi3,ii),gel(borne,ii))));
      }
    }
    if (gcmp(gabs(gcoeff(m, i, j), prec), strtor("0.0001", prec)) > 0){
      /*do nothing if the coordinate is small */
      if (gcmpgs(gcoeff(m, i, j), 0) > 0){
	borneg = gdiv(borneg, gcoeff(m, i, j));
	borned = gdiv(borned, gcoeff(m, i, j));
      }
      else{
	temp = borneg;
	borneg = gdiv(borned, gcoeff(m, i, j));
	borned = gdiv(temp, gcoeff(m, i, j));
      }
      b1 = gmax(borneg,b1);
      b2 = gmin(borned,b2); 
    }
  }
  *t = malloc(2*sizeof(GEN));
  (*t)[0] = gfloor(b1);
  (*t)[1] = gceil(b2);
  gerepileall(ltop,2,&((*t)[0]),&((*t)[1]));
}


/* tests if q is a nonzero non-unit algebraic integer */
/* should be useless with the small norms computed    */
long est_entier_non_nul_ni_unite(GEN q, long prec){
  if (gequal0(gmodgs(gtovec(minpoly(q, -1)), 1)))
    return (gcmpgs(gabs(gnorm(q), prec), 1) > 0);
  return 0;
}



/* treatment applied in the following recursive map */
void traitement_vecteur(GEN v, GEN* mini,GEN point, GEN base_entiers,long prec){
  GEN nor ;
  pari_sp ltop = avma;
  nor = calcul_norme(point,v,base_entiers,prec);
  if(gcmp(nor,*mini)<0){
    *mini = gerepilecopy(ltop,nor);
  } 
  else{
    avma=ltop;
  }
}

/* recursive map to compute the bounds   */
void parcours_rapide(int l, GEN vv,GEN borne, GEN xi2, GEN gamma_k, GEN k,GEN m, GEN ll, int n, GEN point, GEN base_entiers, GEN* mini, long prec){
  GEN *a;
  if(l > n){
    traitement_vecteur(vv,mini,point,base_entiers,prec);
  }
  else{
    bornes_fines(l,  borne, xi2, vv,  gamma_k, k, m, ll,  n, prec,&a);
    for(gel(vv,l) = a[0] ; gcmp(gel(vv,l), a[1])<=0 ; gel(vv,l)= gaddgs(gel(vv,l),1))
      parcours_rapide(l+1,vv,borne,xi2,gamma_k,k,m,ll,n,point,base_entiers,mini,prec);
  }
}




/* application of the recursive map  */
GEN calcul_minimum_un_point_orbite_dim_qcq(GEN point, GEN borne_g, GEN gam, GEN base_entiers, GEN k, GEN m, GEN ll, int n,long prec){
  GEN v,nor2,xi2,res;
  GEN* mini;
  pari_sp ltop,lbot;
  ltop=avma;
  mini=malloc(sizeof(GEN));
  k = checkbnf(k);
  v = gneg(borne_g);
  nor2 = calcul_norme(point, v, base_entiers, prec);
  *mini = nor2;
  xi2 = gmul(ginv(m), calcul_plongement(point, nf_get_r1(bnf_get_nf(k)), nf_get_r2(bnf_get_nf(k)), n, m, ll, prec));
  for( ; gcmp(gel(v, 1), gel(borne_g, 1)) <= 0 ;  gel(v,1)=gaddgs(gel(v,1),1))
    parcours_rapide(2,v,borne_g, xi2,  gam, k,m, ll, n, point, base_entiers, mini, prec);
  lbot=avma;
  res = gcopy(*mini);
  free(mini);
  return gerepile(ltop,lbot,res);
}


/* computes the minimum of a point in the orbit */
GEN calcul_minimum_un_point_orbite(GEN point, GEN borne_g, GEN gam, long dim, GEN base_entiers, GEN k, GEN m, GEN ll, long prec){
  k = checkbnf(k);
  if (est_entier_non_nul_ni_unite(ginv(point), prec)){
    if(NIV_AFF > 0){
      if(LANGUAGE){
	fprintf(stdout,"Inverse of an integer (non-zero, non-unit)");
      }
      else{
	fprintf(stdout,"Inverse d'un entier (ni nul, ni unité)");
      }
    }
    return(gabs(gnorm(point),prec));
  }
  return calcul_minimum_un_point_orbite_dim_qcq(point, borne_g, gam, base_entiers, k, m, ll, dim,prec);
}



/* general function to compute the minimum */
GEN calcul_minimum(GEN point, GEN k, GEN borne_g, GEN gam, long dim, GEN base_entiers, GEN m, GEN ll, long prec){
  long l;
  pari_sp ltop,lbot;
  GEN s,v,nor, nor2,p1;
  ltop=avma;
  k = checkbnf(k);
  s = orbite3(point, k, ll);
  if (NIV_AFF >0){
    if (LANGUAGE){
      pari_printf("bounds: ");
    }
    else{
      pari_printf("bornes: ");
    }
    pari_printf("%Ps\n", borne_g);
  }
  /*
  pari_printf("orbite :");
  pari_printf("%Ps\n", s);
  */
  l = lg(s)-1;
  /*
  pari_printf("(longueur ");
  pari_printf("%ld", l);
  pari_printf(")\n");
  */
  v = calcul_rapide(k, gtos(powis(gen_2, 2+degree(nf_get_pol(bnf_get_nf(k))))), s);
  if (!gequal0(gel(v, 1)))
  {
    lbot=avma;
    p1 = cgetg(4, t_VEC);
    gel(p1, 1) = gcopy(s);
    gel(p1, 2) = stoi(l);
    gel(p1, 3) = gcopy(gel(v, 2));
    return gerepile(ltop,lbot,p1);
  }
  else
  {
    nor = calcul_minimum_un_point_orbite(gmodulo(geval(gel(s, 1)), nf_get_pol(bnf_get_nf(k))), borne_g, gam, dim, base_entiers, k, m, ll, prec);
    {
      long i;
      for (i = 2; i <= l; ++i)
	{
	  nor2 = calcul_minimum_un_point_orbite(gmodulo(geval(gel(s, i)), nf_get_pol(bnf_get_nf(k))), borne_g, gam, dim, base_entiers, k, m, ll, prec);
	  if (gcmpgs(gsub(nor2, nor), 0) < 0)
	    nor = nor2;
	}
    }
    lbot=avma;
    p1 = cgetg(4, t_VEC);
    gel(p1, 1) = gcopy(s);
    gel(p1, 2) = stoi(l);
    gel(p1, 3) = gcopy(nor);
    return gerepile(ltop,lbot,p1);
  }
}


/* minimum of a circuit */
GEN calcul_min_cycle(GEN ent, GEN unite, GEN borne_g, GEN gam, long dim, GEN base_entiers, GEN k, GEN m, GEN ll, GEN pts_deja_vus, long prec)
{
  GEN point,res;
  long i;
  pari_sp ltop,lbot;
  ltop=avma;
  k = checkbnf(k);
  point = calcul_pt_pb(k, ent, unite, base_entiers);
  if(NIV_AFF > 0)
    pari_printf("%Ps\n", point);
  i = setsearch(pts_deja_vus, lift(point), 0);
  lbot=avma;
  if (i > 0){
    res = cgetg(3, t_VEC);
    gel(res, 1) = strtoGENstr("déjà vu");
    gel(res, 2) = stoi(i);
    return gerepile(ltop,lbot,res);
  }
  return gerepile(ltop,lbot,calcul_minimum(point, k, borne_g, gam, dim, base_entiers, m, ll, prec));
}


/* minimum of a point                               */
/* not used in general but useful to check examples */
GEN calcul_min_point(GEN p, GEN k, GEN mini, long prec){
  GEN a, base_ent,gam,borne_g;
  pari_sp ltop,lbot;
  ltop=avma;
  k = checkbnf(k);
  a = matrice_plongement(k);
  base_ent = gtrans(gmul(nf_get_zk(bnf_get_nf(k)), gel(a, 2)));
  gam = calcul_gamma(k, mini, gel(a, 1), gel(a, 2), prec);
  borne_g = calcul_bornes_grossieres(gam, gel(a, 1), prec);
  lbot=avma;
  return gerepile(ltop,lbot,calcul_minimum(p, k, borne_g, gam, glength(nf_get_pol(bnf_get_nf(k))) - 1, base_ent, gel(a, 1), gel(a, 2), prec));
}


/* non-trivial ideal of smallest norm */
GEN ideal_plus_petit(GEN k, long petite_n){
  long i;
  GEN l, r,res;
  pari_sp ltop,lbot;
  ltop = avma;
  k = checkbnf(k);
  l = ideallist0(k, gtos(powis(gen_2, nf_get_r1(bnf_get_nf(k)) + (2*nf_get_r2(bnf_get_nf(k))))), 4);
  for(i=2 ; glength(gel(l, i)) == 0 ; i++) ;
  r = gel(gel(l, i), 1);
  if (petite_n == 0){
    lbot=avma;
    res = cgetg(3, t_VEC);
    gel(res, 1) = stoi(i);
    gel(res, 2) = matid(glength(nf_get_pol(bnf_get_nf(k))) - 1);
    return gerepile(ltop,lbot,res);
  }
  lbot=avma;
  res = cgetg(3, t_VEC);
  gel(res, 1) = stoi(i);
  gel(res, 2) = gcopy(r);
  return gerepile(ltop,lbot,res);
}


/* computation of the minimum of a list of circuits */
GEN calcul_min_liste_cycles(GEN cyc, GEN unite, GEN k, GEN mini, long petite_n, long prec){
  GEN a,id,base_ent,v,w,gam,m,pts_deja_vus,borne_g;
  pari_sp ltop, lbot;
  long i;
  ltop=avma;
  a=matrice_plongement(k);
  if (NIV_AFF >0){
    if (LANGUAGE){
      pari_printf("We use the unit %Ps.\n",unite);
    }
    else{
      pari_printf("On raisonne avec l'unité %Ps.\n",unite);
    }
  }
  id = ideal_plus_petit(k, petite_n);
  if(petite_n==0){
    gel(id,1) = gen_1;
    /*fprintf(stdout,"On raisonne avec les entiers.\n");*/
  }
  else{
    if (NIV_AFF >0){
      if (LANGUAGE){
	pari_printf("We use the ideal %Ps of minimal norm %Ps.\n",gel(id,2),gel(id,1));
      }
      else{
	pari_printf("On utilise l'idéal %Ps de norme minimale %Ps.\n",gel(id,2),gel(id,1));
      }
    }
  }
  gel(a, 1) = gmul(gmul(gmul(gel(a, 1), ginv(gel(a, 2))), gel(id, 2)), gel(a, 2));
  gel(a, 2) = gmul(gel(id, 2), gel(a, 2));
  base_ent = gtrans(gmul(nf_get_zk(bnf_get_nf(k)), gel(a, 2)));
  gam = calcul_gamma(k, gmul(mini, gel(id, 1)), gel(a, 1), gel(a, 2), prec);
  borne_g = calcul_bornes_grossieres(gam, gel(a, 1), prec);
  v=cgetg(lg(cyc),t_VEC);
  pts_deja_vus = cgetg(1,t_VEC);
  m=gen_0;
  for(i=1 ; i < lg(cyc) ; i++){
    if (NIV_AFF >0)
      fprintf(stdout,"circuit %ld\n",i);
    gel(v,i) = calcul_min_cycle(gel(cyc, i), unite, borne_g, gam, nf_get_r1(bnf_get_nf(k)) + (2*nf_get_r2(bnf_get_nf(k))), base_ent, k, gel(a, 1), gel(a, 2), pts_deja_vus, prec);
    if (NIV_AFF >0)
      pari_printf("%Ps\n",gel(v,i));
    if (glength(gel(v, i)) > 2){
        gel(gel(v, i), 3) = gdiv(gel(gel(v, i), 3), gel(id, 1));
	m = gmax(m, gel(gel(v, i), 3));
        pts_deja_vus = setunion(pts_deja_vus, gel(gel(v, i), 1));
    }
  }
  w = gsub(gmul(mini, denom(m)), numer(m));
  if(gcmpgs(w,0) < 0){
    if (NIV_AFF >0){
      if (LANGUAGE){
	pari_printf("New computation (m=%Ps, mini =%Ps).\n",m,mini);
      }
      else{
	pari_printf("Nouveau calcul (m=%Ps, mini =%Ps).\n",m,mini);
      }
    }
    lbot = avma;
    return gerepile(ltop,lbot,calcul_min_liste_cycles(cyc, unite, k, gdiv(gadd(numer(m), strtor("0.0001", prec)), denom(m)), petite_n, prec));
  }
  if (NIV_AFF >0){
    pari_printf("m=%Ps <= mini = %Ps\n",m,mini);
  }
  lbot= avma;
  return gerepile(ltop,lbot,gcopy(v));
}


/* computation of data of the nf given by the polynomial P */
/* petite_n tells if we consider ideal classes             */
/*                                                         */
/* note: stupid function which returns r1,r2,n ...         */ 
/* as GEN type objects; with care, we should do with       */
/* integers directly.                                      */
GEN calculs_corps_de_nombres(GEN P, long petite_n, long prec){
  GEN k,a,id,res,unites,unit,rac1,rac_un,rac_un2,be,me;
  pari_sp ltop,lbot;
  long r1,r2,n,r,j,i,ordre_rac_un;
  ltop=avma;
  k=Buchall(P,nf_FORCE,prec);
  a=matrice_plongement(k);
  id=ideal_plus_petit(k,petite_n);
  if (NIV_AFF >1){
    if (LANGUAGE){
      pari_printf("Smallest norm used: %Ps, of matrix %Ps\n", gel(id, 1), gel(id, 2));
    }
    else{
      pari_printf("Plus petite norme utilisée : c'est %Ps, de matrice %Ps\n", gel(id, 1), gel(id, 2));
    }
  }
  gel(a, 1) = gmul(gmul(gmul(gel(a, 1), ginv(gel(a, 2))), gel(id, 2)), gel(a, 2));
  gel(a, 2) = gmul(gel(id, 2), gel(a, 2));
  r1 = nf_get_r1(bnf_get_nf(k));
  r2 = nf_get_r2(bnf_get_nf(k));
  n = r1 + (2*r2);
  r = (r1 + r2) - 1;
  if(r1+r2 < 2){
    lbot = avma;
    res = cgetg(5,t_VEC);
    gel(res,1) = stoi(n);
    gel(res,2) = stoi(r1);
    gel(res,3) = stoi(r2);
    gel(res,4) = P;
    return gerepile(ltop,lbot,res);
  }
  unit = conjvec(gmodulo(lift(bnf_get_fu(k)), P),prec);
  unites = cgetg(r+1, t_MAT);
  for(j=1 ; j<= r ; j++){
    gel(unites,j) = cgetg(n+1,t_COL);
    for(i=1 ; i<= r1 ; i++)
      gcoeff(unites, i, j) = gcopy(gcoeff(unit, i, j));
    for(i=r1+1 ; i < n ; i=i+2){
      gcoeff(unites,r1+(i-r1-1)/2+1,j) = greal(gcoeff(unit,i,j));
      gcoeff(unites,r1+(i-r1-1)/2+1+r2,j) = gimag(gcoeff(unit,i+1,j));
    }
  }
  rac1 = member_tu(k);
  ordre_rac_un = gtos(gel(rac1, 1));
  rac_un2 = conjvec(gmodulo(lift(gel(rac1, 2)), P), prec);
  rac_un = cgetg(n+1, t_VEC);
  for(j=1 ; j<=r1 ; j++)
    gel(rac_un,j) = gel(rac_un2,j);
  for(j=r1+1 ; j<n ; j++){
    gel(rac_un,r1+(j-r1-1)/2+1) = greal(gel(rac_un2,j));
    gel(rac_un,r1+(j-r1-1)/2+1+r2) = gimag(gel(rac_un2,j+1));
  }
  be = lift(nf_get_zk(bnf_get_nf(k)));
  me = cgetg(n+1,t_MAT);
  for(i=1 ; i<= n ; i++)
    gel(me,i) = cgetg(n+1,t_COL);
  for(i=1 ; i<= n ; i++){
    for(j=1 ; j<= n ;j++){
      gcoeff(me,i,j) = polcoeff0(gel(be,i),j-1,-1);
    }
  }
  lbot=avma;
  res = cgetg(14,t_VEC);
  gel(res,1) = stoi(n);
  gel(res,2) = stoi(r1);
  gel(res,3) = stoi(r2);
  gel(res,4) = gcopy(gel(a,1));
  gel(res,5) = ginv(gel(a,1));
  gel(res,6) = gtrans(unites);
  gel(res,7) = gcopy(rac_un);
  gel(res, 8) = stoi(ordre_rac_un);
  gel(res, 9) = gcopy(gel(id, 1));
  gel(res, 10) = icopy(gel(bnf_get_clgp(k), 1));
  gel(res, 11) = gmul(me, gel(a, 2));
  gel(res, 12) = gtovec(P);
  gel(res, 13) = gcopy(k);
  return gerepile(ltop,lbot,res);
}



/* squarefree part of n  */
GEN sqfp(GEN n){
  GEN v,res;
  long i,l;
  pari_sp ltop,lbot;
  ltop = avma;
  v=Z_factor(n);
  res = gen_1;
  l=lg(gtrans(v));
  for(i=1 ; i < l ; i++){
    if (gequal1(gmodgs(gcoeff(v, i, 2), 2)))
      res = gmul(res,gcoeff(v,i,1));
  }
  lbot = avma;
  return gerepile(ltop,lbot,gcopy(res));  
}


/* Euclidean minimum of an imaginary quadratic nf given by the polynomial Q */
GEN minimum_quad_im(GEN Q){
  GEN a,b,m,m2,res;
  pari_sp ltop, lbot;
  ltop = avma;
  a=polcoeff0(Q,1,-1);
  b=polcoeff0(Q,0,-1);
  m=gsub(gmulgs(b,4),gsqr(a));
  m2=sqfp(m);
  res = gen_0;
  if(itos(gmodgs(m2,4)) == 3)
    res = gdiv(gsqr(gaddgs(m2,1)),gmulgs(b,4));
  else
    res = gdivgs(gaddgs(m2,1),4);
  lbot=avma;
  return gerepile(ltop,lbot,gcopy(res));
}


/* speed-up: often, the minimum is 1/smthg, where   */
/* smthg is the absolute value of some norm;        */
/* -> computations of inverse of small norms to     */
/* recognize them                                   */
GEN liste_inv_entiers_pas_unites(GEN k, long n){
  GEN l,u,v;
  long i,j;
  pari_sp ltop,lbot;
  ltop=avma;
  k=checkbnf(k);
  l=cgetg(n+1,t_VEC);
  gel(l,1)=cgetg(1,t_VEC);
  fflush(stdout);
  for(i=2 ; i<= n ; i++){
    u=bnfisintnorm(k,stoi(i));
    v=bnfisintnorm(k,stoi(-i));
    gel(l,i)=cgetg(glength(u)+glength(v)+1,t_VEC);
    for(j=1 ; j <= glength(u) ; j++)
      gel(gel(l,i),j) = ginv(gel(u,j));
    for(j=1 ; j<= glength(v) ; j++)
      gel(gel(l,i),j+glength(u)) = ginv(gel(v,j));
  }
  lbot=avma;
  return gerepile(ltop,lbot,gcopy(l));
}


/* given s, tells if s is the inverse of a non-unit element */
/* of k whose absolute value of norm is at most n           */
GEN calcul_rapide(GEN k, long n, GEN s){
  GEN l,res;
  pari_sp ltop,lbot;
  long i,ii,j;
  ltop=avma;
  k=checkbnf(k);
  res=cgetg(3,t_VEC);
  gel(res,1)=gen_0;
  gel(res,2)=gen_0;
  l = liste_inv_entiers_pas_unites(k, n);
  for(i=1 ; i<= n ;i++){
    for(ii=1 ; ii<= glength(gel(l,i)) ; ii++){
      for(j=1 ; j<= glength(s); j++){
	if (gequal0(gfrac(gtovec(charpoly0(gmodulo(lift(gsub(geval(gel(s, j)), gel(gel(l, i), ii))), nf_get_pol(bnf_get_nf(k))), -1, 3))))){
	  if (NIV_AFF >0){
	    if (LANGUAGE){
	      fprintf(stdout,"Fast computation successful.\n");
	    }
	    else{
	      fprintf(stdout,"Calcul rapide réussi.\n");
	    }
	  }
	  gel(res,1)=gen_1;
	  gel(res,2)=ginv(stoi(i));
	  lbot=avma;
	  return gerepile(ltop,lbot,gcopy(res));
	}
      }
    }
  }
  lbot=avma;
  return gerepile(ltop,lbot,gcopy(res));
}


/* basic function to compute data on k defined by P*/
#ifdef SQLITE
GEN calcul_pari2(GEN P, long prec, numberfield* k , int c_i, double *K2, double *K4,sqlite3 *db){
  char *poly; 
#else /*SQLITE*/
GEN calcul_pari2(GEN P, long prec, numberfield* k , int c_i, double *K2, double *K4){
#endif /*SQLITE*/
  GEN calc;
  int dim,r1,r2,ordre_rac_unit,norme_ideal,h_k;
  double **sigma2,**sigma2_inv,**unit2 ;
  double * gen_rac_unit2;
  calc  = calculs_corps_de_nombres(P,c_i,prec);
  /*pariprintf("%Ps\n",calc);*/
#ifdef SQLITE  
  /*test pour voir si le polynôme est connu */
  poly= GENtostr(polredabs0(P,0));
  recherche_pol(db, poly);
#endif /*SQLITE*/
  /* suite du programme */
  dim = itos(gel(calc,1));
  r1 = itos(gel(calc,2));
  r2 = itos(gel(calc,3));
  /*fprintf(stdout,"n=%d,r1=%d,r2=%d\n",*dim,*r1,*r2);*/
  if(r1+r2>1){
    sigma2 = allocMat(dim,dim);
    sigma2_inv = allocMat(dim,dim);
    unit2 = allocMat(r1 + r2 -1, dim);
    gen_rac_unit2 = allocVec(dim);
    GENtoMat(&sigma2,gel(calc,4));
    GENtoMat(&sigma2_inv,gel(calc,5));
    GENtoMat(&unit2,gel(calc,6));
    affGENtoVec(gen_rac_unit2,gel(calc,7));
    ordre_rac_unit = itos(gel(calc,8));
    norme_ideal =itos( gel(calc,9));
    h_k = itos(gel(calc,10));
    *K2 = 1. / (norme_ideal);
    if(!c_i){
      norme_ideal = 1;
    }
    /*pariprintf("%Ps\n",gel(calc,10));*/
    init_numberfield(k,dim,r1,r2,sigma2,sigma2_inv,unit2,ordre_rac_unit,gen_rac_unit2,norme_ideal,h_k);
  }
  else{
    k->r1 = r1;
    k->r2 = r2;
    k->dim = dim;
  }
  return calc;
}

/*  Launch of the pari library         */
#ifdef SQLITE
GEN lancement_pari(int argc, char **argv, numberfield *k, int c_i, double *K2, double *K4,sqlite3 *db){
#else /*SQLITE*/
GEN lancement_pari(int argc, char **argv, numberfield *k, int c_i, double *K2, double *K4){
#endif /*SQLITE*/
  GEN P;
  pari_init(PARI_STACK_SIZE, 500000);
  if (argc==1){
    if (LANGUAGE){
      fprintf(stdout,"Please enter a polynomial to define the number field: ");
    }
    else{
      fprintf(stdout,"Entrez un polynôme définissant le corps de nombres : ");
    }
    P = gp_read_stream(stdin);
  }
  else{
    P = gp_read_str(argv[1]);
    if (typ(P) == t_VEC) P = gtopoly(P, 0);
  }
  if(typ(P) != t_POL){
    if(LANGUAGE){
      fprintf(stdout,"Input is not a polynomial.\n");
      exit(-1);
    }
    else{
      fprintf(stdout,"L'entrée n'est pas un polynôme.\n");
      exit(-1);
    }
  }
#ifdef SQLITE
  return calcul_pari2(P,DEFAULTPREC, k, c_i, K2,K4,db);
#else /*SQLITE*/
  return calcul_pari2(P,DEFAULTPREC, k, c_i, K2,K4);
#endif /*SQLITE*/
  /*pari_close();*/
}

