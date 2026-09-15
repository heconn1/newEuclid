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

/* This files contains basic functions to work with sqlite3 */

#include "header.h"



char* string_numberfield(int n, int r1, int r2, char *DISC, char *DISC2, int clnb, char *POL, char *M1, int t1, char * C1, float approx, int euclidean){
  char N[30];
  char R1[30];
  char R2[30];
  char CLNB[30];
  char T1[30];
  char APPROX[30];
  char EUCLIDEAN[30];
  char *result=malloc(20000*sizeof(char));
  snprintf(N, 29, "%d", n);
  snprintf(R1, 29, "%d", r1);
  snprintf(R2, 29, "%d", r2);
  snprintf(CLNB, 29, "%d", clnb);
  snprintf(T1, 29, "%d", t1);
  snprintf(APPROX, 29, "%f", approx);
  snprintf(EUCLIDEAN,29, "%d",euclidean);
  strcpy(result,"insert or replace into "); 
  if(PRINCIPAL<2){ 
    strcat(result,"euc");
  }
  else{
    strcat(result,"euc_class");
  }
  strcat(result," (deg,r1,r2,disc,app_disc,clnb,pol,M1,T1,C1,approxim,euclidean) values (");
  strcat(result,N);
  strcat(result,",");
  strcat(result,R1);
  strcat(result,",");
  strcat(result,R2);
  strcat(result,",\"");
  strcat(result,DISC);
  strcat(result,"\",");
  strcat(result,DISC2);
  strcat(result,",");
  strcat(result,CLNB);
  strcat(result,",\"");
  strcat(result,POL);
  strcat(result,"\",\"");
  strcat(result,M1);
  strcat(result,"\",");
  strcat(result,T1);
  strcat(result,",\"");
  strcat(result,C1);
  strcat(result,"\",");
  strcat(result,APPROX);
  strcat(result,",");
  strcat(result,EUCLIDEAN);
  strcat(result,");");
  return result;
}

char* string_numberfield2(int n, int r1, int r2, char *DISC, char *DISC2, int clnb, char *POL, int euclidean){
  /* quand on obtient seulement l'euclidianité*/
  char N[30];
  char R1[30];
  char R2[30];
  char CLNB[30];
  char EUCLIDEAN[30];
  char *result=malloc(10000*sizeof(char));
  snprintf(N, 29, "%d", n);
  snprintf(R1, 29, "%d", r1);
  snprintf(R2, 29, "%d", r2);
  snprintf(CLNB, 29, "%d", clnb);
  snprintf(EUCLIDEAN,29, "%d",euclidean);
  strcpy(result,"insert into ");
  if(PRINCIPAL<2){
    strcat(result,"euc");
  }
  else{
    strcat(result,"euc_class");
  }
  strcat(result," (deg,r1,r2,disc,app_disc,clnb,pol,euclidean) values (");
  strcat(result,N);
  strcat(result,",");
  strcat(result,R1);
  strcat(result,",");
  strcat(result,R2);
  strcat(result,",\"");
  strcat(result,DISC);
  strcat(result,"\",");
  strcat(result,DISC2);
  strcat(result,",");
  strcat(result,CLNB);
  strcat(result,",\"");
  strcat(result,POL);
  strcat(result,"\",");
  strcat(result,EUCLIDEAN);
  strcat(result,");");
  return result;
}



static int callback(void *NotUsed, int argc, char **argv, char **azColName){
  int i;
  NotUsed=0;
  if(NIV_AFF>0){
    for(i=0; i<argc; i++){
      printf("%d. %s = %s\n", i, azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
  }
  if(argv[7]){
    if(strlen(argv[7]) >0 ){
      if(LANGUAGE){
	fprintf(stdout,"This number field is already known.\n");
	fprintf(stdout,"If we define it by the polynomial %s, ",argv[6]);
	fprintf(stdout,"its Euclidean minimum is %s, reached at %s critical point(s):%s.\n",argv[7],argv[8],argv[9]);
      }
      else{
	fprintf(stdout,"Corps déjà connu.\n");
	fprintf(stdout,"Si on le définit par le polynôme %s, ",argv[6]);
	fprintf(stdout,"son minimum euclidien est %s, atteint en %s  point(s) critique(s) :%s.\n",argv[7],argv[8],argv[9]);
      }
      exit(1);
    }
  }
  if(argv[11]){
    if(strlen(argv[11]) > 0){
      if(MINIMAL_VALUE_K == 1){
	if(LANGUAGE){
	  fprintf(stdout,"This field is known to be norm-Euclidean.\n");
	}
	else{
	  fprintf(stdout,"On connaît déjà l'euclidianité pour la norme.\n");
	}
	exit(1);
      }
    }
  }
  return 0;
}

void execution(sqlite3 *db, char *st){
  char *zErrMsg = 0;
  if(sqlite3_exec(db,st,callback,0,&zErrMsg)!= SQLITE_OK){
    fprintf(stderr, "SQL error : %s\n",zErrMsg);
  }
}

void creation_table_numberfield(sqlite3 *db){
  char *s = malloc(10000*sizeof(char));
  strcpy(s, "create table if not exists ");
  if(PRINCIPAL<2){
    strcat(s,"euc");
  }
  else{
    strcat(s,"euc_class");
  }
  strcat(s,"(deg smallint, r1 smallint, r2 smallint, disc text, app_disc int,clnb smallint, pol text primary key, M1 text, T1 text, C1 text, approxim real, euclidean smallint);");
  execution(db, s);
} 

void ajout_numberfield(sqlite3 *db, int n, int r1, int r2, char *disc, char *disc2, int clnb, char *pol, char *M1, int T1, char * C1, double approx, long euclidean){
  execution(db,string_numberfield(n,r1,r2,disc,disc2,clnb,pol,M1,T1,C1,approx,euclidean));
}


void ajout_numberfield2(sqlite3 *db, int n, int r1, int r2, char *disc, char *disc2, int clnb, char *pol, long euclidean){
  execution(db,string_numberfield2(n,r1,r2,disc,disc2,clnb,pol,euclidean));
}


void recherche_pol(sqlite3 *db, char *pol){
  char *s;
  s = malloc(400);
  strcpy(s,"select * from ");
  if(NIV_AFF>0){
    if (LANGUAGE){
      fprintf(stdout,"Looking up in the database.\n");
    }
    else{
      fprintf(stdout,"Recherche dans la base de données.\n");
    }
  }
  if(PRINCIPAL < 2){
    strcat(s,"euc");
  }
  else{
    strcat(s,"euc_class");
  }
  strcat(s," where pol=\"");
  strcat(s,pol);
  strcat(s,"\"");
  if(NIV_AFF >0){
    if(LANGUAGE){
      fprintf(stdout,"sqlite request: %s\n",s);
    }
    else{
      fprintf(stdout,"Requête sqlite : %s\n",s);
    }
  }
  execution(db,s);
  if(NIV_AFF >0){
    if(LANGUAGE){
      fprintf(stdout,"unknown field or incomplete data, running the program.\n");
    }
    else{
      fprintf(stdout,"corps inconnu ou données incomplètes, retour au programme\n");
    }
  }
}

