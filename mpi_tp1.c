#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"


int main(int argc, char** argv) 
{
	double tiempoTotalInicia;
	double tiempoTotalFinal;
    srand(time(NULL)); 
    int n;				/*N  */
	int  p;            /*  Numero de procesos que corren */
	int  proc_id; 
    int i;
	int j;
	int y;
	
	double tiempoEjecucionInicia;
	double tiempoEjecucionFinal;

	
	int *M;
    int *MLocal;
    int *V;
    int *B;// 
	int *BLocal;
	int  *sendcounts;  /* nuevo para Scatterv  para 10 procesos */
    int  *displs;      /* nuevo para Scatterv*/
	int *Q;				//Vector resultado de M*V
	int *MyQ;//
	int *MyP;
	int *P;
	int tp;
	int MyTp;

	tiempoTotalInicia = MPI_Wtime(); 
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &p); 
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

   if (proc_id==0)
   {
	  //Generacion de Matriz M y Vector V
      printf("Ingrese la dimension de la matriz: \n");
      scanf("%d", &n);
	  tiempoEjecucionInicia = MPI_Wtime(); 
      printf("n es %d\n",n);
      sendcounts=malloc(p*sizeof(int));
	displs=malloc(p*sizeof(int));
	B= malloc(n * n * sizeof(int));
      M=malloc(n * n * sizeof(int));
      V = malloc(n * sizeof(int));
	  Q= malloc(n *  sizeof(int));
	  P= malloc(n *  sizeof(int));
	  tp=0;
	
       for (i = 0; i < p; i++)
		 {
			 if (i==0 || i==(p-1))
			 {
				 if(p==1)
				 {
					 sendcounts[i]=n*n;
				 }	
				 else
				 {
				 sendcounts[i]=((n*n)/p)+n;
				 }
			 }
			 else
			 {
				 sendcounts[i]=((n*n)/p)+ 2*n; 
			 }
			if (i==0)
			 {
				displs[0]=0;
			 }
			 else
			 {
				 displs[i]=((n/p)*n*i)-n;
			 }			 
		 }

      //Revisar que se haya asignado memoria
      if (M == NULL || V == NULL){
        printf("Problemas reservando memoria\n");
      }

      for (i=0;i<n;i++)
		 {
			 for (j=0;j<n;j++)
			{
				//M[i * ncolumns + j]=rand() % 10; 
				 M[i * n + j]=rand() % 10; //random del 0 al 9
			}
			V[i] = rand()%6;
		 }
   }

   MPI_Barrier(MPI_COMM_WORLD);
    
   MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);


    MLocal= malloc(n * n * sizeof(int));/*La parte de la matriz M que recibe cada proceso*/
   	BLocal= malloc(n * n * sizeof(int));
     MyQ=malloc(n * sizeof(int));
	 MyP=malloc(n * sizeof(int));
	MyTp =0;

	 if (proc_id > 0){
		 V= malloc(n *  sizeof(int));	
    }

   //MPI_Scatter(M,(n*n)/p, MPI_INT,MLocal,(n*n)/p, MPI_INT,0,MPI_COMM_WORLD);
    MPI_Scatterv(M, sendcounts,displs,MPI_INT,MLocal, (n*n),MPI_INT ,0,MPI_COMM_WORLD);
	MPI_Bcast(V,n, MPI_INT,0,MPI_COMM_WORLD);

     int k=0;
	 int h=(n/p);
    if (proc_id > 0){
        k=1;
		h=h+1;
    }

//Todos los proceso calculan MyQ[i]
for (i=k;i<h;i++)
    {
	    for (j=0;j<n;j++)
	    {
			//printf("estoy en el proceso %d, %d * %d = %d \n",my_rank,MLocal[i*n+j],V[j],(MLocal[i*n+j]*V[j]));
			//printf("estoy en proceso %d, i vale %d y j vale %d \n", my_rank, i, j);
			MyQ[i]+=(MLocal[i*n+j]*V[j]);
		}
	   
    }
	
	//Esto que es?

	if(proc_id>0)
{
	for(i=0;i<n;i++)
	{
		MyQ[i]=MyQ[i+1];
	}
}
MPI_Gather(MyQ, (n/p), MPI_INT, Q, (n/p), MPI_INT, 0, MPI_COMM_WORLD);

//Calculo de P[i]

for (i=0;i<n;i++){
        MyP[i] = 0;
      }

for (i=k;i<h;i++) //filas
    {
	    for (j=0;j<n;j++) //columnas
	    {	
            if (MLocal[i*n+j]==2 || MLocal[i*n+j]==3  || MLocal[i*n+j]==5  || MLocal[i*n+j]==7){
                MyP[j]++;
				MyTp++;
			}
		}
    }

MPI_Reduce(MyP,P,((n*n)/p),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Reduce(&MyTp,&tp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

///Calculo cruz de Matriz B
	for (i=k;i<h;i++)
    {
	    for (j=0;j<n;j++)
	    {	
			BLocal[i*n+j]=0;
			BLocal[i*n+j]+=MLocal[i*n+j];
			if (i>=0){
				BLocal[i*n+j]+=MLocal[(i-1)*n+j];
			}
			if(j>=0){
				BLocal[i*n+j]+=MLocal[i*n+(j-1)];
			}
			if (i<n){
				BLocal[i*n+j]+=MLocal[(i+1)*n+j];
			}
			if(j<n){
				BLocal[i*n+j]+=MLocal[i*n+(j+1)];
			}
	    }
    }

	if(proc_id>0)
		{
		for(i=0;i<(n*n)/p;i++)
			{
		BLocal[i]=BLocal[i+n];
		}		
	}

	MPI_Gather(BLocal, ((n*n)/p), MPI_INT, B, ((n*n)/p), MPI_INT, 0, MPI_COMM_WORLD);

 MPI_Barrier(MPI_COMM_WORLD);


//TODO
   // free toda la memoria asignada
   // tiempo
 MPI_Finalize(); 

if (proc_id == 0){
	if (n < 100){
		printf("La Matriz M es: \n");
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					//printf("%d",M[i * ncolumns + j]); 
					printf(" %d ",M[i * n + j]); 
				}
				printf("\n");
			}
		printf("\n");

		printf("El Vector V es: \n");
		for (i=0;i<n;i++){
			printf(" %d ", V[i]);
		}
		printf("\n\n");

		printf("El resultado de Q es: \n");
		for (i=0;i<n;i++){
			printf(" %d ",Q[i]);
		}
		printf("\n\n");

		printf("El resultado de P es: \n");
		for (i=0;i<n;i++){
			printf(" %d ",P[i]);
		}
		printf("\n\n");
		printf("TP es = %d",tp);
		printf("\n\n");
		printf("El resultado de B: \n");

			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					printf(" %d ",B[i * n + j]);
				}
				printf("\n");
			}
	}
	else{ //n > 100 guardar en un archivo
		FILE *fp;
	fp = fopen("/resultados.txt", "w+");
	fprintf(fp,"La Matriz M es: \n");
			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					//printf("%d",M[i * ncolumns + j]); 
					fprintf(fp," %d ",M[i * n + j]); 
				}
				fprintf(fp,"\n");
			}
		fprintf(fp,"\n");

		fprintf(fp,"El Vector V es: \n");
		for (i=0;i<n;i++){
			fprintf(fp," %d ", V[i]);
		}
		fprintf(fp,"\n\n");

		fprintf(fp,"El resultado de Q es: \n");
		for (i=0;i<n;i++){
			fprintf(fp," %d ",Q[i]);
		}
		fprintf(fp,"\n\n");

		fprintf(fp,"El resultado de P es: \n");
		for (i=0;i<n;i++){
			fprintf(fp," %d ",P[i]);
		}
		fprintf(fp,"\n\n");
		fprintf(fp,"TP es = %d",tp);
		fprintf(fp,"\n\n");
		fprintf(fp,"El resultado de B: \n");

			for (i=0;i<n;i++)
			{
				for (j=0;j<n;j++)
				{
					fprintf(fp," %d ",B[i * n + j]);
				}
				fprintf(fp,"\n");
			}
	fclose(fp);
	}
	tiempoEjecucionFinal = MPI_Wtime(); 
tiempoTotalFinal = MPI_Wtime(); 
printf("El tiempo total que tardo el programa es de %f \n", tiempoTotalInicia-tiempoTotalFinal);
printf("El tiempo que tardo el programa desde que ingreso los valores es de %f \n", tiempoEjecucionInicia-tiempoEjecucionFinal);
}

   /* Se termina el ambiente MPI */

  
    return (EXIT_SUCCESS);
}


