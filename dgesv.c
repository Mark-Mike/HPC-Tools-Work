#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <openblas/lapacke.h>
//#include <mkl_lapacke.h>

double *generate_matrix(int size)
{
  int i;
  double *matrix = (double *) malloc(sizeof(double) * size * size);

  srand(1);

  for (i = 0; i < size * size; i++) {
    matrix[i] = rand() % 100;
  }

  return matrix;
}

int is_nearly_equal(double x, double y)
{
  const double epsilon = 1e-5 /* some small number */;
  return abs(x - y) <= epsilon * abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

int check_result(double *bref, double *b, int size)
{
  int i;

  for(i = 0; i < size*size; i++) {
    if (!is_nearly_equal(bref[i], b[i]))
      return 0;
  }

  return 1;
}

/*int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{

  //Replace next line to use your own DGESV implementation
  LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);

}*/

/*int my_dgesv(int n, double *a, int *ipiv, double *b)
{
  printf("Hi, the size is : %d\n", n);
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      //printf("pow(*(a + i*n + j)-*(b + i*n + j),2) = %d\n", (int)sqrt(pow(*(a + i*n + j)-*(b + i*n + j),2)));
      *(ipiv + i*n + j) = (int)sqrt(pow(*(a + i*n + j)-*(b + i*n + j),2));
      //printf("ipiv[%d][%d] = %d\n", i, j, *(ipiv + i*n + j));
      if (*(ipiv + i*n + j) >= 1)
      {
        return 0;
      }
    }
  }
  return 1;
}*/

int my_dgesv(int n, double *a, int *ipiv, double *b)
{
  int i,j,k;
  float ratio;
  for(i=0;i<n;i++)
  {
	
	  printf("a[%d][%d] = %.10f\n", i, j, *(a + i*n + j));
	  printf("a[%d][%d] = %.10f\n", i, j, a[i*n+j]);
	  if (*(a + i*n + i) == 0)
	  {
		printf("Mathematical Error!");
		exit(0);
	  }
	  for(j=i+1;j<n;j++)
	  {
		ratio = *(a + i*n + j)/(*(a + i*n + i));
		for(k=0;k<n+1;k++)
		{
			*(a + i*n + k) = *(a + j*n + k) - ratio*(*(a + i*n + k));
		}
	  }
  }
  printf("float[%d] = %f\n",n-1, *(b + (n-1)*n + (n)) / (*(a + (n-1)*n + (n))));
  printf("float[%d] = %f\n",0, b[0] / a[0]);
  printf("float[%d] = %f\n",1, b[1*n] / a[1*n]);
  //ipiv[n-1] = *(b + (n-1)*n + (n-1)) / (*(a + (n-1)*n + (n-1)));
  ipiv[n-1] = *(b + (n-1)*n + (n)) / (*(a + (n-1)*n + (n)));
  printf("ipiv2[%d] = %d\n", n-1, ipiv[n-1]);

  for(i=n-2;i>=0;i--)
  {
	ipiv[i] = *(a + i*n + n);
	for (j=i+1;j<n;j++)
	{
	  ipiv[i] = ipiv[i] - *(a + i*n + j) * ipiv[j];
	}
	ipiv[i] = ipiv[i] / *(a + i*n +i);
  }
  printf("new ipiv2 : \n");
  print_matrix_1d_int(ipiv,n);
  return 1;
}

void print_matrix_1d(double *matrix, int n)
{
	for (int i=0; i<n; i++)
	{
		printf("matrix[%d] = %f\n",i,matrix[i]);
	}
}

void print_matrix_1d_int(int *matrix, int n)
{
	for (int i=0; i<n; i++)
	{
		printf("matrix[%d] = %d\n",i,matrix[i]);
	}
}



void main(int argc, char *argv[])
{
  int size = atoi(argv[1]);

  double *a, *aref;
  double *b, *bref;

  a = generate_matrix(size);
  aref = generate_matrix(size);
  b = generate_matrix(size);
  bref = generate_matrix(size);

  printf("b = \n");
  print_matrix_1d(b,size);
  printf("bref = \n");
  print_matrix_1d(bref,size);

  // Using LAPACK dgesv OpenBLAS implementation to solve the system
  int n = size, nrhs = size, lda = size, ldb = size, info;
  int *ipiv = (int *) malloc(sizeof(int) * size);

  clock_t tStart = clock();
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
  printf("info = %d\n",info);
  printf("new bref = \n");
  print_matrix_1d(bref,size);
  printf("Time taken by OpenBLAS LAPACK: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);
  printf("ipiv 1 = %d\n", ipiv[n-1]);
  int *ipiv2 = (int *) malloc(sizeof(int) * size);

  tStart = clock();
  //my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
  printf("result : %d\n", my_dgesv(size, a, ipiv2, b));
  printf("Time taken by my implementation: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

  if (check_result(bref, b, size) == 1)
    printf("Result is ok!\n");
  else
    printf("Result is wrong!\n");
}
