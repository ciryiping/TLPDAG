#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b) ((a) < (b) ?  (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define SIGN(a) ((a) < 0 ? (-1) : (1))
//#define DEBUG



const int maxIter = 1000;
int m = 0;
int n = 0;
double lambda;
double tau;
int nonzero;
double mu = 0;
double opts_tol;
using namespace std;
extern "C" {
    /*
    	check whether adding a directed edge from node a to node b induces cycles in graph G.
    	ij-th element of G indicates an edge j --> i
    */
    int Cycle(int node, int *G, int a, int b)
    {
        if(a == b)  return 1;
        int i, j, lo, nBot = 0, nTop = 0, SLeng = 1, bCycle = 0;
        int *color, *S;
        color = (int *) Calloc(node, int);
        S = (int *) Calloc(node, int);
        color[a] = 1;
        S[0] = a;
        while(SLeng > 0)
        {
            i = S[nBot];
            SLeng--;
            nBot++;
            for(j = 0; j < node; j++)
            {
                lo = j * node;
                if(G[lo + i] == 1)
                {
                    if(j == b)
                    {
                        bCycle = 1;
                        break;
                    }
                    else if(color[j] == 0)
                    {
                        nTop++;
                        S[nTop] = j;
                        SLeng++;
                        color[j] = 1;
                    }
                }
            }
            if(bCycle == 1) break;
        }
        Free(color);
        Free(S);
        return bCycle;
    }


    /*************************************************************************/
    /*********************************************************************/
    // Arguments
    // X: n by m data matrix
    // A: m by m adjacency matrix
    // m: the number of variables
    // n: sample size
    // lambda = (*Lambda): penalty (\mu in the paper)

    // tau: (\tau in the paper)
    // A_NZ: m by m  matrix (0 indicates the coresponding element of A being 0,
    //       1 indicates nonzero. -1 indicates a suggested direction, which will be
    //       checked outside this function. If this suggested direction does not violate
    //       acyclicity, -1 will be converted to 1.
    // NZ: the number of nonzero elements in A. ( the number of 1's in A_NZ)
    // sigma: m by 1 variances vector (\sigma in the paper). In this version, they are fixed at 1's.
    // tol: tolerance level
    // obj: objective function value
    // XTX: m by m matrix X^T %*% X
    // XTX_inv: inverse matrix of XTX

    //  Other variables:
    // B: m by m matrix (B in the paper)
    // LLambda: m by m dual variable matrix (\lambda in the paper)
    /*********************************************************************/

    void uncon_dag_(double *X, double *A, int *mm, int *nn, double *Lambda,
                   double *tautau, int *A_NZ, int *NZ, double *sigma, double *tol, double *obj,
                   double *XTX, double *XTX_inv, double *rhorho, double *rhorho2)
    {

        m = (*mm);
        int m_A = m - 1;
        n = (*nn);
        lambda = (*Lambda); //changed &
        tau = (*tautau);
        nonzero = (*NZ);
        opts_tol = (*tol);
        mu = lambda / tau;
        //int count_proj = 0;

        int i = 0;
        int j = 0;
        int k = 0;
        int pos = 0;
        int iterDC = 0;


        double valueDC = 0.0;
        double valueDCp = 0.0;


        double Bab_temp = 0.0;
        double rho = (*rhorho);
        double rho2 = (*rhorho2);

        // define: B, U, Bab,
        double *B = new double[m * m];
        if (B == 0)
            Rprintf("Error: memory could not be allocated for B");
        double *Bab = new double[m * m]; // absolute value of B
        if (Bab == 0)
            Rprintf("Error: memory could not be allocated for Bab");
        double *U = new double[m * m]; // U in the paper
        if (U == 0)
            Rprintf("Error: memory could not be allocated for U");
        double *XAT = new double[n * m]; // X%*%A^T
        if (XAT == 0)
            Rprintf("Error: memory could not be allocated for XAT");
        double *V_A = new double[(m - 1)];
        double *A_row = new double[(m - 1)];
        double obj_cur = 0.0;
        double obj_pre = 0.0;

        /*FOR OTHER DIRECTION */
        int iter = 0;

        double *LLambda = new double[(m * m)]; //lambda
        //initiate LLambda, B, Bab, U
        for(i = 0; i < (m * m); i++)
        {
            B[i] = 0;
            Bab[i] = 0;
            U[i] = 0;
            LLambda[i] = 1;
        }
        double *Mat_step = new double[(m * m)];

        // initiate Mat_step
        for(i = 0; i < m; i++)
        {
            for(j = 0; j < m; j++)
            {
                if(j == 0)
                {
                    Mat_step[j * m + i] = 1;
                }
                else if(i == 0)
                {
                    Mat_step[j * m + i] = 0;
                }
                else if(i == j)
                {
                    Mat_step[j * m + i] = 2 / ((double)m);
                }
                else
                {
                    Mat_step[j * m + i] = 1 / ((double)m);
                }
            }
        }

        double *V_step = new double[m * m];
        for(i = 0; i < (m * m); i++)
        {
            V_step[i] = 0;
        }
        //double AA_temp=0;
        //for blas
        double one = 1.0;
        double zero = 0.0;
        int one_int = 1;
        char transa = 'N';
        char transb = 'N';
        char transT = 'T' ;
        // primal residual
        double res = 0;
        double res_temp = 0;
        double ***Xi_3array;
        // Allocate memory
        Xi_3array = new double **[m]; //\Xi in the paper
        for ( i = 0; i < m; ++i)
        {
            Xi_3array[i] = new double*[m];

            for ( j = 0; j < m; ++j)
                Xi_3array[i][j] = new double[m];
        }


        double ***y_3array; // y in the paper
        // Allocate memory
        y_3array = new double **[m];
        for ( i = 0; i < m; ++i)
        {
            y_3array[i] = new double*[m];

            for (j = 0; j < m; ++j)
                y_3array[i][j] = new double[m];
        }



        // initiate Xi_3array, y_3array;
        for ( i = 0; i < m; ++i)
        {
            for ( j = 0; j < m; ++j)
            {
                for ( k = 0; k < m; ++k)
                {
                    Xi_3array[i][j][k] = 0;
                    y_3array[i][j][k] = 0;
                }
            }
        }
        // calculate objective function value: norm((X - X*A'), 'fro')-Mu*sum(sum(abs(A)))

        F77_CALL(dgemm)(&transa, &transT, &n, &m, &m, &one, X, &n, A, &m, &zero, XAT, &n);
        obj_cur = 0.0;
        for(pos = 0; pos < (n * m); pos++)
        {
            obj_cur += (X[pos] - XAT[pos]) * (X[pos] - XAT[pos]) / 2 ;
        }
        for(pos = 0; pos < (m * m); pos++)
        {
            if(A_NZ[pos] != 1)
            {
                obj_cur += mu * ABS(A[pos]);
            }
        }

        obj_pre = obj_cur;



        for(iterDC = 1; iterDC <= 6; iterDC ++) //DC loop
        {


            for(iter = 0; iter < maxIter; iter++)
            {
                // A direction

                for(k = 0; k < m; k++)
                {
                    //calculate V_A
                    for(i = 0; i < (m - 1); i++)
                    {
                        if(i >= k) j = i + 1;
                        else j = i;
                        V_A[i] = XTX[j * m + k] + rho2 * (B[j * m + k] - U[j * m + k]);
                    }
                    //calcluate A
                    F77_CALL(dgemv)( &transa, &m_A, &m_A, &one, &XTX_inv[k * (m - 1) * (m - 1)], &m_A, V_A, &one_int,
                                     &zero, A_row, &one_int );
                    for(i = 0; i < (m - 1); i++)
                    {
                        if(i >= k) j = i + 1;
                        else j = i;
                        A[j * m + k] = A_row[i];
                    }
                }
         
                //B direction
                for(pos = 0; pos < (m * m); ++pos)
                {
                    j = pos / m;
                    i = pos % m;
                    Bab_temp = 0.0;
                    for(k = 0; k < m; k++)
                    {
                        Bab_temp = Bab_temp + LLambda[k * m + i] - LLambda[k * m + j] -
                                   Xi_3array[i][j][k] - y_3array[i][j][k];

                    }
                    Bab_temp = Bab_temp + (m - 1) * tau;
                    if(A_NZ[pos] == 1)
                    {
                        Bab[pos] = tau;
                        B[pos] = A[pos] + U[pos];
                    }
                    else
                    {
                        Bab[pos] = MAX(0, ((rho2 * ABS(A[pos] + U[pos]) + rho * Bab_temp - mu) / (rho2 + m * rho)));
                        B[pos] = Bab[pos] * (A[pos] + U[pos] >= 0 ? 1 : -1);
                    }

                }
                //U direction
                for(pos = 0; pos < (m * m); ++pos)
                {
                    U[pos] = U[pos] + A[pos] - B[pos];
                }

                //initiate V_step;
                for ( i = 0; i < m; ++i)
                {
                    for ( k = 0; k < m; ++k)
                    {
                        V_step[k * m + i] = 0;
                        for(j = 0; j < m; ++j)
                        {
                            V_step[k * m + i] += Bab[j * m + i] + Xi_3array[i][j][k] + y_3array[i][j][k];
                        }
                    }
                }
                for ( j = 0; j < m; ++j)
                {
                    for ( k = 0; k < m; ++k)
                    {
                        for(i = 0; i < m; ++i)
                        {
                            V_step[k * m + j] -= (Bab[j * m + i] + Xi_3array[i][j][k] + y_3array[i][j][k]);
                        }
                    }
                }
                for ( i = 1; i < m; ++i)
                {
                    for ( k = 0; k < m; ++k)
                    {
                        if(i == k)
                        {
                            V_step[k * m + i] = (-(m - 1) * tau + V_step[k * m + i]) / 2;
                        }
                        else
                        {
                            V_step[k * m + i] = (tau + V_step[k * m + i]) / 2;
                        }
                    }
                }
                //V_step(1,:) = ones(1,p);
                for ( k = 0; k < m; ++k)
                {
                    V_step[k * m] = 1;
                }

                //update LLambda; Mat_step * V_step;
                res = 0;

                F77_CALL(dgemm)(&transa, &transb, &m, &m, &m, &one, Mat_step, &m, V_step, &m, &zero, LLambda, &m);
                //update Xi_3array and y_3array
                for ( i = 0; i < m; ++i)
                {
                    for ( j = 0; j < m; ++j)
                    {
                        if(i != j)
                        {
                            for ( k = 0; k < m; ++k)
                            {
                                if(k == j)
                                {
                                    res_temp = Bab[j * m + i] + Xi_3array[i][j][k] - LLambda[k * m + i] + LLambda[k * m + j];
                                    Xi_3array[i][j][k] = MAX(0, (LLambda[k * m + i] - LLambda[k * m + j] - Bab[j * m + i] - y_3array[i][j][k]));
                                    y_3array[i][j][k] = y_3array[i][j][k] + res_temp;
                                    //res+=res_temp*res_temp;
                                    res = MAX(res, (res_temp * res_temp));
                                }
                                else
                                {
                                    res_temp = Bab[j * m + i] + Xi_3array[i][j][k] - LLambda[k * m + i] - tau + LLambda[k * m + j];
                                    Xi_3array[i][j][k] = MAX(0, (LLambda[k * m + i] + tau - LLambda[k * m + j] - Bab[j * m + i] - y_3array[i][j][k]));
                                    y_3array[i][j][k] = y_3array[i][j][k] + res_temp;
                                    //res+=res_temp*res_temp;
                                    res = MAX(res, (res_temp * res_temp));
                                }
                            }
                        }
                    }
                }
                obj_pre = obj_cur;
                F77_CALL(dgemm)(&transa, &transT, &n, &m, &m, &one, X, &n, A, &m, &zero, XAT, &n);
                obj_cur = 0.0;
                for(pos = 0; pos < (n * m); pos++)
                {
                    obj_cur += (X[pos] - XAT[pos]) * (X[pos] - XAT[pos]) / 2 ;
                }
                for(pos = 0; pos < (m * m); pos++)
                {
                    if(A_NZ[pos] != 1)
                    {
                        obj_cur += mu * ABS(A[pos]);
                    }
                }



                if(ABS((obj_cur - obj_pre)) < opts_tol && (res <= tau * tau * 1e-1)) break;
                //if((res <= tau*1e-2)&&(iter>=3)&&(ABS((ADMM_value-ADMM_valuep))<=opts_tol*tau)) break;
            } // end of ADMM
#ifdef DEBUG
            Rprintf("iter=%d\n", iter);
#endif /* DEBUG */
            valueDCp = valueDC;
            valueDC = obj_cur;
            nonzero = 0;
            for(i = 1; i < m; i++)
            {
                for(j = 0; j < i; j++)
                {
                    if((ABS(A[j * m + i]) + ABS(A[i * m + j])) >= tau * 0.99)
                    {

                        if(ABS(A[j * m + i]) > tau * 0.5)
                        {
                            if( Cycle(m, A_NZ, j, i) == 0)
                            {

                                A_NZ[j * m + i] = 1;
                                A_NZ[i * m + j] = 0;
                                nonzero++;
                            }
                            else        // if (Cycle(m,A_NZ, i,j) == 0) ...
                            {
                                A_NZ[j * m + i] = 0;
                                A_NZ[i * m + j] = 0;
                            }
                        }
                        else
                        {
                            if( Cycle(m, A_NZ, i, j) == 0)
                            {

                                A_NZ[j * m + i] = 0;
                                A_NZ[i * m + j] = 1;
                                nonzero++;
                            }
                            else
                            {
                                A_NZ[j * m + i] = 0;
                                A_NZ[i * m + j] = 0;
                            }


                        }

                    }
                }
            }
            

            if(ABS((valueDC - valueDCp)) <= opts_tol && iterDC >= 2 ) break; 



        }//end of DC loop
#ifdef DEBUG
        Rprintf("iterDC=%d\n", iterDC);
        Rprintf("valueDC=%f\n", valueDC);              
#endif /* DEBUG */
        *obj = valueDC;
         
        delete [] B;
        delete [] Bab;
        delete [] U;
        delete [] V_A;
        delete [] A_row;
        delete [] XAT;
        delete [] LLambda;
        delete [] Mat_step;
        delete [] V_step;


        // De-Allocate memory to prevent memory leak
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m; ++j)
                delete [] Xi_3array[i][j];

            delete [] Xi_3array[i];
        }
        delete [] Xi_3array;
        // De-Allocate memory to prevent memory leak
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m; ++j)
                delete [] y_3array[i][j];

            delete [] y_3array[i];
        }
        delete [] y_3array;

    }

}
