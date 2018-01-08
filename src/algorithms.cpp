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
//R CMD SHLIB ADMMDAG_cycle.cpp

const int maxIter = 1000;
const double Large = 100000;
int m=0;
int n=0;
double lambda;
double tau;
int nonzero;
double mu=0;//changed &
double opts_tol;
using namespace std;

//for DAG()
bool equal_sigma = false; 


double K;
double z=0;

extern "C"{
 
	/****************************************/
	//check loop in C
	//check whether adding a directed edge from node a to node b induces cycles in graph G.
	// ij-th element of G indicate j --> i
	int Cycle(int node,int *G, int a, int b)
	{
		if(a==b)  return 1;
		int i,j,lo,nBot=0,nTop=0,SLeng=1,bCycle=0;
		int *color,*S;
		color=(int *) Calloc(node,int);
		S=(int *) Calloc(node,int);
		color[a]=1;
		S[0]=a;
		while(SLeng>0)
		{
			i=S[nBot];
			SLeng--;
			nBot++;
			for(j=0;j<node;j++)
			{
				lo= j*node;			 
				if(G[lo+i]==1)
				{
					if(j==b)
					{
						bCycle=1;
						break;
					}
					else if(color[j]==0)
					{
						nTop++;
						S[nTop]=j;
						SLeng++;
						color[j]=1;
					}
				}
			}
			if(bCycle==1) break;
		}
		Free(color);
		Free(S);
		return bCycle;
	}
	/****************************************/
	
void proj(double *A, double *v, int *B);


void DAG(double* X,double *A, int *mm, int *nn, double* KK,
		 double *tautau, int *B, int *NZ,double *sigma, double *tol, double *obj, int *equalsigma)
{
//
//
    equal_sigma = (*equalsigma);
	m = (*mm);
	n = (*nn);
	K = (*KK);
	tau = (*tautau);
	nonzero = (*NZ);
	opts_tol = (*tol);
	z = 0.0;
	
	int i = 0;
	int j = 0;
	int k = 0;
	int iterDC = 0;
	double L = 32;
	double alphap = 0;
	double alpha = 0;
	double beta = 0;
	int bFlag = 0;
	double valueF_s = 0;
	double valueF_search = 1000000;
	double valueF_searchp = 1000000;
	double valueDC = 0;
	double valueDCp = 0;
	double r_sum = 0;
	double l_sum = 0;
	double sse = 0.0;
	
	double *Ap= new double[m*m];
  	if (Ap == 0)
  		Rprintf("Error: memory could not be allocated for Ap");
  	double *s = new  double[m*m];
  	if (s == 0)
		Rprintf("Error: memory could not be allocated for s");
  	double *g = new double[m*m];
  	if (g == 0)
		Rprintf("Error: memory could not be allocated for g");
  	double *v = new  double[m*m];
  	if (v == 0)
		Rprintf("Error: memory could not be allocated for v");  
		double *Av = new double[n*m];
	if (Av == 0)
		Rprintf("Error: memory could not be allocated for Av");	  	
  	double *XXs = new double[n*m];
	if (XXs == 0)
		Rprintf("Error: memory could not be allocated for XXs");
	for(i = 0; i< (m*m);i++){
		Ap[i] = 0;
		s[i]=0;
		g[i]=0;
		v[i]=0;		
	} 
	for(i=0;i<(n*m);i++){
		XXs[i]=0;
	}
	
// use vectors for valueL, valueF, n_steps

	z = (K - nonzero)*tau;  // type conversion
	
	for(iterDC = 1; iterDC <=5;iterDC ++){ //DC loop
#ifdef DEBUG
    Rprintf("iterDC=%d\n", iterDC);  
#endif /* DEBUG */
		L = 16;
		for(i = 0; i< (m*m);i++){
			Ap[i] = A[i];			
		} //Ap = A
		alphap = 0;
		alpha =1;
		
		for(int iterStep = 1; iterStep <= maxIter; iterStep++)
		{ // gradient projection loop
 
			bFlag = 0;
			beta = (alphap - 1)/alpha;
			//s = A + beta * (A - Ap);
			for(i = 0; i<(m*m);i++){
				s[i] = A[i] + beta * (A[i] - Ap[i]);
			}
			//XXs =X - X%*%sT
			sse = 0;
			for(i = 0; i<n; i++){
				for(j = 0; j<m;j++){
					XXs[n*j+i] = X[n*j+i];
					for(k = 0; k<m; k++){
						XXs[n*j+i] -= X[n*k+i]*s[m*k + j];
					}
					sse += XXs[n*j+i]*XXs[n*j+i];
				}
			}
			if (equal_sigma){
			    valueF_s = sse;
			}else{
                //sigma = colMeans((XXs)^2) 
                for(j = 0; j<m; j++){
                    sigma[j]=0;
                    for(i=0;i<n;i++){
                        sigma[j] += XXs[n*j+i]*XXs[n*j+i];
                    }
                    sigma[j] /= n;
                }
                        
                valueF_s = 0;
                for(i=0;i<m;i++){
                    valueF_s += log(sigma[i]);    //???
                }   //valueF_s = sum(log(sigma))
			}
			
			// calculate the gradient g at s
			for(i=0;i<m;i++){
				for(j=0;j<m;j++){
					if(i==j) g[j*m+i]=0;
					else{
						g[j*m+i] = 0;
						for(k=0;k<n;k++){
							g[j*m+i] -= XXs[n*i+k]*X[n*j+k];
						}
						g[j*m+i] = g[j*m+i]*2/n/sigma[i];
					}
				}
			}
//#ifdef DEBUG
//Rprintf("sigma=%f %f %f \n", sigma[0],sigma[1],sigma[2]); 
//Rprintf("g=%f %f %f %f %f %f %f %f %f \n", g[0],g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8]); 	
//#endif /* DEBUG */	
			
			for(i = 0; i< (m*m);i++){
				Ap[i] = A[i];			
			} //Ap = A
			
			
			
			
			
			while(1){
				for(i = 0; i< (m*m);i++){
					v[i] = s[i] - g[i]/L;			
				} //v = s - g/L
			
				
				if(z<0){
					Rprintf("nonzero=%d \n", nonzero); 
					Rprintf("z<0 happened. Return.");
					return;
				}
				//projection AT = proj_mat(v,m,tau,z,B)
				proj(A,v,B);
				
				
				for(i = 0; i< (m*m);i++){
					v[i] = A[i]-s[i];			
				} //v = A-s 
				sse = 0;
				//XXs =X - X%*%AT  (might need to define XXA?)
				for(i = 0; i<n; i++){
					for(j = 0; j<m;j++){
						XXs[n*j+i] = X[n*j+i];
						for(k = 0; k<m; k++){
							XXs[n*j+i] -= X[n*k+i]*A[m*k + j];
						}
						sse += XXs[n*j+i]*XXs[n*j+i];
					}
				}
                

                if (equal_sigma){
                    valueF_search = sse;
                }else{
                    //sigma = colMeans((XXs)^2) 
                    for(j = 0; j<m; j++){
                        sigma[j]=0;
                        for(i=0;i<n;i++){
                            sigma[j] += XXs[n*j+i]*XXs[n*j+i];
                        }
                        sigma[j] /= n;
                    }
                        
                    valueF_search = 0;
                    for(i=0;i<m;i++){
                        valueF_search += log(sigma[i]);    //???
                    }   //valueF_search = sum(log(sigma))
                }
			    
			    
 
				
				// update Xv;
				for(i = 0; i<n; i++){
					for(j = 0; j<m;j++){
						Av[n*j+i] = 0;
						for(k = 0; k<m; k++){
							Av[n*j+i] += X[n*k+i]*v[m*j + k];
						}
					}
				}
				
								// l_sum=Av'*Av;
								l_sum = 0;
								for(i=0;i<(n*m);i++){
									l_sum += Av[i]*Av[i];
								}
				r_sum = 0;
				for(i=0;i<(m*m);i++){
					r_sum += v[i]*v[i];  //r_sum=sum(v^2)
										
				}
				if(r_sum <= 1e-20){      //???
					bFlag = 1;//this shows that the gradient step makes litt
				}
				if(l_sum <= r_sum*L/2) break;    //???le improvement
				else L= MAX( 2*L,l_sum/r_sum);  

			} // end of line search while loop
			
			alphap = alpha;
			alpha = (1 + sqrt(4*alpha*alpha + 1))/2;
			
			if(bFlag){
				Rprintf("The program terminates as the gradient step changes the solution very small.");
				bFlag = 0;
				break;
			}
			
						if(ABS((valueF_search - valueF_searchp)) <= opts_tol&& iterStep>=3){  //???
				break;
			}
#ifdef DEBUG
Rprintf("L=%f\n", L); 
Rprintf("valueF_search - valueF_searchp=%f\n",valueF_search - valueF_searchp);
Rprintf("valueF_search=%f\n",valueF_search);
Rprintf("valueF_s=%f\n",valueF_search);
#endif /* DEBUG */				
			valueF_searchp = valueF_search;

		} //end of projection loop
		valueDC = valueF_search;
		nonzero = 0;
/*		for(i=0;i<(m*m);i++){
			if(ABS((A[i]))>=tau){            
				B[i] = 1;
				nonzero++;
			}else B[i]=0;
		} */ //change to the bellow
		for(i=1;i<m;i++){
		    for(j=0;j<i;j++){
				if(ABS(A[j*m+i])+ABS(A[i*m+j]) >= tau){
				    
				    
				    if(ABS(A[j*m+i])>=ABS(A[i*m+j])){
				    //added 05192014
				    	if( Cycle(m,B, j, i)==0){
							B[j*m+i]=1;
							B[i*m+j]=0;
							nonzero++; 
						}else if(Cycle(m,B, i, j )==0){
							B[j*m+i]=0;
							B[i*m+j]=1;
							nonzero++; 
						}else {
							B[j*m+i]=0;
							B[i*m+j]=0;
						}
				    }
				    else{
				    	if(Cycle(m,B, i, j)==0){
							B[j*m+i]=0;
							B[i*m+j]=1;
							nonzero++; 
						} else if( Cycle(m,B, j, i)==0){
							B[j*m+i]=1;
							B[i*m+j]=0;
							nonzero++; 
						}else {
							B[j*m+i]=0;
							B[i*m+j]=0;
						}					 
				    }
				}else{
					B[j*m+i]=0;
					B[i*m+j]=0;
				}
		    }
		}
		z = (K-nonzero)*tau;
		
		if(ABS((valueDC - valueDCp)) <= opts_tol && iterDC>=2) break;  //???
		
		valueDCp=valueDC;
		
	}//end of DC loop
	*obj = valueDC;
	delete [] Ap;
	delete [] s;
	delete [] g;
	delete [] v;
	delete [] XXs;
	delete [] Av;

}//end of DAG



double func(const double mu, double *u, int *B){	
	
	double res = 0.0;
	int ii,jj;
	double max12;
	for(ii=1; ii < m;ii++){
		for(jj=0;jj<ii;jj++){
			if((B[jj*m+ii]!=1) && (B[ii*m+jj]!=1)){
				max12 = MAX((u[jj*m+ii] - mu),0) + MAX((u[ii*m+jj] - mu), 0);
				res += MIN(max12,tau);
			}
		}		
	}
		return (res - z);
} //end of func

void proj(double *A, double *v, int *B){
	double *u = new double[(m*m)];
	int i,j,iter;
	double mu = 0;	
	double u_max = 0.0; //the maximal absolute value in u
	const int NITER = 500;
	double diff;
	for(i=0;i<(m*m);i++){
		u[i] = ABS(v[i]);
		if (u[i] > u_max)
        	u_max=u[i];  
	}
	if(func(0,u,B)<=0) mu = 0;   //??
	else{
		double mu_1, mu_2,f_1,f_2,f_mu;
		mu_1=0; mu_2=u_max;
    	f_1=func(mu_1,u,B);
    
    	f_2=func(mu_2,u,B);
    	if(f_1*f_2 >0){
    		Rprintf("\n Endpoints should have opposite signs!");
			return;
    	} 
    
    	iter = 1;
    	while(iter<=NITER){
    		mu = (mu_1 + mu_2)/2;
    		f_mu = func(mu,u,B);
    		if(f_mu == 0.0 || (mu_2 - mu_1)/2 < opts_tol) break;
    		iter++;
    		if(f_mu * f_1 >= 0){ mu_1 = mu;f_1=f_mu;}
    		else {mu_2 = mu;f_2 = f_mu;}
    	}
    	// check iter==NITER
    	
	
	}//bisection
	for(i=0;i<(m*m);i++){
    	A[i]=0;
    }
    
	for(i=1;i<m;i++){
		for(j=0;j<i;j++){
			if(B[j*m+i]==1){
				A[j*m+i]=v[j*m+i];
				A[i*m+j]=0;
			}
			else if(B[i*m+j]==1){
				A[j*m+i]=0;
				A[i*m+j]=v[i*m+j];
			}
			else{
				if((u[j*m+i]>=mu) && (u[i*m+j]>=mu)){   //???
					diff = u[j*m+i] - u[i*m+j];
					if(diff >= tau) A[j*m+i] = tau*SIGN(v[j*m+i]);  //???
					else if(diff <= -tau) A[i*m+j]=tau*SIGN(v[i*m+j]);
					else{
						A[j*m+i]=MIN((tau+diff)/2,u[j*m+i]-mu)*SIGN(v[j*m+i]);
						A[i*m+j]=MIN((tau-diff)/2,u[i*m+j]-mu)*SIGN(v[i*m+j]);
					}
				}
				else if(u[j*m+i]>=mu){
					A[j*m+i] = MIN(tau,u[j*m+i]-mu) * SIGN(v[j*m+i]);
				}
				else if(u[i*m+j]>=mu){
					A[i*m+j] = MIN(tau,u[i*m+j]-mu) * SIGN(v[i*m+j]);
				}
			
			}
		}
	}
	
	

	delete [] u;
} //end of proj

void DAG_uncon(double* X,double *A, int *mm, int *nn, double* Lambda,//changed &
		 double *tautau, int *A_NZ, int *NZ,double *sigma, double *tol, double *obj,
		 double* XTX, double* XTX_inv,double *rhorho,double *rhorho2)  			
{
//
//
	m = (*mm);
	int m_A = m-1;
	n = (*nn);
	lambda = (*Lambda); //changed &
	tau = (*tautau);
	nonzero = (*NZ);
	opts_tol = (*tol);
	mu = lambda/tau;  
	int count_proj =0;
	
	int i = 0;
	int j = 0;
	int k = 0;
	int pos = 0;
	int iterDC = 0;
	
	double valueF_s = 0.0;  //need to work on these
	double valueF_search = 0.0;
	double valueF_searchp = 0.0;
	double valueDC = 0.0;
	double valueDCp = 0.0;
	double r_sum = 0.0;
	double l_sum = 0.0;
	double mu_temp=0.0;
	double ADMM_value = 0.0;
	double ADMM_valuep = 0.0;
	
	double Bab_temp = 0.0;
	double rho = (*rhorho); //changed
	double rho2 = (*rhorho2);

	// define: B, U, Bab,
	double *B= new double[m*m];
  	if (B == 0)
  		Rprintf("Error: memory could not be allocated for B");
	double *Bab= new double[m*m];  // absolute value of B
  	if (Bab == 0)
  		Rprintf("Error: memory could not be allocated for Bab");
	double *U= new double[m*m];  // U in the paper
  	if (U == 0)
  		Rprintf("Error: memory could not be allocated for U");
	double *XAT = new double[n*m];  // X%*%A^T
	if (XAT == 0)
		Rprintf("Error: memory could not be allocated for XAT");
	double *V_A = new double[(m-1)];
	double *A_row = new double[(m-1)];
	double obj_cur = 0.0;
	double obj_pre = 0.0;
	
  	
	
	
	
/****************** FOR OTHER DIRECTION *********************************/
	int iter=0;
// new variable for projection

	
	double *LLambda = new double[(m*m)]; //lambda
	//initiate LLambda, B, Bab, U
	for(i=0;i<(m*m);i++){
	  B[i]=0;
	  Bab[i]=0;
	  U[i]=0;
	  LLambda[i]=1;
	}
	double *Mat_step = new double[(m*m)];
	
	// initiate Mat_step
	for(i = 0; i<m;i++){
	    for(j =0; j<m; j++){
		if(j==0) {
		  Mat_step[j*m + i] = 1;
		}else if(i==0){
		  Mat_step[j*m + i] = 0;
		}else if(i==j){
		   Mat_step[j*m + i] = 2/((double)m);
		}else{
		   Mat_step[j*m + i] = 1/((double)m);
		}	      
	    }	  
	}
	
	
	double *V_step = new double[m*m];
	for(i=0;i<(m*m);i++){
		V_step[i]=0;
	}
		double AA_temp=0;
		//for blas
			double one=1.0; double zero =0.0; int one_int =1;
			char transa='N'; char transb='N'; char transT = 'T' ;
			// primal residual
			double res=0;
			double res_temp=0;
			double ***Xi_3array;
			// Allocate memory
			Xi_3array = new double**[m];//\Xi in the paper
			for ( i = 0; i < m; ++i) {
				Xi_3array[i] = new double*[m];
				
				for ( j = 0; j < m; ++j)
					Xi_3array[i][j] = new double[m];
			}
						
							    
			double ***y_3array; // y in the paper
			// Allocate memory
			y_3array = new double**[m];
				for ( i = 0; i < m; ++i) {
					y_3array[i] = new double*[m];
									
					for (j = 0; j < m; ++j)
						y_3array[i][j] = new double[m];
				}
											
												
													
			// initiate Xi_3array, y_3array;  
			for ( i = 0; i < m; ++i) {
				for ( j = 0; j < m; ++j){
					for ( k = 0; k < m; ++k){
						Xi_3array[i][j][k] = 0;
						y_3array[i][j][k] = 0;
					}
				}
			}
// calculate objective function value: norm((X - X*A'), 'fro')-Mu*sum(sum(abs(A)))
			
			F77_CALL(dgemm)(&transa,&transT,&n,&m,&m,&one,X,&n,A,&m,&zero,XAT,&n);
			obj_cur = 0.0;
			for(pos=0;pos<(n*m);pos++){
				obj_cur += (X[pos] - XAT[pos])*(X[pos] - XAT[pos])/2 ;
			}
			for(pos=0;pos<(m*m);pos++){
				if(A_NZ[pos]!=1){
					obj_cur += mu*ABS(A[pos]);
				}
			}
			 
			obj_pre = obj_cur;
/****************** FOR OTHER DIRECTION *********************************/															
// use vectors for valueL, valueF, n_steps

	
	
	for(iterDC = 1; iterDC <=6;iterDC ++){ //DC loop
 
		
		for(iter=0;iter<maxIter;iter++){
		// A direction
		
		for(k=0;k<m;k++){
			//calculate V_A
			for(i=0;i<(m-1);i++){
				if(i>=k) j=i+1;else j=i;
				V_A[i]=XTX[j*m+k]+rho2*(B[j*m+k]-U[j*m+k]);
			}
			//calcluate A			
			F77_CALL(dgemv)( &transa, &m_A,&m_A, &one, &XTX_inv[k*(m-1)*(m-1)], &m_A, V_A, &one_int,
							   &zero, A_row, &one_int );
			for(i=0;i<(m-1);i++){
				if(i>=k) j=i+1;else j=i;
				A[j*m+k]=A_row[i];
			}
		}
			//ADMM_valuep = ADMM_value;
			//ADMM_value = valueF_search;
			// update AA
			//for(i = 0 ; i<m*m; i++){
			//	if(B[i]==1){
			//		AA[i]=tau;
			//	}else{
			//		AA[i]=ABS(A[i]);
			//	}
			//}
		//B direction
		for(pos = 0; pos<(m*m);++pos){
			j = pos/m;
			i = pos%m;
			Bab_temp = 0.0;
			for(k = 0;k<m;k++){
				Bab_temp = Bab_temp + LLambda[k*m + i] - LLambda[k*m+j] -
				Xi_3array[i][j][k] - y_3array[i][j][k];
				
			}
			Bab_temp = Bab_temp + (m-1)*tau;
			if(A_NZ[pos]==1){
				Bab[pos] = tau;
				B[pos] = A[pos]+U[pos];
			}else{
				Bab[pos] = MAX(0,((rho2 * ABS(A[pos]+U[pos]) + rho*Bab_temp - mu)/(rho2+m*rho)));
				B[pos] = Bab[pos]*(A[pos]+U[pos]>=0?1:-1);
			}
			
		}
		//U direction
		for(pos = 0; pos<(m*m);++pos){
			U[pos] = U[pos] + A[pos] - B[pos];
		}
		
			//initiate V_step;
			for ( i = 0; i < m; ++i) {
			  for ( k = 0; k < m; ++k){
			    V_step[k*m+i]=0;	   
			    for(j=0; j<m; ++j){
			      V_step[k*m+i] += Bab[j*m+i] + Xi_3array[i][j][k] + y_3array[i][j][k];		      
			    }
			  }
			}
			for ( j = 0; j < m; ++j) {
			  for ( k = 0; k < m; ++k){	   
			    for(i=0; i<m; ++i){
			      V_step[k*m+j] -= (Bab[j*m+i] + Xi_3array[i][j][k] + y_3array[i][j][k]);		      
			    }
			  }
			}
			for ( i = 1; i < m; ++i) {
			  for ( k = 0; k < m; ++k){	  	   
			    if(i==k){
			      V_step[k*m+i] = (-(m-1)*tau + V_step[k*m+i])/2;	
			    }else{
			      V_step[k*m+i] = (tau + V_step[k*m+i])/2;		      
			    }
			  }
			}
			//V_step(1,:) = ones(1,p);
			for ( k = 0; k < m; ++k){
			  V_step[k*m] = 1;
			}

			//update LLambda; Mat_step * V_step;
			res=0;
			
			F77_CALL(dgemm)(&transa,&transb,&m,&m,&m,&one,Mat_step,&m,V_step,&m,&zero,LLambda,&m);
			//update Xi_3array and y_3array
			for ( i = 0; i < m; ++i) {
			  for ( j = 0; j < m; ++j){
			    if(i!=j){
			    for ( k = 0; k < m; ++k){	      
			      if(k==j){
				res_temp= Bab[j*m+i]+Xi_3array[i][j][k] -LLambda[k*m+i]+LLambda[k*m+j];
				Xi_3array[i][j][k] = MAX(0,(LLambda[k*m+i]-LLambda[k*m+j]-Bab[j*m+i]-y_3array[i][j][k]));
				y_3array[i][j][k] =y_3array[i][j][k] + res_temp;
				//res+=res_temp*res_temp;
				res = MAX(res,(res_temp*res_temp));
			      }else{
				res_temp = Bab[j*m+i]+Xi_3array[i][j][k] -LLambda[k*m+i]-tau+LLambda[k*m+j];
				Xi_3array[i][j][k] = MAX(0,(LLambda[k*m+i]+tau-LLambda[k*m+j]-Bab[j*m+i]-y_3array[i][j][k]));
				y_3array[i][j][k] = y_3array[i][j][k] + res_temp;
				//res+=res_temp*res_temp;
				res = MAX(res,(res_temp*res_temp));
			      }
			    }
			    }
			  }
			}	
			obj_pre = obj_cur;
			F77_CALL(dgemm)(&transa,&transT,&n,&m,&m,&one,X,&n,A,&m,&zero,XAT,&n);
			obj_cur = 0.0;
			for(pos=0;pos<(n*m);pos++){
				obj_cur += (X[pos] - XAT[pos])*(X[pos] - XAT[pos])/2 ;
			}
			for(pos=0;pos<(m*m);pos++){
				if(A_NZ[pos]!=1){
					obj_cur += mu*ABS(A[pos]);
				}
			}
			 
			
			
						if(ABS((obj_cur-obj_pre))<opts_tol &&(res <= tau*tau*1e-1)) break;
			//if((res <= tau*1e-2)&&(iter>=3)&&(ABS((ADMM_value-ADMM_valuep))<=opts_tol*tau)) break;
		} // end of ADMM	
#ifdef DEBUG
			Rprintf("iter=%d\n", iter); 
#endif /* DEBUG */	
		valueDCp=valueDC;		
		valueDC = obj_cur;
		nonzero = 0;
/*		for(i=0;i<(m*m);i++){
			if(ABS((A[i]))>=tau*(1.0-0.2)){            
				A_NZ[i] = 1;
				nonzero++;
			}else A_NZ[i]=0;
		}  //change to the bellow */
		for(i=1;i<m;i++){
		    for(j=0;j<i;j++){
				if((ABS(A[j*m+i])+ABS(A[i*m+j])) >= tau*0.99){
				    
					if(ABS(A[j*m+i])>tau*0.5){
						if( Cycle(m,A_NZ, j,i)==0){
							
							A_NZ[j*m+i]=1;
							A_NZ[i*m+j]=0;
							nonzero++;
						}else{      // if (Cycle(m,A_NZ, i,j) == 0) ...
							A_NZ[j*m+i]=0;
							A_NZ[i*m+j]=0;
						}
				    }else{
						if( Cycle(m,A_NZ, i,j)==0){
							
							A_NZ[j*m+i]=0;
							A_NZ[i*m+j]=1;
							nonzero++;
						}else {
							A_NZ[j*m+i]=0;
							A_NZ[i*m+j]=0;
						}

					
					}
		    
				}
			}
		}
		// changed & z = (K-nonzero)*tau; */
		
		if(ABS((valueDC - valueDCp)) <= opts_tol && iterDC>=2 ) break;  //???
		
		
		
	}//end of DC loop
		#ifdef DEBUG
		Rprintf("valueDC=%f\n", valueDC); 
		#endif /* DEBUG */				
	*obj = valueDC;
	Rprintf("iterDC=%d\n", iterDC);
	Rprintf("count_proj=%d\n", count_proj);
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
	for (int i = 0; i < m; ++i) {
	  for (int j = 0; j < m; ++j)
	    delete [] Xi_3array[i][j];

	  delete [] Xi_3array[i];
	}
	delete [] Xi_3array;
	// De-Allocate memory to prevent memory leak
	for (int i = 0; i < m; ++i) {
	  for (int j = 0; j < m; ++j)
	    delete [] y_3array[i][j];

	  delete [] y_3array[i];
	}
	delete [] y_3array;
	
}//end of DAG


/*
double func(const double mu, double *u, int *B){	
	
	double res = 0.0;
	int ii,jj;
	double max12;
	for(ii=1; ii < m;ii++){
		for(jj=0;jj<ii;jj++){
			if((B[jj*m+ii]!=1) && (B[ii*m+jj]!=1)){
				max12 = MAX((u[jj*m+ii] - mu),0) + MAX((u[ii*m+jj] - mu), 0);
				res += MIN(max12,tau);
			}
		}		
	}
		return (res - z);
} //end of func
*/  //changed &


}//end of extern
