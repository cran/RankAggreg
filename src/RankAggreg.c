#include <R.h>
//#include <Rinternals.h>
#include <Rmath.h>

void cumsum(double* x, double* csum, int n)
{
	for(int i=0; i< n; i++)
		csum[i]=0;
	csum[0] = x[0];
    for (int i=1; i < n; i++)
		csum[i] = csum[i-1]+x[i];
}

double absd(double x){
	if(x > 0)
		return(x);
	else
		return(-x);
}

int absi(int x){
	if(x > 0)
		return(x);
	else
		return(-x);
}

void sampling(int* NN, double* v, int* nn, int* kk, int* res){
	int i,j,t,m,z;
	int N = *NN;
	int n = *nn;
	int k = *kk;
	double runif;
	int excluded[k];
	int skip = 0;
	double sum = 0;
	
	GetRNGstate();
	for(i=0; i<N; i++){
		for(t=0; t<k; t++)
			excluded[t] = -1;
		for(j=0; j<k; j++){
			runif = unif_rand();
			double probs[n-j];
			double csum[n-j];
			int cands[n-j];
		
			z = 0; //select jth element for the ith candidate list
			for(t=0; t<n; t++){
				skip=0;
				for(m=0; m<j; m++){
					if((t+1)==excluded[m]){
						skip=1;
						break;
					}
				}
				if(skip==0){
					probs[z] = v[n*j+t];
					cands[z++] = t+1;
				}
			}
			
			sum=0;
			for(t=0; t<(n-j); t++)
				sum += probs[t];
			for(t=0; t<(n-j); t++)
				probs[t] = probs[t]/sum;
			
			cumsum(probs, csum, n-j);
			for(t=0; t<(n-j); t++)
				if(runif < csum[t])
					break; //get the t
			
			res[j*N+i] = cands[t];
			excluded[j] = cands[t];
		}
	}
	PutRNGstate();
}	

double spearman(int*x, int*y, int k){
	int i,j;
	double sum = 0;

	for(i=0; i<k; i++){ //go through x first
		for(j=0; j<k; j++)
			if(x[i]==y[j]){
				sum += abs(i-j);
				break;
			}
		if(j==k)
			sum += abs(i-k);
	}
	
	for(i=0; i<k; i++){ //go through second
	//pick only y's not found in x
		for(j=0; j<k; j++)
			if(y[i]==x[j])
				break;
		if(j==k)
			sum += abs(i-k);
	}
	return(sum);
}

double spearmanW(int*x, int*y, int k, double* weights){
	int i,j;
	double sum = 0;
	
	for(i=0; i<k; i++){ //go through x first
		for(j=0; j<k; j++)
			if(x[i]==y[j]){
				sum += absi(i-j)*absd(weights[i]-weights[j]);
				break;
			}
		if(j==k)
			sum += absi(i-k)*weights[i];			
	}
	
	for(i=0; i<k; i++){ //go through second
	//pick only y's not found in x
		for(j=0; j<k; j++)
			if(y[i]==x[j])
				break;
		if(j==k)
			sum += absi(i-k)*weights[i];
	//Rprintf ("sum: %f %d %d %f %f \n", sum, i, j, weights[i], weights[j]);			
	}
	return(sum);
}	
				
void spearmanCands(int* cands, int* x, int* NN, int*kk, 
		int* nn, double* imp, double* res){
	int N = *NN; //number of candidates
	int k = *kk; //top-k list
	int n = *nn; //number of lists to be combined
	double sum = 0;
	int i,j,t;
	int a[k];
	int b[k];
	
	for(i=0; i<N; i++){
		sum = 0;
		for(j=0; j<n; j++){
			for(t=0; t<k; t++){
				a[t] = cands[t*N+i];
				b[t] = x[t*n+j];
			}
			sum += imp[j] * spearman(a, b, k);
		}
		res[i]=sum;
	}
}

void spearmanCandsW(int* cands, int* x, int* NN, int*kk, 
		int* nn, double* imp, double* weights, double* res){
	int N = *NN; //number of candidates
	int k = *kk; //top-k list
	int n = *nn; //number of lists to be combined
	double sum = 0;
	int i,j,t;
	int a[k];
	int b[k];
	double w[k];
	
	for(i=0; i<N; i++){
		sum = 0;
		for(j=0; j<n; j++){
			for(t=0; t<k; t++){
				a[t] = cands[t*N+i];
				b[t] = x[t*n+j];
				w[t] = weights[t*n+j];
			}
			sum += imp[j] * spearmanW(a, b, k, w);
		}
		res[i]=sum;
	}
}

