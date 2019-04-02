/*
Assumptions
1. Even length of filters
*/
/*
Major code borrows:
	1. rosetta fft code
*/
#include<stdio.h>
#include<math.h>
#include"stdio.h"
#define N 16
//#define N 16384
//#define Nfil 256
#define BlkSize 1024
#define pi 3.1415
#define fc 50
#define Acheb 40 //db/decade
#define k 20
#include "math.h"
#include<cufft.h>
#include<stdlib.h>
#include<complex.h>
//#include<fftw.h>

//Defined macros 
double PI;
typedef double complex cplx;
typedef struct{
	double *sig_t;
	cplx *sig_f;
}Filter;

void init_fil(Filter *fil,int arr_len)
{
	fil->sig_t = (double*)calloc(arr_len,sizeof(double));
	fil->sig_f = (cplx*)calloc(arr_len,sizeof(cplx));
}
void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) 
	{
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) 
		{
			cplx t = cexp(-I * PI * i / n) * out[i + step];
			buf[i / 2]     = out[i] + t;
			buf[(i + n)/2] = out[i] - t;
		}
	}
}
 
void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++) out[i] = buf[i];
		_fft(buf, out, n, 1);
}
 
 
void show(const char * s, cplx buf[]) {
	printf("%s", s);
	for (int i = 0; i < 16; i++)
		if (!cimag(buf[i]))
			printf("%g ", creal(buf[i]));
		else
			printf("(%g + j%g) ", creal(buf[i]), cimag(buf[i]));
}
 
void printUnsignedRange(int bytes)
{
    int bits = 8 * bytes;
    unsigned long long to = (1LLU << (bits - 1)) + ((1LL << (bits - 1)) - 1);;    
    printf(" 0 to %llu\n\n", to);
}

void shift_left_cplx(cplx *arr,int len,int shift,int offset)
{
	memcpy(arr+len-offset,arr,shift*sizeof(cplx));
	arr = arr+shift; 
}

void shift_left_fl(float *arr,int len,int shift,int offset)
{
	memcpy(arr+len-offset,arr,shift*sizeof(float));
	arr = arr+shift; 
}
 
unsigned int next_pow2(unsigned int n) 
{ 
    unsigned int p = 1; 
    if (n && !(n & (n - 1))) 
        return n; 
    while (p < n)  
        p <<= 1;      
    return p; 
} 
  

double cheby_poly(double n, double x)
{
		double res;
    if (fabs(x) <= 1)
      res = cos(n*acos(x));
    else
      res = creal(ccosh(n*cacosh(x)));
    return res;
}

void cheby_win(double *out,cplx *fout,int Nfil,double tolerance)
{
		//out = (float*)calloc(Nfil,sizeof(float));
		int _temp = next_pow2(Nfil);
		cplx *temp = (cplx*)calloc(_temp,sizeof(cplx));
    int idx, i;
    //float M, n, sum = 0, max=0;
    double tg = 1/tolerance;//  /* 1/r term [2], 10^gamma [2] */
    double expr = cosh(acosh(tg)/(Nfil-1));
    //M = Nfil/2;
    printf("Nfil:%d x0:%lf _temp:%lu",Nfil,expr,_temp);
    //exit(0);
    for(idx=0; idx<Nfil; idx++)
    {
			temp[idx] = cheby_poly(Nfil-1,expr*cos((pi*idx)/Nfil))*tolerance;//cosf(2.0*n*pi*i/Nfil);
     	fout[idx] = temp[idx];
			printf("Filter:%lf\n",creal(temp[idx]));
    }
		for(idx=0;idx<_temp-Nfil;idx++)
			temp[idx]=0;
		//printf("\nCount:%d\n",idx);
		/*
    for(idx=0; idx<Nfil; idx++)
      out[idx] /= max; // normalise everything 
		*/
		fft(temp,_temp);
		show("\ntempt:",temp);
		exit(0);
		shift_left_cplx(temp,Nfil,Nfil/2,Nfil/2);
		exit(0);
		show("\ntempsht:",temp);		
		//exit(0);
		for(idx=0;idx<Nfil;idx++)
			out[idx] = creal(temp[idx]);
}

float sum_arr_fl(float *arr,int num)
{
	float sum = 0;
	for(int i=0;i<num/4;i++)
		sum += arr[i]+arr[i+int(num/4.0)]+arr[i+int(num/2.0)]+arr[i+int(num*3/4)];
}

float sum_arr_cplx(cplx *arr,int num)
{
	cplx sum = 0;
	for(int i=0;i<num/4;i++)
		sum += arr[i]+arr[i+int(num/4.0)]+arr[i+int(num/2.0)]+arr[i+int(num*3/4)];
}

void normalize_max(cplx *tempx2,int len,int max)
{
	for(int i=0;i<len/2;i++)
	{
		tempx2[i] /= max;
		tempx2[i+len/2] /= max;
	}
}

void ifft_set(cplx *arr,int norm)
{
	normalize_max(arr,norm,norm);
	cplx *arr2;
	arr2 = arr+1;
	cplx tmp = arr[0];
	for(int i=0;i<norm-1;i++)
	{
		cplx tmp = arr[i];
		arr[i] = arr[norm-1-i];
		arr[norm-1-i] = tmp;
	}
}

void init_cplx_from_db(cplx *dest,int len,double *src_real,double *src_imag)
{
	if(src_real==NULL)
	{
		for(int i=0;i<len;i++)
			dest[i] = I*src_imag[i];
	}	
	else if(src_imag==NULL)
	{
		for(int i=0;i<len;i++)
			dest[i] = src_real[i];
	}
	else
	{
		for(int i=0;i<len;i++)
			dest[i] = src_real[i]+I*src_imag[i];
	}
}

void init_db_from_cplx(double *dest,int len,cplx *src)
{
	for(int i=0;i<len;i++)
		dest[i] = creal(src[i]);
}
void conv_window(Filter *fil,int nfil,int bin_size)
{
	//fil.sig_t = (float*)calloc(Nan,sizeof(float));
	//fil.sig_f = (cplx*)calloc(Nan,sizeof(cplx));
	cplx *tempx = (cplx*)calloc(N,sizeof(cplx));
	cplx *tempx2 = (cplx*)calloc(N,sizeof(cplx));
	init_cplx_from_db(tempx,N,fil->sig_t,NULL);
	int width = N - nfil, bias = bin_size/2;
	int x=0;
	show("\nFFT first:\n",tempx);
	exit(0);
	memcpy(tempx+nfil,&x,width*sizeof(cplx));	// Time shifting and padding the signal
	shift_left_cplx(tempx,nfil,nfil/2,0);
	show("\nNo FFT:\n",tempx);
	fft(tempx,N);
	show("\nFFT tempx\n",tempx);
	float max = 0;
	cplx sum = sum_arr_cplx(tempx,nfil);
	for(int i=0;i<N;i++)
	{
		tempx2[(i+bias)%N] = sum;
		int intl = cabs(sum);
		//printf("(%lf,j%lf)",creal(sum),cimag(sum));
		if(max<=intl)
			max = intl;
		sum = sum+tempx[(i+bias)%N]-tempx[i%N];
	} 
	normalize_max(tempx2,nfil,max);
	show("\ntempx2:\n",tempx2);
	cplx mag = 1; 
	cplx step = cexp(-2*pi*I*(nfil/2)/N);
	for(int i=0;i<N;i++)
	{
		tempx2[i] *= mag;
		mag *= step;
	} 
	tempx = tempx2;
	show("\nTempx\n",tempx);
	fft(tempx2,nfil);
	ifft_set(tempx2,nfil);
	init_db_from_cplx(fil->sig_t,nfil,tempx2);
	fil->sig_f = tempx;
	show("\n\nTEMPX\n",fil->sig_f);
	exit(0);
	//free(tempx);
	//free(tempx2);
}
/*
void inner_loop()
{

}
*/
void generate_signal(double x[N],cplx xf[N],int fIdx[k])
{
	unsigned int i,n,idx;
	for(i=0;i<k;i++)
	{
		idx = rand()%N;
		x[idx] = 1;
		fIdx[i] = idx;
		xf[idx] = 1;
	}	
	fft(xf,N);
}

void init_double(double* arr,int num)
{
	for(int i=0;i<num;i++)
		arr[i]=0.0;
}

void init_cplx(cplx* arr,int num)
{
	for(int i=0;i<num;i++)
		arr[i] = 0.0+I*0.0;
}

int main()
{
	// System Constants
	int Bcst_loc = 2;
	int Bcst_est = 2;
	float lobefrac_loc =.003;
	float tolerance_loc = .00000001;
	int b_loc = 68;//Number of bins formed
	int B_loc = 128;
	int B_thresh = 20;
	int loc_loops = 7;
	int threshold_loops = 6;
	float lobefrac_est = .003;
	float tolerance_est = .00000001;
	int b_est = 79;
	int B_est = 128;
	int est_loops = 16;
	int Comb_loops = 10;
	//int k = 10;
	// Secondary Parameters
	int Nfil_loc = (1/pi)*(1.0/lobefrac_loc)*acosh(1.0/tolerance_loc);
	int Nfil_est =  (1/pi)*(1.0/lobefrac_est)*acosh(1.0/tolerance_est);
	if(not(Nfil_loc%2))//Making filter element number odd
		Nfil_loc--;
	if(not(Nfil_est%2))
		Nfil_est--;
	double x[N];		
	double *tfilLoc = (double*)calloc(Nfil_loc,sizeof(double));
	double *tfilEst = (double*)calloc(Nfil_est,sizeof(double));
	cplx xf[N];
	cplx *ffilLoc = (cplx*)calloc(Nfil_loc,sizeof(cplx));
	cplx *ffilEst = (cplx*)calloc(Nfil_est,sizeof(cplx));
	int freqIdx[k],i,j;
	//Nfil_loc = 9;
	//Nfil_est = 9;
	//exit(0);
	//printUnsignedRange(sizeof(unsigned int));
	printf("\nNfil: %d %d",Nfil_loc,Nfil_est);
	init_double(x,N);
	init_double(tfilLoc,Nfil_loc);
	init_double(tfilEst,Nfil_est);
	//exit(0);
	init_cplx(xf,N);
	init_cplx(ffilLoc,Nfil_loc);
	init_cplx(ffilEst,Nfil_est);
	show("\ntest:",xf);
	//x = (float*)calloc(N,sizeof(float));
	//xf = (float*)calloc(N,sizeof(float));
	generate_signal(x,xf,freqIdx);
	
	// Filter Design
	cheby_win(tfilLoc,ffilLoc,Nfil_loc,tolerance_loc);
	//cheby_win(tfilEst,ffilEst,Nfil_est,tolerance_est);
	
	printf("\nFilter:\n");
	for(i=0;i<N;i++)
	{
		if(x[i]!=0.0)
			printf("\n%d:%f",i,x[i]);
	}

	Filter f_loc,f_est;
	init_fil(&f_loc,Nfil_loc);
	init_fil(&f_est,Nfil_est);
	f_loc.sig_f = ffilLoc;
	f_loc.sig_t = tfilLoc;
	f_est.sig_f = ffilLoc;
	f_est.sig_t = tfilLoc;
	conv_window(&f_loc,Nfil_loc,b_loc);
	show("Floc:",f_loc.sig_f);
	//show("\nFFT:\n",tfilLoc);
	//double complex z = 2+3*I;
	//printf("%f+I*%f",creal(z),cimag(z));
	return 0;
}
