#include<stdio.h>
#include<math.h>
#include<complex.h>
#include"stdio.h"
#include"math.h"
#include<fftw3.h>
#define N 8192
#define k 10
#define BcstLoc 2
#define BcstEst 2
#define Combcst 32
#define combLoops 8
#define estLoops 16
#define locLoops 7
#define thresholdLoops 6
#define tolerance_loc 1e-8
#define tolerance_est 1e-8
#define pi 3.1415926535
typedef double complex cplx;

//////////////////////////////////////////////////////
//////////////Current FFT/////////////////////////////
/*
void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) 
	{
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);
 
		for (int i = 0; i < n; i += 2 * step) 
		{
			cplx t = cexp(-I*pi*i/n)*out[i + step];
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
*/
void show(const char * s, cplx buf[]) {
	printf("%s", s);
	for (int i = 0; i < 16; i++)
		if (!cimag(buf[i]))
			printf("%g ", creal(buf[i]));
		else
			printf("(%g + j%g) ", creal(buf[i]), cimag(buf[i]));
}
/*
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
	*/

void fftRoutineC2C(cplx *in,cplx *out,int len)
{
	fftw_plan p;
	p = fftw_plan_dft_1d(len,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void fftRoutineC2R(cplx *in,double *out,int len)
{
	fftw_plan p;
	p = fftw_plan_dft_c2r_1d(len,in,out,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

void fftRoutineR2C(float *in,cplx *out,int len)
{
	fftw_plan p;
	p = fftw_plan_dft_r2c_1d(len,in,out,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}
//////////////////////////////////////////////////////
////////// Filter Structure Description///////////////
typedef struct
{
	double *sigt;
	cplx *sigf;
	int numElements;
}Filter;

/////////////////////////////////////////////////////
////////////////UTILITY FUNCTIONS////////////////////

void generate_signal(double x[N],cplx xf[N],int fIdx[k]) // Generate signal
{
	unsigned int i,n,idx;
	for(i=0;i<k;i++)
	{
		idx = rand()%N;
		x[idx] = 1;
		fIdx[i] = idx;
		xf[idx] = 1;
	}	
	fftRoutineR2C(x,xf,N);
}

float sum_arr_cplx(cplx *arr,int num)
{
	cplx sum = 0;
	for(int i=0;i<num/4;i++)
		sum += arr[i]+arr[i+(int)(num/4.0)]+arr[i+(int)(num/2.0)]+arr[i+(int)(num*3/4.0)];
}

void normalize_max(cplx *tempx2,int len,int max)
{
	for(int i=0;i<len/2;i++)
	{
		tempx2[i] /= max;
		tempx2[i+len/2] /= max;
	}
}

void shift_left_cplx(cplx *arr,int len,int shift,int offset)
{
	memcpy(arr+len-offset,arr,shift*sizeof(cplx));
	arr = arr+shift; 
}

void shift_left_db(double *arr,int len,int shift,int offset)
{
	memcpy(arr+len-offset,arr,shift*sizeof(double));
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

int floor_to_pow2(double x)
{
  unsigned int ans;
  for (ans = 1; ans <= x; ans <<= 1) ;
  return ans / 2;
}

int modInverse(int a,int m)
{
	a = a%m; 
  for (int x=1; x<m; x++)
	{ 
		if ((a*x) % m == 1) 
    	return x;
	} 
}

void sqArray(int *a,int len)
{
	for(int i=0;i<len;i++)
		a[i] = pow(a[i],2); 
}

void findlargeIdx(cplx *arr,int *idxList,int num,int len)
{
	double absArr[len];
	for(int i=0;i<len;i++)
		absArr[i] = cabs(arr[i]);
	int tempList[num];
	int flag = 0;
	while(flag<len)
	{	
		int max = 0;
		int idx = 0;
		int cnt = 0;
		for(int j=0;j<len;j++)
		{
			if(max<=absArr[j])
			{
				max = absArr[j];
				tempList[cnt] = j;
				if(max==absArr[j])
					cnt++;
			}
		}
		for(int i=0;i<cnt;i++)
			idxList[flag+i] = tempList[i];
		flag+=cnt;
	}  
}
//values[a] = nthElement(values[a],location);

void nthElement(int arr[], int target)
{
	int high = loops
	int low = 0;
  if (high - low < target) 
	{
    // A single element is sorted
    return 0;
  }
  int mid = low + ((high - low) / 2);      
  // Put the pivot at the last
  int temp = arr[mid];
  arr[mid] = arr[high];
  arr[high] = temp;
      
  int pivot = arr[high];
  int i = low;
  int j = high - 1;
      
    // The partition
  while (j >= i) 
	{
      if (arr[j] < pivot && arr[i] > pivot) {
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;          
        ++i;
        --j;
  } else if (arr[j] < pivot && arr[i] <= pivot) {
            ++i;
        } else if (arr[j] >= pivot && arr[i] > pivot) {
            --j;
        } else if (arr[j] >= pivot && arr[i] <= pivot) {
            ++i;
            --j;
        }
    }
      
    // Bring the pivot back to its
    // appropriate place
  temp = arr[i];
  arr[i] = arr[high];
  arr[high] = temp;
      
  if (target >= i) {
      quickSelect(arr, i + 1, high, target);
  } else {
      quickSelect(arr, low, i, target);
  }
}
/*
void find_largest_indices(int *output, int num, real_t * samples, int n, real_t * tmp_storage)
{
  assert(n >= num + 1);
  //use num+1 so we can use > cutoff and probably get exactly num.
  //if we get fewer, the second pass uses == cutoff.
  real_t cutoff =
    nth_element_immutable(samples, n, n - num - 1, tmp_storage);

  int count = 0;
  for (int i = 0; i < n; i++)
    if (samples[i] > cutoff)
      output[count++] = i;
  if (count < num)
    {
      for (int i = 0; i < n; i++)
        {
          if (samples[i] == cutoff)
            {
              output[count++] = i;
              if (count >= num)
                break;
            }
        }
      std::sort(output, output + count);
    }
  assert(count == num);
}
*/
///////////////////////////////////////////////////////////
///////////////////Filter Functions////////////////////////

double chebyPoly(double n, double x)
{
	double res;
  if (fabs(x) <= 1)
    res = cos(n*acos(x));
  else
    res = creal(ccosh(n*cacosh(x)));
  return res;
}

void chebyWin(double *out,int Nfil,double tolerance)
{
	int _temp = next_pow2(Nfil);
	//cplx *temp = (cplx*)calloc(_temp,sizeof(cplx));
	double *temp = (double*)calloc(_temp,sizeof(double));
  int idx, i;
  //float M, n, sum = 0, max=0;
  double tg = 1/tolerance;//  /* 1/r term [2], 10^gamma [2] */
  double expr = cosh(acosh(tg)/(Nfil-1));
  //M = Nfil/2;
  printf("Nfil:%d x0:%lf _temp:%lu",Nfil,expr,_temp);
  for(idx=0; idx<Nfil; idx++)
  {
		temp[idx] = chebyPoly(Nfil-1,expr*cos((pi*idx)/Nfil))*tolerance;//cosf(2.0*n*pi*i/Nfil);
 		fout[idx] = temp[idx];
		//printf("Filter:%lf\n",creal(temp[idx]));
  }
	for(idx=0;idx<_temp-Nfil;idx++)
		temp[idx]=0;
		//printf("\nCount:%d\n",idx);
		/*
    for(idx=0; idx<Nfil; idx++)
      out[idx] /= max; // normalise everything 
		*/
	double fsig[Nfil];
	fftRoutineC2R(temp,out,Nfil);
	//show("\ntempt:",temp);
	shift_left_db(out,Nfil,Nfil/2,Nfil/2);
	//show("\ntempsht:",temp);		
	//for(idx=0;idx<Nfil;idx++)
	//	out[idx] = creal(temp[idx]);
}

void initFilter(Filter *fil,int arr_len)
{
	fil->sigt = (double*)calloc(arr_len,sizeof(double));
	fil->sigf = (cplx*)calloc(arr_len,sizeof(cplx));
	fil->numElements = arr_len;
}

void convWindow(Filter *fil,int nfil,int bin_size)
{
	//fil.sig_t = (float*)calloc(Nan,sizeof(float));
	//fil.sig_f = (cplx*)calloc(Nan,sizeof(cplx));
	cplx *tempx = (cplx*)calloc(N,sizeof(cplx));
	cplx *tempx2 = (cplx*)calloc(N,sizeof(cplx));
	init_cplx_from_db(tempx,N,fil->sigt,NULL);
	int width = N - nfil, bias = bin_size/2;
	int x=0;
	show("\nFFT first:\n",tempx);
	memcpy(tempx+nfil,&x,width*sizeof(cplx));	// Time shifting and padding the signal
	shift_left_cplx(tempx,nfil,nfil/2,0);
	show("\nNo FFT:\n",tempx);
	fftRoutineR2C(fil->sigt,tempx,N);
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
	fftRoutineC2R(tempx2,fil->sigt,nfil);
	//ifft_set(tempx2,nfil);
	show("\n\nTEMPX\n",fil->sigf);
	//exit(0);
	//free(tempx);
	//free(tempx2);
}

////////////////////////////////////////////////////////////
////////////////////Loop Functions//////////////////////////

//Locating the points
void innerloopLocate(double *origx,double *sampx,Filter *fil,int *J,int num,int B,int a,int ai,int b)
{
	if(N%B==1)
		printf("Warning N is not divisible by B!");
	else
	{
		sampx = (cplx*)calloc(B,sizeof(float));
		for(int i=0,idx = b;i<fil->numElements;i++)
		{
			sampx[i%B] = origx[idx]*(fil->sigt)[i];
			idx = (idx+ai)%N;
		}
		fftRoutineR2C(origx,sampx,B);//fft procedure
		findlargeIdx(sampx,J,num,B);
	}
}

//loop threshold
//Predict the regions
int innerloopFilter(int *J,int num,int B,int a,int ai,int b,int score[N],int hits[N],int hitsFound)
{
	for(int i=0;i<num;i++)
	{
		int low = ((int)(ceil((J[i]-.5)*N/B)+N))%N;
		int high = ((int)(ceil((J[i]+.5)*N/B)+N))%N;
		int loc = (low*a)%N;
		int j = low;
		while(j!=high)
		{
			score[loc]++;
			if(score[loc]==thresholdLoops)
			{
				hits[hitsFound] = loc;
				hitsFound++;
			}
			loc = (loc+a)%N;
			j = (j+1)%N;
		}
	}
	return hitsFound;
}

//Estimate the Values
void estValues(int *hits,int *hitsFound,int *xsamp,int loops,int *permute,int B2,int B,int N2,Filter *filLoc,Filter *filEst,cplx *ans)
{
	double *values[2];
	for(int i=0;i<2;i++)
		values[i] = (double*)calloc(loops,sizeof(double));
	for(int i=0;i<hitsFound;i++)
	{
		int pos = 0,currB = 0;
		Filter *fil;
		for(int j=0;j<loops;j++)
		{
			if(j<locLoops)
			{
				currB = B;
				fil = filLoc;
			}	
			else if(j>locLoops)
			{
				currB = B2;
				fil = filEst;
			}
			int factor = ((double)N)/currB;
			int pidx = (permute[j]*hits[i])%N;//permute_idx
			int hashed_to = pidx/(factor);
			int dist = pidx%factor;
			if(dist>(factor/2))
			{
				hashed_to = (hashed_to+1)%B;
				dist -= factor;
			}
			dist = (N-dist)%N;
			cplx filValue = fil->sigf[dist];
			values[0][pos] = xsamp[j][hashed_to]/creal(filValue);
			values[1][pos] = xsamp[j][hashed_to]/cimag(filValue);
			pos++;
		}
		location = (loops-1)/2;
		for(int a=0;a<2;a++)
			nthElement(values[a],location);
		double realv = values[0][location];
		double imagv = values[1][location];
		ans[i] = realv+I*imagv;
		idxList[i] = hits[i];
		//ans[hits[i]] = realv+I*imagv;
	}
}
////////////////////////////////////////////////////////////
////////////////////Host Code///////////////////////////////

int main()
{
	int *index,Arr,i,j,a,b,ai,bi;
	printf("N:%d,k:%d,(N*K)/log2(N):%d\n",N,k,(float)(N*k)/log2(N));
	double BBLoc = BcstLoc*sqrt(((float)(N*k))/log2(N));
	double BBEst = BcstEst*sqrt(((float)(N*k))/log2(N));
	double lobefracLoc = .5/BBLoc;
	double lobefracEst = .5/BBEst;
	int bLoc = floor(1.2*1.1*(N/BBLoc));
	int bEst = floor(1.4*1.1*(N/BBEst));
	
	int BLoc = floor_to_pow2(BBLoc);
	int BThresh = 2*k;
	int BEst = floor_to_pow2(BBEst);
	int loops = locLoops+estLoops;
	printf("BBLoc:%lf,BBEst:%lf,lobefracLoc:%lf",BBLoc,BBEst,lobefracLoc);
	//exit(0);
	int nfilLoc = floor((1/pi)*(1/lobefracLoc)*acosh(1.0/tolerance_loc));
	int nfilEst = floor((1/pi)*(1/lobefracEst)*acosh(1.0/tolerance_est));
	printf("%d:%d",nfilLoc,nfilEst);
	//exit(0);
	cplx sigHt[N],*sigDt,sigHf[N],*sigDf;
	Filter filLoc,filEst;
	//init_filter(filLoc);
	//sigHt = (cplx*)calloc(N,sizeof(cplx));
	//sigHf = (cplx*)calloc(N,sizeof(cplx));
	size_t size_sigDt = N*sizeof(cplx);
	cudaMalloc(&sigDt,size_sigDt);
	cudaMalloc(&sigDf,size_sigDt);
	initFilter(&filLoc,nfilLoc);
	initFilter(&filEst,nfilEst);
	chebyWin(filLoc.sigt,nfilLoc,tolerance_loc);
	//exit(0);
	chebyWin(filEst.sigt,nfilEst,tolerance_est);
	int permute[loops],permuteb[loops],score[N],hits[N],hitsFound;
	double *xsamp[locLoops];
	int *J;
	for(i=0;i<loops;i++)
		xsamp[i] = (float*)calloc(BLoc,sizeof(float));
	J = (int*)calloc(BThresh,sizeof(int));
	for(i=0;i<loops;i++)
	{
		a = 7773;
		b = 0;
		ai = modInverse(a,N);
		permute[i] = ai;
		permuteb[i] = b;
		int targetLoc;
		int J[BThresh];
		if(i<locLoops)
		{
			innerloopLocate(sigHt,xsamp[i],&filLoc,J,BThresh,BLoc,a,ai,b);
			hitsFound = innerloopFilter(J,BThresh,BLoc,a,ai,b,score,hits,hitsFound);
		}
		else if(i>=locLoops)
			innerLoopLocate(sigHt,xsamp[i],&filEst,J,BThresh,BLoc,a,ai,b);
	}
	cplx ans[hitsFound],idxList[hitsFound];
	estimateValues(hits,hitsFound,xsamp,loops,permute,BEst,BLoc,filLoc,filEst,ans,idxList);
	
	//Checking the result
	
	return 0;
}
