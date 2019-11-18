//Spi4.2. Use the parameters of the distribution as the input variables.
//Spi4.1. Select 64 sampling cells uniformly in the biofilm.
//Spi4. Use mean and variance as the inputs indtead of parameters of the distribution. Fix a bug in deciding the wave pattern: 
        //if phi is too small to find out 1 capable cell, the code will just use non-capable cell.
//Spi3.1. New algorithm to find brek up: N_firing/N_sample< 0.8;
//Spi3: antomatically catigrize the wave style: 0 died; 2, spiral wave+ bereak up; 3, break up, 5. normal; 8. secondary+Spiral; 9. secondary; 
//binary noise
//Spi2: radial lineage inheritance for any distributions
//_Lognor: use \nu and \sigma as the input arguments of lognormal distribution
//_Gau: parameter, for example tau is log normal distributed
//Spi ->Spi1: Eliminate parameter v_0;
//SpiWave_SS_WT.c is from Fit17Wave_SS_WT_spiral.c. Use square unit, rectangular lattice, reflective boundary in top and both sides, obsorbing in bottom. 
//Fit17: updated method to visulize wavelength.
//_SS: same tructure. In one repeating, all combinations of (a,b) will use the same structure.
//Fit16: replace b to tau. tau=1/b;
//FitClu15: dt_write=0.5; only consider these periods from firing capable cells; T_th=0.1;
//FitClu11->FitClu14: boundary condition changed. FitClu11: 4 absorbing boundaries (dxdt_Tri); FitClu14: obsorbing in the bottom, 3 reflective boundaries in the other three sides (dxdt_Tri_Top).
//FitClu11->FitClu14: Biofilm is by 200, which is located inside a 100 by 200 biofilm. row 0 is triggered, row 1 to row 35 is the observed area. The left rows is used to reduce the influence of bottom obsorbing boundary.
//11Loop_a_1: explore the parameter space of F-N model to explain the trkA, wild type and ktrA.
//FitClu10->11: Directionality or transmission speed as a function of phi; Trigger in the middle: row and column 48 49 50 51 52  ## Delete neighbor in SearchTree_G4
//FitClu9->10: introduce parameter c: 
//FitClu8->9: multi-generation inheritance to distribute the parameter.
//FitClu7->FitClu8: Triangle lattice takes the place of square lattice.
//FitClu6->FitClu7Introduce noise which follows Gaussian distribution or uniform distribution in log space.
//Introduce the parameter gamma: {a,b,gamma}: wt=b*(v-gamma*w);
//Attention: M, the number of columns; N, the number of rows. Size of matrix: N by M; 
//Delete the limit that v>=0;
//Change the model: When v decreased less than a/4, make w=0;
//The bug1: CluGen: Initiation of cluster,neighbor solved; 
//The original code, for all the c files in folder "Fitzhugh"
//aNoise contains the distribution of parameter a;

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define RAND_MAX 2147483647
#define PI 3.1415926
#define EPSILON 0.0000000001
#define EPSILON2 0.00001


//Use the equation dEdt=-delta(-D*deltaE);
//Add noise to explain the double spiking wich  occurses in the experiment of Gurol. Applied to 2D grid. Parallel the Runge Kunta part within while loop.
// Parallel Algorithm, run in cluster: No iteractive jobs.
//Trigger 9 cells in the ceter of cell community, use Runge Kutta Algorithm.
// Reproducing the work in Fig 3d and 3e in Pindle et al's Nature, 2015,"Ion channels enable electrical communication in bacterial communities" 
//May 11, 2016, Purdue; 
//2d N*M cells,Protocol b: Equilibrate for 100 min, then pulse then pulse S=100 uM ?or E=200*1000 uM?

void dxdt_Squ(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void dxdt_Squ_o4(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void dxdt_Tri(double *EVnST,double *deri,double *aNoise,double *tauNoise, double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);   //deri[] is used to store the data of dVdt,dndt......
void dxdt_Tri_Top(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void print1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt);
void print1_1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt);
int print5(int N, int M, double t, double *EVnST, FILE *fpt);
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri);     // addition of N*M matrix;
void Transport(int N,int M,int triy,int trix,int Ln,int Lm,int n,double *source,double *Tprint);
void Ini_Matri(int n,int m, double constant,double *matrix);
void Ini_Matri_I(int n,int m, int constant,int *matrix);
FILE * name(int N,int M,int *coefx, int *coefy,char *filename);
void Island(double *aNoise,double percentage,int N,int M,double lisland,double F);
int InOut(int target,int *source);
int CluGen(double *aNoise,double phi,double psi, int N,int M,double aLow,double DHigh);
int CluGen_G5_1(double *nvth,int height,int width,double phi,double rho,double DLow,double DHigh,
              double alpha, double beta, double a1, double a2, double a3);
int CluRadLin(double *nvth, int height, int width, char mode, double rho, double phi, double DLow, 
             double median, double lambda);
int mothergrandmother(char mode, double *nvth, int HeiMat, int WidMat, int height, double DLow, double DHigh,
           double *possibility, double ponoffoff, double ponoffon, double pononoff, double pononon);
double AssignPhi(double phi,double low,double high);
int SearchTree_G4(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DLow,int *cluster);
int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square);
double minimum(double a,double b);


int equal(double a, double b, double error);
int WheMultiple(double numerator, double denominator, double error);

int MeaLam(double *wave, int height, int width, int TimStep, double threshold);

int wavefront(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
int waveback(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
double wavelength(double *wave, int height, int width, int LenSec, double threshold);
int NFir(double *wave, int height, int width, double threshold);
double TMean(double *Tij, int height, int width, double threshold);

double AssignPhi(double phi,double low,double high);
double DisUni(double mean, double HalWid);
double gaussrand_K(void);

int FindPeaks(double deltat, double TimLen, double dt_write, double *trace, const int TimStep, double v_th, double *PeakTim, int LenPeakTime);
int WheDied(double * wave, int height, int width, int TimStep, double v_th);
int WheNeverDied(double *wave, int height, int width, int TimStep, double v_th);

int main(void){ 

//$$ 1. change the value of parameter a and b,c;
double aglobal, aLow, as[]={0};  ////phi: percentage of responsive cells; rho, clustering rate, aLow, the parameter of unresponsive cell.                                                 
double tau, tau_off, sigma, rho, rhos[]={0, 0.25, 0.5, 0.75};
double cglobal, cLow, cs[]={10};
double gamma=0;
double taus[]={2, 2.56,  3.11,  3.67,  4.22,  4.78,  5.33,  5.89,  6.44,  7}, sigmas[]={ 1, 1.33, 1.67,  2, 2.33, 2.67, 3, 3.33, 3.67,  4};

int  Len_taus=10, Len_sigmas=10, Len_rhos=4, NumRep=1;
int i,j, ia, ic, itau, isigma, irho, ip; //i, i is only used for any short, self-closed loop. 



//$$ 2. Number of row: N, column, M. Print type,1:all cells,5:row 0:5:30. 
    //Noise:0. no noise;1. clustering, bimodal distribution; 2.clustring, uniform distribution in log space; 250. other noise. 
    //Change phi,psi,aLow,tau_off;   $$$$$$$$$$$$$$$$*//
    //Change deltat, TimeLen, deltax, deltay;
int N=100,M=100, Leng=N*M, WinHei=35, print=1,coefx=1,coefy=1;   
 
//int WF_t, WF, WB_t;
//double WL_t, WL;
double v_th=0.6, T_th=0.1;
//double WF_ab[Len_as][Len_taus], WL_ab[Len_as][Len_taus], PhiFir_ab[Len_as][Len_taus], TraEff_ab[Len_as][Len_taus],
       //T_ab[Len_as][Len_taus], Tij[WinHei][M];   


double deltat=0.002,t=0.00,TimLen=400.00,deltax=1, deltay=1, dt_write=2;          //$$ TimLen incresed to 400 min. 
double EVnST[2][N*M],k0[2][N*M],EVnST1[2][N*M],k1[2][N*M],k2[2][N*M],k3[2][N*M],aNoise[N*M],tauNoise[N*M],cNoise[N*M];   //row 0 stores the v, row 1 stores w.

FILE *fptPhase, *fptWave, *fptMap;     //fpt: ThT; fpt2: parameter;
char filename[100], filename2[100];
char  num[100]; 
 time_t tt;
 srand((unsigned) time(&tt));  //Initiating random number generator only once. 

//fpt=fopen("Spi_T_100by100_Phi.45Rho.6.dat","w");
 //fpt1=fopen("Spi_sta_100by100_Phi.45Rho.6.dat","w");
 
     //Ini_Matri(Len_as, Len_taus,0.00, &WF_ab[0][0]);
     //Ini_Matri(Len_as, Len_taus,0.00, &WL_ab[0][0]);
     //Ini_Matri(Len_as, Len_taus,0.00, &PhiFir_ab[0][0]);
     //Ini_Matri(Len_as, Len_taus,0.00, &TraEff_ab[0][0]);
     //Ini_Matri(Len_as, Len_taus,0.00, &T_ab[0][0]);
const int TimStep=(int)(TimLen/dt_write+EPSILON2)+1; 
double wave[N*TimStep][M];  //Used only for test. Please
int PhaseDiagram[Len_taus][Len_sigmas];
fptPhase=fopen("Spi421_LN_Phase.dat","w");
   
  

for(ip=0;ip<NumRep; ip++){
for(irho=0; irho<Len_rhos; irho++){
	
  Ini_Matri_I(Len_taus, Len_sigmas, -1, &PhaseDiagram[0][0]);
for(itau=0;itau<Len_taus;itau++){	
for(isigma=0;isigma<Len_sigmas;isigma++){ 
	   
for(ia=0;ia<1;ia++){
for(ic=0;ic<1;ic++){
	
      
    aglobal=as[ia];
	aLow=as[ia];
	cglobal=cs[ic];
	cLow=cs[ic];	
	
	rho=rhos[irho]; 
	tau_off=5;	//Tau_off is not applied in LogNormal distribution
	tau=taus[itau]; 
	sigma=sigmas[isigma];
    
     
     strcpy(filename, "Spi421_LN"); 
     strcat(filename,"Rho");   sprintf(num, "%d",(int)(rho*100+EPSILON));strcat(filename,num); 
     strcat(filename,"u");   sprintf( num, "%d",(int)(taus[itau]*100+EPSILON) );    strcat(filename,num);
     strcat(filename,"s");   sprintf(num, "%d",(int)(sigmas[isigma]*100+EPSILON) );strcat(filename,num);
     
      //strcat(filename,"u");   sprintf(num, "%d",(int)(tau+EPSILON);strcat(filename,num);     
      //strcat(filename,"s");  sprintf( num, "%d",(int)(sigma) );strcat(filename,num);
      
    strcat(filename,"c");   sprintf(num, "%d",(int)(cglobal));strcat(filename,num);
    strcpy(filename2,filename);
    
    strcat(filename,"_wave.dat");	
	strcat(filename2,"_map.dat");
	
	fptWave=fopen(filename,"w");
	fptMap=fopen(filename2,"w");
	//'s', homogeneous, single value; 'b', binary distribution; 'N', LogNormal; 'U' Loguniform;
	CluRadLin(aNoise,   N, M, 's', -999, -999, -999, aglobal, -999);   //int CluRadLin(double *nvth, int height, int width, char mode, double rho, double phi, double DLow, double median, double lambda);
    CluRadLin(tauNoise, N, M, 'N', rho, -999, -999, tau, sigma); //-999 means that input is not necessary in the function.
    CluRadLin(cNoise,   N, M, 's', -999, -999, -999, cglobal, -999);  
    
    

	t=0;
	Ini_Matri(2,N*M,0.00,&EVnST[0][0]);        //Initiate the EVnST and t. I don't need to Iniate the EVnST1, k0,k1,k2,k3 and aNoise.    	  
 //$$ 4. Change the Trigger Time,length, and number of rows $$$$$$$$$$$$$$$$*//  
  while (t<TimLen) { 
	//if(fmod(t+EPSILON2,5)<10*EPSILON2)     //because t=5.0 or 10.0, the remainder might be 5.0, which might be due to some very small error, because t is a variable. 
		//  printf("time=%lf v0=%lf v1=%lf v9=%lf v99=%lf\n",t,EVnST[0][100],EVnST[0][301],EVnST[0][1909],EVnST[0][19999]);
	//MaxColumn=LargestCluster(tauNoise, N, M, N,M, 0, tau_off, tau, 5);
   
     if(fabs(t-0.00)<EPSILON2){     //Because t is a variable, it might have some little error.   
		for (i=N/2-3; i<=N/2+2; i++){
			for(j=M/2-3; j<=M/2+2; j++)
		       EVnST[0][i*M+j]=1;  	   
       		    }
       		    }     // Trigger: pulse v:1;        
	
    dxdt_Squ_o4(&EVnST[0][0],&k0[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay); 
	#pragma omp parallel for private(i)
	for(j=0;j<N*M;j++){
		for(i=0;i<2;i++)
		EVnST1[i][j]=EVnST[i][j]+deltat*k0[i][j];}    
		
		//if(t>=20-EPSILON2&&t<=20+0.006+EPSILON2)  {
		     //for(i=0;i<5;i++){
			 //for(j=0;j<5;j++)
				 //printf("%lf ",EVnST[0][i*M+j]); 
				 ////printf("\n");}

			     	 
		     //for(i=0;i<5;i++){
			 //for(j=0;j<5;j++)
				 //printf("%lf ",EVnST[1][i*M+j]);
				 //printf("\n");}
	
				 
		     //for(i=0;i<5;i++){
			 //for(j=0;j<5;j++)
				 //printf("%lf ",k0[0][i*M+j]);
				 //printf("\n");}
	
			//for(i=0;i<5;i++){
			 //for(j=0;j<5;j++)
				 //printf("%lf ",k0[1][i*M+j]);
				 //printf("\n");}   
				     //printf("\n");        }  
		
		
        
      dxdt_Squ_o4(&EVnST1[0][0],&k1[0][0],aNoise,tauNoise,cNoise,gamma,coefx,coefy,N,M,Leng,deltax,deltay);  
		#pragma omp parallel for private(i)
		for(j=0;j<N*M;j++){
			for(i=0;i<2;i++)
	       EVnST1[i][j]=EVnST[i][j]+deltat*k1[i][j]; }
     
        
       dxdt_Squ_o4(&EVnST1[0][0],&k2[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay);
        #pragma omp parallel for private(i)
        for(j=0;j<N*M;j++)  {
			for(i=0;i<2;i++)
            EVnST1[i][j]=EVnST[i][j]+2*deltat*k2[i][j];  }   //E1,V1... calculated by k0;       

       dxdt_Squ_o4(&EVnST1[0][0],&k3[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay); 
        #pragma omp parallel for private(i)
        for(j=0;j<N*M;j++)  {
			for(i=0;i<2;i++)
           EVnST[i][j]=EVnST[i][j]+deltat*(k0[i][j]+2*k1[i][j]+2*k2[i][j]+k3[i][j])/3;    }  //E1,V1... calculated by k0      
                
        t+=2*deltat;        // Runge Kutta Algorithm;  
        
        //if(equal(t,2*deltat,deltat/10)){
			//WB_t=waveback(&EVnST[0][0], N, M, 0, M-1, v_th);
			//WF_t=wavefront(&EVnST[0][0], N, M, 0, M-1, v_th);  
			//WL_t=wavelength(&EVnST[0][0], N, M, 25, v_th);
			
			//WF=WF_t;
			//WL=WL_t;
			
			//for (i=1;i<WinHei+1;i++)
				//for(j=0; j<M; j++)
					//Tij[i-1][j]=dt_write*(EVnST[0][i*M+j]> v_th);
					                 //}
        
       //else if(WheMultiple(t, dt_write, deltat/10)){
			//WB_t=waveback(&EVnST[0][0], N, M, 0, M-1, v_th);
			//WF_t=wavefront(&EVnST[0][0], N, M, 0, M-1,v_th);
			//WL_t=wavelength(&EVnST[0][0], N, M, 25, v_th);
			
			//WF=WF*(WF>WF_t)+ WF_t*(WF_t >= WF);
			//WL=WL*(WL>WL_t)+ WL_t*(WL_t>= WL);
			
			//for (i=1;i<WinHei+1;i++)
				//for(j=0; j<M; j++)
					//Tij[i-1][j]+= dt_write*(EVnST[0][i*M+j]> v_th);
			  //}
			  
		 if(ip==NumRep-1 && print==1)
         print1_1(N,M,Leng,t,dt_write,deltat, &EVnST[0][0],fptWave);
			  
	    if(WheMultiple(t, dt_write, deltat/10)|| equal(t,2*deltat,deltat/10) ){
			for (i=0; i<N; i++)
			  for (j=0;j<M;j++)
			   wave[i+ (int)(t/dt_write+EPSILON2) *N][j]=EVnST[0][i*M+j];}
   
             			                      
  }//end of while loop;
     //Find out 16 firing-capable cells
     int k, it, NumSample=64, LenPeakTime=10;
     double theta, trace[TimStep], PeakTim[NumSample][LenPeakTime];     
     double FirSamMean=0, FirSamMean_th=0.6, PeakMean=0, PeakMean_th=1.2;
     
     Ini_Matri(NumSample, LenPeakTime, TimLen+1, &PeakTim[0][0]);  //If no peak, the the peak time is TimLen+1;
       for (k=0; k<NumSample; k++){   //Find out NumSample capable sampling cells
		     
		 theta=k*2*PI/NumSample;
		 j=(int)( N/3*cos(theta)+(M-1)/2.0 +0.5);  
		 i=(int)( N/3*sin(theta)+(N-1)/2.0 +0.5);
				  
		   for (it=0; it< TimStep; it++)
		     trace[it]=wave[i+it*N][j];
		     
		   FindPeaks(deltat, TimLen, dt_write, trace, TimStep, v_th, &PeakTim[k][0], LenPeakTime);
		   } 
	 //for(i=0; i<NumSample; i++)
	   //PeakTimMean+=PeakTim[i][0];	   
	 //PeakTimMean/=NumSample;
	  //for(i=0; i<NumSample; i++)
	   //PeakTimDev+=(PeakTim[i][0]-PeakTimMean)*(PeakTim[i][0]-PeakTimMean);	 
	  //PeakTimDev=sqrt( PeakTimDev/ (NumSample-1) );
	  
	  for(i=0; i<NumSample; i++){
	   FirSamMean+= ( PeakTim[i][LenPeakTime-1]>= 1-EPSILON2 );   //Compare two double numbers, reduce the latter one by EPSILON2 to avoid the false negative.
	   PeakMean+=PeakTim[i][LenPeakTime-1];	 }  
	  FirSamMean/=NumSample;
	  PeakMean/=NumSample;
     
      if ( WheDied(&wave[0][0], N, M, TimStep, v_th) )
       PhaseDiagram[itau][isigma]=0;
     else if( FirSamMean< FirSamMean_th){
		 if( WheNeverDied(&wave[0][0], N, M, TimStep, v_th) )
		   PhaseDiagram[itau][isigma]=2;//Spiral-brek up
		 else
		   PhaseDiagram[itau][isigma]=3;}
     else if( PeakMean>= PeakMean_th-EPSILON2 ){
		 if( WheNeverDied(&wave[0][0], N, M, TimStep, v_th) )
		   PhaseDiagram[itau][isigma]=8;  //spiral_secondary
		 else
		   PhaseDiagram[itau][isigma]=9;}  //secondary
	 else
	    PhaseDiagram[itau][isigma]=5;     //normal
       
  
  
  //for (i=1;i<WinHei+1;i++)
				//for(j=0; j<M; j++)
					//Tij[i-1][j]= Tij[i-1][j]* (tauNoise[i*M+j]==tau); 
  
  //WF_ab[ia][itau]+=WF;
  //WL_ab[ia][itau]+=WL;
  //PhiFir_ab[ia][itau]+= (double)NFir(& Tij[0][0], WinHei, M, T_th)/WinHei/M;
  
  //if (NFir(&Tij[0][0], 5, M, T_th)==0)
  //TraEff_ab[ia][itau]+=0;
  //else  
  //TraEff_ab[ia][itau]+= (double)NFir(&Tij[30][0], 5, M, T_th)/(double)NFir(&Tij[0][0], 5, M, T_th);
  
  //T_ab[ia][itau]+= TMean(& Tij[0][0], WinHei, M, T_th);
  
   
  
 if(ip==NumRep-1){
	 
		printf("mean=%lf, variance=%lf, tau=%lf tau_off %lf sigma=%lf rho=%lf  c=%lf, cLow=%lf\n",taus[itau], sigmas[isigma], tau, tau_off, sigma, rho,  cglobal, cLow);
		for (i=0;i<N*M;i++){ 	
				if(i%M==M-1)
		        fprintf(fptMap,"%lf\n",aNoise[i]);	
		        else
				fprintf(fptMap,"%lf  ",aNoise[i]);}
         for (i=0;i<N*M;i++){ 	
				if(i%M==M-1)
		        fprintf(fptMap,"%lf\n",tauNoise[i]);	
		        else
				fprintf(fptMap,"%lf  ",tauNoise[i]);}
		 for (i=0;i<N*M;i++){ 	
				if(i%M==M-1)
		        fprintf(fptMap,"%lf\n",cNoise[i]);	
		        else
				fprintf(fptMap,"%lf  ",cNoise[i]);}
	 fclose(fptWave);
	 fclose(fptMap);
				 
				}
				
   
   
}   //end of loop ic
}   //end of loop ia

}  //end of isigma
} //end of loop itau

//for(ia=0; ia<Len_as; ia++){ 
//for(itau=0;itau<Len_taus;itau++){ 	
				//if(itau==Len_taus-1)
		        //fprintf(fpt1,"%lf\n",WF_ab[ia][itau]/NumRep);	
		        //else
				//fprintf(fpt1,"%lf ",WF_ab[ia][itau]/NumRep);} }
				
//for(ia=0; ia<Len_as; ia++){ 
//for(itau=0;itau<Len_taus;itau++){ 	
				//if(itau==Len_taus-1)
		        //fprintf(fpt1,"%lf\n",WL_ab[ia][itau]/NumRep);	
		        //else
				//fprintf(fpt1,"%lf ",WL_ab[ia][itau]/NumRep);} }

//for(ia=0; ia<Len_as; ia++){ 
//for(itau=0;itau<Len_taus;itau++){ 	
				//if(itau==Len_taus-1)
		        //fprintf(fpt1,"%lf\n",PhiFir_ab[ia][itau]/NumRep);	
		        //else
				//fprintf(fpt1,"%lf ",PhiFir_ab[ia][itau]/NumRep);} }
	
//for(ia=0; ia<Len_as; ia++){ 
//for(itau=0;itau<Len_taus;itau++){ 	
				//if(itau==Len_taus-1)
		        //fprintf(fpt1,"%lf\n",TraEff_ab[ia][itau]/NumRep);	
		        //else
				//fprintf(fpt1,"%lf ",TraEff_ab[ia][itau]/NumRep);} }
				
//for(ia=0; ia<Len_as; ia++){ 
//for(itau=0;itau<Len_taus;itau++){ 	
				//if(itau==Len_taus-1)
		        //fprintf(fpt1,"%lf\n",T_ab[ia][itau]/NumRep);	
		        //else
				//fprintf(fpt1,"%lf ",T_ab[ia][itau]/NumRep);} }
				
 for(i=0; i<Len_taus; i++){ 
 for(j=0;j<Len_sigmas;j++)	
	fprintf(fptPhase,"%d  ",PhaseDiagram[i][j]);
	fprintf(fptPhase, "\n");	
		         }
		         
} //end of loop irho		
}   //end of loop ip

fclose(fptPhase);
				
//fclose(fpt);
//fclose(fpt1);


printf("running this tourine costs %ld seconds of time\n",time(NULL)-tt);
return 0;	
}

//Reflective boundary in the top and two sides, observing in the bottom.
void dxdt_Squ(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
 for(i=0;i<N*M;i++){  
  deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                       //dwdt;   
   deri[i]=cNoise[i]*(  EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i] );   //dvdt: This step also initiates the deri matrix.
 
    
  
  if(i/M==0)
	 deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2);
  else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
  else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
     
   if(i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2);
   else if(i%M==M-1)
      deri[i]+=(EVnST[i-1]-EVnST[i])/pow(deltax,2); 
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)+(EVnST[i-1]-EVnST[i])/pow(deltax,2);  	     
 }     //end of for loop
 } // end of parallel
}

void dxdt_Squ_o4(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
 for(i=0;i<N*M;i++){  
  deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                       //dwdt;   
   deri[i]=cNoise[i]*(  EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i] );   //dvdt: This step also initiates the deri matrix.
 
    
  
  if(i/M==0)
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i+M]-EVnST[i])/pow(deltay,2);
  else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
  else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
     
   if(i%M==0)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)+(EVnST[i+1]-EVnST[i])/pow(deltax,2);
   else if(i%M==M-1)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)+(EVnST[i-1]-EVnST[i])/pow(deltax,2); 
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)+(EVnST[i-1]-EVnST[i])/pow(deltax,2);  	     
 }     //end of for loop
 } // end of parallel
}

void dxdt_Tri(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
for(i=0;i<N*M;i++){  
   deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                      //dwdt;
   
   deri[i]=cNoise[i]*(EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i]);   //dvdt: This step also initiates the deri matrix.
   
   if(i/M==0)    //v_yy
	 deri[i]+=(EVnST[i+M]-2*EVnST[i])/pow(deltay,2);
   else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
   else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2); 	 
    
  if((i%M)%2==0){ 
	       
   if(i/M==N-1&&i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==N-1)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M-1]-EVnST[i])/pow(deltax,2)/2;
            }
            
            
   else if((i%M)%2==1){
	       
   if(i/M==0&&i%M==M-1)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==M-1)
      deri[i]+=(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
            }     
}     //end of for loop
 } // end of parallel
}

void dxdt_Tri_Top(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
for(i=0;i<N*M;i++){     
   deri[Leng+i]=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                      //dwdt;
   
   deri[i]=cNoise[i]*(EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i]);   //dvdt: This step also initiates the deri matrix.
   
   if(i/M==0)
	 deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2);   //Reflective on the top
   else if(i/M==N-1)	  
	 deri[i]+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);  //Obsorbing on the bottom
   else
     deri[i]+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2); 	 
    
  if((i%M)%2==0){ 
	       
   if(i/M==N-1&&i%M==0)   //Last row, first column
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==N-1)  //Last row
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2+(0-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i+M-1]-EVnST[i])/pow(deltax,2)/2;
            }
            
            
   else if((i%M)%2==1){  //Odd columns
	       
   if(i/M==0&&i%M==M-1)   //First row, last column
      deri[i]+=(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2;
   else if(i/M==0)
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2;
   else if(i%M==M-1)
      deri[i]+=(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
   else
      deri[i]+=(EVnST[i+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M+1]-EVnST[i])/pow(deltax,2)/2+(EVnST[i-M-1]-EVnST[i])/pow(deltax,2)/2;
            }     
}     //end of for loop
 } // end of parallel
}

double minimum(double a,double b)
{if (a<=b)
	return a;
else 
	return b;	
}


void print1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt)
{         //print the data of every cell into file.
   
    int j=0; 	      
    	      
    if(WheMultiple(t,tInterval,deltat/10)|| equal(t,0.5,deltat/10)){   
        	for(j=0;j<N*M;j++){
        		if(j%M==M-1)
        		   fprintf(fpt,"%f\n",*(EVnST+j));	
        		else
		           fprintf(fpt,"%f  ",*(EVnST+j));
	        }	          
      	
        } 
         
}  

void print1_1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt)
{         //print the data of every cell into file.
   
    int j=0; 	      
    	      
    if(WheMultiple(t,tInterval,deltat/10)|| equal(t,2*deltat,deltat/10)){   
        	for(j=0;j<N*M;j++){
        		if(j%M==M-1)
        		   fprintf(fpt,"%f\n",*(EVnST+j));	
        		else
		           fprintf(fpt,"%f  ",*(EVnST+j));
	        }	          
      	
        } 
         
} 

int print5(int N, int M, double t, double *EVnST, FILE *fpt)   //print row 1:5:201;
{  int i=0,j=0;
	if(N!=1&&M!=1&&fabs(fmod(t-0.00005,0.5)-0.5)<0.00008){   
        	for(i=0;i<N;i=i+5){
			for(j=0;j<M;j=j+1){
				if(j==M-1)
        		   fprintf(fpt,"%lf\n",*(EVnST+i*M+j));	
        		else
		           fprintf(fpt,"%lf  ",*(EVnST+i*M+j));
	          }}          
    }
return 1;
}

//Define the addition of Matrix.
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri)
{    // addition of N*M matrix;
   int i,j;
   for(i=0;i<n;i++)
       for(j=0;j<m;j++)
           *(Matout+i*Length+j)=*(Matin+i*Length+j)+*(Deri+i*Length+j)*deltat;
    }
//function Transport: write the data in aNoise or EVnST into Tprint
void Transport(int N,int M,int triy,int trix,int Ln,int Lm,int n,double *source,double *Tprint)
{
	int i=0,j=0,m=0;
	for(i=(N-1)/2-(Ln-1)/2;i<=(N-1)/2+(Ln-1)/2;i++)      //the triy  by trix grid in the left side.  Write the data to the Tprint(m,n); source==&EVnST[4][0]or &aNoise[0][0]
		for(j=(M-1)/4-(Lm-1)/2;j<=(M-1)/4+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);	
      for(i=(N-1)/2-(Ln-1)/2;i<=(N-1)/2+(Ln-1)/2;i++)      //the triy  by trix grid in the right side.
		for(j=3*(M-1)/4-(Lm-1)/2;j<=3*(M-1)/4+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);
		for(i=(N-1)/4-(Ln-1)/2;i<=(N-1)/4+(Ln-1)/2;i++)      //the triy  by trix grid in the top side.
	       for(j=(M-1)/2-(Lm-1)/2;j<=(M-1)/2+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);
	  	for(i=3*(N-1)/4-(Ln-1)/2;i<=3*(N-1)/4+(Ln-1)/2;i++)      //the triy  by trix grid in the bottom side.
	       for(j=(M-1)/2-(Lm-1)/2;j<=(M-1)/2+(Lm-1)/2;j++)
		  *(Tprint+400*(m++)+n)=*(source+i*M+j);
}

//Ini_Matri: make every element of the matrix to be a constant.
void Ini_Matri(int n,int m, double constant,double *matrix)
{int i=0,j=0;
   for(i=0;i<n;i++)
	   for(j=0;j<m;j++)
		   *(matrix+i*m+j)=constant;
	
}

void Ini_Matri_I(int n,int m, int constant,int *matrix)
{int i=0,j=0;
   for(i=0;i<n;i++)
	   for(j=0;j<m;j++)
		   *(matrix+i*m+j)=constant;
	
}
//Define the function, name to build a new file.
FILE * name(int N,int M,int *coefx, int *coefy,char *filename)
{
static FILE * fpt;
 if(N==1&&M==1){     
      *coefy=0;
      *coefx=0;
	   fpt=fopen(filename,"w"); 
     }
    else if(N==1&&M!=1) {
      *coefx=1;
	  *coefy=0;
	  fpt=fopen(filename,"w");
      }
      
      else if(N!=1&&M!=1){
      *coefx=1;
      *coefy=1;
	   fpt=fopen(filename,"w");}
   
   	 else if(N!=1&&M==1){
     *coefx=0;
	  *coefy=1;
	  fpt=fopen(filename,"w");
	  }
	  else
  		printf("Wrong");
  		return fpt; 
}

void Island(double *aNoise,double percentage,int N,int M,double lisland,double F)
   {
	//int gk=30; //min^-1
	//double gl=0.2; //min^-1
	//int vk0=-380; //mV
	//int vl0=-156; //mv
	//int sth=40;   //uM
	 int vth=-150; //mV
	//int alpha0=2;  //min^-1
	//double beta=1.3;  //min^-1
	//int mglobal=1;
	//int Fglobal=5600;  //uM/mV
	//double sigma=0.2;  //mV
	//double deltak=0.001; //mV/uM;
	//double deltal=0.008; //mV/uM;
	//double gammas=0.1;  //min^-1
	//int gammae=10;   //min^-1
	//int gammat=4;  //min^-1
	//int alphas=1;  //uM/min/mV
	//int alphat=1;   //uM/min/mV
	//int Dglobal=82800;   //1380*60um^2/min
	int mark[201][201],irdom=0,jrdom=0,badposition=0,Nrdom=0,i=0,j=0,Ntotal=0;
	double Temvth,cos,sin,width,wij,Rij0,kwidth=0.2;
	time_t tij;                  
	srand((unsigned) time(&tij));
	width=lisland/3;
	for(i=0;i<N;i++)
      for(j=0;j<M;j++)
        mark[i][j]=0;		  
	
	
	while(Nrdom<percentage*N*M&&Ntotal<=1000000){   // avoid infinite loop
		Ntotal=Ntotal+1;
		 
		if(Ntotal%100000==0) 
		printf("%d   %d   %lf\n",Ntotal,Nrdom,lisland);   
		
		
        irdom=rand()%(N-(int)lisland/3*4)+(int)lisland/3*2;   //Make sure there is engouth space to form a island with center in this point.
        jrdom=rand()%(M-(int)lisland/3*4)+(int)lisland/3*2;   //If lisland is even, then there are some small probblems: the distribution of irdom and jrdom will not be symmetric.
		Rij0=sqrt(pow(irdom-(N-1)/2,2)+pow(jrdom-(M-1)/2,2));
		
		if(Rij0>=lisland/2){    //Dont't want the island cross the center;		
		cos=(irdom-(N-1)/2)/Rij0;
		sin=(jrdom-(M-1)/2)/Rij0;
		badposition=0;
		for(i=irdom-(int)lisland/3*2;i<=irdom+(int)lisland/3*2;i++){
			for(j=jrdom-(int)lisland/3*2;j<=jrdom+(int)lisland/3*2;j++){
			wij=kwidth*((i-irdom)*cos+(j-jrdom)*sin)+width;   // write wij as a function of radial distance between two points.
			if(mark[i][j]==1&&fabs((i-irdom)*cos+(j-jrdom)*sin)<=(lisland)/2.0+0.00001&&fabs((i-irdom)*sin-(j-jrdom)*cos)<=(wij)/2+0.00001){
			badposition=1;
            break;	
			}   
			}
		    if(badposition==1)
		    break;
			}
		
		if(badposition==0){
			Temvth=-exp((double)rand()/2147483647.0*2*log(F)-log(F)+log(-vth)); 
			for(i=irdom-(int)lisland/3*2;i<=irdom+(int)lisland/3*2;i++){
			for(j=jrdom-(int)lisland/3*2;j<=jrdom+(int)lisland/3*2;j++){
			wij=kwidth*((i-irdom)*cos+(j-jrdom)*sin)+width;   // write wij as a function of R
			if(fabs((i-irdom)*cos+(j-jrdom)*sin)<=(lisland)/2+0.00001&&fabs((i-irdom)*sin-(j-jrdom)*cos)<=(wij)/2+0.00001){
			mark[i][j]=1;
            Nrdom++;
            *(aNoise+i*M+j)=Temvth;	
			}	
			}}
		}  //end of giving value to island.
	   }
	}  //end of while loop.
	
	
}

int CluGen(double *aNoise,double phi,double psi, int N,int M,double aLow,double DHigh)
{int cluster[N*M],neighbor[N*M],NClu,NNei,NNeiE,NHigh=0,i;
	for(i=0;i<N*M;i++)  cluster[i]=-1;       //Initiate it before the first use.
	for(i=0;i<N*M;i++)  neighbor[i]=-1;
	
  while(NHigh<phi*N*M){         //compare two doubles: Be careful for possible calculation errors of computer.
	 i=0;while(cluster[i]!=-1&&i<N*M)	  //Initiate the cluster. -1 is the end of effective part. -2 means that position is concealed.
		  cluster[i++]=-1;
     i=0; while(neighbor[i]!=-1&&i<N*M)
		neighbor[i++]=-1;
		 NClu=0;
		 NNei=0;       //The length of array neighbor;
		 NNeiE=0;       //The number of neighbor cells.
	 
	  do{
		cluster[NClu]=(rand()%N)*M+rand()%M;  
	  } while(*(aNoise+cluster[NClu])>aLow);   
	  
		*(aNoise+cluster[NClu++])=DHigh;
	     NHigh+=1;		
		
		while((double)rand()/(RAND_MAX+0.1)<psi&&NHigh<phi*N*M){     
		  for(i=NClu-1;i<NClu;i++){
			  if(*(aNoise+cluster[i]-1)<=aLow&&cluster[i]%M!=0&&InOut(cluster[i]-1,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]-1;
			  NNeiE++;}
			  if(*(aNoise+cluster[i]+1)<=aLow&&cluster[i]%M!=M-1&&InOut(cluster[i]+1,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]+1;
			  NNeiE++;}
			  if(*(aNoise+cluster[i]-M)<=aLow&&cluster[i]/M!=0&&InOut(cluster[i]-M,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]-M;
		       NNeiE++;}
			  if(*(aNoise+cluster[i]+M)<=aLow&&cluster[i]/M!=N-1&&InOut(cluster[i]+M,&neighbor[0]))
			  {neighbor[NNei++]=cluster[i]+M;
		        NNeiE++;}			  
			  }   
		  
		   if(NHigh>=N*M-3)  printf("%d %d %d %d\n",NHigh,NClu,NNei,NNeiE);//$$BreakPoint     
		  
		       if(NNeiE==0)
		        break;
			  
 			     i=0;
		         while(1){
				   if(neighbor[i]!=-2&&(double)rand()/(RAND_MAX+0.1)<1.0/(double)NNeiE){
				    *(aNoise+neighbor[i])=DHigh;
					cluster[NClu++]=neighbor[i];
				    NHigh++;
					neighbor[i]=-2;      //canceal this selected point from neighbor.
					NNeiE=NNeiE-1;
					break;}
					 if(neighbor[++i]==-1)
				     i=0;}
				     
		 
			    
		}  //end of while((double)rand()/2147483647.0<=psi&&space==1);	end of generating one cluster.  
  }  //end of while(NHigh<phi*N*M);  completing all clusters.
  return 1;
}

int CluGen_G5_1(double *nvth,int height,int width,double phi,double rho,double DLow,double DHigh,
              double alpha, double beta, double a1, double a2, double a3){

 int i,j,m, HeiMat=height+50, WidMat=width;
double pon, poff, ponon,ponoff, poffon, poffoff, pononon, pononoff, ponoffon, ponoffoff, possibility[width],Prenvth[HeiMat*WidMat];

    pon=phi;
    poff=1-pon;
	ponon=phi/(1-rho+rho*phi);
	poffon=1-ponon;	
	ponoff=ponon*(1-rho);
	poffoff=1-ponoff;
	
	pononoff=ponon/(alpha*ponon+ponoff*poff/pon);
	pononon=alpha*pononoff;
	ponoffoff=ponoff/(beta*poffon*pon/poff+poffoff);
	ponoffon=beta*ponoffoff;

  for(i=0;i<HeiMat;i++){
	if(i==0){
		for(j=0;j<width;j++)
		Prenvth[j]=AssignPhi(phi,DLow,DHigh);
	 }
				
	else if(i==1){			
		for(j=0;j<width;j++){
			m=i*WidMat+j;
			if(Prenvth[m-WidMat]==DHigh)   
			 Prenvth[m]=AssignPhi(ponon,DLow,DHigh);				
			else
			 Prenvth[m]=AssignPhi(ponoff,DLow,DHigh);    
		      }
	}
	
	else{
		mothergrandmother('O', Prenvth, HeiMat, WidMat, i,DLow, DHigh, possibility, ponoffoff, ponoffon, pononoff, pononon);
		
		for(j=1;j<width;j=j+2){
			m=i*WidMat+j;
			if(j==0)
			Prenvth[m]=AssignPhi((a1+a2)*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			else if(j==width-1)
			Prenvth[m]=AssignPhi((a2+a3)*possibility[j]+a1*possibility[j-1],DLow, DHigh);
			else
			Prenvth[m]=AssignPhi(a1*possibility[j-1]+a2*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			}
			
		mothergrandmother('E', Prenvth, HeiMat, WidMat, i,DLow, DHigh, possibility, ponoffoff, ponoffon, pononoff, pononon);
		for(j=0;j<width;j=j+2){
			m=i*WidMat+j;
			if(j==0)
			Prenvth[m]=AssignPhi((a1+a2)*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			else if(j==width-1)
			Prenvth[m]=AssignPhi((a2+a3)*possibility[j]+a1*possibility[j-1],DLow, DHigh);
			else
			Prenvth[m]=AssignPhi(a1*possibility[j-1]+a2*possibility[j]+a3*possibility[j+1],DLow, DHigh);
			}			
		}
	}
	
	
	for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			m=i*width+j;
	        nvth[m]=Prenvth[(i+HeiMat-height)*WidMat+j]; } }
		
	return 0;
 }
 
int CluRadLin(double *nvth, int height, int width, char mode, double rho, double phi, double DLow, 
             double median, double lambda){
				 
int cluster[height*width], state[height*width], NClu=0, NClu0, NEffClu=0, NEffClu0, NMother=0, daughter,mother, StoLoc, ABLR,m, i,j;
	Ini_Matri_I(height, width, -999, &cluster[0]);
	Ini_Matri_I(height, width, -999, &state[0]);
	Ini_Matri(height, width, -999, nvth);				 

	
 if(mode=='s'){
	 for (m=0; m<height*width; m++)
	    nvth[m]=median;          
	          }
 else{
	 
	//Initiation: seeding with 6*6 grid of cells, same as the trigger, drawn independently from p(m).
	for (i=height/2-3; i<=height/2+2; i++){
			for(j=width/2-3; j<=width/2+2; j++){
	mother=i*width+j;
	if (mode=='b')
	  nvth[mother]=AssignPhi(phi, DLow, median);
	else if(mode=='N')
	  nvth[mother]=exp( gaussrand_K()*lambda+ median );
	else if(mode=='U')
	  nvth[mother]=exp( DisUni( median, lambda ) );
	 cluster[NClu]=mother;
	 state[NClu++]=0;
	 NEffClu++;
           }
    }
    NClu0=NClu;
	NEffClu0=NEffClu;

	do{  //produce the 2nd, 3rd... generations until the height by width box is full.
		
		while(NMother<NEffClu0)  //produce the next generation from the mothers.  
		{ 
			do{StoLoc=rand()%NClu0;  //Randomly select one mother which has unoccupied neighbors.
			  } while(state[StoLoc]!=0);
			mother=cluster[StoLoc];
			NMother++;
			//Decide whether this mother have unoccupied neighbors to give birth to one daughter
			if( (mother/width!=0&& nvth[mother-width]==-999) || (mother/width!=(height-1) && nvth[mother+width]==-999) || 
			    (mother%width!=0&& nvth[mother-1]==-999) || (mother%width!=(width-1) && nvth[mother+1]==-999) ){
					 do{ABLR=rand()%4;
					    daughter=(mother-width)*(ABLR==0)+(mother+width)*(ABLR==1)+(mother-1)*(ABLR==2)+(mother+1)*(ABLR==3);
					   }  while( ((mother/width==0)&&(ABLR==0))|| ((mother/width==(height-1))&&(ABLR==1)) || ((mother%width==0)&&(ABLR==2)) 
					             || ((mother%width==(width-1))&&(ABLR==3)) || (nvth[daughter]!=-999) );
					 //Assign a number to the daughter: has probability of rho to be the same with the mother; 1-rho to be uncorrelated.
				if(rand()/(RAND_MAX+0.1)<rho)
				   nvth[daughter]=nvth[mother];
				else{ 	
					 if (mode=='b')
						  nvth[daughter]=AssignPhi(phi, DLow, median);
						else if(mode=='N')
						  nvth[daughter]=exp( gaussrand_K()*lambda+ median );
						else if(mode=='U')
						  nvth[daughter]=exp( DisUni( median,lambda ) ); 
				    }
				     //Put daughter into cluster, increase the NClu and NEffClu	    
			    cluster[NClu]=daughter;
			    state[NClu++]=0;
			    NEffClu++;
			    state[StoLoc]=1;  //mark that this mother has generated daughter 
					}
			else  {        //Mark the mother as unable to give birth to daughters.
				state[StoLoc]=-1;
		        NEffClu--;       }   
		}
		//printf("NClu %d, NClu0 %d, NEffClu %d, NEffClu0 %d, NMother %d\n",NClu, NClu0, NEffClu, NEffClu0, NMother);	
		//for (m=0;m<NClu;m++)		  
		//  printf("Storage location %d, mother %d, value of mother %lf, state of mother %d\n",m,cluster[m], nvth[cluster[m]], state[m]);
	      //if(NClu==height*width){		
	    	//exit(0);}
		//Recover the variables to be ready for the production of the next generation.
		NMother=0;
		for (m=0;m<NClu0;m++)
		  state[m]=0*(state[m]==1)-1*(state[m]==-1);
		NEffClu0=NEffClu;
		NClu0=NClu;
			
			
		
	  } while (NEffClu>0);    
 } //end of mode 'b','U','N'
	
	return 0;	 
}
 
 
 int mothergrandmother(char mode, double *nvth, int HeiMat, int WidMat, int height, double DLow, double DHigh,
           double *possibility, double ponoffoff, double ponoffon, double pononoff, double pononon){
int j,m;
if (mode=='O'){
	for(j=0;j<WidMat;j++){
	 m=height*WidMat+j;
	 if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=ponoffoff;
	 else if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DHigh)
	   possibility[j]=ponoffon;
	 else if(nvth[m-WidMat]==DHigh&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=pononoff;
	 else 
	   possibility[j]=pononon;
	 }
	 }
	 
else{
  for(j=0;j<WidMat;j=j+2){
	 m=height*WidMat+j;
	 if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=ponoffoff;
	 else if(nvth[m-WidMat]==DLow&& nvth[m-2*WidMat]==DHigh)
	   possibility[j]=ponoffon;
	 else if(nvth[m-WidMat]==DHigh&& nvth[m-2*WidMat]==DLow)
	   possibility[j]=pononoff;
	 else 
	   possibility[j]=pononon;
	 }
	 
  for(j=1;j<WidMat;j=j+2){
	 m=height*WidMat+j;
	 if(nvth[m]==DLow&& nvth[m-WidMat]==DLow)
	   possibility[j]=ponoffoff;
	 else if(nvth[m]==DLow&& nvth[m-WidMat]==DHigh)
	   possibility[j]=ponoffon;
	 else if(nvth[m]==DHigh&& nvth[m-WidMat]==DLow)
	   possibility[j]=pononoff;
	 else 
	   possibility[j]=pononon;
	 }
	 }
	return 0;
}

double AssignPhi(double phi,double low,double high){
	         if((double)rand()/(RAND_MAX+0.1)<phi)
			    return high;
			    else
			    return low;
}

int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square){
	int j,state[HeiMat][WidMat],cluster[HeiMat][WidMat], clustersize[WidMat];
	int MaxColumn=0;
	
	Ini_Matri_I(HeiMat,WidMat, 0, &state[0][0]);
	Ini_Matri_I(HeiMat,WidMat, 0, &cluster[0][0]);
	 for(j=0;j<=width-square;j++){
	  clustersize[j]=SearchTree_G4(nvth, &state[0][0],row,j,height, width, HeiMat, WidMat,
                   DHigh, &cluster[0][0]);
      if (j==0)       
       MaxColumn=0;
       else if(clustersize[j]>clustersize[MaxColumn])
       MaxColumn=j;       
             }
     
    return MaxColumn;              
      
	   }
	   
int SearchTree_G4(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DHigh,int *cluster)  //il is the dimension of effective area of nvth and state;
{  int column=i*WidMat+j,i1,j1,m,p;  //i1: the No. of row; j1, the No. of column; m, the No. in the srorage, m= i1*WidMat+j1;
	if(state[column]==0&& nvth[column]==DHigh){
		  
	     int NClu0=0,NClu=0,NNei=0;  //The length of array neighbor;    //The number of effictive.neighbor cells.       	
		  cluster[NClu++]=column;
		  NNei+=1;
			state[column]=2;     //for(p=0;p<NClu;p++) printf("%d ",cluster[p]); printf("\n");
		     
		   while(1){
		   for(p=NClu0;p<NClu;p++){    //look for the neighbors of these new cells in cluster. 
			   m=cluster[p];
			   i1=m/WidMat;
			   j1=m%WidMat;
			   
		  if(nvth[m-1]==DHigh&& state[m-1]==0&& j1!=0)
		  {cluster[NNei++]=m-1;	
			 state[m-1]=2; }
			 
		  if(nvth[m+1]==DHigh&& state[m+1]==0&& j1!=width-1)
		  {cluster[NNei++]=m+1;	
			 state[m+1]=2; }
			 
		 if(nvth[m-WidMat]==DHigh&& state[m-WidMat]==0&& i1!=0)
		  {cluster[NNei++]=m-WidMat;	
			 state[m-WidMat]=2; }	
			 
		 if(nvth[m+WidMat]==DHigh&& state[m+WidMat]==0&& i1!=height-1)
		  {cluster[NNei++]=m+WidMat;	
			 state[m+WidMat]=2; }	
				   
		 if(m%2==0){
			 if(nvth[m+WidMat-1]==DHigh&& state[m+WidMat-1]==0&& i1!=height-1&& j1!=0)
		     {cluster[NNei++]=m+WidMat-1;	
			 state[m+WidMat-1]=2; }
			 
			 if(nvth[m+WidMat+1]==DHigh&& state[m+WidMat+1]==0&& i1!=height-1&& j1!=width-1)
		     {cluster[NNei++]=m+WidMat+1;	
			 state[m+WidMat+1]=2; }	   }
			 
			else{
			if(nvth[m-WidMat-1]==DHigh&& state[m-WidMat-1]==0&& i1!=0&& j1!=0)
		     {cluster[NNei++]=m-WidMat-1;	
			 state[m-WidMat-1]=2; }
			 
			 if(nvth[m-WidMat+1]==DHigh&& state[m-WidMat+1]==0&& i1!=0&& j1!=width-1)
		     {cluster[NNei++]=m-WidMat+1;	
			 state[m-WidMat+1]=2; }			}	
		                            } 
		                                          		       
		          if(NNei==NClu)			   		     
				  return NClu;
								
 					  
              NClu0=NClu;
              NClu=NNei;
		    }
	}
	
   else if(*(nvth+column)!=DHigh)
        *(state+column)=1;       
    
   return 0;	   
}


int InOut(int target,int *source){   //If target is not in the array source, return 1; otherwise return 0.
	int i=0,sign=1;
	while(source[i]!=-1){
		if(target==source[i++]){
			sign=0;
			break;
		}			
	}
	return sign;		
}

double DisUni(double mean, double HalWid){
	double RanNum;
	RanNum=(double)rand()/(RAND_MAX);
	RanNum=RanNum*2*HalWid-HalWid+mean;
	return RanNum;	
	}
	
int equal(double a, double b, double error){
	if (fabs(a-b)<=error)
	 return 1;
	else 
	 return 0;}  
	 
int WheMultiple(double numerator, double denominator, double error){
	if (fabs(numerator-denominator*(int)((numerator+denominator/2)/denominator))<=error)
	  return 1;
	else
	  return 0;}
	  
int MeaLam(double *wave, int height, int width, int TimStep, double threshold){
  int i,j, it;
  for (i=height-1;i>=0; i--){
	    for(j=0;j<width;j++){
			for(it=0;it<TimStep;it++){
				if( wave[(it*height+i)*width+j]>=threshold )
				  return i+1;}}
	  }
	  
	  return 0;}
	  
int wavefront(double *wave, int height, int width, int ColSta, int ColEnd, double threshold){
  int i,j;
  for (i=height-1;i>=1; i--){
	    for(j=ColSta;j<=ColEnd;j++){
		if( wave[i*width+j]>=threshold )
				  return i;}}
  
	  return 0;}	
	  
int waveback(double *wave, int height, int width, int ColSta, int ColEnd, double threshold){
  int i,j;
  for (i=1;i<height; i++){
	    for(j=ColSta;j<=ColEnd;j++){
		if( wave[i*width+j]>=threshold )
				  return i;}}
  
	  return 99;}

double wavelength(double *wave, int height, int width, int LenSec, double threshold){
int i, WF, WB;
double WL=0;
	for(i=0;i<width;i=i+LenSec){
		WF=wavefront(wave, height, width, i, i+LenSec-1, threshold);
		WB=waveback(wave, height, width, i, i+LenSec-1, threshold);
		if(WF!=0 && WB!=99)
		WL=WL+ WF-WB+1;
		else
		WL=WL+0;               }
		
	WL=WL/(width/LenSec);
	return WL;
		}
	  
int NFir(double *wave, int height, int width, double threshold){
	int i,j, N_firing=0;
	
	for (i=0;i<height; i++)
	    for(j=0;j<width;j++)
		 N_firing+= ( wave[i*width+j]>=threshold );
		 
		 return N_firing;
			}	
			
double TMean(double *Tij, int height, int width, double threshold){
	int i,j, N_firing=0;
	double T_sum=0;
	
	for (i=0;i<height; i++)
	    for(j=0;j<width;j++){
		 N_firing+= ( Tij[i*width+j]>=threshold );
		 T_sum+= Tij[i*width+j]*( Tij[i*width+j]>=threshold );
		   }
		 
		 if(N_firing>=1)
		 return T_sum/N_firing;
		 else 
		 return 0;}  
		 
double gaussrand_K()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}	  
			  
    
int FindPeaks(double deltat, double TimLen, double dt_write, double *trace, const int TimStep, double v_th, double *PeakTim, int LenPeakTime){
	int i=0, NumPeaks=0, state;	
	for (i=0; i<TimStep; i++){
		if (NumPeaks>= LenPeakTime-1)
		  break;
		  
	  state=(i+1)*(i<=1)-(TimStep-i)*(i>=TimStep-2);  //First, second and the last secod last moments.
	  switch(state) {
      case 1 :
        if(trace[i]>= trace[i+1] && trace[i]>= trace[i+2]&& trace[i]>= v_th)
		  PeakTim[NumPeaks++]=2*deltat; 
         break;
      case 2 :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i+1] && trace[i]>= trace[i+2]&& trace[i]>= v_th)			
		  PeakTim[NumPeaks++]=i*dt_write; 
		  break;
      case ( -1 ) :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i-2]&& trace[i]>= v_th)
		  PeakTim[NumPeaks++]=i*dt_write; 
         break;
      case ( -2 ) :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i-2] && trace[i]>= trace[i+1]&& trace[i]>= v_th)			
		  PeakTim[NumPeaks++]=i*dt_write;
         break;
      
      default :
         if(trace[i]>= trace[i-1] && trace[i]>= trace[i-2] && trace[i]>= trace[i+1]&& trace[i]>= trace[i+2]&& trace[i]>= v_th)			
		  PeakTim[NumPeaks++]=i*dt_write;
                  }		
		}
		
	PeakTim[LenPeakTime-1]=NumPeaks;
	
	return 0;	
	}
	
int WheDied(double *wave, int height, int width, int TimStep, double v_th){
	int i,j, it;
   for(i=0; i<=height-1; i=i+height-1){
	for (j=0; j<width; j++){
		for (it=0; it<TimStep; it++){
			if (wave[ (i+it*height) *width +j]>=v_th) 
			  return 0;} 
  } 	}
			  
    for(i=0; i<height; i++){
	for (j=0; j<width; j=j+width-1){
		for (it=0; it<TimStep; it++){
			if (wave[ (i+it*height) *width +j]>=v_th) 
			  return 0;} 
  } 	}
  
  return 1;	
	}
	
int WheNeverDied(double *wave, int height, int width, int TimStep, double v_th){
 int i,j, it=TimStep-1;
	for(i=0; i<height; i++){
	for (j=0; j<width; j++){
		if ( wave[ (i+it*height) *width +j]>=v_th) 
			  return 1;} 
  } 
  return 0;
  	}
