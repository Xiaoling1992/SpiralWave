//Spi11: Decide whether it is still live in solving the PDE. Stop it if it is dead and go to the next iteration to save time.
//Spi10: Follow the coding habit in the "C Primer Plus" for loops and pointers, delete the non-necessary external function. 
        //New definition for spiral and normal. 
//Spi9: Increase phi to 0.98; repeat 1000 times; output phi_simulation;
//Spi8: Parallel the process of categorizing and period. TimLen=300; NumRep=800 to shrink the running time.
//_n6: define anchor size as the part, which is inside the convex hull of the off-cell cluster which has the most cells inside the convex hull.
//_n5: consider all the off-cell clusters. Find out the one which is most likely to be the anchor (The cluster has the most off-cells inside the anchor). Use the cluster size of this off-cell cluster as the anchor size n5.
//n4: exclude all the clusters which are less than 200.
//_n3: calculate both n1 and n2
//Spi7 _n1: anchor size the the average of all te cluster sizes. Each cluster has the same probability to be picked up.
//Spi7: Visulize the structural anchor size. 
//Spi6*b5: period and convex hull, time duration: 50->TimLen.
//Spi6*b5: Convex Hull: typo is verified. consider wave from LocSta->TimLen. if the cells with most peaks are less than 10, then return failure. Anchor size is 0. 
//Spi6*b4: period and convex hull: find out all the cells with the most peaks. (there are too many cells if use the cells with peraks of MCP-1)
//spi6*b3: 1. updated method to categorize transmission: use the ratio of firing cells before the wave reaches the edge to tell break up. 
//         2. period and convex hull: find out all the cells with the most peaks (MCP) and peaks of MCP-1.  
//Spi6*b2: use the area of the convex hull as the anchor size.
//Spi6*b1: fix the bug in ConvexHull: j>2 in finding the convex hull.
//Spi61_b2: anchor is defined as the sum of cluster size containing any off cells inside the Convex Hull
//Spi6 #1. New method to calculate the period: 20 cells, average the 10 cells whose periods are less deviated from the mean; #2. anchor size is defined as the off cells inside the convex hull.
//Spi5. 1. repeat 500 times for each parameter combination; When calculating the period, consider only the wave between 100 to 300.
//Spi4.7. Test the relationsip between period and phi. Calculate the period of each simulaton.
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
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri);     // addition of N*M matrix;
void IniMat_D(int n,int m, double constant,double *matrix);
void IniMat_I(int n,int m, int constant,int *matrix);
int PriMat_D(double *source,int height,int width,int HeiMat,int WidMat,FILE *fpt);
int PriMat_I(int *source,int height,int width,int HeiMat,int WidMat,FILE *fpt);
double SumMat_D(double *source, int height, int width, int HeiMat, int WidMat, double threshold, char mode);
int equal(double a, double b, double error);
int WheMultiple(double numerator, double denominator, double error);

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
                  double DHigh,int *cluster);
int SearchClu_G4_Squ(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DHigh,int *cluster);
int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square);
double minimum(double a,double b);

int MeaLam(double *wave, int height, int width, int TimStep, double threshold);
int wavefront(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
int waveback(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
double wavelength(double *wave, int height, int width, int LenSec, double threshold);
int NFir(double *wave, int height, int width, double threshold);
double TMean(double *Tij, int height, int width, double threshold);

double AssignPhi(double phi,double low,double high);
double DisUni(double mean, double HalWid);
double gaussrand_K(void);

int FindPeaks(double deltat, double TimLen, double dt_write, double *trace, int LocSta, int LocEnd, const int TimStep, double v_th, double *PeakTim, int LenPeakTime);
int WheDied(double * wave, int height, int width, int TimStep, double v_th);
int WheDied_Tim(double *wave, int height, int width, int TimStep, double v_th);   //died, -1; non-died, return the time sequence.
int WheNeverDied(double *wave, int height, int width, int TimStep, double v_th);
int WheLive_t(double *wave, double t, double dt_write, int height, int width, int TimStep, double v_th);

int ConvexHull(double *SetX, int SizeSetX, int length, double *ConHull, int LenConHull);	
int ExchangeD(double * a, double *b);
int WheInConHull(double *ConHull, int NumPoints, int LenConHull, double x, double y);
double CroPro2D(double x1, double y1, double x2, double y2);
double PerCelMaxPea(double *PerCMP, int NumCMP);

int MinArrXZ_I(int *array, int length);
int MaxArrXZ_I(int *array, int length);

struct ConHulPoint{
	double xaxis;
	double yaxis;
	double cos;
	};

int main(void){ 

//$$ 1. change the value of parameter a and b,c;
double aglobal, aLow, as[]={0};  ////phi: percentage of responsive cells; rho, clustering rate, aLow, the parameter of unresponsive cell.                                                 
double tau, tau_off, phi, rho;
double cglobal, cLow, cs[]={10};
double gamma=0;
double rhos[]={0, 0.2, 0.4, 0.6, 0.8, 1.0}, taus[]={9}, phis[]={0.46, 0.5 };

int  NumRep=1, Len_as=sizeof(as)/sizeof(as[0]), Len_cs=sizeof(cs)/sizeof(cs[0]),  Len_rhos=sizeof(rhos)/sizeof(rhos[0]), Len_taus=sizeof(taus)/sizeof(taus[0]), Len_phis=sizeof(phis)/sizeof(phis[0]);
int i,j, ia, ic, itau, iphi, irho, ip; //i, j is only used for any short, self-closed loop. 



//$$ 2. Number of row: N, column, M. Print type,1:all cells,5:row 0:5:30. 
    //Noise:0. no noise;1. clustering, bimodal distribution; 2.clustring, uniform distribution in log space; 250. other noise. 
    //Change phi,psi,aLow,tau_off;   $$$$$$$$$$$$$$$$*//
    //Change deltat, TimeLen, deltax, deltay;
int N=100,M=100, Leng=N*M, coefx=1,coefy=1;   //N, the number of rows, M, the number of columns.

double v_th=0.6;
 
double deltat=0.01,t=0.00,TimLen=400.00,deltax=1, deltay=1, dt_write=2;          //$$ TimLen incresed to 400 min. 
double EVnST[2][N*M],k0[2][N*M],EVnST1[2][N*M],k1[2][N*M],k2[2][N*M],k3[2][N*M],aNoise[N*M],tauNoise[N*M],cNoise[N*M];   //row 0 stores the v, row 1 stores w.


char filename[100], filename2[100];
char  num[100]; 
 time_t tt;
 srand((unsigned) time(&tt));  //Initiating random number generator only once. 


const int TimStep=(int)(TimLen/dt_write+EPSILON2)+1; 
double wave[N*TimStep][M];  //Used only for test. Please

FILE *fptPeriod,  *fptWave, *fptMap, *fptPhi;     //fpt: ThT; fpt2: parameter;
fptPeriod=fopen("Spi11_Period_Phi_Rho95Tau9_delete.dat", "w");
fptPhi=    fopen("Spi11_PhiAct_Phi_Rho95Tau9_delete.dat", "w");
   

for(irho=0; irho<Len_rhos; irho++){
for(itau=0;itau<Len_taus;itau++){
  double period[Len_phis][NumRep], PhiAct[Len_phis][NumRep];  //define the variable in necessary postion to make it easier to understand. 
  IniMat_D(Len_phis, NumRep, -1, &period[0][0]);	
  IniMat_D(Len_phis, NumRep, -1, &PhiAct[0][0]);	//-1 means it is a non-spiral
 
for(iphi=0;iphi<Len_phis;iphi++){ 
	   
for(ia=0;ia<Len_as;ia++){
for(ic=0;ic<Len_cs;ic++){
	
	aglobal=as[ia];
	aLow=as[ia];
	cglobal=cs[ic];
	cLow=cs[ic];	
	
	rho=rhos[irho]; 
	tau_off=5;	//Tau_off is not applied in LogNormal distribution
	tau=taus[itau]; 
	phi=phis[iphi];    
     
     strcpy(filename, "Spi11"); 
     strcat(filename,"Rho");   sprintf(num, "%d",(int)(rho*100+EPSILON));strcat(filename,num); 
     strcat(filename,"Phi");   sprintf(num, "%d",(int)(phis[iphi]*1000+EPSILON) );strcat(filename,num);
     strcat(filename,"t");   sprintf( num, "%d",(int)(taus[itau]+EPSILON) );    strcat(filename,num);     
     strcat(filename,"c");   sprintf(num, "%d",(int)(cglobal));strcat(filename,num);
     strcpy(filename2,filename);
    
     strcat(filename,"_wave_Fig3A.dat");	
    strcat(filename2,"_map_Fig3A.dat");
	
	fptWave=fopen(filename,"w");
	fptMap=fopen(filename2,"w");

int NumSpi=0;	
for(ip=0;ip<NumRep; ip++){         
      //strcat(filename,"u");   sprintf(num, "%d",(int)(tau+EPSILON);strcat(filename,num);     
      //strcat(filename,"s");  sprintf( num, "%d",(int)(sigma) );strcat(filename,num);
      
  
	//'s', homogeneous, single value; 'b', binary distribution; 'N', LogNormal; 'U' Loguniform;
	CluRadLin(aNoise,   N, M, 's', -999, -999, -999, aglobal, -999);   //int CluRadLin(double *nvth, int height, int width, char mode, double rho, double phi, double DLow, double median, double lambda);
    CluRadLin(tauNoise, N, M, 'b', rho, phi, tau_off, tau, -999); //-999 means that input is not necessary in the function.
    CluRadLin(cNoise,   N, M, 's', -999, -999, -999, cglobal, -999);    
     
    int WheSpiral; 
	t=0;
	IniMat_D(2,N*M,0.00,&EVnST[0][0]);        //Initiate the EVnST and t. I don't need to Iniate the EVnST1, k0,k1,k2,k3 and aNoise.    
	for (i=N/2-3; i<=N/2+2; i++){       //Move the trigger out of while loop to reduce the running time.  
		for(j=M/2-3; j<=M/2+2; j++)
		      EVnST[0][i*M+j]=1;  	       // Trigger: pulse v:1;   
    }
   
 //$$ 4. Change the Trigger Time,length, and number of rows $$$$$$$$$$$$$$$$*//  
  while (t<TimLen) {            
	
    dxdt_Squ_o4(&EVnST[0][0],&k0[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay); 
	#pragma omp parallel for private(i)
	for(j=0;j<N*M;j++){
		for(i=0;i<2;i++)
		    EVnST1[i][j]=EVnST[i][j]+deltat*k0[i][j];}    
		      
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
        
			  
	    if(WheMultiple(t, dt_write, deltat/10)|| equal(t,2*deltat,deltat/10) ){
			for (i=0; i<N; i++){
			  for (j=0;j<M;j++)
			    wave[i+ (int)(t/dt_write+EPSILON2) *N][j]=EVnST[0][i*M+j];
			}
			
			if (t> 24-EPSILON2 && (WheSpiral=WheLive_t(&wave[0][0], t, dt_write, N, M, TimStep, v_th) ) ==0 ) //24 is the time when the wave reached the edge.
			    break;			   
		}
   
             			                      
  }//end of while loop;   
 
    //Visualize the structural anchor size: off-cell wighted average cluster size
    //int state[N*M], cluster[N*M], numerator=0, denominator=0, numerator2=0,CluSiz;
	 //IniMat_I(1, N*M, 0, state);
	 //for (i=0; i<N; i++){
		 //for (j=0; j<M; j++){			 
		    //CluSiz=SearchClu_G4_Squ(tauNoise, state, i, j, N, M, N, M, tau_off, cluster);  //Because iterations (i,j) are not independent from each other, we cann't parallel them.
		    //if (CluSiz>=1){     //account all the clusters.
				//denominator+=1;
				//numerator+=CluSiz;
				//numerator2+=CluSiz*CluSiz; }		      
		  //}
	  //}	  
	  //if(denominator >=1)
      //AncSizeStr[iphi][ip]=(double) numerator2/ (double)numerator;   //S=n1^2+ n2^2 n3^3.../n1+n2+n3...
      
     //Visualize the PhiAct
     
     //Calculate the period of the spiral wave: 1. find out all cell with the most peaks (CMP) and peaks of CMP-1; 
     //      2. give up half cells whose periods are more deviated from the mean; 
     //      3. average the periods of the left cells as the period of the spiral Wave.  
     //      The head of the spiral is more frequent than the tail of the spiral.  
              //Categorize the wave   
              
                
  int NumSample=N*M, k, it, LenPeakTime=30;
  double PeakTim[NumSample][LenPeakTime]; 
    
   if( WheSpiral  ){
	   //Print the first spiral wave.	
	  NumSpi++;
	  if(NumSpi<=1){
		  PriMat_D(&wave[0][0], N*TimStep, M, N*TimStep, M, fptWave);
		  
		  PriMat_D(aNoise, N, M, N, M, fptMap);  
		  PriMat_D(tauNoise, N, M, N, M, fptMap); 
		  PriMat_D(cNoise, N, M, N, M, fptMap); 	  
	 }
	 
	 PhiAct[iphi][ip]=SumMat_D(&tauNoise[0], N, M, N, M, (tau+tau_off)/2, 't' )/N/M;        //Visulize the PhiAct only when it is a spiral.	
		
	  int  NoMaxPea=0, MaxPea=0,  peaks_k, count=0, LocSta=(int)(100/dt_write+EPSILON), NumCMP=0;// CelMaxPeak[2][N*M];
	  double PerCMP[N*M];
	  
	//Visulize the period and convex hull. Exclude the staring 50 or 100 mins of chaos. wave duration: LocSta-> TimStep.
	IniMat_D(NumSample, LenPeakTime, TimLen+1, &PeakTim[0][0]);  //If no peak, the the peak time is TimLen+1;  	
  #pragma omp parallel private(j, it)
  {	
	double trace[TimStep];
	#pragma omp for
	for (i=0; i<N; i++){  
	  for (j=0; j<M; j++){ 
			for (it=0; it< TimStep; it++)
		      trace[it]=wave[i+it*N][j];
	   	    FindPeaks(deltat, TimLen, dt_write, trace, LocSta ,TimStep,TimStep, v_th, &PeakTim[i*M+j][0], LenPeakTime); //consider only the trace before the LocEnd
	   } 
	 }
  }
		   
	  //puts("The PeakTim");
	  //for(i=0; i<LenPeakTime; i++)
	     //printf("%lf  ", PeakTim[74*M+71][i]);
	  //puts("\n");
	   
	  //find out the max peaks. 
	  for (k=0; k<N*M; k++){
		  if (PeakTim[k][LenPeakTime-1]> PeakTim[NoMaxPea][LenPeakTime-1])
		    NoMaxPea=k;}	  
	  MaxPea=(int) (PeakTim[NoMaxPea][LenPeakTime-1]+EPSILON);  //Notice the data type. Always use int or double to transfer the data type.
	  //printf("2. MaxPeaks %d\n", MaxPea);
	  
	  //find out 20 cells with the most peaks (CMP)
	  for(k=0; k<N*M; k++){
			  if ((int)(PeakTim[k][LenPeakTime-1]+ EPSILON2)== MaxPea){  //I can't compare a float with an int. 
				  //CelMaxPeak[0][count]= k/M;  //i=k/M, x axis
				  //CelMaxPeak[1][count]= k%M;  //j=k%M, y axis
				  peaks_k=(int)(PeakTim[k][LenPeakTime-1]+ EPSILON2);
				  if (peaks_k >= 2)   //Consider all the possibilities. Make it 
				    PerCMP[count]= (PeakTim[k][ peaks_k-1] -PeakTim[k][0])/(peaks_k-1);
				  else
				    PerCMP[count]=(int)TimLen+1;
				  count++;
			  }
		  }
	  NumCMP=count;
	  if (NumCMP<=1)
	    period[iphi][ip]=PerCMP[0];
	  else  
	    period[iphi][ip]=PerCelMaxPea(PerCMP,  NumCMP); //The order of PerCMP is changed by the deviation.	
      
   }  //End of if it is a spiral
   
  
   
}   //end of loop ip   

        printf("rho %lf tau=%lf phi %lf\n",rho, tau, phi);  
        if(NumSpi<=0){
		 PriMat_D(&wave[0][0], N*TimStep, M, N*TimStep, M, fptWave);	 
		 PriMat_D(aNoise, N, M, N, M, fptMap);  
		 PriMat_D(tauNoise, N, M, N, M, fptMap); 
		 PriMat_D(cNoise, N, M, N, M, fptMap); 
		  }
        
	 fclose(fptWave);
	 fclose(fptMap);

}   //end of loop ic
}   //end of loop ia

}  //end of iphi
		         
PriMat_D(&period[0][0], Len_phis, NumRep, Len_phis, NumRep, fptPeriod); 
PriMat_D(&PhiAct[0][0], Len_phis, NumRep, Len_phis, NumRep, fptPhi); 
		         
		         
//PriMat_D(&AncSizeDyn_n6[0][0], Len_phis, NumRep, Len_phis, NumRep, fptAnc);

} //end of loop itau


				        		         
} //end of loop irho		


fclose(fptPeriod);
fclose(fptPhi);
				
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


//Define the addition of Matrix.
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri)
{    // addition of N*M matrix;
   int i,j;
   for(i=0;i<n;i++)
       for(j=0;j<m;j++)
           *(Matout+i*Length+j)=*(Matin+i*Length+j)+*(Deri+i*Length+j)*deltat;
    }


//IniMat_D: make every element of the matrix to be a constant.
void IniMat_D(int n,int m, double constant,double *matrix)
{int i=0,j=0;
   	
   for(i=0;i<n;i++)
	   for(j=0;j<m;j++)
		   *(matrix+i*m+j)=constant;
	
}

void IniMat_I(int n,int m, int constant,int *matrix)
{int i=0,j=0;
   for(i=0;i<n;i++)
	   for(j=0;j<m;j++)
		   *(matrix+i*m+j)=constant;
	
}

int PriMat_D(double *source,int height,int width,int HeiMat,int WidMat,FILE *fpt){
	int i,j;
	
		for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			  if(j==width-1)
	          fprintf(fpt,"%lf\n",*(source+i*WidMat+j));
	          else
	          fprintf(fpt,"%lf ",*(source+i*WidMat+j));}}
	 return 0;
	
	}
	
double SumMat_D(double *source, int height, int width, int HeiMat, int WidMat, double threshold, char mode){
	int i, j;
	double sum=0;
	if (mode== 't'){  //With threshold
	  for(i=0; i<height; i++){
		  for(j=0; j<width; j++){
			  sum+= (source[i*WidMat+j]>= threshold ); 
		  }
	   }
	 }
	 
	else if (mode== 's'){  //With threshold
	  for(i=0; i<height; i++){
		  for(j=0; j<width; j++){
			  sum+= source[i*WidMat+j]; 
		  }
	   }
	 }
	return sum;	
	}

int PriMat_I(int *source,int height,int width,int HeiMat,int WidMat,FILE *fpt){
	int i,j;
	
		for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			  if(j==width-1)
	          fprintf(fpt,"%d\n",*(source+i*WidMat+j));
	          else
	          fprintf(fpt,"%d ",*(source+i*WidMat+j));}}
	          
	 return 0;
	
	}

//Create a biofilm with correlation strength psi. All cells are off at the beginning. Find out an off cell in the biofilm, turn it on. It has a probability psi 
     //to give birth to a daughter (on): Find out all the empty neighbors, and randomly pick up one as the daughter and turn it on. The the 
     //they have a probability psi to give birth to one more daughter (on)....If They failed to give birth to the daughter, then reselect an off cell and repeat.  
int CluGen(double *aNoise,double phi,double psi, int N,int M,double aLow,double DHigh)
{int cluster[N*M],neighbor[N*M],NClu,NNei,NNeiE,NHigh=0,i;
	for(i=0;i<N*M;i++)  
	  cluster[i]=-1;       //Initiate it before the first use.
	for(i=0;i<N*M;i++)  
	  neighbor[i]=-1;
	
  while(NHigh<phi*N*M){         //compare two doubles: Be careful for possible calculation errors of computer.
	 i=0;
	 while(cluster[i]!=-1&&i<N*M)	  //Initiate the cluster. -1 is the end of effective part. -2 means that position is concealed.
		cluster[i++]=-1;
     i=0; 
     while(neighbor[i]!=-1&&i<N*M)
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
		     i=0;
		 }
				     
		 
			    
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
	IniMat_I(height, width, -999, &cluster[0]);
	IniMat_I(height, width, -999, &state[0]);
	IniMat_D(height, width, -999, nvth);				 

	
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
	
	IniMat_I(HeiMat,WidMat, 0, &state[0][0]);
	IniMat_I(HeiMat,WidMat, 0, &cluster[0][0]);
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

//We need to initiate the matrix state before using this function. 0: unchecked. 1: checked to be non-targe. 2. checked to be the target. 	   
int SearchClu_G4_Squ(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
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
				   
		 if(j1%2==0){
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
				  return i+1;
			}
		}
   }
	  
	  return 0;
}
	  
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
		 return 0;
}  
		 
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
			  
    
int FindPeaks(double deltat, double TimLen, double dt_write, double *trace, int LocSta, int LocEnd, const int TimStep, double v_th, double *PeakTim, int LenPeakTime){
	int i=0, NumPeaks=0, state;	
	for (i=LocSta; i<LocEnd; i++){
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

//Mean the half perid data which is less deviated. 	
double PerCelMaxPea(double *PerCMP, int NumCMP){
	int i=0, j;
	double mean=0, deviation[NumCMP];
	
	for (i=0; i<NumCMP; i++)
		mean+= PerCMP[i];
	mean/= NumCMP;
	
	for (i=0; i<NumCMP; i++)
	    deviation[i]= fabs( PerCMP[i]- mean); 
	    
	for(i=1; i< NumCMP; i++){
		j=i;
		while(j>=1 && deviation[j]< deviation[j-1]){
			ExchangeD(&deviation[j], &deviation[j-1]);
			ExchangeD(&PerCMP[j], &PerCMP[j-1]);
			j--;
			}
	}
	
	mean=0;
	for(i=0; i< NumCMP/2; i++)
	  mean+= PerCMP[i];
	mean/= (NumCMP/2);
	return mean;  
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
	
int WheDied_Tim(double *wave, int height, int width, int TimStep, double v_th){   //died, -1; non-died, return the time sequence. 
	int i,j, it;
  for (it=0; it<TimStep; it++){	
	  
   for(i=0; i<height; i=i+height-1){
	for (j=0; j<width; j++){		
			if (wave[ (i+it*height) *width +j]>=v_th) 
			  return it;} 
   }
  
  for(i=0; i<height; i++){
	for (j=0; j<width; j=j+width-1){
		if (wave[ (i+it*height) *width +j]>=v_th) 
			  return it;
  } 	
  } 	
  
  
  }
			  
    
  
  return -1;	
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
  	
int WheLive_t(double *wave, double t, double dt_write, int height, int width, int TimStep, double v_th){
	int i, j, it= (int)(t/dt_write+ EPSILON2);
	
	for(i=0; i<height; i++){
	for (j=0; j<width; j++){
		    if ( wave[ (i+it*height) *width +j]>= v_th ) 
			  return 1;
    } 
    } 
  return 0;
	
     }
  	
int ConvexHull(double *SetX, int SizeSetX, int length, double *ConHull, int LenConHull){
	struct ConHulPoint SorPoi[SizeSetX];
	int i, j;
	for (i=0; i< SizeSetX; i++){
		SorPoi[i].xaxis= SetX[i];
		SorPoi[i].yaxis= SetX[i+length];
		SorPoi[i].cos= 0;
		}
	//Find out A0: the 1)Lowest 2)most left cell 	
	int A0=0;
	for (i=0; i<SizeSetX; i++){
		if( SorPoi[i].yaxis < SorPoi[A0].yaxis  ||  (SorPoi[i].yaxis == SorPoi[A0].yaxis && SorPoi[i].xaxis < SorPoi[A0].xaxis) )  
	      A0=i;
	 }
	//Put A0 at the beginning of the array. 
	ExchangeD(&SorPoi[0].xaxis, &SorPoi[A0].xaxis);
	ExchangeD(&SorPoi[0].yaxis, &SorPoi[A0].yaxis);
	A0=0;
	
	//Calculate the cosine of angle between A0Ai and x axis, vector (1,0)
	SorPoi[A0].cos=1;
	for (i=1; i<SizeSetX ; i++){
	   SorPoi[i].cos=(SorPoi[i].xaxis - SorPoi[A0].xaxis )
	           /sqrt( pow(SorPoi[i].xaxis -SorPoi[A0].xaxis, 2) + pow(SorPoi[i].yaxis - SorPoi[A0].yaxis, 2) ); }
	           
	//Sort the array SorPoi[i] by deceeding the SorPoi[i].cos
	for (i=1; i<SizeSetX; i++){
		j=i;
		//If angle(A0Aj,x) < angle(A0Aj-1, x) switch. If the angle is the same, but the distance(A0Aj)< distance(A0Aj-1), switch.  
		while (j>=1 &&(SorPoi[j].cos > SorPoi[j-1].cos || (SorPoi[j].cos == SorPoi[j-1].cos && pow(SorPoi[j].xaxis- SorPoi[A0].xaxis, 2)+pow(SorPoi[j].yaxis- SorPoi[A0].yaxis, 2) < pow(SorPoi[j-1].xaxis- SorPoi[A0].xaxis, 2)+pow(SorPoi[j-1].yaxis- SorPoi[A0].yaxis, 2) ) ) ){
		   ExchangeD(&SorPoi[j].xaxis, &SorPoi[j-1].xaxis);
	       ExchangeD(&SorPoi[j].yaxis, &SorPoi[j-1].yaxis);
	       ExchangeD(&SorPoi[j].cos, &SorPoi[j-1].cos);
		   j--;
		   }
		j++;		   
	 }
	 
	 //puts("sorted array");
	 //for (i=0; i<SizeSetX; i++)
	   //printf("%lf %lf\n", SorPoi[i].xaxis, SorPoi[i].yaxis);
	
	 
	 
	 
	//Find out the minmal cells necessary for the convex hull. Colinear cells are discarded. Use the structure of stack: last in, first out. Push and pop.
	ConHull[0]=SorPoi[0].xaxis; ConHull[0+LenConHull]=SorPoi[0].yaxis;
	ConHull[1]=SorPoi[1].xaxis; ConHull[1+LenConHull]=SorPoi[1].yaxis;
	j=2;
	for (i=2; i<SizeSetX; i++){
	  ConHull[j]=SorPoi[i].xaxis; ConHull[j+LenConHull]=SorPoi[i].yaxis;	
	  
	  while (j>=2&& ((ConHull[j-1]-ConHull[j-2])*(ConHull[j+LenConHull]- ConHull[j-1+LenConHull])-(ConHull[j-1+LenConHull]-ConHull[j-2+LenConHull])*(ConHull[j]- ConHull[j-1])<= 0) ){ //turn right or colinear
	    ConHull[j-1]=ConHull[j]; 
	    ConHull[j-1+LenConHull]=ConHull[j+LenConHull] ;
	    j--;
	    
	   // printf("discard one cell %lf %lf \n", ConHull[j], ConHull[j+LenConHull]);
	   }
	   j++; 
	}
	
return j;
}

int WheInConHull(double *ConHull, int NumPoints, int LenConHull, double x, double y){
	int i;
	//On the right or colinear with AN-1A0
	if ( CroPro2D(ConHull[0]-ConHull[NumPoints-1], ConHull[0+LenConHull]-ConHull[NumPoints-1+LenConHull], x- ConHull[0], y- ConHull[0+LenConHull]) <=0 )
		  return 0;
	for (i=1; i<NumPoints; i++){
		if ( (ConHull[i]-ConHull[i-1])*(y- ConHull[i+LenConHull])-(ConHull[i+LenConHull]-ConHull[i-1+LenConHull])*(x- ConHull[i])<= 0)   //On the right or colinear with Ai-1Ai
		  return 0;
	}
	
	return 1;
	}

int ExchangeD(double * a, double *b){
	double TemExc;
	TemExc= *a;
	*a= *b;
	*b=TemExc;
	return 0; }
	
double CroPro2D(double x1, double y1, double x2, double y2){
	return x1*y2-y1*x2;}

int MinArrXZ_I(int *array, int length){
	int k=0, mini;
	for(k=0; k<length; k++){
		//the above: minimal i
		if(array[k]< array[mini]) 
		  mini=k;
	 }
	 
	 return mini;	
	}

int MaxArrXZ_I(int *array, int length){
	int k=0, max=0;
	for(k=0; k<length; k++){
		//the above: minimal i
		if(array[k]> array[max]) 
		  max=k;
	 }
	 
	 return max;	
	}
