//This is used to only output the time trace of very cells. The analysis of dutation at. al is stopped. 
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

#define HEIGHT 100
#define WIDTH 200         //Use as less global variable in main function; don't use global variable in external function except defining a matrix.
#define LENG2MAT HEIGHT*WIDTH 
//Use the equation dEdt=-delta(-D*deltaE);
//Add noise to explain the double spiking wich  occurses in the experiment of Gurol. Applied to 2D grid. Parallel the Runge Kunta part within while loop.
// Parallel Algorithm, run in cluster: No iteractive jobs.
//Trigger 9 cells in the ceter of cell community, use Runge Kutta Algorithm.
// Reproducing the work in Fig 3d and 3e in Pindle et al's Nature, 2015,"Ion channels enable electrical communication in bacterial communities" 
//May 11, 2016, Purdue; 
//2d N*M cells,Protocol b: Equilibrate for 100 min, then pulse then pulse S=100 uM ?or E=200*1000 uM?

void dxdt_Tri(double *EVnST,double *deri,double *aNoise,double *tauNoise, double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);   //deri[] is used to store the data of dVdt,dndt......
void dxdt_Tri_Top(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay);
void print1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt);
void print1_1(int N,int M,int Length, double t, double tInterval,double deltat,double *EVnST, FILE *fpt);
int print5(int N, int M, double t, double *EVnST, FILE *fpt);
void Mplus(int n, int m, int Length, double deltat, double *Matin, double *Matout, double *Deri);     // addition of N*M matrix;
void Transport(int N,int M,int triy,int trix,int Ln,int Lm,int n,double *source,double *Tprint);
void Ini_Matri(int n,int m, double constant,double *matrix);
FILE * name(int N,int M,int *coefx, int *coefy,char *filename);
void Island(double *aNoise,double percentage,int N,int M,double lisland,double F);
int InOut(int target,int *source);
int CluGen(double *aNoise,double phi,double psi, int N,int M,double aLow,double DHigh);
int CluGen_G5_1(double *nvth,int height,int width,double phi,double rho,double DLow,double DHigh,
              double alpha, double beta, double a1, double a2, double a3);
int mothergrandmother(char mode, double *nvth, int HeiMat, int WidMat, int height, double DLow, double DHigh,
           double *possibility, double ponoffoff, double ponoffon, double pononoff, double pononon);
double AssignPhi(double phi,double low,double high);
int SearchTree_G4(double *nvth,int *state,int i,int j,int height,int width,int HeiMat,int WidMat,
                  double DLow,int *cluster);
int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square);
double minimum(double a,double b);
double DisUni(double mean, double HalWid);

int equal(double a, double b, double error);
int WheMultiple(double numerator, double denominator, double error);

int MeaLam(double *wave, int height, int width, int TimStep, double threshold);

int wavefront(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
int waveback(double *wave, int height, int width, int ColSta, int ColEnd, double threshold);
double wavelength(double *wave, int height, int width, int LenSec, double threshold);
int NFir(double *wave, int height, int width, double threshold);
double TMean(double *Tij, int height, int width, double threshold);


long MaxRandom=2147483647;
double EPSILON=0.0000000001;
double EPSILON2=0.00001;




int gk=30; //min^-1
double gl=0.2; //min^-1
int vk0=-380; //mV
int vl0=-156; //mv
int sth=40;   //uM
int vth=-150; //mV
int alpha0=2;  //min^-1
double beta=1.3;  //min^-1
int mglobal=1;
int Fglobal=5600;  //uM/mV
double sigma=0.2;  //mV
double deltak=0.001; //mV/uM;
double deltal=0.008; //mV/uM;
double gammas=0.1;  //min^-1
int gammae=10;   //min^-1
int gammat=4;  //min^-1
int alphas=1;  //uM/min/mV
int alphat=1;   //uM/min/mV
int Dglobal=82800;   //1380*60um^2/min


int main(void){ 

//$$ 1. change the value of parameter a and b,c;
double aglobal,as[]={0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1};  //WT, trkA, sinR, ktrA
double tau, taus[]={300, 200, 100, 75, 50, 25, 20, 15};
double cglobal, cs[]={10};
double gamma=0, phi_real;
double phi,rho, phis[]={0.5};

int  Len_as=9, Len_taus=8, Len_phis=1, NumRep=1;



//$$ 2. Number of row: N, column, M. Print type,1:all cells,5:row 0:5:30. 
    //Noise:0. no noise;1. clustering, bimodal distribution; 2.clustring, uniform distribution in log space; 250. other noise. 
    //Change phi,psi,aLow,tau_off;   $$$$$$$$$$$$$$$$*//
    //Change deltat, TimeLen, deltax, deltay;
int N=HEIGHT,M=WIDTH,WinHei=35, print=1,noise=3,coefx,coefy;   
double aLow, tau_off, cLow;   //phi: percentage of responsive cells; psi, clustering probability, aLow, the parameter of unresponsive cell.
double lambda_a,lambda_b, lambda_c;   //phi: percentage of responsive cells; psi, clustering probability, aLow, the parameter of unresponsive cell.
 
 
int WF_t, WF, WB_t;
double WL_t, WL;
double v_th=0.6, T_th=0.1;
double WF_ab[Len_as][Len_taus], WL_ab[Len_as][Len_taus], PhiFir_ab[Len_as][Len_taus], TraEff_ab[Len_as][Len_taus],
       T_ab[Len_as][Len_taus], Tij[WinHei][M];   
FILE *fpt, *fpt1, *fptWave, *fptMap;     //fpt: ThT; fpt2: parameter;
char filename[100], filename2[100];
int i=0,j=0,ia,itau,ic, iphi, ip;
double deltat=0.01,t=0.00,TimLen=200.00,deltax=1, deltay=2, dt_write=0.5;          //$$ TimLen incresed to 400 min.    
int TimStep=(int)(TimLen/dt_write+1+EPSILON2);

double EVnST[2][LENG2MAT],k0[2][LENG2MAT],EVnST1[2][LENG2MAT],k1[2][LENG2MAT],k2[2][LENG2MAT],k3[2][LENG2MAT],aNoise[LENG2MAT],tauNoise[LENG2MAT],cNoise[LENG2MAT];   //row 0 stores the v, row 1 stores w.
int Leng=LENG2MAT;   //*$$$$$$Must Must change the value of Leng when the  dimensions of EVnST,ko are changed.   
char  num[100]; 
 time_t tt;
 srand((unsigned) time(&tt));  //Initiating random number generator only once. 

fpt=fopen("F17_spiral_T.dat","w");
 fpt1=fopen("F17_spiral_sta.dat","w");
 
     Ini_Matri(Len_as, Len_taus,0.00, &WF_ab[0][0]);
     Ini_Matri(Len_as, Len_taus,0.00, &WL_ab[0][0]);
     Ini_Matri(Len_as, Len_taus,0.00, &PhiFir_ab[0][0]);
     Ini_Matri(Len_as, Len_taus,0.00, &TraEff_ab[0][0]);
     Ini_Matri(Len_as, Len_taus,0.00, &T_ab[0][0]);
//double wave[WinHei*TimStep][M];  //Used only for test. Please

//This is a loop along a verctor, iphi in the only loop variable, Len_phis is the length Lenth of the vector.
//Len_as, Len_taus must be 1;
for(iphi=0;iphi<Len_phis;iphi++){   
	
for(ip=0;ip<NumRep; ip++){
	   
for(ia=0;ia<Len_as;ia++){
for(itau=0;itau<Len_taus;itau++){
for(ic=0;ic<1;ic++){
	phi=phis[iphi];
    rho=0;    
    aglobal=as[ia];
	aLow=as[ia];	
	tau=taus[itau];
	tau_off=5;	
	
	cglobal=cs[ic];
	cLow=cs[ic];	
     
	 
     
    strcpy(filename, "F17"); 
    strcat(filename,"Phi");   sprintf(num, "%d",(int)(phi*100));strcat(filename,num); 
    strcat(filename,"Rho");   sprintf(num, "%d",(int)(rho*100));strcat(filename,num); 
      strcat(filename,"a"); sprintf(num, "%d",(int)(aglobal*1000)); strcat(filename,num); 
      strcat(filename,"t");   sprintf(num, "%d",(int)(tau));strcat(filename,num);
         strcat(filename,"tL");  sprintf(num, "%d",(int)(tau_off));strcat(filename,num);
      strcat(filename,"c");   sprintf(num, "%d",(int)(cglobal));strcat(filename,num);
    strcpy(filename2,filename);
    
    strcat(filename,"_wave.dat");	
	strcat(filename2,"_map.dat");
	
	fptWave=fopen(filename,"w");
	fptMap=fopen(filename2,"w");
	
 
	
		
  
   if (noise==0)           { 
	   for (i=0;i<N*M;i++){          //Initiate the parameters in cell grid.  
	     aNoise[i]=aglobal;
	      tauNoise[i]=tau;
	      cNoise[i]=cglobal;}  }
      
   else if(noise==1){                //Because responsive cell has a=0.2 (DHigh), while unresponsive cell has a=0.3 (aLow), DHigh<aLow. Thus make DHigh=-0.2, aLow=-0.3. Finally reverse the values. 
    do{
     CluGen_G5_1(aNoise, N, M, phi, rho, aLow, aglobal, 1.5873*rho+1, 1.5873*rho+1, 0, 1, 0);
     CluGen_G5_1(tauNoise, N, M, phi, rho, tau_off, tau, 1.5873*rho+1, 1.5873*rho+1, 0, 1, 0);
     CluGen_G5_1(cNoise, N, M, phi, rho, cLow, cglobal, 1.5873*rho+1, 1.5873*rho+1, 0, 1, 0);
    
      phi_real=0;
      for(i=0;i<N*M;i++) 
         phi_real+=(tauNoise[i]==tau); 
           phi_real=phi_real/N/M; 
            printf("ip=%d phi_measure=%lf  phi_design=%lf error=%lf \n",ip, phi_real, phi, fabs(phi_real-phi)); }           
      while  ( fabs(phi_real-phi) > 0.003 );
        }
     
    
   else if(noise==2){
	    
	   	     for (i=0;i<N*M;i++){		                    
		     aNoise[i]=exp(DisUni(log(aglobal),log(lambda_a))); 
		     tauNoise[i]=exp(DisUni(log(tau),log(lambda_b)));
		     cNoise[i]=exp(DisUni(log(cglobal),log(lambda_c))); }                     //Because responsive cell has a=0.2 (DHigh), while unresponsive cell has a=0.3 (aLow), DHigh<aLow. Thus make DHigh=-0.2, aLow=-0.3. Finally reverse the values. 	
		}
		
  else if(noise==250){
	  for (i=0;i<M;i++){		                    
		     aNoise[i]=exp(DisUni(log(aglobal),log(lambda_a)));  
		     tauNoise[i]=exp(DisUni(log(tau),log(lambda_b))); 
		     cNoise[i]=exp(DisUni(log(cglobal),log(lambda_c))); }  
		     
	  for (i=M;i<N*M;i++){		                    
		     aNoise[i]=aNoise[i-M];  
		     tauNoise[i]=tauNoise[i-M];
		     cNoise[i]=cNoise[i-M];} 
	}
	
  else if (noise==3){
//$$ 5. Decide when to generate new structure.
   if (ia+itau+ic==0){
	   
	do{
    CluGen_G5_1(aNoise, N, M, phi, rho, aLow, aglobal, 1.5873*rho+1, 1.5873*rho+1, 0, 1, 0);
     CluGen_G5_1(tauNoise, N, M, phi, rho, tau_off, tau, 1.5873*rho+1, 1.5873*rho+1, 0, 1, 0);
     CluGen_G5_1(cNoise, N, M, phi, rho, cLow, cglobal, 1.5873*rho+1, 1.5873*rho+1, 0, 1, 0);
    
      phi_real=0;
      for(i=0;i<N*M;i++) 
         phi_real+=(tauNoise[i]==tau); 
           phi_real=phi_real/N/M; 
            printf("ip=%d phi_measure=%lf  phi_design=%lf error=%lf \n",ip, phi_real, phi, fabs(phi_real-phi)); }           
      while  ( fabs(phi_real-phi) > 0.003 );
     
       }

  else{	
		for (i=0;i<N*M;i++){
			aNoise[i]=aglobal;
			cNoise[i]=cglobal;
			
			if (tauNoise[i]!=tau_off)
			  tauNoise[i]=tau;
			                }	
	  }	
			          }   //end of NOise==3
	
	
        
	      
	t=0;
	Ini_Matri(2,N*M,0.00,&EVnST[0][0]);        //Initiate the EVnST and t. I don't need to Iniate the EVnST1, k0,k1,k2,k3 and aNoise.    	  
 //$$ 4. Change the Trigger Time,length, and number of rows $$$$$$$$$$$$$$$$*//  
  while (t<TimLen) { 
	//if(fmod(t+EPSILON2,5)<10*EPSILON2)     //because t=5.0 or 10.0, the remainder might be 5.0, which might be due to some very small error, because t is a variable. 
		//  printf("time=%lf v0=%lf v1=%lf v9=%lf v99=%lf\n",t,EVnST[0][100],EVnST[0][301],EVnST[0][1909],EVnST[0][19999]);
	//MaxColumn=LargestCluster(tauNoise, N, M, N,M, 0, tau_off, tau, 5);
   
    if(fabs(t-0.00)<EPSILON2){     //Because t is a variable, it might have some little error.   
		for (i=0;i<M;i++)
		 EVnST[0][i]=1;  	   
       		    }   // Trigger: pulse v:1;        
	
    dxdt_Tri_Top(&EVnST[0][0],&k0[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay); 
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
		
		
        
      dxdt_Tri_Top(&EVnST1[0][0],&k1[0][0],aNoise,tauNoise,cNoise,gamma,coefx,coefy,N,M,Leng,deltax,deltay);  
		#pragma omp parallel for private(i)
		for(j=0;j<N*M;j++){
			for(i=0;i<2;i++)
	       EVnST1[i][j]=EVnST[i][j]+deltat*k1[i][j]; }
     
        
       dxdt_Tri_Top(&EVnST1[0][0],&k2[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay);
        #pragma omp parallel for private(i)
        for(j=0;j<N*M;j++)  {
			for(i=0;i<2;i++)
            EVnST1[i][j]=EVnST[i][j]+2*deltat*k2[i][j];  }   //E1,V1... calculated by k0;       

       dxdt_Tri_Top(&EVnST1[0][0],&k3[0][0],aNoise,tauNoise,cNoise,gamma, coefx,coefy,N,M,Leng, deltax,deltay); 
        #pragma omp parallel for private(i)
        for(j=0;j<N*M;j++)  {
			for(i=0;i<2;i++)
           EVnST[i][j]=EVnST[i][j]+deltat*(k0[i][j]+2*k1[i][j]+2*k2[i][j]+k3[i][j])/3;    }  //E1,V1... calculated by k0      
                
        t+=2*deltat;        // Runge Kutta Algorithm;  
        
        if(equal(t,2*deltat,deltat/10)){
			WB_t=waveback(&EVnST[0][0], N, M, 0, M-1, v_th);
			WF_t=wavefront(&EVnST[0][0], N, M, 0, M-1, v_th);  
			WL_t=wavelength(&EVnST[0][0], N, M, 25, v_th);
			
			WF=WF_t;
			WL=WL_t;
			
			for (i=1;i<WinHei+1;i++)
				for(j=0; j<M; j++)
					Tij[i-1][j]=dt_write*(EVnST[0][i*M+j]> v_th);
					                 }
        
       else if(WheMultiple(t, dt_write, deltat/10)){
			WB_t=waveback(&EVnST[0][0], N, M, 0, M-1, v_th);
			WF_t=wavefront(&EVnST[0][0], N, M, 0, M-1,v_th);
			WL_t=wavelength(&EVnST[0][0], N, M, 25, v_th);
			
			WF=WF*(WF>WF_t)+ WF_t*(WF_t >= WF);
			WL=WL*(WL>WL_t)+ WL_t*(WL_t>= WL);
			
			for (i=1;i<WinHei+1;i++)
				for(j=0; j<M; j++)
					Tij[i-1][j]+= dt_write*(EVnST[0][i*M+j]> v_th);
			  }
			  
		 if(print==1)
         print1_1(N,M,Leng,t,dt_write,deltat, &EVnST[0][0],fptWave);
			  
	   //if(WheMultiple(t, dt_write, deltat/10)|| equal(t,0.5,deltat/10) ){
			//for (i=1; i<1+WinHei; i++)
			  //for (j=0;j<M;j++)
			   //wave[i-1+ (int)(t/dt_write+EPSILON2) *WinHei][j]=EVnST[0][i*M+j];}
   
             			                      
  }//end of while loop;
  
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
  
  for (i=1;i<WinHei+1;i++)
				for(j=0; j<M; j++)
					Tij[i-1][j]= Tij[i-1][j]* (tauNoise[i*M+j]==tau); 
  
  WF_ab[ia][itau]+=WF;
  WL_ab[ia][itau]+=WL;
  PhiFir_ab[ia][itau]+= (double)NFir(& Tij[0][0], WinHei, M, T_th)/WinHei/M;
  
  if (NFir(&Tij[0][0], 5, M, T_th)==0)
  TraEff_ab[ia][itau]+=0;
  else  
  TraEff_ab[ia][itau]+= (double)NFir(&Tij[30][0], 5, M, T_th)/(double)NFir(&Tij[0][0], 5, M, T_th);
  
  T_ab[ia][itau]+= TMean(& Tij[0][0], WinHei, M, T_th);
  
   
  
 if(ip==NumRep-1){
	 
	printf("phi=%lf rho=%lf a=%lf, aLow=%lf, b=%lf, tau_off=%lf, c=%lf, cLow=%lf\n",phi,rho,aglobal,aLow, tau,tau_off, cglobal, cLow);
	  //Print the Tij for the last biofilm
  for(i=0; i< WinHei; i++){ 
      for(j=0; j<M; j++){ 	
				if(j==M-1)
		        fprintf(fpt,"%lf\n",Tij[i][j]);	
		        else
				fprintf(fpt,"%lf ",Tij[i][j]); } }   
				
           
    
   //strcpy(filename, "Fit14.2_test");    
      //strcat(filename,"Phi");   sprintf(num, "%d",(int)(phi*100));strcat(filename,num); 
      //strcat(filename,"Rho");   sprintf(num, "%d",(int)(rho*100));strcat(filename,num); 
      //strcat(filename,"a"); sprintf(num, "%d",(int)(aglobal*1000)); strcat(filename,num); 
      //strcat(filename,"b");   sprintf(num, "%d",(int)(tau*10000));strcat(filename,num);
         //strcat(filename,"bL");  sprintf(num, "%d",(int)(tau_off*10000));strcat(filename,num);
      //strcat(filename,"c");   sprintf(num, "%d",(int)(cglobal*10));strcat(filename,num);      

   //strcpy(filename2, filename);
   //strcat(filename,".dat");
   //strcat(filename2,"Mp.dat");
   
    ////fpt1=fopen(filename,"w");
	////for (i=0;i<WinHei*TimStep;i++){
		////for(j=0;j<M;j++) 					
		        ////fprintf(fpt1,"%lf ",wave[i][j]);
		        ////fprintf(fpt1,"\n");	
		        		        ////}
    ////fclose(fpt1);
     
    //fpt2=fopen(filename2,"w");
  
         //for (i=0;i<N*M;i++){ 	
				//if(i%M==M-1)
		        //fprintf(fpt2,"%lf\n",aNoise[i]);	
		        //else
				//fprintf(fpt2,"%lf  ",aNoise[i]);}
         //for (i=0;i<N*M;i++){ 	
				//if(i%M==M-1)
		        //fprintf(fpt2,"%lf\n",tauNoise[i]);	
		        //else
				//fprintf(fpt2,"%lf  ",tauNoise[i]);}
		 //for (i=0;i<N*M;i++){ 	
				//if(i%M==M-1)
		        //fprintf(fpt2,"%lf\n",cNoise[i]);	
		        //else
				//fprintf(fpt2,"%lf  ",cNoise[i]);}  
				
				//fclose(fpt2); 
				}
				
    fclose(fptWave);
	fclose(fptMap);
   
}   //end of loop ic
}   //end of loop ib
}  //end of ia

} //end of loop repeats

for(ia=0; ia<Len_as; ia++){ 
for(itau=0;itau<Len_taus;itau++){ 	
				if(itau==Len_taus-1)
		        fprintf(fpt1,"%lf\n",WF_ab[ia][itau]/NumRep);	
		        else
				fprintf(fpt1,"%lf ",WF_ab[ia][itau]/NumRep);} }
				
for(ia=0; ia<Len_as; ia++){ 
for(itau=0;itau<Len_taus;itau++){ 	
				if(itau==Len_taus-1)
		        fprintf(fpt1,"%lf\n",WL_ab[ia][itau]/NumRep);	
		        else
				fprintf(fpt1,"%lf ",WL_ab[ia][itau]/NumRep);} }

for(ia=0; ia<Len_as; ia++){ 
for(itau=0;itau<Len_taus;itau++){ 	
				if(itau==Len_taus-1)
		        fprintf(fpt1,"%lf\n",PhiFir_ab[ia][itau]/NumRep);	
		        else
				fprintf(fpt1,"%lf ",PhiFir_ab[ia][itau]/NumRep);} }
	
for(ia=0; ia<Len_as; ia++){ 
for(itau=0;itau<Len_taus;itau++){ 	
				if(itau==Len_taus-1)
		        fprintf(fpt1,"%lf\n",TraEff_ab[ia][itau]/NumRep);	
		        else
				fprintf(fpt1,"%lf ",TraEff_ab[ia][itau]/NumRep);} }
				
for(ia=0; ia<Len_as; ia++){ 
for(itau=0;itau<Len_taus;itau++){ 	
				if(itau==Len_taus-1)
		        fprintf(fpt1,"%lf\n",T_ab[ia][itau]/NumRep);	
		        else
				fprintf(fpt1,"%lf ",T_ab[ia][itau]/NumRep);} }
				
}   //end of loop iphi


				
fclose(fpt);
fclose(fpt1);


printf("running this tourine costs %ld seconds of time\n",time(NULL)-tt);
return 0;	
}


void dxdt_Tri(double *EVnST,double *deri,double *aNoise,double *tauNoise,double *cNoise, double gamma, int coefx,int coefy,int N,int M,int Leng, double deltax, double deltay)
 {  
 int i=0;
 #pragma omp parallel
{
#pragma omp for
for(i=0;i<N*M;i++){  
   *(deri+Leng+i)=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                      //dwdt;
   
   *(deri+i)=cNoise[i]*(EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i]);   //dvdt: This step also initiates the deri matrix.
   
   if(i/M==0)    //v_yy
	 *(deri+i)+=(EVnST[i+M]-2*EVnST[i])/pow(deltay,2);
   else if(i/M==N-1)	  
	 *(deri+i)+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);
   else
     *(deri+i)+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2); 	 
    
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
   *(deri+Leng+i)=1/tauNoise[i]*(EVnST[i]-gamma*EVnST[Leng+i]);                      //dwdt;
   
   *(deri+i)=cNoise[i]*(EVnST[i]*(1-EVnST[i])*(EVnST[i]-aNoise[i])-EVnST[Leng+i]);   //dvdt: This step also initiates the deri matrix.
   
   if(i/M==0)
	 *(deri+i)+=(EVnST[i+M]-EVnST[i])/pow(deltay,2);   //Reflective on the top
   else if(i/M==N-1)	  
	 *(deri+i)+=(0-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2);  //Obsorbing on the bottom
   else
     *(deri+i)+=(EVnST[i+M]-EVnST[i])/pow(deltay,2)+(EVnST[i-M]-EVnST[i])/pow(deltay,2); 	 
    
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
   {int mark[201][201],irdom=0,jrdom=0,badposition=0,Nrdom=0,i=0,j=0,Ntotal=0;
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
{int cluster[LENG2MAT],neighbor[LENG2MAT],NClu,NNei,NNeiE,NHigh=0,i;
	for(i=0;i<N*M;i++)  cluster[i]=-1;       //Initiate it before the first use.
	for(i=0;i<N*M;i++)  neighbor[i]=-1;
	
  while(NHigh<phi*N*M){         //compare two doubles: Be careful for possible calculation errors of computer.
	 i=0;while(cluster[i]!=-1&&i<N*M)	  //Initiate the cluster. -1 means the end of effective part. -2 means that position is concealed.
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
		
		while((double)rand()/(MaxRandom+0.1)<psi&&NHigh<phi*N*M){     
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
				   if(neighbor[i]!=-2&&(double)rand()/(MaxRandom+0.1)<1.0/(double)NNeiE){
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

double pon, poff, ponon,ponoff, poffon, poffoff, pononon, pononoff, ponoffon, ponoffoff, possibility[WIDTH],Prenvth[150*WIDTH];
 int i,j,m, HeiMat=150, WidMat=WIDTH;
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
	         if((double)rand()/(MaxRandom+0.1)<phi)
			    return high;
			    else
			    return low;
}

int LargestCluster(double *nvth, int HeiMat, int WidMat, int height,int width, int row, double DLow, double DHigh, int square){
	int j,state[HEIGHT][WIDTH]={{0}},cluster[HEIGHT][WIDTH]={{0}}, clustersize[WIDTH];
	int MaxColumn=0;
	
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
	RanNum=(double)rand()/(MaxRandom);
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
int i,j, WF, WB;
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
			  
    
