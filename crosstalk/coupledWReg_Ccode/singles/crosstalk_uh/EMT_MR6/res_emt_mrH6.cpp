# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string>
# include <cstring>
# include <mpi.h>
# include <fstream>
# include <iostream>



struct bioLevels{
	 double u,mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh;};
int factorial(int num);
int combination(int k, int n);
double H(double X,double X0,int nX,double lamdaX);
double M(int i,int n,double x,double x0);
double L(double X,double X0,int nX,double liX[] );
double Yu(double X,double X0,double nX,double yuiX[]);
double Ym(double X,double X0,int nX,double ymiX[]);
double CompRmt(double yyX,double gnX,double hX,double h0rmtX,int nhmtX,double AX,double A0rmtX, int namtX);
double CompRn(double g0X,double hX,double h0rnX,int nhnX,double ghnX,double AX,double A0rnX,double ganX, int nanxX);
double R(double R1,double R2);
double uF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double mzF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double ZF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double msF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double u3F(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double SF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double AF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double RmtF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double RnoxF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double hF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
double mhF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh);
void runSimulation(char* finame,char* frname, int numICS,int my_rank, int comm_sz);


bool uhFlag=false;
//Global parameters //!!!
double gz = 100. 	;
double gu = 2100.; 	
double grn = 40.;
double grm =150.;
double gms =90. ;
double gu3 = 1350. ;
double gs = 100.;
double ga = 30. ;
double gn = 0.2;
double g1 = 5. ;
double g2 =  0.2;
double gmh=10.;
double gh=1.5;
double gmz = 11. 	;
double y =8.;

double ks = 0.125;
double ka = 0.2 ;
double kz=0.1 		;
double ku = 0.05 	;
double kmz=0.5		;
double kms =0.5 ;
double ku3 =0.05;
double krm = 5. ;
double krn = 5.;
double kmh=.143;
double kh=1.75;

double ymi[7]={0,0.04, 0.2,1.0,1., 1.,1.};
double yui[7]={0,0.005,0.05,0.5,0.5,0.5,0.5};
double li[7] = {1.,0.6,0.3,0.1,0.05,0.05,0.05};
double ymih[3]={0,0.04, 0.2};
double yuih[3]={0,0.005,0.05};
double lih[3] = {1.,0.6,0.3};

double Z0u = 220000. ;  
double I0m = 50000. ;
double S0ms = 200000.;
double u30 = 10000. ;
double S0u3 = 300000. ;
double h0ha= 250.;
double A0aa= 350.;
double h0hrm = 200. ;
double h0hrn = 250. ;
double A0rn = 150. ;
double A0rm = 150. ;
double A0aR = 350. ;
double R0ra =  100. ;
double Z0u3 = 600000. ;
double Z0m = 25000.	;
double S0u=180000.	;
double S0m=180000.	;
double u0=10000.	;
double h0u = 200.;
double A0u=300.;
double A0m =300.;
double A0ms= 300.;
double h0ms = 200.;
double nu30rn = 10000.;
double nu30rm = 10000.;
double A0ah = 250.;
double R0rh = 300.;
double h0hh = 80. ;

//double lamdahu = 1.5  ;
//double lamda3n = 10.;
//double lamdaAms = 0.4 ; 
//double lamdahms = 7. ;
//double lamda3m= 2.;
//double lamdaAu = 0.6;
//double lamdaAm=0.5;
double lamdahu = 1.  ;
double lamdaAu = 1.;
double lamda3n = 1.;
double lamda3m= 1.;
double lamdaAms = 1. ; 
double lamdahms = 1. ;
double lamdaAm=1.;

double lambdarh = 0.2;
double lambdahh = 0.1;
double lambdaah = 0.1;
double lambdaar = 0.25;
double lambdara = 8.;
double lambdaha= 0.1;
double lambdaaa = 0.2;

double lamdazu=0.1	;
double lamdaSu=0.1	;
double lamdaIm =10.;
double lamdaSms = 0.1;
double lamdazu3 = 0.2;
double lamdaSu3 = 0.1;
double lamdaZm=7.5	;
double lamdaSm=10.	;

int naa = 2;
int nu3 = 2;
int nzu3 = 2;
int nsms = 1;
int nIm = 2;
int nsu3 = 1;
int nra= 4;
int nha= 1;
int nar = 2;
int narm = 4;
int nhrm = 2;
int nhrn = 2;
int nzu =3;
int nsu=2;
int nzm=2;
int nsm=2;
int nu=6;
int narn = 2;
int nah = 1;
int nrh = 4;
int nhh = 4;
int nuh=2;
int nhu = 1;
int nAu = 1;
int nAm=2;
int nAms =2;
int nhms = 2 ;
int n3m=3;
int n3n=2;
//////////

double H(double X,double X0,int nX,double lamdaX){
	return lamdaX+(1.-lamdaX)/(1.+pow((X/X0),nX));}

double M(int i,int n,double x,double x0){
	return pow((x/x0),i)/pow((1+(x/x0)),n);}
double L(double X,double X0,int nX,double liX[]){
	double total=0;
	for(int i=0; i<=nX; i++)
		total= total+ liX[i]*combination(i,nX)*M(i,nX,X,X0);
	return total;
	}
double Ym(double X,double X0,int nX,double ymiX[]){
	double total=0;
	for(int i=0; i<=nX; i++)
		total=total+ ymiX[i]*combination(i,nX)*M(i,nX,X,X0);
	return total;
	}
double Yu(double X,double X0,double nX,double yuiX[]){
	double total=0;
	for(int i=0; i<=nX; i++)
		total= total+i*yuiX[i]*combination(i,nX)*M(i,nX,X,X0);
	return total;}
double CompRmt(double yyX,double gnX,double hX,double h0rmtX,int nhmtX,double AX,double A0rmtX, int namtX){
	return yyX*(gnX+pow((AX/A0rmtX),namtX))/(1.+pow((hX/h0rmtX),nhmtX)+pow((AX/A0rmtX),namtX));}
double CompRn(double g0X,double hX,double h0rnX,int nhnX,double ghnX,double AX,double A0rnX,double ganX, int nanxX){
	return (g0X+ghnX*pow((hX/h0rnX),nhnX)+ganX*pow((AX/A0rnX),nanxX))/(1.+pow((hX/h0rnX),nhnX)+pow((AX/A0rnX),nanxX));

}
double R(double R1,double R2){
	return R1+R2;}
int combination(int k, int n){
	return factorial(n)/factorial(k)/factorial(n-k);
}
int factorial(int num){
	int fact=1;
	for(int i=num;i>0;i--)
		fact = fact*i;
	return fact;
}			
/*
*H(h,h0u,nhu,lamdahu)
*H(A,A0u,nAu,lamdaAu)
*H(A,A0m,nAm,lamdaAm)
*H(A,A0ms,nAms,lamdaAms)
*H(h,h0ms,nhms,lamdahms)
*H(u3,nu30rn,n3n,lamda3n)
*H(u3,nu30rm,n3m,lamda3m)
u->H
*/

double uF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	if (uhFlag){
		return gu*H(Z,Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)-mz*Yu(u,u0,nu,yui)-mh*Yu(u,u0,nuh,yuih)-ku*u;}
	return gu*H(Z,Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)-mz*Yu(u,u0,nu,yui)-mh*Yu(0,u0,nuh,yuih)-ku*u;}
double mzF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return gmz*H(Z,Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)*H(A,A0m,nAm,lamdaAm)-mz*Ym(u,u0,nu,ymi)-kmz*mz;}
double ZF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return gz*mz*L(u,u0,nu,li)-kz*Z;}
double msF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	double I = 50000.;
	return gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)*H(h,h0ms,nhms,lamdahms)*H(A,A0ms,nAms,lamdaAms)-ms*Ym(u3,u30,nu3,ymi)-kms*ms;}
double u3F(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return gu3*H(Z,Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu(u3,u30,nu3,yui)-ku3*u3;}
double SF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return gs*ms*L(u3,u30,nu3,li)-ks*S;}
double AF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A;}
double RmtF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt*H(u3,nu30rm,n3m,lamda3m);}
double RnoxF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	return grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox*H(u3,nu30rn,n3n,lamda3n);}
double hF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	if (uhFlag){
		return gh*mh*L(u,u0,nuh,lih)-kh*h;}
	return gh*mh*L(0,u0,nuh,lih)-kh*h;}
double mhF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){
	if (uhFlag){
		return gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ym(u,u0,nuh,ymih);}
	return gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ym(0,u0,nuh,ymih);}


bioLevels getICS(){
	int mz_range=2000;
	int S_range=250000;
	int ms_range=1000;
	int u3_range=20000;
	int Z_range=700000;
	int u_range=25000;
	int A_range=1000;
	int h_range=1000;
	int Rmt_range=1000;
	int Rnox_range=1000;
	int mh_range=1000;

	bioLevels ics = {rand()%u_range,rand()%mz_range,rand()%Z_range,rand()%ms_range,rand()%u3_range,rand()%S_range,rand()%A_range,rand()%Rmt_range,rand()%Rnox_range,rand()%h_range,rand()%mh_range};
	return ics;
}




bioLevels run(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh){

	int total=10000;

	double dt=0.1;
	int times = total*10;
	for(int i=0; i<times; i++){

		double utmp =     uF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double mztmp =    mzF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double ztmp =     ZF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double mstmp =    msF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double u3tmp =    u3F(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double Stmp =     SF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double Atmp =     AF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double rmttmp =   RmtF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double rnoxtmp =  RnoxF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double htmp =     hF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		double mhtmp =     mhF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		
		u =    dt*utmp +u;
		mz =   dt*mztmp +mz;
		Z =    dt*ztmp +Z;
		ms =   dt*mstmp +ms;
		u3 =   dt*u3tmp +u3;
		S =    dt*Stmp +S;
		A =    dt*Atmp +A;
		Rmt =  dt*rmttmp +Rmt;
		Rnox = dt*rnoxtmp +Rnox;
		h =    dt*htmp +h;
		mh =    dt*mhtmp +mh;
	}
	
	bioLevels res = {u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh};
	return res;
}

void runSimulation(char* finame,char* frname, int numICS,int my_rank, int comm_sz){
	
	MPI_File f_i;
	MPI_File f_r;
	strcat(finame,"_");
	std::string s = std::to_string(numICS);
	strcat(finame,s.c_str());
	strcat(finame,"_ics.txt");
	strcat(frname,"_");
	strcat(frname,s.c_str());
	strcat(frname,"_res.txt");

	MPI_File_open(MPI_COMM_WORLD,frname,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_r);
	MPI_File_open(MPI_COMM_WORLD,finame,MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_i);

	int number;
	if(my_rank==0){
		char tempi_string[1000]="";
		char tempr_string[1000]="";
		snprintf(tempi_string,sizeof(tempi_string),"u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh\n");
		snprintf(tempr_string,sizeof(tempr_string),"u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh\n");
		MPI_File_write_shared(f_i, tempi_string, strlen(tempi_string), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(f_r, tempr_string, strlen(tempr_string), MPI_CHAR, MPI_STATUS_IGNORE);
                number=2;
                for(int new_rank=1;new_rank<comm_sz;new_rank++){
                        MPI_Send(&number,1,MPI_INT,new_rank,0,MPI_COMM_WORLD);}
		fflush(stdout);
        }
        else{
                // wait until the header has been written
                MPI_Recv(&number,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}

        fflush(stdout);
	int timeL=numICS/(comm_sz);

	if(my_rank==comm_sz-1)
		timeL = numICS-numICS/comm_sz*(comm_sz-1);
	
	for(int i=0; i< timeL; i++){
		char tempi_string[1000]="";
		char tempr_string[1000]="";
		bioLevels ics = getICS();
		double u = ics.u;
		double mz = ics.mz;
		double Z = ics.Z;
		double ms = ics.ms;
		double u3 = ics.u3;
		double S = ics.S;
		double A = ics.A;
		double Rmt = ics.Rmt;
		double Rnox = ics.Rnox;
		double h = ics.h;
		double mh = ics.mh;

		bioLevels res = run(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);

		snprintf(tempi_string,sizeof(tempi_string),"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh);
		snprintf(tempr_string,sizeof(tempr_string),"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",res.u,res.mz,res.Z,res.ms,res.u3,res.S,res.A,res.Rmt,res.Rnox,res.h,res.mh);
	
		MPI_File_write_shared(f_i, tempi_string, strlen(tempi_string), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(f_r, tempr_string, strlen(tempr_string), MPI_CHAR, MPI_STATUS_IGNORE);

	}
	fflush(stdout);
	MPI_File_close(&f_i);
	MPI_File_close(&f_r);

}

int  main(void){

	//Initialize everythign necessary for mpi
	int my_rank;
	int comm_sz;
	MPI_Init(NULL,NULL);
	//get num processes
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	//get rank
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	std::ifstream seedFile("seedsForSims.txt");
	int seeds[300];
	int count=0;
	std::string line;
	while(std::getline(seedFile,line)){
	   seeds[count] = atoi(line.c_str());
           count++;
	   }//end while
	seedFile.close();

	/// now run the sims
	lamdahu=1.;
	lamdaAu=1.;
	lamdaAm=1.;
	lamdaAms=1.;
	lamdahms=1.;
	lamda3n=1.;
	lamda3m=1.;
	char finame[100]="";
	char frname[100]="";
	std::string string_addto = "";

	uhFlag=true;
	/// u -> H // double check

	/// u -> H // double check
	for(int i=8; i<11; i++){
	    for(int j=8; j<11; j++){
   	        for(int k=0; k<11; k++){
   	        	for(int l=1; l<11; l++){

				// MR1 : yu=[0,0,0]
				// MR2 : yu=[0,0.01,0.09]
				// MR3/MR6 : yu=[0,0.1,0.2]
				// MR4 : yu=[0,0.5,1.0]
				strcpy(finame,"EMT_MR6_lamdauH_");
				strcpy(frname,"EMT_MR6_lamdauH_");
				std::string s;

				lih[0]=1.;
				ymih[0]=0.;

				yuih[0]=0.;
				yuih[1]=0.1;
				yuih[2]=0.2;
				double tmpi=i*0.1;
				double tmpj=j*0.1;
				double tmpk = k*0.002+0.001;
				double tmpl=l*0.004+0.001;

				if (i==0 && j==0){
					lih[0]=0.;
					lih[1]=0.;
					lih[2]=0.;
				}
				if(tmpi<tmpj){
					continue;
				}
				if(tmpk>tmpl){
					continue;
				}
				
				lih[1]=tmpi;
				lih[2]=tmpj;
				ymih[1]=tmpk;
				ymih[2]=tmpl;
                        
				s = std::to_string(i);
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
				s = "_";
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
				s = std::to_string(j);
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
				s = "_";
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
				s = std::to_string(k);
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
				s = "_";
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
				s = std::to_string(l);
				strcat(finame,s.c_str());
				strcat(frname,s.c_str());
                        
                        
				srand(seeds[my_rank]);
				runSimulation(finame,frname,1000,my_rank,comm_sz);
		    	}//end for l
		    }//end for k
		}//end for j
	}//end for i



	MPI_Finalize();//clean up mpi
	return 0;
}
