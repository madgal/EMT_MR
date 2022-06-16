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
	 double u,mz, Z, ms, u3, S, A, Rmt, Rnox, h, mh,G,O,mg,mo;};
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
double uF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double mzF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double ZF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double msF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double u3F(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double SF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double AF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double RmtF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double RnoxF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double hF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double mhF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double mOF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double mGF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double tF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double OF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
double GF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo);
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
double ymih[3]={0.,0.,0.};//[0,0.04, 0.2}
double yuih[3]={0,0.005,0.05};
double lih[3] = {1.,0.,0.};//[1.,0.6,0.3]

double I = 50000.;
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
//
double o0u=250000.;
double o0z=25000.;
double o0mo=25000.;
double Z0mo=10000.;
double Z0mg=10000;
double g0mz=25000;
/////////////

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

double lou=0.1;
double loz=0.1;
double lomo=0.1;
double lzmo=0.5;
double lzmg=0.5;
double lgmz=0.1;

double rmtprod=1.;
double rnoxprod=1.;

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
//
int nou=1;
int noz=1;
int nomo=2;
int    nzmo=1;
int    nzmg=3;
int    ngmz=1;
//////////

double gg=200;//200 from paper or 70??
double gmg=22.;// 22 from paper or 8??
double gov=200.;//200 from paper  or 70??
double gmo=22.;//?? not listed prob 22 in paper  or 8??
double kov=0.1;//hr-1 
double kmo=0.5;//hr-1 
double kgg=0.1;// 
double kmg=0.5;//hr-1


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
H(h,h0u,nhu,lamdahu)
H(A,A0u,nAu,lamdaAu)
H(A,A0m,nAm,lamdaAm)
H(A,A0ms,nAms,lamdaAms)
H(h,h0ms,nhms,lamdahms)
H(u3,nu30rn,n3n,lamda3n)
H(u3,nu30rm,n3m,lamda3m)
u->H
*/
/*
H(O,o0u,nou,lou)
H(O,o0z,noz,loz)
H(O,o0mo,nomo,lomo)
H(Z,Z0mo,nzmo,lzmo)
H(Z,Z0mg,nzmg,lzmg)
H(G,g0mz,ngmz,lgmz)

H(G,G0mo,nGmo,lGmo)
H(G,G0rm,ngm,lgm)
H(G,G0rn,ngn,lgn)
*/

double o0mt=120;//25000.;
double G0mo=120;//25000.;
double G0rm=120;//25000.;
double G0rn=120;//25000.;

int nGmo = 2;
int ngm=1;
int ngn=1;

double lGmo=1.;
double lgm=1.;
double lgn=1.;

double uF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	if (uhFlag){
		return gu*H(Z,Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)*H(O,o0u,nou,lou)-mz*Yu(u,u0,nu,yui)-mh*Yu(u,u0,nuh,yuih)-ku*u;}
	return gu*H(Z,Z0u,nzu,lamdazu)*H(S,S0u,nsu,lamdaSu)*H(h,h0u,nhu,lamdahu)*H(A,A0u,nAu,lamdaAu)*H(O,o0u,nou,lou)-mz*Yu(u,u0,nu,yui)-mh*Yu(0,u0,nuh,yuih)-ku*u;}
double mzF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gmz*H(Z,Z0m,nzm,lamdaZm)*H(S,S0m,nsm,lamdaSm)*H(A,A0m,nAm,lamdaAm)*H(O,o0z,noz,loz)*H(G,g0mz,ngmz,lgmz)-mz*Ym(u,u0,nu,ymi)-kmz*mz;}
double ZF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gz*mz*L(u,u0,nu,li)-kz*Z;}
double msF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gms*H(S,S0ms,nsms,lamdaSms)*H(I,I0m,nIm,lamdaIm)*H(h,h0ms,nhms,lamdahms)*H(A,A0ms,nAms,lamdaAms)-ms*Ym(u3,u30,nu3,ymi)-kms*ms;}
double u3F(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gu3*H(Z,Z0u3,nzu3,lamdazu3)*H(S,S0u3,nsu3,lamdaSu3)-ms*Yu(u3,u30,nu3,yui)-ku3*u3;}
double SF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gs*ms*L(u3,u30,nu3,li)-ks*S;}
double AF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return ga*H(R(Rmt,Rnox),R0ra,nra,lambdara)*H(h,h0ha,nha,lambdaha)*H(A,A0aa,naa,lambdaaa)-ka*A;}
double RmtF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return rmtprod*grm*H(A,A0aR,nar,lambdaar)*CompRmt(y,gn,h,h0hrm,nhrm,A,A0rm,narm)-krm*Rmt*H(u3,nu30rm,n3m,lamda3m)*H(G,G0rm,ngm,lgm);}
double RnoxF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return rnoxprod*grn*CompRn(gn,h,h0hrn,nhrn,g1,A,A0rn,g2,narn)-krn*Rnox*H(u3,nu30rn,n3n,lamda3n)*H(G,G0rn,ngn,lgn);}
double hF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	if (uhFlag){
		return gh*mh*L(u,u0,nuh,lih)-kh*h;}
	return gh*mh*L(0.,u0,nuh,lih)-kh*h;}
double mhF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	if (uhFlag){
		return gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ym(u,u0,nuh,ymih);}
	return gmh*H(A,A0ah,nah,lambdaah)-kmh*mh*H(h,h0hh,nhh,lambdahh)*H(R(Rmt,Rnox),R0rh,nrh,lambdarh)-mh*Ym(0.,u0,nuh,ymih);}
double mOF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gmo*H(O,o0mo,nomo,lomo)*H(Z,Z0mo,nzmo,lzmo)*H(G,G0mo,nGmo,lGmo)-kmo*mo;}
double mGF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gmg*H(Z,Z0mg,nzmg,lzmg)-kmg*mg;}
double OF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gov*mo-kov*O;}
double GF(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){
	return gg*mg-kgg*G;}


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
	int G_range=50000;
	int O_range=1000;//FIX
	int mg_range=1000;//FIX
	int mo_range=1000;//FIX
	
	double rf_u = rand()%u_range;
	double rf_mz=rand()%mz_range;
	double rf_Z=rand()%Z_range;
	double rf_ms=rand()%ms_range;
	double rf_u3=rand()%u3_range;
	double rf_S=rand()%S_range;
	double rf_A=rand()%A_range;
	double rf_Rmt=rand()%Rmt_range;
	double rf_Rnox=rand()%Rnox_range;
	double rf_h=rand()%h_range;
	double rf_mh=rand()%mh_range;
	double rf_G=rand()%G_range;
	double rf_O=rand()%O_range;
	double rf_mg=rand()%mg_range;
	double rf_mo=rand()%mo_range;

	bioLevels ics = {rf_u,rf_mz,rf_Z,rf_ms,rf_u3,rf_S,rf_A,rf_Rmt,rf_Rnox,rf_h,rf_mh,rf_G,rf_O,rf_mg,rf_mo};
	return ics;
}




bioLevels run(double u,double mz,double Z,double ms,double u3,double S,double A,double Rmt,double Rnox,double h,double mh,double G,double O,double mg,double mo){

	int total=10000;

	double dt=0.1;
	int times = total*10;
	for(int i=0; i<times; i++){

		double utmp =     uF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double mztmp =    mzF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double ztmp =     ZF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double mstmp =    msF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double u3tmp =    u3F(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double Stmp =     SF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double Atmp =     AF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double rmttmp =   RmtF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double rnoxtmp =  RnoxF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double htmp =     hF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double mhtmp =    mhF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double Gtmp =     GF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double Otmp =     OF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double mgtmp =    mGF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		double motmp =    mOF(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		
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
		G =    dt*Gtmp +G;
		O =    dt*Otmp +O;
		mg =    dt*mgtmp +mg;
		mo =    dt*motmp +mo;
	}
	
	bioLevels res = {u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo};
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

	int number=1;
	if(my_rank==0){
		char tempi_string[1000]="";
		char tempr_string[1000]="";
		snprintf(tempi_string,sizeof(tempi_string),"u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo\n");
		snprintf(tempr_string,sizeof(tempr_string),"u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo\n");
		MPI_File_write_shared(f_i, tempi_string, strlen(tempi_string), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(f_r, tempr_string, strlen(tempr_string), MPI_CHAR, MPI_STATUS_IGNORE);
		number=0;
		for(int i=1; i<comm_sz;i++){
	        	MPI_Send(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		fflush(stdout);
	}
	else{
    		MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	int timeL=numICS/(comm_sz);

	if(my_rank==(comm_sz-1))
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
		double G = ics.G;
		double O = ics.O;
		double mg = ics.mg;
		double mo = ics.mo;

		bioLevels res = run(u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);

		snprintf(tempi_string,sizeof(tempi_string),"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh,G,O,mg,mo);
		snprintf(tempr_string,sizeof(tempr_string),"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",res.u,res.mz,res.Z,res.ms,res.u3,res.S,res.A,res.Rmt,res.Rnox,res.h,res.mh,res.G,res.O,res.mg,res.mo);
	
		MPI_File_write_shared(f_i, tempi_string, strlen(tempi_string), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(f_r, tempr_string, strlen(tempr_string), MPI_CHAR, MPI_STATUS_IGNORE);

	}
	fflush(stdout);
	MPI_File_close(&f_i);
	MPI_File_close(&f_r);

}

std::string getCoupledParameters(int a){
    std::string tmp_return="";
    std::string s;        

    uhFlag=false;
    if (a==1){
    	uhFlag=true;

	lih[0]=1.0;
	lih[1]=0.9;
	lih[2]=0.8;

	ymih[0]=0.;
	ymih[1]=0.002;
	ymih[2]=0.01;

	yuih[0]=0.;
	yuih[1]=0.001;
	yuih[2]=0.009;

	tmp_return+="_uh_310";
     }

    int i1,i2,i3,i4,i5,j1,j2,k;
    int i1b=0,i2b=0,j1b=0,j2b=0;
    switch(a){
	case 0://none
		i1=10.;//AS
		i2=10.;//AZ
		i3=10.;//Hu
		i4=10.;//u3m
		i5=10.;//u3n
		j1=1.;//HS
		j2=1.;//Au
		k=5;// (input/10000.)
		break;
	case 1://all
		i1=9.;//AS
		i1b=5.;//AS
		i2=9.;//AZ
		i2b=5.;//AZ
		i3=1.;//Hu
		i4=2.;//u3m
		i5=1.;//u3n
		j1=1.;//HS
		j1b=1.;//HS
		j2=1.;//Au
		j2b=1.;//Au
		k=0;// (input/10000.)
		break;
	case 2://only ROS
		i1=10.;//AS
		i2=10.;//AZ
		i3=10.;//Hu
		i4=0.1;//u3m
		i5=0.1;//u3n
		j1=1.;//HS
		j2=1.;//Au
		k=5;// (input/10000.)
		break;
	case 3://Hu and input
		i1=10.;//AS
		i2=10.;//AZ
		i3=0.;//Hu
		i4=10.;//u3m
		i5=10.;//u3n
		j1=1.;//HS
		j2=1.;//Au
		k=2;// (input/10000.)
		break;
	}//end switch

	lamdaAms=i1*0.1+i1b*0.01;
	if (lamdaAms!=1.){
    	tmp_return="_AS_0_";
	s = std::to_string(i1);
	tmp_return+=s.c_str();
	s = std::to_string(i1b);
	tmp_return+=s.c_str();}

	lamdaAm=i2*0.1+i2b*0.01;
	if (lamdaAm!=1.){
   	tmp_return="_AZ_0_";
	s = std::to_string(i2);
	tmp_return+=s.c_str();
	s = std::to_string(i2b);
	tmp_return+=s.c_str();}

	lamdahu=i3*0.1;
	s = std::to_string(i3);
	if (lamdahu!=1.){
    	tmp_return="_Hu_0_";
	tmp_return+=s.c_str();}

	lamda3m=i4*0.1;
	s = std::to_string(i4);
	if (lamda3m!=1.){
    	tmp_return="_u3m_0_";
	tmp_return+=s.c_str();}

	lamda3n=i5*0.1;
	s = std::to_string(i5);
	if (lamda3n!=1.){
    	tmp_return="_u3n_0_";
	tmp_return+=s.c_str();}

	lamdahms=j1*1.+j1b*0.1;
	if (lamdahms!=1.){
        tmp_return="_HS_";
	s = std::to_string(j1);
	tmp_return+=s.c_str();
	tmp_return+="_";
	s = std::to_string(j1b);
	tmp_return+=s.c_str();}

	lamdaAu=j2*1.+j2b*0.1;
	if (lamdaAu!=1.){
        tmp_return="_Au_";
	s = std::to_string(j2);
	tmp_return+=s.c_str();
	tmp_return+="_";
	s = std::to_string(j2b);
	tmp_return+=s.c_str();}

	I=k*10000.;
	s = std::to_string(k*10000);
	if (k!=5.){
    	tmp_return="_input_";
	tmp_return+=s.c_str();}

	return tmp_return;
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
	/////////////////////////
	o0u=250000.;//
	o0z=25000.;//
	o0mo=25000.;//
	Z0mo=10000.;//
	Z0mg=10000;//
	g0mz=25000;//
	G0mo=25000.;//120;//
	G0rm=25000.;//120;//
	G0rn=25000.;//120;//

	nou=1;//
	noz=1;//
	nomo=2;//
	nzmo=1;//
	nzmg=3;//
	ngmz=1;//
	nGmo=2;
	ngm=1;
	ngn=1;

	lgm=1.;//
	lgn=1.;//

	lamdaIm=16.;
	lou=0.1;
	loz=0.1;
	lomo=0.1;
	lzmo=0.5;
	lzmg=0.5;
	lgmz=0.1;

	lGmo=1.;//.7;

	gg=200;
	gmg=22.;
	gov=200.;
	gmo=22.;
	kov=0.1;
	kmo=0.5;
	kgg=0.1;
	kmg=0.5;


	//////////////////////
	char finame[100]="";
	char frname[100]="";
	std::string string_addto = "";

	for(int a=0;a<1;a++){// set=0 for nothing, 1 for all, 2 for MR, 3 for EMT
	   for(int j=10;j<11;j++){//rmtproduction set =10 for nothing  below 10 for reduction, above for increase
	     for(int k=10;k<11;k++){//rnoxproduction set =10 for nothing below 10 for reduction, above for increase
	      for(int a1=10;a1<11;a1++){// gz I
	       for(int b1=5;b1<6;b1++){// zg  I
	        for(int c1=10;c1<11;c1++){// go I
	         for(int d1=10;d1<11;d1++){// grm 
	          for(int e1=10;e1<11;e1++){// grn 
	           for(int f1=1;f1<2;f1++){// ts 
	            for(int g1=1;g1<2;g1++){// trm
	             for(int h1=1;h1<2;h1++){// trn
	              for(int i1=10;i1<11;i1++){// ut
	               for(int j1=10;j1<11;j1++){// ot
	               for(int k1=10;k1<11;k1++){// oz I
	               for(int a2=10;a2<11;a2++){// zo I
	               for(int b2=10;b2<11;b2++){// oo I
	               for(int c2=10;c2<11;c2++){// ou I

			for(int i=0; i<6;i++){
			strcpy(finame,"EMT_MR_test");
			strcpy(frname,"EMT_MR_test");
			std::string str_cat=getCoupledParameters(a);
			strcat(finame,str_cat.c_str());
			strcat(frname,str_cat.c_str());


			 switch(i){
			   case 0:
				lGmo=0.9;
				lgm=.26;
				lgn=.25;
				break;
			   case 1:
				lGmo=0.8;
				lgm=.26;
				lgn=.25;
				break;
			   case 2:
				lGmo=0.7;
				lgm=.26;
				lgn=.25;
				break;
			   case 3:
				lGmo=0.6;
				lgm=.26;
				lgn=.25;
				break;
			   case 4:
				lGmo=0.5;
				lgm=.26;
				lgn=.25;
				break;
			   case 5:
				lGmo=0.4;
				lgm=.26;
				lgn=.25;
				break;
			}
			str_cat=std::to_string(i);
			strcat(finame,str_cat.c_str());
			strcat(frname,str_cat.c_str());

			srand(seeds[my_rank]);
			runSimulation(finame,frname,100,my_rank,comm_sz);
			}
	}}}//a,j,k
	}}}}}}}}}}}//a1-k1
	}}}//a2-c2

	MPI_Finalize();//clean up mpi
	return 0;
}
