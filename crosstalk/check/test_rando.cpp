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
bioLevels getICS();
//FIX
//seeds = pd.read_csv("seedsForSims.txt",header=None)
//np.random.seed(seeds.values[myrank][0])
//int combination(int k,int n){
//	return math.factorial(n)/math.factorial(k)/math.factorial(n-k);}

int  main(void){

	MPI_Init(NULL,NULL);
	int my_rank,comm_sz;
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	
	std::ifstream seedFile("seedsForSims.txt");
	int seeds[300];

	int count=0;
	std::string line;
	while(std::getline(seedFile,line)){
	   seeds[count] = atoi(line.c_str());
	   if(my_rank==0)
		   std::cout <<seeds[count]<<'\n';
	    count++;
	   }
	seedFile.close();

	count=0;
	for(int i=0;i<1000000;i++){
		count=count*i+i;}

	if(my_rank==1){
	count=0;
	for(int i=0;i<1000000;i++){
		count=count*i+i;}
	}
	if(my_rank==2){
	count=0;
	for(int i=0;i<2000000;i++){
		count=count*i+i;}
	}

	std::cout<<"u,mz,Z,ms,u3,S,A,Rmt,Rnox,h,mh\n";
	for(int i=0;i<10;i++){
		std::cout<<i<<','<<my_rank<<","<<my_rank+i*3<<","<<seeds[my_rank+i*3]<<'\n';
		srand(seeds[my_rank+i*3]);
		bioLevels ics=getICS();
		for (int k=0;k<10;k++){
			std::cout<<ics.u<<","<<ics.mz<<","<<ics.Z<<","<<ics.ms<<","<<ics.u3<<","<<ics.S<<","<<ics.A<<","<<ics.Rmt<<","<<ics.Rnox<<","<<ics.h<<","<<ics.mh<<"\n";
		}
	}


	//std::cout<< factorial(5)<<std::endl;
	//std::cout<< combination(2,5)<<std::endl;

	MPI_Finalize();

	return 0;
}


int factorial(int num){
	int fact=1;
	for(int i=num;i>0;i--)
		fact = fact*i;
	return fact;
}			

int combination(int k, int n){
	return factorial(n)/factorial(k)/factorial(n-k);
}

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
