/*###############################################################################
#                                                                               #
#    These mpi runner are developed for parallelized run of CORSIKA code        #
#    based on separation of secondary showers by energy threshold               #
#                                                                               #
#    v0.1  01-06-11   Sushant Sharma                                            #
#          -Implementing MPI communication with bookkeeping system              #
#          -Developing Book-Keeping and statistic system                        #
#    v0.2  02-10-11   Tanguy Pierog, Gevorg Poghosyan                           #
#          -new syntax to run mpi_runner with CORSIKA steering file only        #
#    v0.3  22-11-11   Gevorg Poghosyan                                          #
#          -Optimising the read of steering file only in master                 #
#          -Using debuger switch to reduce printouts                            #
#    v0.31 24-11-11   Gevorg Poghosyan                                          #
#          -Optimising the outputs                                              #
#          -Bug fix: master is broadcasting global variables                    #
#    v0.32 05-12-11   Gevorg Poghosyan                                          #
#          -Introducing debugger switch                                         #
#    v0.33 23-07-12   Tanguy Pierog, Gevorg Poghosyan                           #
#          -new syntax to have same names  for cut and list files               #
#    v0.33 15-04-13   Gevorg Poghosyan                                          #
#          -reducing debugging outputs                                          #
#    v0.34 07-06-13   Gevorg Poghosyan                                          #
#          -fixing the bookkeeping bugs                                         #
#          -fixing the MPI communication bugs                                   #
#          -optimized debugging system and outputs                              #
#    v0.35 22-07-13    Arashdeep Kaur                                           #
#          -sum and join the .long files into one file                          #
#          -join DATnnnnnn files by removing some data                          #
#    v0.36 18-10-13    Gevorg Poghosyan, Juergen Oehlschlaeger                  #
#          -fixing bugs in DAT and Long files joiners                           #
#    v0.37 18-10-13    Gevorg Poghosyan, Juergen Oehlschlaeger                  #
#          -advancing step size for longitudinal distribution file              #
#          -limiting generation of list output files                            #
#    v0.38 15-11-13    Gevorg Poghosyan                                         #
#          -Bugfix and optimisation for longitudinal file joiner                #
#    v0.39 28-11-13    Gevorg Poghosyan, Juergen Oehlschlaeger                  #
#          -making parallel runner compatible for use with THINNING             #
#    v0.40 23-01-14    Gevorg Poghosyan                                         #
#          -Bugfix:eliminating writeout of empty particle Buffer as a file      #
#    v0.41 07-02-14    Gevorg Poghosyan                                         #
#          -Implementation of new parallel algorithm for joint resultfiles      #
#          -Bugfix: last seeds used when not 6 of them is given in intput file  #           
#    v0.42 09-04-14    Gevorg Poghosyan                                         #
#          -Implementation of communication and compatibility with CoREAS       #
#    v0.43 11-06-14    Gevorg Poghosyan, Pranav Bisht                           #
#          -generating "corsika_timetable" per parallel slave                   #
#          -bugfix: elimination of memory leaks                                 #
#          -Barrierfree communication for logifile generation                   #
#    v0.44 11-06-14    Gevorg Poghosyan,                                        #
#          -unified structure for MPI communication                             #
#          -Barriered communication for longifile generation                    #
#    v0.5  11-07-14    Gevorg Poghosyan,Tim Huege                               #
#          -parallel algorithm of joining CoREAS data                           #
#                                                                               #
#    ToDo:                                                                      #
#         -Advanced command line options recognition                            #
#         -Optimize usage of global variables between slaves and corsika        #
#         -Reduce the usage of printouts for productive runs                    #
#                                                                               #
# Acknowledgement:                                                              #
#   Special Thanks to Tanguy Pierog, Dieter Heck and Ralph Engel for            #
#   contributing their ideas used for development of this algorithms            #
#                                                                               #
# Contact:                                                                      #
#    Gevorg.Poghosyan@kit.edu SimLab E&A Particle Physics of SCC,KIT            #
#                                                                               #
# Reference (please cite when using this program):                              #
#   Heck,D. et.al; Computational Methods in Science and Engineering (2011)      #
#                   pp 87-90, KIT Scientific Publishing                         #
#################################################################################
*/
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "config.h"

/////////////////////initializing parameters/////////////////// 
//the number of data type block in the MPI message
#define BLOCK_COUNT_MESSAGE 2
//rank of the master ---(DONT CHANGE IT)          
#define MASTER 0

#define MXLIST 200001
#define MAX_GROUP_SIZE 201
#if __MULTITHIN__
#define PARTICLE_INFO_COUNT 47
#else
#define PARTICLE_INFO_COUNT 19
#endif
#define NSTEP_MAX 15000
#define LISTFILE_MAX 1

////////////////////flags(TYPES OF MPI MESSAGES)//////////////////////
#define START 1000
#define REQUEST 1001
#define FINISH 1002
#define STOP 1003
#define VALN 2120
#define STRP 2121
#define VALP 2122
#define STRE 2123
#define VALE 2124
#define STOREVNH 3131
#define STORRUNH 3232

//////////////////error codes///////////////////////////////
#define NA -1
#define NOT_FOUND -2

//1 GByte Memory of 268435456 floats occupying each 4 bytes (use 1/10);
#define MAXBLOCKS 26830440 // 4680 records * 5733 = 4095 records * 6552; 

double *particle_stack_or_cut_file,*particles_from_stack;
int job_count, runnumber, join_counts=1;
int II1;			//Pointer of current particle in group
int NSTEP1=0;			//number of vertical steps for generating longfile
int dbg_infos=0;	//Debugging switch to show details about MPI communication
int lcout;          // LCOUT(BOOLEAN)- tells weather to generate CUTFILES or not (run always done in parallel).
int MAXBUFSLAVE=0;	// Current length of subblock of particles
char statdir[255];  // output directory name for bookkeeping: information about running status, finished free slaves etc.
float Block[MAXBLOCKS]; 
float Event_header[312], Run_header[312];
int particle_block_counts=0, maxblock_switch=0, switch_first_runH=1, switch_first_evnH=1;
struct message //structure of MPI message.
{
	int id1;        // I1CUTPAR (INT) : FIRST INDEX OF PARTICLE TO READ FROM CFILINP
	int id2;        // I2CUTPAR (INT) : LAST INDEX OF PARTICLE TO READ FROM CFILINP
	int lprim;      // (BOOLEAN) : TRUE(1), IF SHOWER IS PRIMARY || FALSE(0), IF SUBSHOWER IS SECONDARY
	int rnnum;      // RUN NUMBER
	int seed;       // SEED USED BY THE SUB SHOWER
	int uid;        // THE UNIQUE ID OF THE MPI JOB
	int parent_id;  // THE PARENT ID OF THE MPI JOB(USED TO RECONSTRUCT THE CUTFILE_NAME)
	double stack[MAX_GROUP_SIZE*(PARTICLE_INFO_COUNT+1)];
};

struct message_queue //structure of MPI message.
{
	int id1;        // I1CUTPAR (INT) : FIRST INDEX OF PARTICLE TO READ FROM CFILINP
	int id2;        // I2CUTPAR (INT) : LAST INDEX OF PARTICLE TO READ FROM CFILINP
	int lprim;      // (BOOLEAN) : TRUE(1), IF SHOWER IS PRIMARY || FALSE(0), IF SUBSHOWER IS SECONDARY
	int rnnum;      // RUN NUMBER
	int seed;       // SEED USED BY THE SUB SHOWER
	int uid;        // THE UNIQUE ID OF THE MPI JOB
	int parent_id;  // THE PARENT ID OF THE MPI JOB(USED TO RECONSTRUCT THE CUTFILE_NAME)
	double* stack;
};

struct book   //Keeps the information about whether a processor is "Busy(1)" or "Free(0)"
{
	int busy;
	int id; 
};

struct node   //An entry for the QUEUE
{
	struct message_queue queue_data;
	struct node* next;
};

struct process //structure for longfile
{
	char pstr[300];//text for number of particles
	char estr[300];//text for energy of particles
	double pbuffer[NSTEP_MAX][10];//buffer for containing number of particles
	double ebuffer[NSTEP_MAX][10];//buffer for containing energy of particles
};

struct process lngfile,pro; //object for containing the structure of longfile for each parallel task
/*
C  CORSIKA as a function 
C
C  Linked as EXTERNAL library generated from CORSIKA Fortran PROGRAM with Parallel PARAMETERS.
C  ARGUMENTS:
C     LPRIM    (Integer4LOG of Fortran) : PRIMARY INTERACTION (T) OR NOT (F)
C     DECTCUT  (double4DBL) : Threshold energy for sub-showers in gev. All particles with energy above dectcut will have a seed from the 6th sequence of random numbers. All particles above dectcut will be used as a group sub-shower
C     DECTMAX  (double4DBL) : Maximum energy in gev for a complete parallel separate sub-showers. All particles above dectmax will be considered as a primary particles for a new sub-shower.
C     I1CUTPAR (Integer4INT) : First id of particle to read from cutfile or stack buffer
C     I2CUTPAR (Integer4INT) : Last id of particle to read from cutfile or stack buffer
C     CFILINP  (Character4CHA) : Cutfile name to read particles from and fill 2nd stack
C     CFILOUT  (Character4CHA) : Standard output file name
C     CFILSTE  (Character4CHA) : Steering-input file name to read parameters
C     IDMPI    (Integer4INT) : ID OF THE MPI TASK
*/
void corsika_(int* ,double* ,double* ,int* ,int* ,char* ,char* ,char* ,int*);
#if __COAST__
void mpinodedone_ ();
#endif

/*function called by master to write down the resultant Particle Files per parallel task*/
void write_block_fort(int, float *);

/*function called by master to write down the resultant longfile after adding the content of all longfiles received from slaves.
This function is called inside the "finalize_simulation()" function*/
void write_joined_longfile(char *joinedlongfile)
{   
if (dbg_infos)
	{	int rank;
		MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		char statfile[255],numstr[9];
		strcpy(statfile,statdir);
		strcat(statfile,"mpiid-");
		sprintf(numstr, "%d", rank);
		strcat(statfile,numstr);
		strcat(statfile,"-runprotocol.txt");
		FILE* file =fopen(statfile,"a");
		fprintf(file,"\n[%f] LONGFILE:The master processor %d is writing down the joined longfile with %d steps",MPI_Wtime(),MASTER,NSTEP1);
		fclose(file);
	}
	int i,j;
	FILE*longfile;
//	char filename[255];
//	strcpy(filename,joinedlongfile);
	longfile=fopen(joinedlongfile,"a");
	
	fprintf(longfile,"%s",lngfile.pstr);
	for(i=0;i<NSTEP1;i++)
	{	fprintf(longfile,"%6.1f",lngfile.pbuffer[i][0]);
		for(j=1;j<10;j++)
		{  if ( lngfile.pbuffer[i][j] < 0. ) lngfile.pbuffer[i][j]=1.e-12;//normalize values to zero
			fprintf(longfile,"%12.5e",lngfile.pbuffer[i][j]);
//fprintf(longfile,"%6.1f%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\n",lngfile.pbuffer[i][0],lngfile.pbuffer[i+1][1],lngfile.pbuffer[i+2][2],lngfile.pbuffer[i+3][3],lngfile.pbuffer[i+4][4],lngfile.pbuffer[i+5][5],lngfile.pbuffer[i+6][6],lngfile.pbuffer[i+7][7],lngfile.pbuffer[i+8][8],lngfile.pbuffer[i+9][9]);
		}
        fprintf(longfile,"\n");
	}
	
    fprintf(longfile,"%s",lngfile.estr);
	for(i=0;i<NSTEP1;i++)
	{	fprintf(longfile,"%6.1f",lngfile.ebuffer[i][0]);
		for(j=1;j<10;j++)
		{  if ( lngfile.ebuffer[i][j] < 0. ) lngfile.ebuffer[i][j]=1.e-12;//normalize values to zero
			fprintf(longfile,"%12.5e",lngfile.ebuffer[i][j]);
//fprintf(longfile,"%6.1f%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\n",lngfile.ebuffer[i][0],lngfile.ebuffer[i+1][1],lngfile.ebuffer[i+2][2],lngfile.ebuffer[i+3][3],lngfile.ebuffer[i+4][4],lngfile.ebuffer[i+5][5],lngfile.ebuffer[i+6][6],lngfile.ebuffer[i+7][7],lngfile.ebuffer[i+8][8],lngfile.ebuffer[i+9][9]);
		}
        fprintf(longfile,"\n");
	}

    fclose(longfile);
}
/*
function finalize_simulations()
Called from main() by the MASTER to finalize the slave environment.
It broadcast messages to each slave giving order to FINALIZE environment.
It also receives longfiles generated by each slave and add their content to generate the final longfile.
*/
void finalize_simulations(int num_process,double start_time,double stop_time) 
{
 MPI_Status status;
 int i,j;
 int dummy=0,unused=0;
 int mpierr;
 char statfile[255],statfiles[255],numstr[9];
 
 strcpy(statfile,statdir);
 strcat(statfile,"time.txt");
 FILE* file = fopen(statfile,"w");
 fprintf(file,"    START TIME\t\t STOP TIME\t  TIME (min)\n%f   %f   %f\n",start_time,stop_time,(stop_time-start_time)/60.);
 fclose(file);
 
 //Stopping all waiting slaves
 for(i=1;i<num_process;i++) 
	{	 
	 if (dbg_infos) 
		{
		strcpy(statfile,statdir);
		strcat(statfile,"mpiid-");
		sprintf(numstr, "%d", MASTER);
		strcat(statfile,numstr);
		strcat(statfile,"-runprotocol.txt");
		file =fopen(statfile,"a");
		fprintf(file,"\n[%f] STOP: Sending STOP request to slave= %d",MPI_Wtime(),i);
		fclose(file);
		}
	MPI_Send(&dummy,1,MPI_SHORT,i,STOP,MPI_COMM_WORLD);
	}

  for(i=1;i<num_process;i++) 
	{	 
	MPI_Recv(&dummy,1,MPI_SHORT,i,VALN,MPI_COMM_WORLD,&status);
	if (NSTEP1<dummy) NSTEP1=dummy;
	if (dummy==0) unused++;
	if (dbg_infos) 
		{
		file =fopen(statfile,"a");
		fprintf(file,"\n[%f]: RECEIVE #%d from slave %d maximal vertical steps in long files as %d and will use %d. Unused tasks are %d",MPI_Wtime(), i, status.MPI_SOURCE, dummy, NSTEP1,unused);
		fclose(file);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

// int dummybuffer; //collecting the information about maximal vertical steps in long files
// MPI_Reduce (&dummybuffer, &NSTEP1, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
 
 //MASTER receiving the longfile from each slave and adding their content to generate a single longfile.
 for(i=1;i<num_process;i++) 
	{
	if (dbg_infos) 
		{
		file =fopen(statfile,"a");
		fprintf(file,"\n[%f]: RECEIVE #%d MASTER waiting to receiving the string, related to number of particles,from each slave. Unused tasks are %d. Latest communication were with slave %d",MPI_Wtime(),i,unused,status.MPI_SOURCE);
		fclose(file);
		}
	mpierr=MPI_Recv(lngfile.pstr,300,MPI_CHAR,i,STRP,MPI_COMM_WORLD,&status);
	//strcpy(lngfile.pstr,pro.pstr);
	if (dbg_infos)
		{
			strcpy(statfiles,statdir);
			strcat(statfiles,"mpiid-");
			sprintf(numstr, "%d", status.MPI_SOURCE);
			strcat(statfiles,numstr);
			strcat(statfiles,"-recv.long");
			file = fopen(statfiles,"a");
			if(mpierr==MPI_SUCCESS)
			{ 
				fprintf(file,"[%f]:pstr from slave %d is successfully received:\n",MPI_Wtime(),status.MPI_SOURCE);
				fprintf(file,"%s",pro.pstr);
			}
			else
				fprintf(file,"[%f]:pstr from slave %d is not successfully received\n",MPI_Wtime(),status.MPI_SOURCE);
			fclose(file);
		}
	}
 if (dbg_infos)
	{
		file = fopen(statfile,"a"); 
		fprintf(file,"\n[%f] Barrier: waiting to other slaves. Error code of MPI_Recv=%d",MPI_Wtime(),mpierr);
		fclose(file);
	}
 MPI_Barrier(MPI_COMM_WORLD);
 
 for(i=1;i<num_process;i++) 
	{	int ii;
		if (dbg_infos) 
			{
			file =fopen(statfile,"a");
			fprintf(file,"\n[%f]: RECEIVE #%d MASTER waiting to receiving the buffer containing the number of particles from each slave. Unused tasks are %d. Latest communication were with slave %d",MPI_Wtime(), i,unused,status.MPI_SOURCE);
			fclose(file);
			}  
		mpierr=MPI_Recv(&(pro.pbuffer[0][0]),10*NSTEP1,MPI_DOUBLE,i,VALP,MPI_COMM_WORLD,&status);
		if (dbg_infos)
			{
				strcpy(statfiles,statdir);
				strcat(statfiles,"mpiid-");
				sprintf(numstr, "%d", status.MPI_SOURCE);
				strcat(statfiles,numstr);
				strcat(statfiles,"-recv.long");
				file = fopen(statfiles,"a");
				if(mpierr==MPI_SUCCESS)
				{
					fprintf(file,"\n[%f]:pro.pbuffer[%d][10] from slave %d is successfully received\n",MPI_Wtime(),NSTEP1,status.MPI_SOURCE);
					for(ii=0;ii<NSTEP1;ii++)
					{
						fprintf(file,"\n%6.1f ",pro.pbuffer[ii][0]);	 
						for(j=1;j<10;j++)
							fprintf(file,"%12.5e",pro.pbuffer[ii][j]);	 
					}
				}
				else
					fprintf(file,"\n[%f]:pro[].pbuffer[%d][10] from slave %d is not successfully received",MPI_Wtime(),NSTEP1,status.MPI_SOURCE);
				fclose(file);
			}
		//MASTER adding the content of particle buffer received from each slave
		for(ii=0;ii<NSTEP1;ii++)
			{
				lngfile.pbuffer[ii][0]=pro.pbuffer[ii][0];
				for(j=1;j<10;j++)
					lngfile.pbuffer[ii][j]=lngfile.pbuffer[ii][j]+pro.pbuffer[ii][j];
			}
	}
	if (dbg_infos)
	 {
		file = fopen(statfile,"a"); 
		fprintf(file,"\n[%f] Barrier: waiting to other slaves. Error code of MPI_Recv=%d",MPI_Wtime(),mpierr);
		fclose(file);
	 }
	 MPI_Barrier(MPI_COMM_WORLD);
	 
 
for(i=1;i<num_process;i++) 
	{
		if (dbg_infos) 
			{
			file =fopen(statfile,"a");
			fprintf(file,"\n[%f]: RECEIVE #%d MASTER waiting to receiving the string related to energy of particles from each slave. Unused tasks are %d. Latest communication were with slave %d",MPI_Wtime(), i,unused,status.MPI_SOURCE);
			fclose(file);
			} 
		mpierr=MPI_Recv(pro.estr,300,MPI_CHAR,i,STRE,MPI_COMM_WORLD,&status); 
		strcpy(lngfile.estr,pro.estr);
		if (dbg_infos)
		{
			strcpy(statfiles,statdir);
			strcat(statfiles,"mpiid-");
			sprintf(numstr, "%d", status.MPI_SOURCE);
			strcat(statfiles,numstr);
			strcat(statfiles,"-recv.long");
			file = fopen(statfiles,"a");
			if(mpierr==MPI_SUCCESS)
			{ 
				fprintf(file,"\n[%f]:estr from slave %d is successfully received",MPI_Wtime(),status.MPI_SOURCE);
				fprintf(file,"\n%s",pro.estr);
			}
			else
				fprintf(file,"\n[%f]:estr from slave %d is not successfully received",MPI_Wtime(),status.MPI_SOURCE);
			fclose(file);
		}
	}
 if (dbg_infos)
 {
	file = fopen(statfile,"a"); 
	fprintf(file,"\n[%f] Barrier: waiting to other slaves. Error code of MPI_Recv=%d",MPI_Wtime(),mpierr);
	fclose(file);
 }
 MPI_Barrier(MPI_COMM_WORLD);
 
 
for(i=1;i<num_process;i++) 
 {int ii;
	if (dbg_infos) 
		{
		file =fopen(statfile,"a");
		fprintf(file,"\n[%f]: RECEIVE #%d MASTER waiting to receiving the buffer containing the energy of particles from each slave. Unused tasks are %d. Latest communication were with slave %d",MPI_Wtime(), i,unused,status.MPI_SOURCE);
		fclose(file);
		} 
	mpierr=MPI_Recv(&(pro.ebuffer[0][0]),10*NSTEP1,MPI_DOUBLE,i,VALE,MPI_COMM_WORLD,&status);
	if (dbg_infos)
	{
		strcpy(statfiles,statdir);
		strcat(statfiles,"mpiid-");
		sprintf(numstr, "%d", status.MPI_SOURCE);
		strcat(statfiles,numstr);
		strcat(statfiles,"-recv.long");
		file = fopen(statfiles,"a");
		if(mpierr==MPI_SUCCESS)
		{
			fprintf(file,"\n[%f]:pro[].ebuffer[%d][10] from slave %d successfully received",MPI_Wtime(),NSTEP1,status.MPI_SOURCE);
			for(ii=0;ii<NSTEP1;ii++)
			{
				fprintf(file,"\n%6.1f ",pro.ebuffer[ii][0]);	 
				for(j=1;j<10;j++)
					fprintf(file,"%12.5e",pro.ebuffer[ii][j]);	 
			}
//		 for(i=0;i<NSTEP1;i++)
//		 fprintf(file,"%6.1f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",pro.ebuffer[i][0],pro.ebuffer[i+1][1],pro.ebuffer[i+2][2],pro.ebuffer[i+3][3],pro.ebuffer[i+4][4],pro.ebuffer[i+5][5],pro.ebuffer[i+6][6],pro.ebuffer[i+7][7],pro.ebuffer[i+8][8],pro.ebuffer[i+9][9]);
		}
		else
			fprintf(file,"\n[%f]:pro[].ebuffer[%d][10] from slave %d not successfully received",MPI_Wtime(),NSTEP1,status.MPI_SOURCE);
		fclose(file);
	}
	//MASTER adding the content of energy buffer received from each slave  
	for(ii=0;ii<NSTEP1;ii++)
		{
		   lngfile.ebuffer[ii][0]=pro.ebuffer[ii][0];
		   for(j=1;j<10;j++)
		   lngfile.ebuffer[ii][j]=lngfile.ebuffer[ii][j]+pro.ebuffer[ii][j];
		}
 }
 if (dbg_infos)
 {
	file = fopen(statfile,"a"); 
	fprintf(file,"\n[%f] Barrier: waiting to other slaves. Error code of MPI_Recv= %d",MPI_Wtime(),mpierr);
	fclose(file);
 }
 MPI_Barrier(MPI_COMM_WORLD);
}
void finalize_simulation_slave()
{			int y1,y2,y3,y4,rank;
			int dummybuffer; //Dummy buffer for receive variable in collective communication
			char statfile[255],numstr[9];
			FILE *file;
			MPI_Comm_rank (MPI_COMM_WORLD, &rank);
			if (dbg_infos) 
			{
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d", rank);
			strcat(statfile,numstr);
			strcat(statfile,"-runprotocol.txt");
			file = fopen(statfile,"a"); 
			fprintf(file,"\n[%f] FINALIZE: Going to finalize sending %d for number of particles in long file",MPI_Wtime(),NSTEP1);
			fclose(file);
			}
//			MPI_Reduce (&NSTEP1,&dummybuffer, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
			y1=MPI_Send(&NSTEP1,1,MPI_SHORT,MASTER,VALN,MPI_COMM_WORLD);
			if (dbg_infos)
			{
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] BARRIER: to wait for other processes completing submission of number of particles: MPI communication error code %d",MPI_Wtime(),y1);
				fclose(file);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			
			if (NSTEP1)
				{
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] SEND: Going to send to MASTER the string related to particles generated inside longfile \n %s",MPI_Wtime(),lngfile.pstr);
					fclose(file);
					}
				y1=MPI_Send(lngfile.pstr,300,MPI_CHAR,MASTER,STRP,MPI_COMM_WORLD);
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					if(y1==MPI_SUCCESS)  
						fprintf(file,"\n[%f]:pstr by slave %d is successfully sent",MPI_Wtime(),rank);
					else
						fprintf(file,"\n[%f]:pstr by slave %d is not successfully sent",MPI_Wtime(),rank);
					fclose(file);
					}
				if (dbg_infos)
				{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] BARRIER: to wait for other processes completing submission of string about particles: MPI communication error code %d",MPI_Wtime(),y1);
					fclose(file);
				}
				MPI_Barrier(MPI_COMM_WORLD);
				
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] SEND: Going to send to MASTER the buffer containing counts of particles generated inside longfile",MPI_Wtime());
					fclose(file);
					}
				y2=MPI_Send(&(lngfile.pbuffer[0][0]),10*NSTEP1,MPI_DOUBLE,MASTER,VALP,MPI_COMM_WORLD);
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					if(y2==MPI_SUCCESS)
						{int i;
							fprintf(file,"\n[%f]:pbuffer by slave %d is successfully sent",MPI_Wtime(),rank);
							for(i=0;i<NSTEP1;i++)
								{
									fprintf(file,"\n%6.1f ",lngfile.pbuffer[i][0]);	 
									for(dummybuffer=1;dummybuffer<10;dummybuffer++)
										fprintf(file,"%12.5e",lngfile.pbuffer[i][dummybuffer]);	 
								}
						}
					else
						fprintf(file,"\n[%f]:pbuffer by slave %d is not successfully sent",MPI_Wtime(),rank);
					fclose(file);
					}
				if (dbg_infos)
				{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] BARRIER: to wait for other processes completing submission of count of particles: MPI communication error code %d",MPI_Wtime(),y1);
					fclose(file);
				}
				MPI_Barrier(MPI_COMM_WORLD);
				
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] SEND: Going to send to MASTER the string related to energy of particles generated inside longfile",MPI_Wtime());
					fclose(file);
					}
				y3=MPI_Send(lngfile.estr,300,MPI_CHAR,MASTER,STRE,MPI_COMM_WORLD);
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					if(y3==MPI_SUCCESS)  
					{
						fprintf(file,"\n[%f]:estr by slave %d is successfully sent",MPI_Wtime(),rank);
						fprintf(file,"\n%s",lngfile.estr);
					}
					else
						fprintf(file,"\n[%f]:estr by slave %d is not successfully sent",MPI_Wtime(),rank);
					fclose(file);
					}
				if (dbg_infos)
				{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] BARRIER: to wait for other processes completing submission of string about energy: MPI communication error code %d",MPI_Wtime(),y1);
					fclose(file);
				}
				MPI_Barrier(MPI_COMM_WORLD);

				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] SEND: Going to send to MASTER the buffer containing energy of particles generated inside longfile",MPI_Wtime());	
					fclose(file);
					}
				
				y4=MPI_Send(&(lngfile.ebuffer[0][0]),10*NSTEP1,MPI_DOUBLE,MASTER,VALE,MPI_COMM_WORLD);
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					if(y4==MPI_SUCCESS)
						{int i;
							fprintf(file,"\n[%f]:ebuffer by slave %d is successfully sent",MPI_Wtime(),rank);
							for(i=0;i<NSTEP1;i++)
								{
									fprintf(file,"\n%6.1f ",lngfile.ebuffer[i][0]);	 
									for(dummybuffer=1;dummybuffer<10;dummybuffer++)
										fprintf(file,"%12.5e",lngfile.ebuffer[i][dummybuffer]);	 
								}
						}
						else
							fprintf(file,"\n[%f]:ebuffer by slave %d is not successfully sent",MPI_Wtime(),rank);
					fclose(file);
					}
				if (dbg_infos)
				{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] BARRIER: to wait for other processes completing submission of buffer of energies: MPI communication error code %d",MPI_Wtime(),y1);
					fclose(file);
				}
				MPI_Barrier(MPI_COMM_WORLD);
				
				}
			else
				{
				if (dbg_infos)
					{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] UNUSED: Slave %d were not participating in simulations or longfile is empty %d",MPI_Wtime(),rank,NSTEP1);
					fclose(file);
					}
				}
}
/*
function check_for_finish()
Called from main() by the MASTER to check if the simulation is over.
Returns TRUE if no processor is busy(SIMULATION OVER), else returns FALSE
*/
int check_for_finish(struct book* book,int num_process)
{
 int i,check=1;
 for(i=1;i<num_process;i++)
 {
 if(book[i].busy==1)
  {
  check=0;
  break;
  }
 }
 return check;
}

/*
function get_free_rank()
Called from main() by the MASTER to get a free processor rank.
Returns the RANK if found else returns the error code NOT_FOUND
*/
int get_free_rank(struct book* book,int num_process)
{
 int i,check=0;
 for(i=1;i<num_process;i++)
 {
  if(book[i].busy==0)
   {
   check=1;
   break;
   }
 }
 if(check==1)
 return i;
 else return NOT_FOUND;
}

/*
function read_input()
Called from main.
This function reads the MPI input file for the parameters.
*/
void read_input(double* dectcut,double* dectmax,int* lprim,char* cutfilename,int* id1,int* id2,int* rnnum,int* seed,char * input_fname)
{
 int x,dum;
 char* end;
 char str1[10],str2[180],str3[32],str4[32],str5[1],str[255];
 FILE *file = fopen(input_fname,"r");
 if(file==NULL)
  {
  printf("Error opening the file for mpi input %s\n",input_fname);
  MPI_Abort(MPI_COMM_WORLD,19);
  }
 else
  {
    lcout=0;
    *id1=0;
    *id2=0;
    *seed=6;
    *rnnum=1;
    *dectcut=1E4;
    *dectmax=1E7;
    *lprim=1;
    int seeder=1;
 while(!feof(file)) 
 {
  strcpy(str1, " ");
  strcpy(str2, " ");
  strcpy(str3, " ");
  strcpy(str4, " ");
  strcpy(str5, " ");
  strcpy(str, " ");
  fgets (str , 255 , file);
  x= sscanf(str,"%s %s %s %s %s\n",str1,str2,str3,str4,str5);
  //puts(str);
  if(x>0)
	{
	if (dbg_infos) 
		{
		printf("\ntag %s value %s %s %s %s",str1,str2,str3,str4,str5);
		fflush(stdout);
		}
    if(!strcmp(str1,"RUNNR")||!strcmp(str1,"runnr"))
       {
       *rnnum= atoi(str2);
		}
    else if(!strcmp(str1,"SEED")||!strcmp(str1,"seed"))
			{       
			//if(seeder) 
				//{if (dbg_infos) printf("\n First Seed would be read and used for first output file names. The rest Seeds would be ignored");
				*seed= atoi(str2);
				//seeder= 0;
				//}
			}
		 else if(!strcmp(str1,"PARALLEL")||!strcmp(str1,"parallel"))
				{
				*dectcut= strtod(str2,&end);
				*dectmax= strtod(str3,&end);
				dum= atoi(str4);
				if(!strcmp(str5,"T")||!strcmp(str5,"t")) lcout= 1;
				}
			  else if(!strcmp(str1,"DIRECT")||!strcmp(str1,"direct"))
						{
						strcpy(statdir,str2);
						#include <sys/stat.h>
						struct stat sb;
						if (stat(statdir, &sb) == 0 && S_ISDIR(sb.st_mode))
						{
							if (dbg_infos) 
							{
							printf("\ntag %s as %s will be used for outputs",str1,statdir);
							fflush(stdout);
							}
						}
						else
						{
							if (dbg_infos) 
							{
							printf("\nDirectory %s dose not exist. Will be created and used for outputs",statdir);
							fflush(stdout);
							}
							int result = mkdir(statdir, 0700);
							if (result) 
							{
							printf("\nFailed to create directory %s using mkdir. Will abort with error %d",statdir, result);
							MPI_Abort(MPI_COMM_WORLD,result);
							}
						}
						}
					else if(!strcmp(str1,"CUTFILE")||!strcmp(str1,"cutfile"))
								{
								*lprim= 0;
								strcpy(cutfilename,str2);
								*id1= atoi(str3);
								*id2= atoi(str4);
								}
	}
 }
 if (dbg_infos) 
	{
	printf("Successfully red the input file %s\n",input_fname);
	fflush(stdout);
	}
 fclose(file);
}
}

/*
function queue_is_empty()
Called from main() by the MASTER
returns TRUE if the QUEUE is empty, else returns FALSE.
*/
int queue_is_empty(struct node** head,struct node** tail)
{
if(*head==NULL)
 return 1;
else return 0;
}

/*
function add_new_particle_to_queue()
Called from main() by the MASTER
adds a new item into the QUEUE.
*/
void add_new_particle_to_queue(struct node** head, struct node** tail, struct message item)
{
int i,j;
FILE* file;
char statfile[255],numstr[9];
if (dbg_infos) 
{
 strcpy(statfile,statdir);
 strcat(statfile,"mpiid-");
 sprintf(numstr, "%d", MASTER);
 strcat(statfile,numstr);
 strcat(statfile,"-runprotocol.txt");
 file =fopen(statfile,"a");
 fprintf(file,"\n[%f] QUEUE: Adding request %d %d %d DAT%06d-%09d-%09d.cut %d %d %p",MPI_Wtime(),item.uid,item.parent_id,item.lprim,item.rnnum,item.seed,item.uid,item.id1,item.id2,item.stack);
 fclose(file);
}
struct node* temp= (struct node*)malloc(sizeof(struct node));
temp->queue_data.stack=(double*)malloc(sizeof(double)*(item.id2-item.id1+1)*(PARTICLE_INFO_COUNT+1));
temp->queue_data.id1 = item.id1;
temp->queue_data.id2 = item.id2;
temp->queue_data.lprim = item.lprim;
temp->queue_data.rnnum = item.rnnum;
temp->queue_data.seed = item.seed;
temp->queue_data.uid = item.uid;
temp->queue_data.parent_id = item.parent_id;

for(i=0;i<(item.id2-item.id1+1);i++)
 for(j=0;j<(PARTICLE_INFO_COUNT+1);j++)
  temp->queue_data.stack[(PARTICLE_INFO_COUNT+1)*i+j]=item.stack[(PARTICLE_INFO_COUNT+1)*i+j];

temp->next = NULL;

if(*head==NULL)
 {
 *head = temp;
 *tail = temp;
 }
else
 {
 (*tail)->next = temp;
 *tail = temp;
 }
 
if (dbg_infos) 
{
 file =fopen(statfile,"a");
 fprintf(file,"\n[%f] QUEUE: Added %d %d %d DAT%06d-%09d-%09d.cut %d %d %p",MPI_Wtime(),temp->queue_data.uid,temp->queue_data.parent_id,temp->queue_data.lprim,temp->queue_data.rnnum,temp->queue_data.seed,temp->queue_data.uid,temp->queue_data.id1,temp->queue_data.id2,temp->queue_data.stack);
 fclose(file);
}
return;
}

/*
function queue_del()
Called from main() by the MASTER
Identifies if there is a job waiting
Deletes an element from the QUEUE and is returns it to MASTER for further run.
*/
struct message queue_del(struct node** head, struct node** tail)
{FILE* file;
 char statfile[255],numstr[9];
 if (dbg_infos) 
 {
	strcpy(statfile,statdir);
	strcat(statfile,"mpiid-");
	sprintf(numstr, "%d", MASTER);
	strcat(statfile,numstr);
	strcat(statfile,"-runprotocol.txt");
	file =fopen(statfile,"a");
	fprintf(file,"\n[%f] QUEUE: Deleting Request",MPI_Wtime());
	fclose(file);
 }

 struct message a;
 if(queue_is_empty(head,tail))
  {
	if (dbg_infos) 
	{
		file =fopen(statfile,"a");
		fprintf(file,"\n[%f] QUEUE: QUEUE_UNDERFLOW_ERROR, QUEUE IS EMPTY",MPI_Wtime());
		fclose(file);
	}
	else 
	{
	printf("\n[%f] QUEUE: QUEUE_UNDERFLOW_ERROR, QUEUE IS EMPTY",MPI_Wtime());
	fflush(stdout);
	}
  }
 else
  {
	  int i,j;
	  a.id1 = (*head)->queue_data.id1;
	  a.id2 = (*head)->queue_data.id2;
	  a.lprim = (*head)->queue_data.lprim;
	  a.rnnum = (*head)->queue_data.rnnum;
	  a.seed = (*head)->queue_data.seed;
	  a.uid = (*head)->queue_data.uid;
	  a.parent_id = (*head)->queue_data.parent_id;

	  for(i=0;i<(a.id2-a.id1+1);i++)
	   for(j=0;j<(PARTICLE_INFO_COUNT+1);j++)
		a.stack[(PARTICLE_INFO_COUNT+1)*i+j]=(*head)->queue_data.stack[(PARTICLE_INFO_COUNT+1)*i+j];

	  struct node* temp = *head;
	  if((*head)->next==NULL)
	   {
	   *head=*tail=NULL;
	   }
	  else
	   {
	   *head = (*head)->next;
	   }
	  free(temp->queue_data.stack);//free the memory of stack of deleted particle
	  free(temp);
	  if (dbg_infos) 
	  {
		file =fopen(statfile,"a");
		fprintf(file,"\n[%f] QUEUE: DELETED %d %d %d DAT%06d-%09d-%09d.cut %d %d %p",MPI_Wtime(),a.uid,a.parent_id,a.lprim,a.rnnum,a.seed,a.uid,a.id1,a.id2,a.stack);
		fclose(file);
	  }
  }
 return(a);
}

/*
function clear_buffer()
Called from main() by the MASTER/SLAVE
Called from new_particle()
Clears the contents of the buffer to avoid ambiguity.
*/
void clear_buffer(struct message* buff)
{
buff->id1=0;
buff->id2=0;
buff->lprim=-1;
buff->rnnum=-1;
buff->seed=-1;
buff->uid=-1;
buff->parent_id=-1;
return;
}

/*
function entry_update_start()
Called from main() by the MASTER
Writes the status about the jobs started into a file called "status_start"
*/
void entry_update_start(struct message* buff,int rank,int* no_running,int* no_queued,int* no_finished,int* no_lost)
{
 char statfile[255];
 strcpy(statfile,statdir);
 strcat(statfile,"status_start");
 FILE* file = fopen(statfile,"a");
 fprintf(file,"%d %d %d %d %d %d %d %d %d %f %d %d %d %d\n",buff->uid,buff->parent_id,buff->lprim,buff->rnnum,buff->seed,buff->uid,buff->id1,buff->id2,rank,MPI_Wtime(),*no_running,*no_queued,*no_finished,*no_lost);
 fclose(file);
 
 strcpy(statfile,statdir);
 strcat(statfile,"relation");
 file = fopen(statfile,"a");
 fprintf(file,"%d %d\n",buff->uid,buff->parent_id);
 fclose(file);
}

/*
function entry_update_queue()
Called from main() by the MASTER
Writes the status about the jobs queued into a file called "queue_add"
*/
void entry_update_queue(struct message* buff,int* no_running,int* no_queued,int* no_finished,int* no_lost)
{
 char statfile[255];
 strcpy(statfile,statdir);
 strcat(statfile,"queue_add");
 FILE* file = fopen(statfile,"a");
 fprintf(file,"%d %d %d %d %d %d %d %d %f %d %d %d %d\n",buff->uid,buff->parent_id,buff->lprim,buff->rnnum,buff->seed,buff->uid,buff->id1,buff->id2,MPI_Wtime(),*no_running,*no_queued,*no_finished,*no_lost);
 fclose(file);
}


/*
function entry_update_lost()
Called from main() by the MASTER
Writes the status about the jobs lost into a file called "lost_jobs"
*/
void entry_update_lost(struct message* buff,int* no_running,int* no_queued,int* no_finished,int* no_lost)
{
   (*no_lost)++;

   char statfile[255];
   strcpy(statfile,statdir);
   strcat(statfile,"lost_jobs");
   FILE* file = fopen(statfile,"a");
   fprintf(file,"%d %d %d %d %d %d %d %d %d %d %d %d\n",buff->uid,buff->parent_id,buff->lprim,buff->rnnum,buff->seed,buff->uid,buff->id1,buff->id2,*no_running,*no_queued,*no_finished,*no_lost);
   fclose(file);
}

/*
function entry_update_resumed()
Called from main() by the MASTER
Writes the status about the jobs resumed into a file called "status_start"
*/
void entry_update_resume(struct message* buff,int rank,int* no_running,int* no_queued,int* no_finished,int* no_lost)
{
   char statfile[255];
 
   strcpy(statfile,statdir);
   strcat(statfile,"status_start");
   FILE* file = fopen(statfile,"a");
   fprintf(file,"%d %d %d %d %d %d %d %d %d %f %d %d %d %d\n",buff->uid,buff->parent_id,buff->lprim,buff->rnnum,buff->seed,buff->uid,buff->id1,buff->id2,rank,MPI_Wtime(),*no_running,*no_queued,*no_finished,*no_lost);
   fclose(file);
 
   strcpy(statfile,statdir);
   strcat(statfile,"relation");
   file =fopen(statfile,"a");
   fprintf(file,"%d %d\n",buff->uid,buff->parent_id);
   fclose(file);
}  

/*
function entry_update_finished()
Called from main() by the MASTER
Writes the status about the jobs started into a file called "status_finish"
*/
void entry_update_finish(int uid,int rank,double mpi_time1,double mpi_time2,int* no_running,int* no_queued,int* no_finished,int* no_lost)
{
   (*no_finished)++;
   (*no_running)--;
   char statfile[255];
 
   strcpy(statfile,statdir);
   strcat(statfile,"status_finish");
   FILE* file = fopen(statfile,"a");
   fprintf(file,"%d %d %f %f %f %d %d %d %d\n",uid,rank,MPI_Wtime(),mpi_time1,mpi_time2,*no_running,*no_queued,*no_finished,*no_lost);
   fclose(file);
}       

/*
function reconstruct_cfname()
Called from main() by the MASTER and SLAVE both
Reconstructs the CUTFILE name from the run number,seed and unique ID of parallel simulated particle.
*/
void reconstruct_cfname( char* cutfilename,int rnnum,int seed,int uid)
{int rank;
 MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 sprintf(cutfilename,"DAT%06d-%09d-%09d.cut",rnnum,seed,uid);
 if (dbg_infos) 
   {
	 char statfile[255],numstr[9];
	 strcpy(statfile,statdir);
	 strcat(statfile,"mpiid-");
	 sprintf(numstr, "%d", rank);
	 strcat(statfile,numstr);
	 strcat(statfile,"-runprotocol.txt");
	 FILE* file =fopen(statfile,"a");
	 fprintf(file,"\n[%f] CUTFILENAME: reconstruct of cut file name as %s with rnnum = %d  seed = %d uid = %d",MPI_Wtime(),cutfilename,rnnum,seed,uid);
	 fclose(file);
	}
 return;
}

/*
function reconstruct_ofname()
Called from main() by the SLAVE
Reconstructs the OUTPUTFILE name from the run number,seed, MPI ID and directory.
*/
void reconstruct_ofname( char* outfname, int rnnum, int seed, int uid)
{	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	if (rank>LISTFILE_MAX) 
		sprintf(outfname,"/dev/null");
	else
		sprintf(outfname,"%sDAT%06d-%09d-%09d.lst",statdir,rnnum,seed,uid);
	if (dbg_infos) 
	{sprintf(outfname,"%sDAT%06d-%09d-%09d.lst",statdir,rnnum,seed,uid);
	 char statfile[255],numstr[9];
	 strcpy(statfile,statdir);
	 strcat(statfile,"mpiid-");
	 sprintf(numstr, "%d", rank);
	 strcat(statfile,numstr);
	 strcat(statfile,"-runprotocol.txt");
	 FILE* file =fopen(statfile,"a");
	 fprintf(file,"\n[%f] OUTPUTFILENAME: reconstruct of output file name as %s with rnnum = %d  seed = %d uid = %d",MPI_Wtime(),outfname,rnnum,seed,uid);
	 fclose(file);
	}
   return;
}

/*
function parse_cfname()
Called from main() by the MASTER to parse the CUTFILENAME read from the command line/input file into the runnum,seed.mpiid
*/
void parse_cfname( char* cutfilename,int* rnnum,int* seed,int* uid)
{
   sscanf(cutfilename,"DAT%06d-%09d-%09d.cut",rnnum,seed,uid);
   if (dbg_infos) 
	{
	char statfile[255],numstr[9];
	strcpy(statfile,statdir);
	strcat(statfile,"mpiid-");
    sprintf(numstr, "%d", *uid);
    strcat(statfile,numstr);
	strcat(statfile,"-runprotocol.txt");
	FILE* file =fopen(statfile,"a");
	fprintf(file,"\n[%f] CUTFILENAME: parse of cut file name as %s with rnnum = %d  seed = %d uid = %d",MPI_Wtime(),cutfilename,*rnnum,*seed,*uid);
	fclose(file);
	}
   return;
}

/*
function main()
Functions as a SLAVE or MASTER depending on the rank of the processor.
*/
int main(int argc, char *argv[])
{
	int no_running,no_queued,no_finished,no_lost; //Number of jobs running,queued,finished and lost respectively.
	no_running=no_queued=no_finished=no_lost=0;
	double dectcut,dectmax;      //DECTCUT and DECTMAX are the energy thresholds
	double start_time,stop_time;     //start and stop time of the complete simulation
	char input_fname[255];       //input file name.(To be read from input file/command line)
	int lprim=1;         //LPRIM(BOOLEAN)- tells weather CUTFILE should be read or not:primary(if 1) or secondary(if 0) particle simulation
	int id1,id2,rnnum,seed;
	char cutfilename[255];       //CUTFILE name
	char statfile[255],statfile2[255],numstr[9];

	FILE* file,*comm,*check;

	int rank, num_process, mpierr , i,j;              //rank of the processor,number of processors.
	mpierr = MPI_Init (&argc, &argv);                //initialize  mpi environment ....
	if(mpierr!= MPI_SUCCESS)
	{
		printf("ERROR: ERROR INITIATING THE MPI ENVIRONMENT, THE PROGRAM WILL ABORT\n");
		MPI_Abort(MPI_COMM_WORLD,mpierr);
	}

	//initializing variables for message
	MPI_Datatype struct_type,oldtypes[BLOCK_COUNT_MESSAGE]={MPI_INT,MPI_DOUBLE};
	int blockcounts[BLOCK_COUNT_MESSAGE]={7,MAX_GROUP_SIZE*(PARTICLE_INFO_COUNT+1)};
	MPI_Aint offsets[BLOCK_COUNT_MESSAGE],extent;
	//setup description for the fields in the message structure 
	int ofst=0;
	for(i=0;i<BLOCK_COUNT_MESSAGE;i++)
	{
		offsets[i]=ofst;
		MPI_Type_extent(oldtypes[i], &extent);
		ofst = ofst + blockcounts[i]*extent;
	} 
   
	// Now define structured type and commit it.
	MPI_Type_create_struct(BLOCK_COUNT_MESSAGE, blockcounts, offsets, oldtypes, &struct_type);
	MPI_Type_commit(&struct_type); 
	//fetching the values of NUMBER OF PROCESSORS and CURRENT RANK.
	MPI_Comm_size (MPI_COMM_WORLD, &num_process);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	//initialization of pbuffer and ebuffer//////////
	for(i=0;i<NSTEP_MAX;i++)
	{
	  for(j=0;j<10;j++)
	  {
	    lngfile.pbuffer[i][j]=0.0;
		lngfile.ebuffer[i][j]=0.0;
		}
	}
	for(i=0;i<MAXBLOCKS;i++) Block[i]=0.0; //Initialise the Particle stack block-array
	//////////////processing by MASTER////////////////
	if(rank==MASTER)
	{
		//Reading the Command line Arguments                  
		if(argc<2)
		{
			// MPI_Abort(MPI_COMM_WORLD,rc);
			//Give default values to the parameters
			strcpy(input_fname,"parallel-inputs");
			printf("WARNING: input file is not defined. Default parameters and file \"%s\" will be used in non debugging mode.\n",input_fname);
			dbg_infos=0;  
			//Read parameters from input file.
			read_input(&dectcut,&dectmax,&lprim,cutfilename,&id1,&id2,&rnnum,&seed,input_fname);
		}
		else
		{
			int i;
			printf("Input: arguments read from the command line are ");
			for(i=1;i<argc;i++) printf("%s ",argv[i]);
			strcpy(input_fname,argv[1]);
			printf("\nFor inputs the FILE \"%s\" will be used\n",input_fname);

			if(argc>=3) 
				if(!strcmp(argv[2],"T")||!strcmp(argv[2],"t")) 
				{ 
					dbg_infos=1; 
					printf("Debugging mode will be used to protocol steps done by MPI-slaves. Check the \"mpiid-$ID$-runprotocol.txt\" , where $ID$ is the unique identification of parallel running slave rank\n");
				}
				else
				{
					dbg_infos=0;  
					printf("Debugging switch \"%d\". No additional information about parallel run will be produced\n",dbg_infos);
				}

			//Read parameters from input file. 
			read_input(&dectcut,&dectmax,&lprim,cutfilename,&id1,&id2,&rnnum,&seed,input_fname);
		}
         //constructing output file name for resultant joinedlongfile generated by MASTER
		 char joinedlongfile[255],joineddatfile[255];
		 sprintf(joinedlongfile,"%sDAT%06d-999999999.long",statdir,rnnum);
         
		if (dbg_infos) 
		{
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d", MASTER);
			strcat(statfile,numstr);
			strcat(statfile,"-runprotocol.txt");
			check = fopen(statfile,"a");
			fprintf(check,"\n[%f]MASTER: Going to broadcast to slaves of MPI_World of size %d. The input file name \"%s\" output folder \"%s\" emax \"%f\" emin \"%f\" CUTFILE switch \"%d\" and Debug switch \"%d\"",MPI_Wtime(),num_process,input_fname,statdir,dectmax,dectcut,lcout,dbg_infos);
			fflush(check);
		}
		// Broadcasting the name of steering file to all slaves
		MPI_Bcast(input_fname, 255, MPI_CHAR, MASTER, MPI_COMM_WORLD); 
		// Broadcasting the name of output directory to all slaves
		MPI_Bcast(statdir, 255, MPI_CHAR, MASTER, MPI_COMM_WORLD); 
		// Broadcasting the Energy threshold for single particle jobs to all slaves
		MPI_Bcast(&dectmax, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD); 
		// Broadcasting the Energy threshold for group of particle jobs to all slaves
		MPI_Bcast(&dectcut, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD); 
		// Broadcasting the Cutfile generation switch to all slaves
		MPI_Bcast(&lcout, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
		// Broadcasting the Debugger switch to all slaves
		MPI_Bcast(&dbg_infos, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
		// Broadcasting the Debugger switch to all slaves
		MPI_Bcast(&rnnum, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
		if (dbg_infos) 
		{
			fprintf(check,"\n[%f] MASTER: Broadcasted to slaves the input file name \"%s\" output folder \"%s\" emax \"%f\" emin \"%f\" CUTFILE switch \"%d\" Debug switch \"%d\" Run number \"%d\"",MPI_Wtime(),input_fname,statdir,dectmax,dectcut,lcout,dbg_infos,rnnum);
			fflush(check);
		}

		//buffer for MPI communications
		struct message buff;
		int i;

		// unique id for the subshowers 
		int uid=0;

		//Initializing the array(book) to keep track of processors status: free/working
		struct book* book = (struct book*) malloc(num_process*sizeof(struct book));
		struct node *head,*tail; //head and tail of the queue used for storing the waiting processes.

		//initialize the queue       
		head=tail=NULL;

		//initialize all ranks as free
		for(i=0;i<num_process;i++) book[i].busy=0;

		//starting the primary particle and allotting a unique id
		++uid;

		if(lprim==0)
		{//setting buffer for secondary particle
			buff.lprim=lprim;
			buff.id1=id1;
			buff.id2=id2;
			parse_cfname(cutfilename,&buff.rnnum,&buff.seed,&buff.parent_id);
			buff.uid=uid;
			if (dbg_infos)
			{
				fprintf(check,"\n[%f] CUTFILENAME:%s PARSED from RUNNUMBER= %d, SEED= %d, MPIID= %d",MPI_Wtime(),cutfilename,buff.rnnum,buff.seed,buff.parent_id);
				fprintf(check,"\n[%f] RUN_INFO: SECONDARY RUN, LPRIM = %d I1CUTPAR = %d I2CUTPAR = %d RUNNUMBER = %d SEED = %d ID = %d",MPI_Wtime(),buff.lprim,buff.id1,buff.id2,buff.rnnum,buff.seed,buff.uid);
			}
		}
		else
		{//setting buffer for primary particle
			buff.lprim=lprim;
			buff.id1=0;
			buff.id2=0;
			buff.rnnum=rnnum;
			buff.seed=seed;
			buff.uid=uid;
			buff.parent_id=MASTER;
			if (dbg_infos) fprintf(check,"\n[%f] RUN_INFO: PRIMARY RUN, LPRIM = %d LCOUT = %d I1CUTPAR = %d I2CUTPAR = %d RUNNUMBER = %d SEED = %d ID = %d",MPI_Wtime(),buff.lprim,lcout,buff.id1,buff.id2,buff.rnnum,buff.seed,buff.uid); 
		}
		if (dbg_infos) fclose(check);
	  
		//allotting number of working/running slaves
		no_running++;
		//Printing status in statistic files
		if (dbg_infos) entry_update_start(&buff,1,&no_running,&no_queued,&no_finished,&no_lost);

		//Recording the start time of the SIMULATION
		start_time=MPI_Wtime();
		//Send the order to start the parallel CORSIKA task on slave with id=buff.uid
		MPI_Send(&buff,1,struct_type,1,START,MPI_COMM_WORLD);

		if (dbg_infos) 
		{  //Record status of Master to Slave communication in statistic file
			strcpy(statfile,statdir);
			strcat(statfile,"Master2SlaveOrder");
			comm = fopen(statfile,"w");
			fprintf(comm,"%d %d %d %d %d %d %d %d %p\n",buff.uid,buff.parent_id,buff.lprim,buff.rnnum,buff.seed,1,buff.id1,buff.id2,buff.stack);
			fclose(comm);
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d", MASTER);
			strcat(statfile,numstr);
			strcat(statfile,"-runprotocol.txt");
			file = fopen(statfile,"a");
			fprintf(file,"\n[%f] START: The initial shower with id = %d and parent id = %d will run on the the rank = %d",start_time,buff.uid,buff.parent_id,1);
			fclose(file); 
		}

		//updating status for bookkeeping
		book[1].busy=1;
		book[1].id=uid;

master_back:; 
		//check if slaves have finished their jobs and finalize simulation ??may be the queue of waiting jobs need to be checked again??
		if(check_for_finish(book,num_process))
		{
			free(book);
			//recording the stop time of the SIMULATION
			stop_time=MPI_Wtime();                
			//Finish the simulation
			if (dbg_infos) 
			{  //Record status of Master to Slave communication in statistic file
				strcpy(statfile,statdir);
				strcat(statfile,"mpiid-");
				sprintf(numstr, "%d", MASTER);
				strcat(statfile,numstr);
				strcat(statfile,"-runprotocol.txt");
				file = fopen(statfile,"a");
				fprintf(file,"\n[%f] FINALIZE: The shower simulation took %f minutes started at %f is finished at %f",MPI_Wtime(),(stop_time-start_time)/60.,start_time,stop_time);
				fclose(file); 
			}
			finalize_simulations(num_process,start_time,stop_time);
			write_joined_longfile(joinedlongfile);
#if __COAST__
			 unsigned int numAntennas,length,mpierr,dummy,k;
			 mpierr=MPI_Reduce(&dummy,&numAntennas,1,MPI_UNSIGNED,MPI_SUM,MASTER,MPI_COMM_WORLD); 
			 if ( numAntennas%(num_process-1)||mpierr!= MPI_SUCCESS )
				{
				printf("ERROR: In MPI communication or the number of antennas=%d are different in parallel simulations\n",numAntennas);
				MPI_Abort(MPI_COMM_WORLD,mpierr);
				}
			if (dbg_infos) 
				{
				file =fopen(statfile,"a");
				fprintf(file,"\n[%f]: Reduced number of antennas=%d",MPI_Wtime(), numAntennas);
				fclose(file);
				} 
			 mpierr=MPI_Reduce(&dummy,&length,1,MPI_UNSIGNED,MPI_SUM,MASTER,MPI_COMM_WORLD); 
			 if ( length%(num_process-1)||mpierr!= MPI_SUCCESS )
				{
				printf("ERROR: In MPI communication or the length of antenna data=%d are different in parallel simulations\n",length);
				MPI_Abort(MPI_COMM_WORLD,mpierr);
				}
			if (dbg_infos) 
				{
				file =fopen(statfile,"a");
				fprintf(file,"\n[%f]: Reduced length of antennas=%d",MPI_Wtime(), length);
				fclose(file);
				} 
			 numAntennas/=(num_process-1);
			 length/=(num_process-1);
			 
			 MPI_Status status;
/*			 MPI_Datatype struct_type_antenna,oldtypes[5]={MPI_CHAR,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
			 int blockcounts[5]={300,length,length,length,length};
*/			 const long bufferlength = 1000000; // hard-coded buffer depth correlated to same size TPlotter.cc
			 
/*			 MPI_Aint offsets[5],extent;
			 //setup description for the fields in the message structure 
			 int ofst=0;
			 for(i=0;i<5;i++)
				{
					offsets[i]=ofst;
					MPI_Type_get_extent(oldtypes[i], ,&extent);
					ofst = ofst + blockcounts[i]*extent;
				} 
			   			
			double *antbuff;
			antbuff = (double*)calloc(length,sizeof(double));
*/
			typedef struct antenna_s
				{
				char filename[300];//text for file name of antenna-identifier
				double *t;//time buffer 
				double *x;//X coordinate buffer
				double *y;//Y coordinate buffer
				double *z;//Z coordinate buffer
				} antennastruc;
			 antennastruc antennabuffer;
			 antennabuffer.t = (double*)calloc(length,sizeof(double));
			 antennabuffer.x = (double*)calloc(length,sizeof(double));
			 antennabuffer.y = (double*)calloc(length,sizeof(double));
			 antennabuffer.z = (double*)calloc(length,sizeof(double));	

/*			 MPI_Aint offsets[5] = {
			   offsetof(antennastruc, filename),
			   offsetof(antennastruc, t),
			   offsetof(antennastruc, x),
			   offsetof(antennastruc, y),
			   offsetof(antennastruc, z),
									};
			 // Now define structured type and commit it.
			 MPI_Type_create_struct(5, blockcounts, offsets, oldtypes, &struct_type_antenna);
			 MPI_Type_commit(&struct_type_antenna); 
*/			
			 antennastruc * antennas=(antennastruc *)malloc((numAntennas+1)*sizeof(antennastruc));//array of structures
				
			 for(j=1; j<=numAntennas; j++)
				{
				antennas[j].t = (double*)calloc(length,sizeof(double));
				antennas[j].x = (double*)calloc(length,sizeof(double));
				antennas[j].y = (double*)calloc(length,sizeof(double));
				antennas[j].z = (double*)calloc(length,sizeof(double));	
				}
				double antenna[bufferlength*4]; // for t, x, y, z  
			 for(i=1;i<num_process;i++) 
				for(j=1; j<=numAntennas; j++)
				{	if (dbg_infos) 
						{
						file =fopen(statfile,"a");
						fprintf(file,"\n[%f]: WAITING to receive the #%d antenna data block {%d} from parallel process #%d", MPI_Wtime(), j, length, i);
						fclose(file);
						} 
//					mpierr=MPI_Recv(&antennabuffer.filename[0],1,struct_type_antenna,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
					mpierr=MPI_Recv(antennabuffer.filename, 300, MPI_CHAR, i, 0, MPI_COMM_WORLD,&status);
					status.MPI_TAG=j;
					if(mpierr!= MPI_SUCCESS)
						{
						printf("ERROR: ERROR in MPI communication, THE PROGRAM WILL ABORT\n");
						MPI_Abort(MPI_COMM_WORLD,mpierr);
						}
					if (status.MPI_TAG!=j)
						{
						printf("ERROR: receiving data of radio-antenna #%d from corresponding parallel process instead the data from #%d task received. WILL ABORT\n",status.MPI_TAG,j);
						MPI_Abort(MPI_COMM_WORLD,mpierr);
						}
					//TODO: check content of actual structure to be sure that filename fits
					strcpy(antennas[j].filename,antennabuffer.filename);
					if (dbg_infos) 
						{
						file =fopen(statfile,"a");
						fprintf(file,"\n[%f]: RECEIVED antenna #%d, expected #%d, from parallel process #%d - waited from #%d, the buffer containing file name of antenna data. File to be written %s name received %s",MPI_Wtime(),status.MPI_TAG,j,status.MPI_SOURCE,i,antennas[j].filename,antennabuffer.filename);
						fclose(file);
						} 
					mpierr=MPI_Recv(antenna, bufferlength*4, MPI_DOUBLE, i, 1, MPI_COMM_WORLD,&status);
					status.MPI_TAG=j;
/*					for(k=0;k<length;k++)
						{
						if (dbg_infos) 
						{
						file =fopen(statfile,"a");
						fprintf(file,"\n[%f]: ANTENNA #%d of length #%d time =%f",MPI_Wtime(),status.MPI_TAG,k,antbuff[k]);
						fclose(file);
						}
						antennas[status.MPI_TAG].t[k] = antbuff[k];
						}
					mpierr=MPI_Recv(antennabuffer.x, length, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
					mpierr=MPI_Recv(antennabuffer.y, length, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
					mpierr=MPI_Recv(antennabuffer.z, length, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD,&status);
*/					
					for(k=0;k<length;k++)
						{
						if (dbg_infos) 
						{
						file =fopen(statfile,"a");
						fprintf(file,"\n[%f]: RECEIVED antenna #%d, #%d expected, from slave #%d - waited from #%d, at length #%d with data time=%.12e x=%.12e y=%.12e z=%.12e",MPI_Wtime(),status.MPI_TAG,j,status.MPI_SOURCE,i,k,antenna[k],antenna[k+bufferlength],antenna[k+bufferlength*2],antenna[k+bufferlength*3]);
						fclose(file);
						}
						}
					for(k=0;k<length;k++)
						{
						antennas[status.MPI_TAG].t[k] = antenna[k];
						//TODO: cross-check content of actual structure to be sure that t are fitting
						antennas[status.MPI_TAG].x[k] +=  antenna[k+bufferlength];
						antennas[status.MPI_TAG].y[k] +=  antenna[k+bufferlength*2];
						antennas[status.MPI_TAG].z[k] +=  antenna[k+bufferlength*3];
						}
				}
			 if (dbg_infos)
			 {
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] Barrier: Master syncronisation with other slaves to be sure receiving all antenna data before writing out. Error code of MPI_Recv= %d",MPI_Wtime(),mpierr);
				fclose(file);
			 }
			 MPI_Barrier(MPI_COMM_WORLD);       

			  for(j=1; j<=numAntennas; j++)
			  {FILE *antfile;
			   antfile=fopen(antennas[j].filename,"w");
			   for(k=0;k<length;k++)
				fprintf(antfile,"%.12e\t%.12e\t%.12e\t%.12e\n",antennas[j].t[k],antennas[j].x[k],antennas[j].y[k],antennas[j].z[k]);
			   fclose(antfile);
			  }
			 free(antennabuffer.t);
			 free(antennabuffer.x);
			 free(antennabuffer.y);
			 free(antennabuffer.z);
			 for(j=1; j<=numAntennas; j++)
				{
				free(antennas[j].t);
				free(antennas[j].x);
				free(antennas[j].y);
				free(antennas[j].z);
				}
//			 MPI_Type_free(&struct_type_antenna);
#endif
		} 
		else
		{
			MPI_Status stat;
			//Clearing the buffer to avoid any ambiguity
			clear_buffer(&buff);

			if (dbg_infos) 
			{
				strcpy(statfile,statdir);
				strcat(statfile,"mpiid-");
				sprintf(numstr, "%d", MASTER);
				strcat(statfile,numstr);
				strcat(statfile,"-runprotocol.txt");
				file = fopen(statfile,"a");
				fprintf(file,"\n[%f] WAITING: Waiting to receive ANY message from ANY slave",MPI_Wtime());
				fclose(file); 
			}
			MPI_Recv(&buff,1,struct_type,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
		  
			switch(stat.MPI_TAG)
			{
			case REQUEST:
			//Slave sent a 'REQUEST' informing to master about being ready for a new job
			{
				int i,j;
				int parent_rank=stat.MPI_SOURCE;
				//get the parent id (requesting id) from the bookkeeping
				int parent_id=book[parent_rank].id;

				if (dbg_infos) 
				{
					//Record status of Slave to Master communication in statistic files
					strcpy(statfile,statdir);
					strcat(statfile,"Slave2MasterRecv");
					FILE* comm = fopen(statfile,"a");
					fprintf(comm,"%d %d %d %d %d %d %d %p\n",buff.uid,buff.lprim,buff.rnnum,buff.seed,buff.uid,buff.id1,buff.id2,buff.stack);
					fclose(comm); 

					strcpy(statfile,statdir);
					strcat(statfile,"mpiid-");
					sprintf(numstr, "%d", MASTER);
					strcat(statfile,numstr);
					strcat(statfile,"-runprotocol.txt");

					file = fopen(statfile,"a");
					fprintf(file,"\n[%f] REQUEST: received request from slave with rank %d running subshower with id %d asking to run another parallel task for %d particles with following parameters",MPI_Wtime(),parent_rank,parent_id,buff.id2-buff.id1+1);
					for(i=0;i<(buff.id2-buff.id1+1);i++)
					{
					fprintf(file,"\n");
					for(j=0;j<(PARTICLE_INFO_COUNT+1);j++) fprintf(file,"%f ",buff.stack[i*(PARTICLE_INFO_COUNT+1)+j]);
					}
//					fprintf(file,"\n[%f] For each above mentioned particle/group a new parallel CORSIKA run will be started",MPI_Wtime());
					fclose(file); 
				}

/*				//Reconstruct the CUT FILENAME to be used in case of secondary shower simulations ?? why it is needed here??
				reconstruct_cfname(cutfilename, buff.rnnum, buff.seed, parent_id);
*/
				//allot a new unique ID to the particle
				int child_id= ++uid;
				buff.uid = child_id;
				buff.parent_id = parent_id;
				//fetch a free rank
				int child_rank= get_free_rank(book,num_process);

				//if all slaves are busy, insert the task into a queue
				if(child_rank==NOT_FOUND)
				{
					add_new_particle_to_queue(&head,&tail,buff);
					no_queued++;
					if (dbg_infos) entry_update_queue(&buff,&no_running,&no_queued,&no_finished,&no_lost);
goto master_back;
				}
				else  //Functioning if a free rank is found
				{
					no_running++;//allotting the number of working/running slaves

					if (dbg_infos)
					{//MPI World tagging in statistic files
					  strcpy(statfile,statdir);
					  strcat(statfile,"mpiid-");
					  sprintf(numstr, "%d", MASTER);
					  strcat(statfile,numstr);
					  strcat(statfile,"-runprotocol.txt");
					  file = fopen(statfile,"a");
					  fprintf(file,"\n[%f] START: The task with unique id = %d and parent id = %d will run on the the rank = %d",MPI_Wtime(),child_id,parent_id,child_rank);
					  fclose(file); 
					  entry_update_start(&buff,child_rank,&no_running,&no_queued,&no_finished,&no_lost);
					}
				 
					MPI_Send(&buff,1,struct_type,child_rank,START,MPI_COMM_WORLD);     //send the order to start the task

					if (dbg_infos)
					 {
					 //Record status of Master to Slave communication in statistic files
					 strcpy(statfile,statdir);
					 strcat(statfile,"Master2SlaveOrder");
					 comm = fopen(statfile,"a");              
					 fprintf(comm,"%d %d %d %d %d %d %d %d %p\n",buff.uid,buff.parent_id,buff.lprim,buff.rnnum,buff.seed,child_rank,buff.id1,buff.id2,buff.stack);
					 fclose(comm); 
					 }

					//updating status for bookkeeping
					book[child_rank].busy=1;
					book[child_rank].id=child_id;

				}
goto master_back;
			}
			break;
			
			case STORRUNH:
			{ 	if (switch_first_runH) 
					{
					for (i=0;i<buff.rnnum;i++) Run_header[i] = buff.stack[i];
					switch_first_runH=0;
					}
				for (i=0;i<buff.rnnum;i++) buff.stack[i] = Run_header[i];
				i=MPI_Send(&buff,1,struct_type,stat.MPI_SOURCE,STORRUNH,MPI_COMM_WORLD);
goto master_back;
			}
			break;
			
			case STOREVNH:
			{ 	if (switch_first_evnH) 
					{
					for (i=0;i<buff.rnnum;i++) Event_header[i] = buff.stack[i];
					switch_first_evnH=0;
					}
				for (i=0;i<buff.rnnum;i++) buff.stack[i] = Event_header[i];
				i=MPI_Send(&buff,1,struct_type,stat.MPI_SOURCE,STOREVNH,MPI_COMM_WORLD);
goto master_back;
			}
			break;
			
			case FINISH:
			//Steps of MASTER when a 'FINISH' message is received from a SLAVE
			{
				int free_rank=stat.MPI_SOURCE;
				int finished_id=book[free_rank].id;

				//Type cast the buffer as the message received now contains timing information instead of the particle information.
				double* mpi_time = (double*)(&buff);           

				//Record status of parallel runs in statistic files
				no_finished++;
				no_running--;
				if (dbg_infos)
				{
					char statfile[255];
					strcpy(statfile,statdir);
					strcat(statfile,"status_finish");
					file = fopen(statfile,"a");
					fprintf(file,"%d %d %f %f %f %d %d %d %d\n",finished_id,free_rank,MPI_Wtime(),mpi_time[0],mpi_time[1],no_running,no_queued,no_finished,no_lost);
					fclose(file);
				}
				
				//updating status for bookkeeping
				book[free_rank].busy=0;               
				book[free_rank].id=NA;
			   
				//MPI World tagging in statistic files
				if (dbg_infos)
				{
					strcpy(statfile,statdir);
					strcat(statfile,"mpiid-");
					sprintf(numstr, "%d", MASTER);
					strcat(statfile,numstr);
					strcat(statfile,"-runprotocol.txt");
					file = fopen(statfile,"a");
					fprintf(file,"\n[%f] FINISH: The rank = %d running unique id = %d is now set as free in bookkeeping",MPI_Wtime(),free_rank,finished_id);
					fclose(file);
				}
				//Check the queue to see if any jobs is waiting 
				if(!queue_is_empty(&head,&tail))              
				{	int mpierr;
					//read a pending task and remove it from the queue 
					struct message temp= queue_del(&head,&tail);
					no_running++;
					no_queued--;
					if (dbg_infos) 
					{
						file = fopen(statfile,"a");
						fprintf(file,"\n[%f] FINISH: The task with unique id = %d previously waiting in queue will run on rank = %d",MPI_Wtime(),temp.uid,free_rank); 
						fclose(file);
						entry_update_resume(&temp,free_rank,&no_running,&no_queued,&no_finished,&no_lost);
					}

					//Sending the order to start the task on free rank
					mpierr=MPI_Send(&temp,1,struct_type,free_rank,START,MPI_COMM_WORLD);

					if (dbg_infos) 
					{
						file = fopen(statfile,"a");
						fprintf(file,"\n[%f] FINISH: Sent of request to slave %d returned error code %d",MPI_Wtime(),free_rank,mpierr);
						fclose(file);
						//Record status of Master to Slave communication in statistic files
						strcpy(statfile,statdir);
						strcat(statfile,"Master2SlaveOrder");
						comm = fopen(statfile,"a");              
						fprintf(comm,"%d %d %d %d %d %d %d %d %p\n",temp.uid,temp.parent_id,temp.lprim,temp.rnnum,temp.seed,free_rank,temp.id1,temp.id2,temp.stack);
						fclose(comm);
					}
					//updating status for bookkeeping
					book[free_rank].busy=1;
					book[free_rank].id=temp.uid;
				}
goto master_back;
			}
			break; 

			default:
			{
				if (dbg_infos)
				{
				  strcpy(statfile,statdir);
				  strcat(statfile,"mpiid-");
				  sprintf(numstr, "%d", MASTER);
				  strcat(statfile,numstr);
				  strcat(statfile,"-runprotocol.txt");
				  file = fopen(statfile,"a"); 
				  fprintf(file,"\n[%f] ERROR: Unknown message type \n",MPI_Wtime());
				  fclose(file);
				}
goto master_back;
			} 
			break;
			} //End of MPI status switches
		} //End of bookkeeping 
	} //End of IF on Rank is Master 
 
	///////////////processing by slaves//////////////
	else                    
	{
		int i,j,length,k;
		MPI_Status stat;
		struct message buff;	//MPI Buffer for Communication
		char outfname[255];		//Output file name for CORSIKA
		double mpi_time[2];		//To record time for CORSIKA

		// Receiving the name of steering file from Master via Broadcast
		MPI_Bcast(input_fname, 255, MPI_CHAR, MASTER, MPI_COMM_WORLD); 
		// Receiving the name of output directory From Master for all slaves
		MPI_Bcast(statdir, 255, MPI_CHAR, MASTER, MPI_COMM_WORLD);
		// Receiving the Energy threshold for single particle as a slave
		MPI_Bcast(&dectmax, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD); 
		// Receiving the Energy threshold for group of particle as a slave
		MPI_Bcast(&dectcut, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD); 
		// Receiving Cutfile switch as a slave
		MPI_Bcast(&lcout, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
		// Receiving the Debugger switch as a slave
		MPI_Bcast(&dbg_infos, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
		// Receiving the Run number as a slave
		MPI_Bcast(&rnnum, 1, MPI_INT, MASTER, MPI_COMM_WORLD); 
		runnumber=rnnum;
		job_count= 0;			//Counter of initial particle-stacks managed by slave

		if (dbg_infos) 
		{
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d", rank);
			strcat(statfile,numstr);
			strcat(statfile,"-runprotocol.txt");
			file = fopen(statfile,"a"); 
			fprintf(file,"\n[%f] SLAVE: Receiving the input file name \"%s\" output folder \"%s\" emax \"%f\" emin \"%f\" CUTFILE switch \"%d\" Debug switch \"%d\" Run number \"%d\"",MPI_Wtime(),input_fname,statdir,dectmax,dectcut,lcout,dbg_infos,rnnum);
			fclose(file);
		}
		
slave_back:;
		clear_buffer(&buff);	//Clearing the Buffer to avoid ambiguity.
		II1=-1;					//Resetting the pointer for current particle in group
		particle_stack_or_cut_file = (double*)calloc((PARTICLE_INFO_COUNT+1)*MXLIST,sizeof(double)); //Stack for possible new secondary particles
		
		//Waiting for ANY type of message(new particle) from the MASTER.
		if (dbg_infos) 
		{
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d", rank);
			strcat(statfile,numstr);
			strcat(statfile,"-runprotocol.txt");
			file = fopen(statfile,"a"); 
			fprintf(file,"\n[%f] WAITING: Waiting for a new message with STOP/START request from master",MPI_Wtime());
			fclose(file);
		}

		MPI_Recv(&buff,1,struct_type,0,MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
		
		//////Function of slave if the a 'STOP' message is recieved
		if(stat.MPI_TAG==STOP)
		{
			if (dbg_infos)
			{
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] STOP: Received a stop command from master the rank %d is exiting ",MPI_Wtime(),rank);
				fprintf(file,"\n[%f] STOP: Going to send the information about maximal vertical steps in long files %d ",MPI_Wtime(),NSTEP1);
				fclose(file);
			}
			finalize_simulation_slave();
			if (MAXBUFSLAVE)
			{
			write_block_fort(MAXBUFSLAVE,Block); //Write out particle block, MAXBUFSLAVE = maxbuf after 1st run
			maxblock_switch = 32000;//Finalising write out of particle blocks, 32000 > 23456 which is used in write_block_fort to check for last run
			write_block_fort(MAXBUFSLAVE,Block);//write out Event end and Run end blocks
			}
			else
			{
			if (dbg_infos)
			{
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] DATA: No data in buffer of particles Buffer is %d ",MPI_Wtime(),MAXBUFSLAVE);
				fclose(file);
			}
			}
			
#if __COAST__
			if (dbg_infos)
			{
			file = fopen(statfile,"a"); 
			fprintf(file,"\n[%f] CoREAS: Slave asking COAST interface to send the particular outputs",MPI_Wtime());
			fclose(file);
			}
			mpinodedone_();
			if (dbg_infos)
			{
			file = fopen(statfile,"a"); 
			fprintf(file,"\n[%f] CoREAS: The particular outputs sent to master successfuly",MPI_Wtime());
			fclose(file);
			}
#endif
			free(particle_stack_or_cut_file);
		}
		//////FUNCTION OF SLAVE IF THE A 'START' MESSAGE IS RECIEVED (REFER THE DOCUMENTAION)
		if(stat.MPI_TAG==START)
		{

			if (dbg_infos) 
			{
				job_count++;		//Counting initial particle-stacks received by slave	  
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] START: Received a request with id \"%d\" generated by parent particle \"%d\" with seed \"%d\" and primary status \"%d\" and have to run particles from \"%d\" till \"%d\" from stack/cutfile %p",MPI_Wtime(),buff.uid,buff.parent_id,buff.seed,buff.lprim,buff.id1,buff.id2,buff.stack);
				fclose(file);
			}

			particles_from_stack = (double*)calloc((PARTICLE_INFO_COUNT+1)*(buff.id2-buff.id1+1)*MAX_GROUP_SIZE,sizeof(double)); //Stack of new secondary particles - stand alone or group
/*				if (dbg_infos) 
				{
					file = fopen(statfile,"a");
					fprintf(file,"\n[%f] Initialized stack ",MPI_Wtime());					
					for(i=0;i<(PARTICLE_INFO_COUNT+1)*(buff.id2-buff.id1+1)*MAX_GROUP_SIZE;i++)
					fprintf(file,"particles_from_stack[%d]=%f ", i, particles_from_stack[i] );
					fclose(file);
				}
*/
			//Storing current particle stack into intermediate buffer
			for(i=0;i<(buff.id2-buff.id1+1);i++)
				{
				if (dbg_infos) 
				{
					file = fopen(statfile,"a"); 
					fprintf(file,"\n[%f] START: Data of received particle nr. %d :",MPI_Wtime(), i);
					fclose(file);
				}
				for(j=0;j<(PARTICLE_INFO_COUNT+1);j++)
					{
					particles_from_stack[(PARTICLE_INFO_COUNT+1)*i+j]=buff.stack[(PARTICLE_INFO_COUNT+1)*i+j];
					if (dbg_infos) 
					{
						file = fopen(statfile,"a"); 
						fprintf(file,"%f -> %f ", buff.stack[(PARTICLE_INFO_COUNT+1)*i+j],particles_from_stack[(PARTICLE_INFO_COUNT+1)*i+j]);
						fclose(file);
					}
					}
				}	
			II1 = buff.id1;
			
			//Reconstruct the CUTFILE NAME to pass to CORSIKA
			reconstruct_cfname(cutfilename, buff.rnnum, buff.seed, buff.parent_id);

			//Reconstruct the OUTPUT FILE NAME to pass to CORSIKA
			if (dbg_infos) 
			{	
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] START: Seed =%f will be used for output filename ",MPI_Wtime(), particles_from_stack[(PARTICLE_INFO_COUNT+1)*(buff.id2-II1) + 18]);
				fclose(file);
			}
			double pseed=particles_from_stack[(PARTICLE_INFO_COUNT+1)*(buff.id2-II1) + 18];
			if ( pseed == 0. ) pseed=buff.seed;
			reconstruct_ofname(outfname,buff.rnnum,pseed,buff.uid);
		 
			if (dbg_infos)
			{
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] START: will start %dth parallel run for sub-shower with ID =%d requested by ParentID =%d with seed=%d using options lprim=%d with particle %d till %d, CUTFILE=%s, outputfile=%s, steering file=%s",MPI_Wtime(),job_count,buff.uid,buff.parent_id,buff.seed,buff.lprim,buff.id1,buff.id2,cutfilename,outfname,input_fname);
				fclose(file);
				//Record status of "Master to " communication in statistic files
				strcpy(statfile2,statdir);
				strcat(statfile2,"Master2SlaveRecv");
				comm = fopen(statfile2,"a");              
				fprintf(comm,"%d %d %d %d %d %d %d %d %p\n",buff.uid,buff.parent_id,buff.lprim,buff.rnnum,buff.seed,buff.uid,buff.id1,buff.id2,buff.stack);
				fclose(comm);
			}
			//recording start time for CORSIKA
			mpi_time[0]=MPI_Wtime();
		 
			//CALL of CORSIKA SUBROUTINE
			corsika_(&(buff.lprim),&dectcut,&dectmax,&(buff.id1),&(buff.id2),cutfilename,outfname,input_fname,&(buff.uid));

			//recording stop time for CORSIKA
			mpi_time[1]=MPI_Wtime();
			
			if (dbg_infos)
			{
				file = fopen(statfile,"a"); 
				fprintf(file,"\n[%f] TIME: The starting and ending time of corsika instance for particle with unique id= %d is (%f,%f)",MPI_Wtime(),buff.uid,mpi_time[0],mpi_time[1]);
				fclose(file);
			}
			///////////////starting time and end time for each corsika call///////////////////////////
			strcpy(statfile2,statdir);
			strcat(statfile2,"corsika_timetable-");
			sprintf(numstr, "%06d", rank);
			strcat(statfile2,numstr);
			file = fopen(statfile2,"a");
			fprintf(file,"%02d %02d %18.6f %18.6f %18.6f\n",rank,buff.uid,mpi_time[0],mpi_time[1],mpi_time[1]-mpi_time[0]);
			fclose(file);
			   
			//send message to the master reporting about the finished simulation
			//sending doubles instead of mpi_structure type. Be careful in the master to type cast the buffer
			MPI_Send(mpi_time,2,MPI_DOUBLE,MASTER,FINISH,MPI_COMM_WORLD);   
			free(particles_from_stack);
			free(particle_stack_or_cut_file);
			//////////////////going back to waiting state//////////////////////
goto slave_back;
		}
	}
	MPI_Type_free(&struct_type);
	mpierr = MPI_Finalize();
	if(mpierr!= MPI_SUCCESS)
	{
		printf("ERROR: ERROR FINALIZING THE MPI ENVIRONMENT, THE PROGRAM WILL ABORT\n");
		MPI_Abort(MPI_COMM_WORLD,mpierr);
	}
	else
	printf("SUCCESSFULLY FINALIZED\n");
  return 0;
}

void readparticle_(double* cutpar,int* index)
{
 int j;
 for(j=0;j<PARTICLE_INFO_COUNT;j++)
  {
  cutpar[j] = particles_from_stack[(PARTICLE_INFO_COUNT+1)*(*index-II1) + j] ;
  }
}

void writeparticle_(double* cutpar,int* uid,int* i)
{
 int j, rank;
 char statfile[255];
 FILE* file;
 MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 
 if (dbg_infos)
 {
  char numstr[9];
  strcpy(statfile,statdir);
  strcat(statfile,"mpiid-");
  sprintf(numstr, "%d", rank);
  strcat(statfile,numstr);
  strcat(statfile,"-runprotocol.txt");
  file =fopen(statfile,"a");
  fprintf(file,"\n[%f] CORSIKA: subshower with ID=%d have sent additional particles that must be stored in global particle stack of slave. Current number of particle in group is %d and parameters are:  ",MPI_Wtime(),*uid,*i);
  for(j=0;j<PARTICLE_INFO_COUNT;j++) fprintf(file,"%f ",cutpar[j]);
  fclose(file);
 } 
 for(j=0;j<PARTICLE_INFO_COUNT;j++) particle_stack_or_cut_file[(*i-1)*(PARTICLE_INFO_COUNT+1) + j] = cutpar[j];
 if (dbg_infos)
 {
  file =fopen(statfile,"a");
  fprintf(file,"\n[%f] MEMORY: Additional particles have been stored in global particle stack of slave. Current number of particle in group is %d and parameters are:  ",MPI_Wtime(),*i);
  for(j=0;j<PARTICLE_INFO_COUNT;j++) fprintf(file,"%f ",particle_stack_or_cut_file[(*i-1)*(PARTICLE_INFO_COUNT+1) + j]);
  fclose(file);
 }
}

/*
function new_particle()
Called from function endoffile() to send a REQUEST message to the MASTER. 
*/
void new_particle(int id1,int id2,int rnnum,int seed,int uid)
{
 FILE *file,*comm;
 char statfile[255],numstr[9];
 int i,rank;

 MPI_Datatype struct_type_new_particle,oldtypes[BLOCK_COUNT_MESSAGE]={MPI_INT,MPI_DOUBLE};
 int blockcounts[BLOCK_COUNT_MESSAGE]={7,MAX_GROUP_SIZE*(PARTICLE_INFO_COUNT+1)};//(id2-id1+1)
 MPI_Aint offsets[BLOCK_COUNT_MESSAGE],extent;
 //setup description for the fields in the message structure 
 
 int ofst=0;
 for(i=0;i<BLOCK_COUNT_MESSAGE;i++)
   {
    offsets[i]=ofst;
    MPI_Type_extent(oldtypes[i], &extent);
    ofst = ofst + blockcounts[i]*extent;
   } 
 //Defining structured type and committing it.
 MPI_Type_create_struct(BLOCK_COUNT_MESSAGE, blockcounts, offsets, oldtypes, &struct_type_new_particle);
 MPI_Type_commit(&struct_type_new_particle);
 
 //MPI Buffer for communication
 struct message buff;
 MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 if (dbg_infos)
 {
	strcpy(statfile,statdir);
	strcat(statfile,"mpiid-");
	sprintf(numstr, "%d", rank);
	strcat(statfile,numstr);
	strcat(statfile,"-runprotocol.txt");
	file = fopen(statfile,"a"); 
	fprintf(file,"\n[%f] NEW PARTICLES: Sending request for a simulation of a new group consisting of %d particles from subshower with id= %d from simulation under run number %d and seed %d: ",MPI_Wtime(),id2-id1+1,uid,rnnum,seed);
	fclose(file); 
 }
 int j,c=0;
      for(i=id1;i<=id2;i++)
      {
      for(j=0;j<PARTICLE_INFO_COUNT;j++) 
        {
        buff.stack[c*(PARTICLE_INFO_COUNT+1)+j] = particle_stack_or_cut_file[((i-1)*(PARTICLE_INFO_COUNT+1))+j];
/*		if (dbg_infos)
			{
				file = fopen(statfile,"a"); 
				fprintf(file,"[%d] %f <- [%d] %f ",c*(PARTICLE_INFO_COUNT+1)+j,buff.stack[c*(PARTICLE_INFO_COUNT+1)+j],((i-1)*(PARTICLE_INFO_COUNT+1))+j ,particle_stack_or_cut_file[((i-1)*(PARTICLE_INFO_COUNT+1))+j]);
				fclose(file); 
			}
*/      }
      c++;
      }
 buff.id1=id1;
 buff.id2=id2;   
 buff.lprim=0; // 0-False as of secondary particle
 buff.rnnum=rnnum;
 buff.seed=seed;
 buff.uid=uid;
 buff.parent_id=NA;//Parent ID is available only in Booking database of Master
 
 //Send the REQUEST message to the master
 int mpierr=MPI_Send(&buff,1,struct_type_new_particle,MASTER,REQUEST,MPI_COMM_WORLD);
 	if(mpierr!= MPI_SUCCESS)
	{
	printf("FATAL ERROR: When sending Message from %d to Master with rank %d\n",rank, MASTER);
	MPI_Abort(MPI_COMM_WORLD,mpierr);
	}
 if (dbg_infos) 
   { 
	 file = fopen(statfile,"a"); 
	 fprintf(file,"\n[%f] NEW PARTICLES: New group with %d particles sent to Master from subshower with id= %d from simulation under run number %d and seed %d from stack %p",MPI_Wtime(),id2-id1+1,uid,rnnum,seed,buff.stack);
	 fclose(file); 
	 strcpy(statfile,statdir);
	 strcat(statfile,"Slave2MasterRequest");
	 comm = fopen(statfile,"a");              
	 fprintf(comm,"\n\n[%f] %d %d %d %d %d %d %d %d %p",MPI_Wtime(),buff.uid,buff.parent_id,buff.lprim,buff.rnnum,buff.seed,buff.uid,buff.id1,buff.id2,buff.stack);
	 fclose(comm); 
	}
 MPI_Type_free(&struct_type_new_particle);
}

/*
function endoffile()
Called from Subroutine CORSIKA() with the CUTFILE and PARICLES INFORMATION. 
FOR every particle/Group(separate job) it calls the function new_particle to send a REQUEST message
*/
void endoffile_(int* array,int* size,int* rnnum,int* seed,int* uid)
{ 
 int i,j,k, rank;
 char statfile[255],numstr[9];
 FILE *file,*cf;
 MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 if (dbg_infos)
 {
	strcpy(statfile,statdir);
	strcat(statfile,"mpiid-");
	sprintf(numstr, "%d", rank);
	strcat(statfile,numstr);
	strcat(statfile,"-runprotocol.txt");
	file = fopen(statfile,"a"); 
	fprintf(file,"\n[%f] CORSIKA: There are new %d particle groups from subshower with id= %d from simulation under run number %d and seed %d \n",MPI_Wtime(),*size-1,*uid,*rnnum,*seed);
	fclose(file); 
 }
 if(lcout==1)
	{//CUT file -checkpoint containing the new particles from CORSIKA stack will be generated
	 char cutfile[255];
	 char cutfilename[255];
	 strcpy(cutfile,statdir);
	 reconstruct_cfname(cutfilename,*rnnum,*seed,*uid);
	 strcat(cutfile,cutfilename);
	 cf=fopen(cutfile,"w");
	 for(i=0;i<=(*size)-2;i++)
		{
		if(array[i]==array[i+1]-1)
		{//one particle stack
			fprintf(cf,"%d ",array[i]);
			for(j=0;j<(PARTICLE_INFO_COUNT+1);j++) fprintf(cf,"%E ", particle_stack_or_cut_file[((array[i]-1)*(PARTICLE_INFO_COUNT+1))+j]  );
			fprintf(cf,"\n");
		}
		if(array[i]<array[i+1]-1)
		{//group of particles in stack
			for(k=array[i];k<=array[i+1]-1;k++)
			{
				fprintf(cf,"%d ",k);
				for(j=0;j<(PARTICLE_INFO_COUNT+1);j++) fprintf(cf,"%E ", particle_stack_or_cut_file[((k-1)*(PARTICLE_INFO_COUNT+1))+j]  );
				fprintf(cf,"\n");
			}
		}
		}
	 fclose(cf);
	  if (dbg_infos)
		{
			file = fopen(statfile,"a"); 
			fprintf(file,"\n[%f] CORSIKA: A cutfile containing requests of new sub-shower simulations have been generated",MPI_Wtime()); 
			fclose(file); 
		}
	}
 for(i=0;i<=(*size)-2;i++)
 {
	new_particle(array[i],array[i+1]-1,*rnnum,*seed,*uid);
	if (dbg_infos)
	{
		file = fopen(statfile,"a"); 
		fprintf(file,"\n[%f] CORSIKA: Requests of new sub-shower simulations have been sent. First particle in group has a number %d in stack",MPI_Wtime(),array[i]); 
		fclose(file); 
	}
	array[i]=0;
 }
}

/*
function printstatusstart()
Called from CORSIKA, for debugging purpose
*/
void printstatusstart_(int* uid,double dectcut,double dectmax,int* id1,int* id2,int* lprim,char cfname[255],char input_fname[255],int* lcout)
{
 FILE* file;
 char statfile[255];
 strcpy(statfile,statdir);
 strcat(statfile,"corsika_status_start");
 if(*uid==1)
  {
  file = fopen(statfile,"w"); 
  fprintf(file,"#Unique ID of shower; ID of first particle; ID of last particle; Primary/secondary switch; cutfilename/stack; when using input %s; outputting into %s; Energy thresholds max %f and cut %f; Cutfile generation switch %d; and debugging switch %d\n", input_fname, statdir, dectcut, dectmax, *lcout, dbg_infos); 
  fprintf(file,"%d %d %d %d %s\n",*uid,*id1,*id2,*lprim,"Initial Particle");
  }
 else
  {
  file = fopen(statfile,"a"); 
  fprintf(file,"%d %d %d %d %s\n",*uid,*id1,*id2,*lprim,cfname);
  }
 fclose(file);
}

/*
function printstatusfinish()
Called from CORSIKA, for debugging purpose
*/
void printstatusfinish_(int* uid,double dectcut,double dectmax,int* id1,int* id2,int* lprim,char cfname[255],char input_fname[255],int* lcout)
{
 FILE* file;
 char statfile[255];
 strcpy(statfile,statdir);
 strcat(statfile,"corsika_status_finish");
 if(*uid==1)
  {
  file = fopen(statfile,"w"); 
  fprintf(file,"#Unique ID of shower; ID of first particle; ID of last particle; Primary/secondary switch cutfilename/stack; when using input %s; outputting into %s; Energy thresholds max %f and cut %f; Cutfile generation switch %d; and debugging switch %d\n", input_fname, statdir, dectcut, dectmax, *lcout, dbg_infos); 
  fprintf(file,"%d %d %d %d %s\n",*uid,*id1,*id2,*lprim,"Initial Particle");
  }
 else
  {
  file = fopen(statfile,"a");
  fprintf(file,"%d %d %d %d %s\n",*uid,*id1,*id2,*lprim,cfname);
  }
 fclose(file);
}

//function storetext() called from CORSIKA, to store text for longfiles.
void storetext_(int*type,int*nstep,double *thstep,int *nrrun)
{	NSTEP1 = *nstep;
	char statfile[255],numstr[9];
	FILE* comm;
	if (dbg_infos) 
	{	int rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		strcpy(statfile,statdir);
		strcat(statfile,"mpiid-");
		sprintf(numstr, "%d", rank);
		strcat(statfile,numstr);
		strcat(statfile,".long");
		comm = fopen(statfile,"a");
		fprintf(comm,"[%f]:Text of type=%d sent by CORSIKA\n",MPI_Wtime(),(*type));
		fclose(comm);
	}
	
	switch(*type)
	{
		case 0:
		{
			sprintf(lngfile.pstr, " LONGITUDINAL DISTRIBUTION IN %5d VERTICAL STEPS OF%6.1f G/CM**2 FOR SHOWER  %06d\n DEPTH     GAMMAS   POSITRONS   ELECTRONS         MU+         MU-     HADRONS     CHARGED      NUCLEI   CHERENKOV\n",NSTEP1,*thstep,runnumber);
			if (dbg_infos) 
			{
				FILE* comm = fopen(statfile,"a");
				fprintf(comm,"pstr=%s\n",lngfile.pstr);
				fclose(comm);
			}
		} break;
		case 1:
		{
			if (dbg_infos) 
			{  
				FILE* comm = fopen(statfile,"a");
				fprintf(comm,"pstr=%s\n",lngfile.pstr);
				fclose(comm);
			}
			
		} break;
		case 2:
		{  	
   sprintf(lngfile.estr, " LONGITUDINAL ENERGY DEPOSIT IN %5d VERTICAL STEPS OF%6.1f G/CM**2 FOR SHOWER  %06d\n DEPTH    GAMMA     EM IONIZ     EM CUT     MU IONIZ     MU CUT     HADR IONIZ   HADR CUT   NEUTRINO      SUM\n",NSTEP1,*thstep,runnumber);
			if (dbg_infos) 
				{  
				FILE* comm = fopen(statfile,"a");
				fprintf(comm,"estr=%s\n",lngfile.estr);
				fclose(comm);
				}
			
		} break;
		case 3:
		{
		sprintf(lngfile.estr, " LONGITUDINAL ENERGY DEPOSIT IN %5d SLANT  STEPS OF%6.1f  G/CM**2 FOR SHOWER %06d \n DEPTH    GAMMA     EM IONIZ     EM CUT     MU IONIZ     MU CUT     HADR IONIZ   HADR CUT   NEUTRINO      SUM\n",NSTEP1,*thstep,runnumber);
		if (dbg_infos) 
			{  
			FILE* comm = fopen(statfile,"a");
			fprintf(comm,"estr=%s\n",lngfile.estr);
			fclose(comm);
			}
		} break; 

	}
}
 
//function joinmatrix() called inside CORSIKA by each slave to add the content of longfiles generated. 
void joinmatrix_(double longmatrix[NSTEP_MAX*10],int*type)
{ 
	int j;
	char statfile[255],numstr[9];
	FILE *comm;
	if (dbg_infos) 
		{
			int rank; 
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);  
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d", rank);
			strcat(statfile,numstr);
			strcat(statfile,".long");
			comm = fopen(statfile,"a");
			fprintf(comm,"[%f]:longmatrix of type=%d sent by CORSIKA \n",MPI_Wtime(),(*type));
			fclose(comm);
		}
	for(j=0;j<NSTEP1*10;j++)
	{
		if (dbg_infos) 
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"Buffer[%d,0-9] and corresponding matrix",j);
			fclose(comm);
		}
		switch (*type) 
		{
			case 0:
			{int i = (int)(j/10);
				if (j%10)
					lngfile.ebuffer[i][j%10]=lngfile.ebuffer[i][j%10]+longmatrix[j];
				else			
					lngfile.ebuffer[i][j%10]=longmatrix[j];
				if (dbg_infos) 
					{  
						comm = fopen(statfile,"a");
						fprintf(comm,"[%d] %12.5f -> %12.5f \n",j%10,longmatrix[j],lngfile.ebuffer[i][j%10]);
						fclose(comm);
					}
			} break;
			case 1:
		   {int i = (int)(j/10);
				if(j%10)
					lngfile.pbuffer[i][j%10]=lngfile.pbuffer[i][j%10]+longmatrix[j];
				else			
					lngfile.pbuffer[i][j%10]=longmatrix[j];
				if (dbg_infos) 
					{  
						comm = fopen(statfile,"a");
						fprintf(comm,"[%d] %12.5f -> %12.5f \n",j%10,longmatrix[j],lngfile.pbuffer[i][j%10]);
						fclose(comm);
					}
			} break;
		}
	}
}

//Function called from corsika to store the particle data in buffer of each parallel task
void joindat_(int *maxbuf, int *nsubbl, float outvect[6552] )
{	int i,j,rank;
	MAXBUFSLAVE=*maxbuf;
	char statfile[255],numstr[9];
	FILE* comm;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (dbg_infos) 
    {		
			strcpy(statfile,statdir);
			strcat(statfile,"mpiid-");
			sprintf(numstr, "%d",rank);
			strcat(statfile,numstr);
			strcat(statfile,"-outvect.txt");
			comm = fopen(statfile,"a");
			fprintf(comm,"\n[%f]: slave=%d calling function joindat() %dth times storing in local memory the record from array outvect[%d] generated inside CORSIKA\n",MPI_Wtime(),rank,join_counts++,(*maxbuf * *nsubbl));
			fclose(comm);
	}

	for(j=0;j<(*nsubbl);j++)
	{
		if (dbg_infos&&(rank<LISTFILE_MAX))
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"\n Subblock Nr.%d outvect[%d]=%f\n",j,j*(*maxbuf),outvect[j*(*maxbuf)]);
			fclose(comm);
		}
		if (switch_first_runH || switch_first_evnH)
		{
		MPI_Datatype struct_type_event,oldtypes[BLOCK_COUNT_MESSAGE]={MPI_INT,MPI_DOUBLE};
		int blockcounts[BLOCK_COUNT_MESSAGE]={7,MAX_GROUP_SIZE*(PARTICLE_INFO_COUNT+1)};
		MPI_Aint offsets[BLOCK_COUNT_MESSAGE],extent;
		//setup description for the fields in the message structure 
		int ofst=0;
		for(i=0;i<BLOCK_COUNT_MESSAGE;i++)
		{
		offsets[i]=ofst;
		MPI_Type_extent(oldtypes[i], &extent);
		ofst = ofst + blockcounts[i]*extent;
		} 
		//Defining structured type and committing it.
		MPI_Type_create_struct(BLOCK_COUNT_MESSAGE, blockcounts, offsets, oldtypes, &struct_type_event);
		MPI_Type_commit(&struct_type_event);
		struct message buff;
		MPI_Status status;
		
		if ( 211285.1 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)] < 211285.4 ) 
			{//========== RUNH subblock found =========================================
			if (dbg_infos&&(rank<LISTFILE_MAX))
				{  
						comm = fopen(statfile,"a");
						fprintf(comm,"\n[%f]: Run Header subblock",MPI_Wtime());
						fclose(comm);
					}
			for(i=0;i<(*maxbuf);i++) 
				{
						Run_header[i]=outvect[j*(*maxbuf)+i];
						if (dbg_infos&&(rank<LISTFILE_MAX)) 
						{  
							comm = fopen(statfile,"a");
							fprintf(comm,"outvect[%d]=%f\t",j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
							fclose(comm);
						}
						buff.stack[i]=Run_header[i];
				}
			buff.rnnum=(*maxbuf);
			i=MPI_Send(&buff,1,struct_type_event,MASTER,STORRUNH,MPI_COMM_WORLD);
			i=MPI_Recv(&buff,1,struct_type_event,MASTER,STORRUNH,MPI_COMM_WORLD,&status);
			if(i!= MPI_SUCCESS)
				{
				printf("FATAL ERROR: When sending Message from %d to Master with rank %d\n",rank, MASTER);
				MPI_Abort(MPI_COMM_WORLD,i);
				}
			for(i=0;i<(*maxbuf);i++) Run_header[i]=buff.stack[i];
			switch_first_runH=0;
			}
//========== EVTH subblock found =========================================
			else if ( 217433.0 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)] < 217433.2 )  
				{
				if (dbg_infos&&(rank<LISTFILE_MAX)) 
						{  
							comm = fopen(statfile,"a");
							fprintf(comm,"\n[%f]: Event Header subblock",MPI_Wtime());
							fclose(comm);
						}
				for(i=0;i<(*maxbuf);i++) 
					{
						Event_header[i]=outvect[j*(*maxbuf)+i];
						if (dbg_infos&&(rank<LISTFILE_MAX)) 
						{  
							comm = fopen(statfile,"a");
							fprintf(comm,"outvect[%d]=%f\t",j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
							fclose(comm);
						}
						buff.stack[i]=Event_header[i];
					}
				buff.rnnum=(*maxbuf);
				i=MPI_Send(&buff,1,struct_type_event,MASTER,STOREVNH,MPI_COMM_WORLD);
				i=MPI_Recv(&buff,1,struct_type_event,MASTER,STOREVNH,MPI_COMM_WORLD,&status);
				if(i!= MPI_SUCCESS)
					{
					printf("FATAL ERROR: When sending Message from %d to Master with rank %d\n",rank, MASTER);
					MPI_Abort(MPI_COMM_WORLD,i);
					}
				for(i=0;i<(*maxbuf);i++) Event_header[i]=buff.stack[i];
				switch_first_evnH=0;
				}
		MPI_Type_free(&struct_type_event);
		}
//========== LONG subblock found =========================================
		if ( 52815.2 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)] < 52815.3 ) 
					{
					if (dbg_infos&&(rank<LISTFILE_MAX)) 
					{  
						comm = fopen(statfile,"a");
						fprintf(comm,"\n[%f]: Long subblock outvect[%d]=%f",MPI_Wtime(),j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
						fclose(comm);
					}
					}
//========== particle data subblock ======================================
				else if ( ( 1000.0 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)] < 3301.3 ) ||
                ( 3301.4 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)]   < 3397.3 ) ||
                ( 3397.4 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)]  < 52815.2 ) ||
               ( 52815.3 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)]  < 77000.0 ) ||
              ( 116000.0 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)] < 195000.0 )) 
					{
						if (dbg_infos&&(rank<LISTFILE_MAX)) 
							{  
								comm = fopen(statfile,"a");
								fprintf(comm,"\n[%f]: Particle subblock",MPI_Wtime());
								fclose(comm);
							}
						for(i=0;i<*maxbuf;i++)
						{
							Block[particle_block_counts++]=outvect[j*(*maxbuf)+i];
							if (dbg_infos&&(rank<LISTFILE_MAX)) 
							{  
								comm = fopen(statfile,"a");
								fprintf(comm,"outvect[%d]=%f\t",j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
								fclose(comm);
							}
						}
					}
//========== EVTE subblock found =========================================
					else if ( 3397.3 < outvect[j*(*maxbuf)] && outvect[j*(*maxbuf)] <  3397.5 ) 
						   {
							if (dbg_infos&&(rank<LISTFILE_MAX)) 
								{  
								comm = fopen(statfile,"a");
								fprintf(comm,"\n[%f]: Event End subblock outvect[%d]=%f",MPI_Wtime(),j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
								fclose(comm);
								}
							}
//========== RUNE subblock found =========================================
						else if ( 3301.2 < outvect[j*(*maxbuf)]&& outvect[j*(*maxbuf)] < 3301.5 )
							   {
								if (dbg_infos&&(rank<LISTFILE_MAX)) 
									{  
									comm = fopen(statfile,"a");
									fprintf(comm,"\n[%f]: Run End subblock outvect[%d]=%f",MPI_Wtime(),j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
									fclose(comm);
									}
								}
							else if ( outvect[j*(*maxbuf)] == 0.000000 ) 
									if (dbg_infos&&(rank<LISTFILE_MAX)) 
									{  
									comm = fopen(statfile,"a");
									fprintf(comm,"\n[%f]: Empty subblock outvect[%d]=%f",MPI_Wtime(),j*(*maxbuf)+i,outvect[j*(*maxbuf)+i]);
									fclose(comm);
									}
	}
	if( particle_block_counts>=(MAXBLOCKS-21*(*maxbuf)) ) //Write down to escape the stack overflow
	{ 	
		write_block_fort((*maxbuf),Block);
		maxblock_switch++;
		if (dbg_infos) 
		{   
			strcpy(statfile,statdir);
			strcat(statfile,"write_1GB_block.txt");
			FILE* file = fopen(statfile,"a");	  
			fprintf(file,"[%f]: slave with rank %d has written out the block containing %d elements due to limitation of memory. maxblock_switch=%d\n",MPI_Wtime(),rank,particle_block_counts,maxblock_switch);
			fclose(file);
		}
		particle_block_counts=0;
	}
} 
	 
void write_block_fort(int maxbuf, float *Block) 
{
	char statfile[250],numstr[9];
	int rank,i,a,mreclength=21*maxbuf;//nsubbl*maxbuf = 6552 (if thinning) or 5733 (no thinning)
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	sprintf(statfile,"%sDAT%06.0f-%06d",statdir,Run_header[1],rank);//ERhdr[1] contains the run number
	FILE *comm, *file = fopen(statfile,"a");
	if (dbg_infos) 
    {  
		strcpy(statfile,statdir);
		strcat(statfile,"mpiid-");
		sprintf(numstr, "%d",rank);
		strcat(statfile,numstr);
		strcat(statfile,"-datablocks.txt");
		comm = fopen(statfile,"a");
		fprintf(comm,"[%f]: Generating file %sDAT%06.0f-%06d\n",MPI_Wtime(),statdir,Run_header[1],rank);
		fclose(comm);
	}
	a=mreclength*sizeof(float);//Record length = size of 6552 floats
	if ( maxblock_switch < 1 ) //First time recording the DAT output
	{	int nz,ins;
		nz=(particle_block_counts+maxbuf+maxbuf)%mreclength;
		ins=mreclength-nz;
		for(i=0;i<ins;i++) Block[particle_block_counts++]=0.0;//filling rest of particle buffer with 0 to fit the mreclength size blocks (including the first 2 headers)
		if (dbg_infos) 
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"[%f]: Check the new size of Block array %d \n",MPI_Wtime(),particle_block_counts%mreclength);
			fclose(comm);
		}
		// Write first record
		fwrite(&a,sizeof(int),1,file);//Writing record length
		fwrite(Run_header,sizeof(float),maxbuf,file); //sizeof(float) is 4 bytes
		fwrite(Event_header,sizeof(float),maxbuf,file); //sizeof(float) is 4 bytes
		fwrite(Block,sizeof(float),mreclength-2*maxbuf,file);//complete the first record with particles
		fwrite(&a,sizeof(int),1,file);//Writing record length
		if (dbg_infos) 
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"[%f]: The first write_block() call wrote Run %f and Event %f headers before particle blocks\n",MPI_Wtime(),Run_header[0],Event_header[0]);
			fclose(comm);
		}
	        // Write the rest of the particle subblocks from second to the last record
		for(i=0;i<particle_block_counts/mreclength;i++) 
		{
			fwrite(&a,sizeof(int),1,file);
			fwrite(&(Block[mreclength-2*maxbuf+i*mreclength]),sizeof(float),mreclength,file);
			fwrite(&a,sizeof(int),1,file);
		}
		
		if (dbg_infos) 
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"[%f]: The first write_block() call wrote next %d particle blocks\n",MPI_Wtime(),particle_block_counts/mreclength);
			fclose(comm);
		}
 	}
	if ( maxblock_switch >=1  && maxblock_switch < 23456 ) //Additional recordings into output
	{int nz,ins;
		nz=particle_block_counts%mreclength;
		ins=mreclength-nz;
		for(i=0;i<ins;i++) Block[particle_block_counts++]=0.0;//filling rest with 0 to fit the mreclength size blocks.
		for(i=0;i<particle_block_counts/mreclength;i++)
			{
				fwrite(&a,sizeof(int),1,file);
				fwrite(&(Block[i*mreclength]),sizeof(float),mreclength,file);
				fwrite(&a,sizeof(int),1,file);
			}
		if (dbg_infos) 
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"[%f]: The additional write_block() call wrote next %d particle blocks\n",MPI_Wtime(),particle_block_counts/mreclength);
			fclose(comm);
		}
	}
	if ( maxblock_switch > 23456 ) //Last recording into output before finalize
	{float Event_and_Run_end[mreclength];
	 int numtasks;
		MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
		// RUNH = 211285.2812500000000000;
		// EVTH = 217433.0781250000000000;
		// LONG =  52815.2968750000000000;
		// EVTE =   3397.3918457031250000;
		// RUNE =   3301.3325195312500000;

		for( i=0; i<mreclength; i++ ) Event_and_Run_end[i] = 0.0; //Initialising with 0

		Event_and_Run_end[0] = 3397.3918457031250000;
		Event_and_Run_end[1] = 1.*Run_header[1]; //ERhdr[1] contains the run number
		for( i=2; i<7; i++ ) Event_and_Run_end[i] = 1.*particle_block_counts/7; //Estimation (<=) of the number of particles

		Event_and_Run_end[maxbuf] = 3301.3325195312500000;
		Event_and_Run_end[maxbuf+1] = 1.*Run_header[1];
		Event_and_Run_end[maxbuf+2] = 1.*Run_header[1];
		Event_and_Run_end[maxbuf+3] = 1.*numtasks;
		
		fwrite(&a,sizeof(int),1,file);				
		fwrite(Event_and_Run_end,sizeof(float),mreclength,file);
		fwrite(&a,sizeof(int),1,file);
		if (dbg_infos) 
		{  
			comm = fopen(statfile,"a");
			fprintf(comm,"[%f]: The write_block() call wrote %f Run and %f Event end before final close\n",MPI_Wtime(),Event_and_Run_end[0],Event_and_Run_end[maxbuf]);
			fclose(comm);
		}

	}
	
	fclose(file);
}  
