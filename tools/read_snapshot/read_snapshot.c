#include"read_snapshot.h"
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
#include "io.h"

using namespace std;

struct Header header;

long *id;		//particle id's, runs from 0 to N
float *x;		//particle positions, of size x(3,N)
				//units will be in kiloparsecs (1 kpc = 3.26e3 lightyears) 
float *v;		//particle velocities, of size v(3,N)
				//units will be in km/s 

float *mass;		//particle masses, of size mass(N)
				//units will be in solar masses, of size mbuffer(N_without_m)

float *mbuffer;		//a buffer for the particle masses
float *u;		//gas internal energy, of size u(npart[0])
float *rho;		//gas densities, of size rho(npart[0])

float *Ne;		//gas electron densities, of size Ne(part[0])
float *Nh;		//gas hydrogen densities, of size Nh(part[0])
float *hsml;		//gas smoothing lengths, of size hsml(part[0])
float *sfr;		//gas star formation rates, of size sfr(npart[0])
float *HeI;		//stellar ages, of size age(npart[4])
float *HeII;		//gas and stellar metallicities, of size z(npart[0]+npart[4])
float *mass_gas;
float *x_gas;
float *v_gas;

unsigned int N;     		//total number of particles in snapshot
unsigned int Ntot;     		//total number of particles
void read_snapshots(char fname_base[200], int nsnaps)
{

	int isnap;

	FILE *fp_snap;  //file pointer for the fp_snap file

	struct Header header1;
	int i,j,k;   //loop variables

	int dummy; //dummy in to read in field separators


	int N_without_m=0;	//number of particles whose masses vary from
				//particle to particle

	int n_offset;		//a dummy counter for use with reading in mass array

	float float_buf;	//buffer to read in floats with
	float int_buf;		//buffer to read in integers with

	float *xb;
	unsigned int *ib;
	float *mb;
	float *ub;


	char fname[200];
	char test[8];

	unsigned int Ngas = 0;
	unsigned int Np = 0;
	unsigned int Npg = 0;
	unsigned int Nb = 0;
	unsigned int Ngas_snap;
	//read in all snapshots first to determine total #
	//of particles of each type
	printf("Reading in %d snapshots.\n",nsnaps);
	Ntot=0;
	for(isnap=0;isnap<nsnaps;isnap++)
	{
		sprintf(fname,"%s.%d",fname_base,isnap);
		printf("Reading snapshot %s\n",fname);

		if(!(fp_snap=fopen(fname,"r")))
		{
			printf("Error opening %s.\n",fname);
			exit(-1);
		}

		fread(&dummy,sizeof(int),1,fp_snap);
		//printf("dummy = %d\n",dummy);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		fread(&dummy,sizeof(int),1,fp_snap);
		//printf("dummy = %d\n",dummy);
		//printf("sizeof Struct Header = %ld\n",sizeof(struct Header));
		fread(&header1,sizeof(struct Header),1,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		//printf("dummy = %d\n",dummy);

		N = 0;
		Np  = 0;
		unsigned int Nm = 0;
		if(isnap==0)
		{
			copy_header(header1,&header);	
			for(i=0;i<6;i++)
			{
				Ntot += header1.Nall[i];
				header.Npart[i] = header1.Nall[i];
			}
			Ngas = header1.Nall[0];

			printf("Total particles in all snapshots = %ld\n",Ntot);
			printf("Total gas particles in all snapshots = %ld\n",Ngas);
			AllocateMemory(Ntot,Ngas); //allocate all the particle arrays

			// PrintHeader(header);
		}
		printf("********************\n");
		
		
		for(i=0;i<6;i++)
			N += header1.Npart[i];

		Nm = header1.Npart[0] + header1.Npart[4];
		Ngas_snap = header1.Npart[0];
		printf("Total particles in snapshot %d     = %ld\n",isnap,N);
		printf("Total gas particles in snapshot %d = %ld\n",isnap,Ngas_snap);

		//read in positions
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,sizeof(double));
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,8*sizeof(char));
		// printf("test = %s\n",test);
		
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,3*N*sizeof(float));
		fread(&x[3*Np],sizeof(float),3*N,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		
		//Copy gas positions only
		for ( int k=0; k<3*Ngas_snap; k++){
			x_gas[3*Npg + k] = x[3*Np + k];
		}

		// printf("Done with x\n");
		
		//read in velocities
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,3*N*sizeof(float));
		fread(&v[3*Np],sizeof(float),3*N,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		
		
		//Copy gas velocities only
		for ( int k=0; k<3*Ngas_snap; k++){
			v_gas[3*Npg + k] = v[3*Np + k];
		}

		// printf("Done with v\n");
		//read in IDs
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,N*sizeof(long));
		fread(&id[Np],sizeof(long),N,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,N*sizeof(long));

		// printf("Done with id\n");

		// printf("Nm = %ld (Np = %ld)\n",Nm,Np);
		//read in particle masses	
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		// printf("dummy %d %ld\n",dummy,Nm*sizeof(float));

		xb = (float *) malloc(Nm*sizeof(float));

		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("dummy %d %ld\n",dummy,Nm*sizeof(float));
		fread(xb,sizeof(float),Nm,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with reading xb (Np = %ld)\n",Np);
		
		//Copy gas velocities only
		for ( int k=0; k<Ngas_snap; k++){
			mass_gas[Npg + k] = xb[k];
		}


		//save particle masses	
		for(j=0;j<header1.Npart[0];j++)
			mass[Np+j] = xb[j];
		for(j=0;j<header1.Npart[1];j++)
			mass[Np+header1.Npart[0] + j] = header1.Massarr[1];
		for(j=0;j<header1.Npart[4];j++)
			mass[Np+header1.Npart[0] +header1.Npart[1]+ j] = xb[header1.Npart[0]+j];

		free(xb);
		// printf("Done with mass\n");


		//read in gas u
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&u[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with u\n");
		//read in gas rho
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		// fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&rho[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with rho\n");
		//read in gas Ne 
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&Ne[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with Ne\n");
		//read in gas Nh 
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&Nh[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with Nh\n");
		//read in gas hsml
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&hsml[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with hsml\n");
		//read in gas sfr
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&sfr[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with sfr\n");
		//read in gas HeI
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&HeI[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with HeI\n");
		//read in gas HeII
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&test,sizeof(char),8,fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);
		// printf("test = %s\n",test);
		fread(&dummy,sizeof(int),1,fp_snap);
		fread(&HeII[Npg],sizeof(float),header1.Npart[0],fp_snap);
		fread(&dummy,sizeof(int),1,fp_snap);

		// printf("Done with HeII\n");
		Npg += header1.Npart[0];
		for(i=0;i<6;i++)
			Np += header1.Npart[i];


		
		fclose(fp_snap);
		
		//PrintHeader(header);

	}
	PrintHeader(header);
	
	
	
	string data_dir, output_dir; 
	data_dir = "/data/groups/comp-astro/bruno/";
	output_dir  = data_dir + "cosmo_sims/ewald_512/";
	
	int n_snapshot = 12;
	
	ostringstream  out_file_name;
	string field_name;
	
	//Save to hdf5 file  
	out_file_name.str("");
	out_file_name.clear();
	out_file_name  << "snapshot_" << n_snapshot << ".h5";
	
	printf("Saving File: %s\n",  (output_dir + out_file_name.str()).c_str() );
	
	
	// Create a new file using default properties.
	hid_t   file_id; /* file identifier */
	herr_t  status;
	file_id = H5Fcreate( (output_dir + out_file_name.str()).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); 
	
	// Write Header  
	hsize_t   attr_dims;
	hid_t     attribute_id, dataspace_id;
	int       int_data[3];
	
	
	
	// Single attributes first
	attr_dims = 1;
	dataspace_id = H5Screate_simple(1, &attr_dims, NULL);
	attribute_id = H5Acreate(file_id, "current_z", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
	status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &header.Z);
	status = H5Aclose(attribute_id);

	attribute_id = H5Acreate(file_id, "BoxSize", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
	status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &header.BoxSize);
	status = H5Aclose(attribute_id);

	attribute_id = H5Acreate(file_id, "Omega_L", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
	status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &header.Omega_L);
	status = H5Aclose(attribute_id);

	attribute_id = H5Acreate(file_id, "Omega_M", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
	status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &header.Omega_0);
	status = H5Aclose(attribute_id);

	attribute_id = H5Acreate(file_id, "h", H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT); 
	status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &header.h);
	status = H5Aclose(attribute_id);

	// Close the dataspace
	status = H5Sclose(dataspace_id);
	
	field_name = "mass";
	Write_field_to_file( field_name, mass_gas, Ngas, file_id  );
	
	field_name = "x";
	Write_field_to_file( field_name, x_gas, 3*Ngas, file_id  );
	
	field_name = "v";
	Write_field_to_file( field_name, v_gas, 3*Ngas, file_id  );
	
	field_name = "rho";
	Write_field_to_file( field_name, rho, Ngas, file_id  );
	
	field_name = "u";
	Write_field_to_file( field_name, u, Ngas, file_id  );
	
	field_name = "hsml";
	Write_field_to_file( field_name, hsml, Ngas, file_id  );
	
	field_name = "Ne";
	Write_field_to_file( field_name, Ne, Ngas, file_id  );
	
	field_name = "Nh";
	Write_field_to_file( field_name, Nh, Ngas, file_id  );
	
	field_name = "HeI";
	Write_field_to_file( field_name, HeI, Ngas, file_id  );
	
	field_name = "HeII";
	Write_field_to_file( field_name, HeII, Ngas, file_id  );
	
	// close the file
	status = H5Fclose(file_id);
	if (status < 0) {printf("File write failed.\n"); exit(-1); }
	
	
	printf("Saved File: %s\n", (output_dir + out_file_name.str()).c_str());



	// printf("Saving Files\n" );
	// 
	// int N_gas_w, N_tot_w;
	// N_gas_w = Ngas;
	// N_tot_w = Ntot;
	// 
	// ofstream myfile;
	
	// printf(" Saving: pos   %d  particles \n", N_tot_w );
  // myfile.open ("snap_12/x.dat");
  // float x_val, y_val, z_val;
	// for ( int i=0; i<N_tot_w; i++ ){
	// 	x_val = x[i*3 + 0];
	// 	y_val = x[i*3 + 1];
	// 	z_val = x[i*3 + 2];
	// 	// printf("%f  %f  %f \n", x_val, y_val, z_val );
	// 	myfile << x_val << " " << y_val << " " << z_val << endl;
	// }
	// myfile.close();
	// 
	// 
	// printf(" Saving: vel   %d  particles \n", N_tot_w );
  // myfile.open ("snap_12/v.dat");
  // float vx_val, vy_val, vz_val;
	// for ( int i=0; i<N_tot_w; i++ ){
	// 	vx_val = v[i*3 + 0];
	// 	vy_val = v[i*3 + 1];
	// 	vz_val = v[i*3 + 2];
	// 	myfile << vx_val << " " << vy_val << " " << vz_val << endl;
	// }
	// myfile.close();
	// 
	// 
	// printf(" Saving: mass   %d  particles \n", N_tot_w );
	// // myfile.open ("snap_12/mass.dat");
	// int zero_counter = 0;
	// float p_mass;
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	p_mass = mass_gas[i];
	// 	if ( i < N_gas_w  && p_mass == 0 ) zero_counter += 1;
	// 
	// 	// myfile << mass[i] << endl;
	// }
	// // myfile.close();
	// printf( "Gas Particles with zero mass: %d   fraction: %f\n", zero_counter, float(zero_counter)/N_gas_w);
	// 
	// 
	
	
	// 
	// printf(" Saving: u   %d  particles \n", N_gas_w );
	// myfile.open ("snap_12/u.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << u[i] << endl;
	// }
	// myfile.close();
	// 
	// printf(" Saving: rho   %d  particles \n", N_gas_w);
	// myfile.open ("snap_12/rho.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << rho[i] << endl;
	// }
	// myfile.close();
	// 
	// 
	// printf(" Saving: Ne   %d  particles \n", N_gas_w );
	// myfile.open ("snap_12/Ne.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << Ne[i] << endl;
	// }
	// myfile.close();
	// 
	// 
	// printf(" Saving: Nh   %d  particles \n", N_gas_w );
	// myfile.open ("snap_12/Nh.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << Nh[i] << endl;
	// }
	// myfile.close();
	// 
	// printf(" Saving: hsml   %d  particles \n", N_gas_w );
	// myfile.open ("snap_12/hsml.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << hsml[i] << endl;
	// }
	// myfile.close();
	
	
	// printf(" Saving: HeI   %d  particles \n", N_gas_w );
	// myfile.open ("snap_11/HeI.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << HeI[i] << endl;
	// }
	// myfile.close();
	// 
	// 
	// printf(" Saving: HeII   %d  particles \n", N_gas_w );
	// myfile.open ("snap_11/HeII.dat");
	// for ( int i=0; i<N_gas_w; i++ ){
	// 	myfile << HeII[i] << endl;
	// }
	// myfile.close();
	

	
	
	
	
}
void AllocateMemory(unsigned int Ntot, unsigned int Ngas)
{
	id = (long *) malloc(Ntot*sizeof(long));
	x = (float *) malloc(3*Ntot*sizeof(float));
	v = (float *) malloc(3*Ntot*sizeof(float));
	mass = (float *) malloc(Ntot*sizeof(float));
	u  = (float *) malloc(Ngas*sizeof(float));
	rho = (float *) malloc(Ngas*sizeof(float));
	Ne = (float *) malloc(Ngas*sizeof(float));
	Nh = (float *) malloc(Ngas*sizeof(float));
	hsml = (float *) malloc(Ngas*sizeof(float));
	sfr = (float *) malloc(Ngas*sizeof(float));
	HeI = (float *) malloc(Ngas*sizeof(float));
	HeII = (float *) malloc(Ngas*sizeof(float));
	x_gas = (float *) malloc(3*Ngas*sizeof(float));
	v_gas = (float *) malloc(3*Ngas*sizeof(float));
	mass_gas = (float *) malloc(Ngas*sizeof(float));
	
}
void PrintHeader(struct Header header)
{

        int i;
        for(i=0;i<6;i++)
                fprintf(stdout,"Header.Npart[%d] = %d\n",i,header.Npart[i]);
        for(i=0;i<6;i++)
                fprintf(stdout,"Header.Massar[%d] = %f\n",i,header.Massarr[i]);
        fprintf(stdout,"Header.A = %f\n",header.A);
        fprintf(stdout,"Header.Z = %f\n",header.Z);
        fprintf(stdout,"Header.FlagSfr = %d\n",header.FlagSfr);
        fprintf(stdout,"Header.FlagFeedback = %d\n",header.FlagFeedback);
        for(i=0;i<6;i++)
                fprintf(stdout,"Header.Nall[%d] = %d\n",i,header.Nall[i]);
        fprintf(stdout,"Header.FlagCooling = %d\n",header.FlagCooling);
        fprintf(stdout,"Header.NumFiles = %d\n",header.NumFiles);
        fprintf(stdout,"Header.BoxSize = %f\n",header.BoxSize);
        fprintf(stdout,"Header.Omega_0 = %f\n",header.Omega_0);
        fprintf(stdout,"Header.Omega_L = %f\n",header.Omega_L);
        fprintf(stdout,"Header.h       = %f\n",header.h);
        fprintf(stdout,"Header.FlagAge = %d\n",header.FlagAge);
        fprintf(stdout,"Header.FlagMetals = %d\n",header.FlagMetals);
        for(i=0;i<6;i++)
                fprintf(stdout,"Header.NallHW[%d] = %d\n",i,header.NallHW[i]);
        fprintf(stdout,"Header.flag_entr_ics = %d\n",header.flag_entr_ics);

        printf("Total number of particles: N = %d\n",Ntot);

}

void copy_header(struct Header A, struct Header *B)
{

    for(int i=0;i<6;i++)
    {
        B->Npart[i]    = A.Npart[i];
        B->Nall[i]     = A.Nall[i];
        B->NallHW[i]     = A.NallHW[i];
        B->Massarr[i]  = A.Massarr[i];
    }

    for(int i=0;i<56;i++)
    {
        B->fill[i] = A.fill[i];
    }

    B->A = A.A;
    B->Z = A.Z;
    B->FlagSfr      = A.FlagSfr;
    B->FlagFeedback = A.FlagFeedback;
    B->FlagCooling  = A.FlagCooling;
    B->NumFiles     = A.NumFiles;
    B->BoxSize      = A.BoxSize;
    B->Omega_0      = A.Omega_0;
    B->Omega_L      = A.Omega_L;
    B->h            = A.h;
    B->flag_entr_ics = A.flag_entr_ics;
    B->FlagAge   = A.FlagAge;
    B->FlagMetals = A.FlagMetals;
}
