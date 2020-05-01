#ifndef BRANT_READ_SNAPSHOT
#define BRANT_READ_SNAPSHOT

void read_snapshots(char fname_base[200], int nsnaps);
struct Header
{
	unsigned int Npart[6];
	double Massarr[6];
	double A;
	double Z;
	int FlagSfr;
	int FlagFeedback;
	unsigned int Nall[6];
	int FlagCooling;
	int NumFiles;
	double BoxSize;
	double Omega_0;
	double Omega_L;
	double h;
	int FlagAge;
	int FlagMetals;
	unsigned int NallHW[6];
	int flag_entr_ics;
 	char fill[56];	
};

extern struct Header header;
extern long *id;        //particle id's, runs from 0 to N
extern float *x;        //particle positions, of size x(3,N)
                //units will be in kiloparsecs (1 kpc = 3.26e3 lightyears)
extern float *v;        //particle velocities, of size v(3,N)
                //                //units will be in km/s
                //
extern float *mass;        //particle masses, of size mass(N)
			//units will be in solar masses, of size mbuffer(N_without_m)
extern float *mbuffer;        //a buffer for the particle masses
                //
extern float *u;        //gas internal energy, of size u(npart[0])
                //
extern float *rho;        //gas densities, of size rho(npart[0])
                //
extern float *Ne;        //gas electron densities, of size Ne(part[0])
extern float *Nh;        //gas hydrogen densities, of size Nh(part[0])
extern float *hsml;        //gas smoothing lengths, of size hsml(part[0])
extern float *sfr;        //gas star formation rates, of size sfr(npart[0])
extern float *HeI;        //stellar ages, of size age(npart[4])
extern float *HeII;        //gas and stellar metallicities, of size z(npart[0]+npart[4])
extern float *x_gas;        //stellar ages, of size age(npart[4])
extern float *v_gas;        //stellar ages, of size age(npart[4])
extern float *mass_gas;        //gas and stellar metallicities, of size z(npart[0]+npart[4])
//
extern unsigned int N;             //total number of particles in snapshot
extern unsigned int Ntot;             //total number of particles

void PrintHeader(struct Header header);
void copy_header(struct Header A, struct Header *B);
void AllocateMemory(unsigned int Ntot, unsigned int Ngas);
#endif //

