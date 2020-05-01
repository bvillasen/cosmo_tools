#include<stdio.h>
#include<stdlib.h>
#include"read_snapshot.h"


int main(int argc, char **argv)
{
	char fname[200];
	int nsnaps;

	if(argc!=3)
	{
		printf("./ReadSnapshot basename nsnaps\n");
		fflush(stdout);
		exit(-1);
	}
	sprintf(fname,"%s",argv[1]);
	nsnaps = atoi(argv[2]);
	
	read_snapshots(fname,nsnaps);

	return 0;
}
