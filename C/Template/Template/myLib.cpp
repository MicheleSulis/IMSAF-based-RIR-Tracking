#include "stdafx.h"
#include "myLib.h"

void write_dat(char *name, double *data, int dim, char *save_name) {

	char name_a[MAX_FILE_NAME_LENGTH];
	memset(name_a,0,MAX_FILE_NAME_LENGTH*sizeof(char));
	strcpy(name_a,save_name);
	strcat(name_a,name);

	char *c;
	c=(char *)(void *)data;

	fstream File;
	File.open(name_a,ios::out | ios::binary);
	File.write(c,dim*sizeof(double));
	File.close();
}

void read_dat(char *name, double *data, int dim) {
	char *c;
	c=(char *)(void *)data;

	fstream File;
	File.open(name,ios::in | ios::binary);
	File.read(c,dim*sizeof(double));
	File.close();
}

void init_vector(Ipp64f*& vector, int size) {
	if (vector == 0) {
		vector = ippsMalloc_64f(size);
		ippsZero_64f(vector, size);
	}
}

void init_vector(Ipp64fc*& vector, int size) {
	if (vector == 0) {
		vector = ippsMalloc_64fc(size);
		ippsZero_64fc(vector, size);
	}
}

void destroy_vector(Ipp64f*& vector) {
	if (vector != 0) {
		ippsFree(vector);
		vector = 0;
	}
}

void destroy_vector(Ipp64fc*& vector) {
	if (vector != 0) {
		ippsFree(vector);
		vector = 0;
	}
}