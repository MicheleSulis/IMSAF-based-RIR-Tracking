#include <math.h>
//#include "Commdlg.h"
#include "io.h"
//#include "string.h"
//#include "stdlib.h"
//#include "stdio.h"
#include <fstream>
#include <iostream>
#include "LEEffect.h"
#include "ipp.h"
using namespace std;

#define MAX_FILE_NAME_LENGTH 255

void read_dat(char *name, double *data, int dim);
void write_dat(char *name, double *data, int dim, char*save_dir);

void init_vector(Ipp64f*& vector, int size);
void init_vector(Ipp64fc*& vector, int size);

void destroy_vector(Ipp64f*& vector);
void destroy_vector(Ipp64fc*& vector);