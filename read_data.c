#include <stdio.h>
#include "global_matrices.h"

/*---------------------------------------------------------*/

void read_data_file(char file_name[30], double data_matrix[289][700])
{
     int col, row;
     FILE   *xfile;
     char   full_file_name[100];
#ifndef DATA_ROOT
#define DATA_ROOT "/home/rossi/src/eccsvm/"
#endif

     sprintf(full_file_name,"%s%s",DATA_ROOT,file_name);
     xfile=fopen(full_file_name,"r");
  
     for (row=0; row<289; ++row)
     {
       for (col=0; col<700; ++col)  
	{
         fscanf(xfile,"%lf",&data_matrix[row][col]);
	}
     }

     fclose(xfile);
}

/*---------------------------------------------------------*/
void read_data (void)
{
  read_data_file("Dom1.dat", Md1);
  read_data_file("Dom2.dat", Md2);
  read_data_file("Dom3.dat", Md3);
  read_data_file("Dom4.dat", Md4);
  read_data_file("Dom5.dat", Md5);
  read_data_file("Dom6.dat", Md6);
  read_data_file("Dom7.dat", Md7);
  read_data_file("Dom8.dat", Md8);
  read_data_file("Dom9.dat", Md9);
  read_data_file("Dom10.dat", Md10);
  read_data_file("Dom11.dat", Md11);
  read_data_file("Dom12.dat", Md12);
  read_data_file("Dom13.dat", Md13);
}
