#ifndef IO_H
#define IO_H

#include <iostream>
#include <string>
#include <stdio.h>
#include<hdf5.h>

using namespace std;


void Write_field_to_file( const string &field_name, float *data_field, int n_data, hid_t file_id  );


#endif //IO_H