#include "io.h"
#include<stdlib.h>
#include <unistd.h>



void Write_field_to_file( const string &field_name, float *data_field, int n_data, hid_t file_id  ){
  
  // printf(" Writing Field: %s\n", field_name.c_str() );
  int       k;
  hid_t     dataset_id, dataspace_id; 
  float      *dataset_buffer;
  herr_t    status;
  hsize_t   dims[1];
  
  // Create the data space for the datasets
  dims[0] = n_data;
  dataspace_id = H5Screate_simple(1, dims, NULL);
  
  
  dataset_buffer  = (float *) malloc(n_data*sizeof(float));
  for (k=0; k<n_data; k++){ 
    dataset_buffer[k] = data_field[k];
  }
  // 
  // Create a dataset id for density
  dataset_id = H5Dcreate(file_id, field_name.c_str(), H5T_IEEE_F32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Write the density array to file  // NOTE: NEED TO FIX FOR FLOAT float!!!
  status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_buffer);
  // Free the dataset id
  status = H5Dclose(dataset_id);
  
  // Free the dataspace id
  status = H5Sclose(dataspace_id);
  
  
  
  
  
  
  
  
  
  
  free(dataset_buffer);
  
}
