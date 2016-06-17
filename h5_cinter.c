#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <hdf5.h>
#include <errno.h>
#include <sys/stat.h>
#include "h5_cinter.h"
//#include "iscl/os/os.h"
//#include "iscl/log/log.h"


/*! 
 * @brief Tests if dirnm is a directory 
 * 
 * @param[in] dirnm    name of directory to test 
 *
 * @result true -> dirnm is an existing directory 
 *         false -> dirnm is not a directory
 * 
 * @author Ben Baker, ISTI
 *
 */
bool os_path_isdir(const char *dirnm)
{
    struct stat s;
    int err;
    if (dirnm == NULL){return false;}
    if (strlen(dirnm) == 0){return false;}
    err = stat(dirnm,&s);
    if (err == -1){
        if (ENOENT == errno) { //Doesn't exist
            return false;
        }else{ // It exists
            return true;
        }
    }else{ //Exists
        if (S_ISDIR(s.st_mode)){
            return true;
        }else{ //Exists but its a file
            return false;
        }
    }   
}
//============================================================================//
/*! 
 * @brief Tests if filenm is a file
 * 
 * @param[in] filenm    name of file to test 
 * 
 * @result  true  -> filenm is an existing file
 *          false -> filenm is not a file
 *
 * @author Ben Baker, ISTI
 *
 */
bool os_path_isfile(const char *filenm)
{
    struct stat info;
    if (filenm == NULL){return false;}
    if (strlen(filenm) == 0){return false;}
    if (stat(filenm, &info) ==-1){ // Doesn't exist
        return false;
    }else{// Exists
        if (S_ISREG(info.st_mode)){
            return true;
        }else{
            return false;
        }
    }   
}
//============================================================================//
/*!
 * @brief Opens an HDF5 file and returns handle for reading only
 *
 * @param[in] flname   name of HDF5 file to open (NULL terminated)
 *
 * @result file handle
 *
 * @author Ben Baker, ISTI
 *
 */
hid_t h5_open_rdonly(const char *flname)
{
    const char *fcnm = "h5_open_rdonly\0";
    hid_t file_id;
    if (!os_path_isfile(flname))
    {
        printf("%s: HDF5 file %s does not exist!\n", fcnm, flname);
    }
    file_id = H5Fopen(flname, H5F_ACC_RDONLY, H5P_DEFAULT);
    return file_id; 
}
//============================================================================//
/*!
 * @brief Opens an HDF5 file and returns handle for reading and writing
 *
 * @param[in] flname   name of HDF5 file to open (NULL terminated)
 *
 * @result file handle
 *
 * @author Ben Baker, ISTI
 *
 */
hid_t h5_open_rdwt(const char *flname)
{
    hid_t file_id;
    file_id = H5Fopen(flname, H5F_ACC_RDWR, H5P_DEFAULT);
    return file_id;
}
//============================================================================//
/*!
 * @brief Closes an HDF5 file
 *
 * @param[in] file_id    HDF5 file handle
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int h5_close(hid_t *file_id_in)
{
    hid_t ierr, file_id;
    file_id = *file_id_in;
    ierr = H5Fclose(file_id);
    return ierr; 
}
//============================================================================//
/*!
 * @brief Writes a float array to HDF5
 *
 * @param[in] dset_name   name of dataset to write
 * @param[in] file_id     HDF5 file handle
 * @param[in] n           size of dataset to write
 * @param[in] x           dataset to write
 * 
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_write_array__float(const char *dset_name, hid_t *file_id_in,
                          int *nin, const float *x)
{
    const char *fcnm = "h5_write_array__float\0";
    char *citem = (char *)calloc(strlen(dset_name)+1, sizeof(char));
    hid_t flt_dataspace_id, flt_dataset_id;
    hsize_t dims[1];
    herr_t status;
    hid_t file_id = *file_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    //  
    // Copy file handle and name
    strcpy(citem,dset_name);
    // Create dataspace
    dims[0] = n;
    flt_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create dataset
    flt_dataset_id = H5Dcreate2(file_id, citem, H5T_NATIVE_FLOAT,
                                flt_dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write data
    status = H5Dwrite(flt_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, x);
    if (status != 0){
        printf(" %s: Write error\n", fcnm);
        return -1;
    }
    // Close the dataspace
    status  = H5Sclose(flt_dataspace_id);
    status += H5Dclose(flt_dataset_id);
    if (status != 0){
        printf(" %s: Close error\n", fcnm);
    }
    free(citem);
    return status;
}
//============================================================================//
/*!
 * @brief Writes a double array to HDF5
 *
 * @param[in] dset_name   name of dataset to write
 * @param[in] file_id     HDF5 file handle
 * @param[in] n           size of dataset to write
 * @param[in] x           dataset to write
 * 
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_write_array__double(const char *dset_name, hid_t *file_id_in,
                           int *nin, const double *x)
{
    const char *fcnm = "h5_write_array__double\0";
    char *citem = (char *)calloc(strlen(dset_name)+1, sizeof(char));
    hid_t dbl_dataspace_id, dbl_dataset_id;
    hsize_t dims[1];
    herr_t status;
    hid_t file_id = *file_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    //
    // Copy file handle and name
    strcpy(citem, dset_name);
    // Create dataspace
    dims[0] = n;
    dbl_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create dataset
    dbl_dataset_id = H5Dcreate2(file_id, citem, H5T_NATIVE_DOUBLE,
                                dbl_dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write data
    status = H5Dwrite(dbl_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, x);
    if (status != 0){
        printf(" %s: Write error\n", fcnm);
        return -1;
    }
    // Close the dataspace
    status  = H5Sclose(dbl_dataspace_id);
    status += H5Dclose(dbl_dataset_id);
    if (status != 0){
        printf(" %s: Close error\n", fcnm);
    }
    free(citem);
    return status;
}
//============================================================================//
/*!
 * @brief Writes an int array to HDF5
 *
 * @param[in] dset_name   name of dataset to write
 * @param[in] file_id     HDF5 file handle
 * @param[in] n           size of dataset to write
 * @param[in] x           dataset to write
 * 
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_write_array__int(const char *dset_name, hid_t *file_id_in,
                        int *nin, const int *x)
{
    const char *fcnm = "h5_write_array__int\0";
    char *citem = (char *)calloc(strlen(dset_name)+1, sizeof(char));
    hid_t int_dataspace_id, int_dataset_id;
    hsize_t dims[1];
    herr_t status;
    hid_t file_id = *file_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    // Copy file handle and name
    strcpy(citem, dset_name);
    // Create dataspace
    dims[0] = n;
    int_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create dataset
    int_dataset_id = H5Dcreate2(file_id, citem, H5T_NATIVE_INT,
                                int_dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write data
    status = H5Dwrite(int_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, x);
    if (status != 0){
        printf(" %s: Write error\n", fcnm);
        return -1;
    }
    // Close the dataspace
    status  = H5Sclose(int_dataspace_id);
    status += H5Dclose(int_dataset_id);
    if (status != 0){
        printf(" %s: Close error\n", fcnm);
    }
    free(citem);
    return status;
}
//============================================================================//
/*!
 * @brief Writes a char ** array to HDF5
 *
 * @param[in] citem_char  name of char** dataset
 * @param[in] file_id     HDF5 file handle
 * @param[in] citem_char  name of char** dataset
 * @param[in] n           number of items in c
 * @param[in] c           char** dataset to write [n]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int h5_write_array__chars(const char *citem_chr, hid_t *file_id_in,
                          int *nin, const char **c)
{
    const char *fcnm = "h5_write_array__chars\0";
    char **cout, *citem_hdf5;
    int len_item = strlen(citem_chr);
    hid_t chr_dataset_id, chr_dataspace_id, cftype, cmtype;
    hsize_t dims[1];
    herr_t status;
    int i, lens;
    hid_t file_id = *file_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    //  
    // Set the name of the attribute while remembering a null terminator 
    citem_hdf5 = (char *)calloc(len_item+1, sizeof(char));
    memset(citem_hdf5, 0, len_item+1);
    strncpy(citem_hdf5, citem_chr, len_item);
    // Create dataspace
    dims[0] = n;
    chr_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create file and memory types
    cftype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(cftype, H5T_VARIABLE);
    if (status < 0){
        printf("%s: Error setting space \n", fcnm);
        return -1;
    }
    cmtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(cmtype, H5T_VARIABLE);
    if (status < 0){
        printf("%s: Error setting memory space\n", fcnm);
        return -1;
    }
    // Create the dataset
    chr_dataset_id = H5Dcreate2(file_id, citem_hdf5, cftype,
                                chr_dataspace_id,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Create output
    cout = (char **)calloc(n, sizeof(char *));
    for (i=0; i<n; i++){
        lens = strlen(c[i]);
        cout[i] = (char *)calloc(lens+1, sizeof(char));
        memset(cout[i], 0, lens+1);
        strcpy(cout[i], c[i]);
    }
    // Write the data
    status = H5Dwrite(chr_dataset_id, cmtype, H5S_ALL,H5S_ALL,
                      H5P_DEFAULT, cout);
    if (status != 0){
        printf("%s: Write error\n", fcnm);
        return -1;
    }
    status += H5Tclose(cftype);
    status += H5Tclose(cmtype);
    status += H5Sclose(chr_dataspace_id);
    status += H5Dclose(chr_dataset_id);
    if (status != 0){
        printf("%s: Close errors\n", fcnm);
    } 
    // Free space
    for (i=0; i<n; i++){
        free(cout[i]);
    }
    free(cout);
    free(citem_hdf5);
    return status;
}

//============================================================================//
/*!
 * @brief Reads a double array from HDF5
 *
 * @param[in] dset_name   name of dataset to read 
 * @param[in] file_id     HDF5 file handle
 * @param[in] nref        size of dataset to read 
 *
 * @param[out] x          dataset read from disk 
 * 
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_read_array__double(const char *dset_name, hid_t *file_id_in,
                          int *nref_in, double *x)
{
    const char *fcnm = "h5_read_array__double\0";
    char *citem = (char *)calloc(strlen(dset_name)+1, sizeof(char));
    hid_t memspace, dbl_dataspace, dbl_dataset, cparms;
    hsize_t *dims;
    herr_t status;
    int nwork, rank, i;
    hid_t file_id = *file_id_in;
    int nref = *nref_in;
    //------------------------------------------------------------------------//
    //  
    // Copy file hand and name
    strcpy(citem,dset_name);
    // Open dataset
    dbl_dataset = H5Dopen(file_id, citem, H5P_DEFAULT);
    // Create dataspace 
    dbl_dataspace = H5Dget_space(dbl_dataset);
    // Get size of dimensions
    rank = H5Sget_simple_extent_ndims(dbl_dataspace);
    dims = (hsize_t *)calloc(rank, sizeof(hsize_t));
    status = H5Sget_simple_extent_dims(dbl_dataspace, dims, NULL);
    nwork = 1;
    for (i=0; i<rank; i++){
        nwork = nwork*dims[i];
    }
    if (nwork > nref){
        printf("%s: Insufficient space!\n", fcnm);
        return -1;
    }
    // Get properties handle
    cparms = H5Dget_create_plist(dbl_dataset);
    // Define memory space
    memspace = H5Screate_simple(1, dims, NULL);
    // Load data
    status = H5Dread(dbl_dataset, H5T_NATIVE_DOUBLE, memspace, dbl_dataspace,
                     H5P_DEFAULT, x);
    if (status != 0){
        printf("%s: Error loading data\n", fcnm);
        return -1;
    }
    // Close it up
    status  = H5Pclose(cparms);
    status += H5Sclose(dbl_dataspace);
    status += H5Sclose(memspace);
    status += H5Dclose(dbl_dataset);
    if (status != 0){
        printf("%s: Error closing space\n", fcnm);
    }
    free(citem);
    free(dims);
    return status;
}
//============================================================================//
/*!
 * @brief Reads a float array from HDF5
 *
 * @param[in] dset_name   name of dataset to read 
 * @param[in] file_id     HDF5 file handle
 * @param[in] nref        size of dataset to read 
 *
 * @param[out] x          dataset read from disk 
 * 
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_read_array__float(const char *dset_name, hid_t *file_id_in,
                         int *nref_in, float *x)
{
    const char *fcnm = "h5_read_array__float\0";
    char *citem = (char *)calloc(strlen(dset_name)+1, sizeof(char));
    hid_t memspace, flt_dataspace, flt_dataset, cparms;
    hsize_t *dims;
    herr_t status;
    int nwork, rank, i;
    hid_t file_id = *file_id_in;
    int nref = *nref_in;
    //------------------------------------------------------------------------//
    //  
    // Copy file hand and name
    strcpy(citem,dset_name);
    // Open dataset
    flt_dataset = H5Dopen(file_id,citem, H5P_DEFAULT);
    // Create dataspace 
    flt_dataspace = H5Dget_space(flt_dataset);
    // Get size of dimensions
    rank = H5Sget_simple_extent_ndims(flt_dataspace);
    dims = (hsize_t *)calloc(rank, sizeof(hsize_t));
    status = H5Sget_simple_extent_dims(flt_dataspace, dims, NULL);
    nwork = 1;
    for (i=0; i<rank; i++){
        nwork = nwork*dims[i];
    }
    if (nwork > nref){
        printf("%s: Insufficient space!\n", fcnm);
        return -1;
    }
    // Get properties handle
    cparms = H5Dget_create_plist(flt_dataset);
    // Define memory space
    memspace = H5Screate_simple(1, dims, NULL);
    // Load data
    status = H5Dread(flt_dataset, H5T_NATIVE_FLOAT, memspace, flt_dataspace,
                     H5P_DEFAULT, x);
    if (status != 0){
        printf("%s: Error loading data\n", fcnm);
        return -1;
    }
    // Close it up
    status  = H5Pclose(cparms);
    status += H5Sclose(flt_dataspace);
    status += H5Sclose(memspace);
    status += H5Dclose(flt_dataset);
    if (status != 0){
        printf("%s: Error closing space\n", fcnm);
    }
    free(citem);
    free(dims);
    return status;
}
//============================================================================//
/*!
 * @brief Reads an integer array from HDF5
 *
 * @param[in] dset_name   name of dataset to read 
 * @param[in] file_id     HDF5 file handle
 * @param[in] nref        size of dataset to read 
 *
 * @param[out] x          dataset read from disk 
 * 
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_read_array__int(const char *dset_name, hid_t *file_id_in,
                       int *nref_in, int *x)
{
    const char *fcnm = "h5_read_array__int\0";
    char *citem = (char *)calloc(strlen(dset_name)+1,sizeof(char));
    hid_t memspace, int_dataspace, int_dataset, cparms;
    hsize_t *dims;
    herr_t status;
    int nwork, rank, i;
    hid_t file_id = *file_id_in;
    int nref = *nref_in;
    //------------------------------------------------------------------------//
    //
    // Copy file hand and name
    strcpy(citem,dset_name);
    // Open dataset
    int_dataset = H5Dopen(file_id,citem,H5P_DEFAULT);
    // Create dataspace 
    int_dataspace = H5Dget_space(int_dataset);
    // Get size of dimensions
    rank = H5Sget_simple_extent_ndims(int_dataspace);
    dims = (hsize_t *)calloc(rank,sizeof(hsize_t));
    status = H5Sget_simple_extent_dims(int_dataspace, dims, NULL);
    nwork = 1;
    for (i=0; i<rank; i++){
        nwork = nwork*dims[i];
    }
    if (nwork > nref){
        printf("%s: Insufficient space!\n", fcnm);
        return -1;
    }
    // Get properties handle
    cparms = H5Dget_create_plist(int_dataset);
    // Define memory space
    memspace = H5Screate_simple(1, dims, NULL);
    // Load data
    status = H5Dread(int_dataset, H5T_NATIVE_INT, memspace, int_dataspace,
                     H5P_DEFAULT, x);
    if (status != 0){
        printf("%s: Error reading data\n", fcnm);
        return -1;
    }
    // Close it up
    status  = H5Pclose(cparms);
    status += H5Sclose(int_dataspace);
    status += H5Sclose(memspace);
    status += H5Dclose(int_dataset);
    if (status != 0){
        printf("%s: Error closing h5\n", fcnm);
    }
    free(citem);
    free(dims);
    return status;
}
//============================================================================//
/*!
 * @brief Writes double attributes to an HDF5 dataset
 *
 * @param[in] citem        attribute name 
 * @param[in] hdf5_id      HDF5 dataset handle
 * @param[in] n            number of attributes
 * @param[in] attr_data    double attributes [n]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int h5_write_attribute__double(const char *citem, hid_t *hdf5_id_in,
                               int *nin, const double *attr_data)
{
    const char *fcnm = "h5_write_attribute__double\0";
    char *citem_hdf5;
    int len_item = strlen(citem);
    hid_t attr_dataspace_id, attribute_id;
    hsize_t dims[1];
    herr_t status;
    hid_t hdf5_id = *hdf5_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    //
    // Set item name
    citem_hdf5 = (char *)calloc(len_item+1, sizeof(char));
    strncpy(citem_hdf5, citem, len_item);
    // Create a dataspace for the attribute
    dims[0] = n;
    attr_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create dataset id for the attribute
    attribute_id = H5Acreate2(hdf5_id, citem_hdf5, H5T_NATIVE_DOUBLE,
                              attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    // Write the attribute data
    status  = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, attr_data);
    status += H5Aclose(attribute_id);
    status += H5Sclose(attr_dataspace_id);
    if (status != 0){
        printf("%s: Error writing attributes\n", fcnm);
    }
    // Free space 
    free(citem_hdf5);
    citem_hdf5 = NULL;
    return status;
}
//============================================================================//
/*!
 * @brief Writes integer attributes to an HDF5 dataset
 *
 * @param[in] citem        attribute name 
 * @param[in] hdf5_id      HDF5 dataset handle
 * @param[in] n            number of attributes
 * @param[in] attr_data    integer attributes [n]
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int h5_write_attribute__int(const char *citem, hid_t *hdf5_id_in,
                            int *nin, const int *attr_data)
{
    const char *fcnm = "h5_write_attribute__int\0";
    char *citem_hdf5;
    int len_item = strlen(citem);
    hid_t attr_dataspace_id, attribute_id;
    hsize_t dims[1];
    herr_t status;
    hid_t hdf5_id = *hdf5_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    //  
    // Set item name
    citem_hdf5 = (char *)calloc(len_item+1, sizeof(char));
    strncpy(citem_hdf5, citem, len_item);
    // Create a dataspace for the attribute
    dims[0] = n;
    attr_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create dataset id for the attribute
    attribute_id = H5Acreate2(hdf5_id, citem_hdf5, H5T_NATIVE_INT,
                              attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    // Write the attribute data
    status  = H5Awrite(attribute_id, H5T_NATIVE_INT, attr_data);
    status += H5Aclose(attribute_id);
    status += H5Sclose(attr_dataspace_id);
    if (status != 0){ 
        printf("%s: Error writing attributes\n", fcnm);
    }   
    // Free space 
    free(citem_hdf5);
    citem_hdf5 = NULL;
    return status;
}
//============================================================================//
/*!
 * @brief Writes character attributes to an HDF5 dataset
 *
 * @param[in] citem       attribute name 
 * @param[in] hdf5_id     HDF5 dataset handle
 * @param[in] n           number of attributes
 * @param[in] cattr       charater attributes [n]
 *
 * @result 0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
int h5_write_attribute__char(const char *citem, hid_t *hdf5_id_in,
                             int *nin, const char **cattr)
{
    const char *fcnm = "h5_write_attribute__char\0";
    char **cout, *citem_hdf5;
    int len_item = strlen(citem);
    int i, lens;
    hid_t attr_dataspace_id, attribute_id, cftype, cmtype;
    hsize_t dims[1];
    herr_t status;
    hid_t hdf5_id = *hdf5_id_in;
    int n = *nin;
    //------------------------------------------------------------------------//
    //
    // Set item name
    citem_hdf5 = (char *)calloc(len_item+1,sizeof(char));
    strncpy(citem_hdf5, citem, len_item);
    // Create a dataspace 
    dims[0] = n;
    attr_dataspace_id = H5Screate_simple(1, dims, NULL);
    // Create the datatypes 
    cftype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(cftype, H5T_VARIABLE);
    if (status < 0){
        printf("%s: Error setting space!\n",fcnm);
        return -1;
    }
    cmtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(cmtype, H5T_VARIABLE);
    if (status < 0){ 
        printf("%s: Error setting memory space\n", fcnm);
        return -1; 
    }
    // Create dataset for attribute
    attribute_id = H5Acreate2(hdf5_id, citem_hdf5, cftype,
                              attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    // Create output
    cout = (char **)calloc(n, sizeof(char *));
    for (i=0; i<n; i++){
        lens = strlen(cattr[i]);
        cout[i] = (char *)calloc(lens+1, sizeof(char));
        strncpy(cout[i], cattr[i], lens);
    } 
    // Write the attributes
    status  = H5Awrite(attribute_id, cmtype, cout);
    // Close it up
    status += H5Aclose(attribute_id);
    status += H5Tclose(cftype);
    status += H5Tclose(cmtype);
    status += H5Sclose(attr_dataspace_id);
    if (status != 0){
        printf("%s: Error writing attribute!\n", fcnm);
    }
    // Free space
    for (i=0; i<n; i++){
        free(cout[i]);
        cout[i] = NULL;
    }
    free(cout);
    free(citem_hdf5);
    citem_hdf5 = NULL;
    cout = NULL;
    return status;
} 
//============================================================================//
/*!
 * @brief Gets the number of members in a group
 *
 * @param[in] group_name    name of HDF5 group
 * @param[in] file_id       HDF5 file handle
 *
 * @result number of members in HDF5 group (0 is empty group)
 *
 * @author Ben Baker, ISTI
 *
 */
int h5_n_group_members(const char *group_name, hid_t *file_id_in)
{
    int nmember;
    char *citem = (char *)calloc(strlen(group_name)+1, sizeof(char));
    hid_t group_id;
    hsize_t nobj;
    herr_t status;
    hid_t file_id = *file_id_in;
    //------------------------------------------------------------------------//
    //
    // copy h5 handle and name
    strcpy(citem,group_name);
    // Check groups exists
    status = H5Gget_objinfo(file_id, citem, 0, NULL);
    if (status == 0){
        // Open the group and get the number of objects
        group_id = H5Gopen2(file_id, citem, H5P_DEFAULT);
        status = H5Gget_num_objs(group_id, &nobj);
        // close it up
        status = H5Gclose(group_id);
        nmember = (int) nobj;
    }else{
        nmember = 0;
    }
    // Free space
    free(citem);
    return nmember;
}
//============================================================================//
/*!
 * @brief Gets the size of an array with name citem
 *
 * @param[in] citem         null terminated name of item of which to 
 *                          identify size
 * @param[in] file_id       HDF5 file handle 
 *
 * @result length of citem
 *
 * @author Ben Baker, ISTI
 *
 */
int h5_get_array_size(const char *citem, hid_t *file_id_in)
{
    const char *fcnm = "h5_get_array_size\0";
    char *citem_hdf5;
    int nwork_out;
    hid_t dataspace_id, dataset_id;
    hsize_t nwork;
    herr_t status;
    hid_t file_id = *file_id_in;
    //------------------------------------------------------------------------//
    //  
    // Copy h5 handle and item name
    citem_hdf5 = (char *)calloc(strlen(citem)+1, sizeof(char));
    strcpy(citem_hdf5,citem);
    // Open dataset
    dataset_id = H5Dopen(file_id, citem_hdf5, H5P_DEFAULT);
    // Open dataspace 
    dataspace_id = H5Dget_space(dataset_id);
    // Get the size of the dataspace 
    nwork = H5Sget_simple_extent_npoints(dataspace_id);
    // Close it up 
    status  = H5Sclose(dataspace_id);
    status += H5Dclose(dataset_id);
    if (status != 0){ 
        printf("%s: Error closing h5\n", fcnm);
        return nwork;
    }   
    free(citem_hdf5);
    nwork_out = (int) nwork;
    return nwork_out;
}
//============================================================================//
/*!
 * @brief Determines if an item exists
 *
 * @param[in] file_id       HDF5 file handle
 * @param[in] citem_in      name of item to find
 * 
 * @result 1 indicates item exists; 0 if the item does not exist
 *
 * @author Ben Baker, ISTI
 *
 */
bool h5_item_exists(const char *citem_in, hid_t *file_id_in)
{
    char *citem = (char *)calloc(strlen(citem_in)+1, sizeof(char));
    int lexist;
    bool exist;
    hid_t file_id = *file_id_in;
    strcpy(citem, citem_in);
    lexist = H5Lexists(file_id, citem, H5P_DEFAULT);
    if (lexist == 1){
        exist = true;
    }else{
        exist = false;
    }
    free(citem);
    return exist;
}
//============================================================================//
/*! 
 * @brief Creates a group in HDF5 with name cgroup
 * 
 * @param[in] cgroup     name of group to create
 * @param[in] file_id    HDF5 file handle 
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
hid_t h5_create_group(hid_t *file_id_in, const char *cgroup)
{
    const char *fcnm = "h5_create_group\0";
    char *cgroup_hdf5;
    int len_item = strlen(cgroup);
    hid_t group_id;
    herr_t status;
    int lexist;
    hid_t file_id = *file_id_in;
    //------------------------------------------------------------------------//
    //  
    // Set group name 
    cgroup_hdf5 = (char *)calloc(len_item+1, sizeof(char));
    memset(cgroup_hdf5, 0, len_item+1);
    strncpy(cgroup_hdf5, cgroup, len_item);
    // Make sure group does not exist 
    lexist = H5Lexists(file_id, cgroup, H5P_DEFAULT);
    if (lexist == 1){
        printf("%s: Warning group exists; skipping\n", fcnm);
        free(cgroup_hdf5);
        return 0;
    }
    // Create group
    group_id = H5Gcreate2(file_id, cgroup_hdf5, H5P_DEFAULT, H5P_DEFAULT,
                          H5P_DEFAULT);
    // Close group
    status = H5Gclose(group_id);
    if (status < 0){
        printf("%s: Error closing group\n", fcnm);
    }
    // Free space
    free(cgroup_hdf5);
    return status;
}
