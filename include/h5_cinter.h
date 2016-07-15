#include <stdbool.h>
#include <hdf5.h>
#ifndef __H5_CINTER_H__
#define __H5_CINTER_H__
#ifdef __cplusplus
extern "C" 
{
#endif

hid_t h5_open_rdonly(const char *flname);
hid_t h5_open_rdwt(const char *flname);
int h5_close(hid_t *file_id);
int h5_write_array__float(const char *dset_name, hid_t *file_id,
                          int *n, const float *x);
int h5_write_array__double(const char *dset_name, hid_t *file_id,
                           int *n, const double *x);
int h5_write_array__int(const char *dset_name, hid_t *file_id,
                        int *n, const int *x);
int h5_write_array__chars(const char *citem_chr, hid_t *file_id,
                          int *n, const char **c);
int h5_read_array__double(const char *dset_name, hid_t *file_id,
                          int *nref, double *x);
int h5_read_array__float(const char *dset_name, hid_t *file_id,
                         int *nref, float *x);
int h5_read_array__int(const char *dset_name, hid_t *file_id,
                       int *nref, int *x);
int h5_write_attribute__double(const char *citem, hid_t *hdf5_id,
                               int *n, const double *attr_data);
int h5_write_attribute__int(const char *citem, hid_t *hdf5_id,
                            int *n, const int *attr_data);
int h5_write_attribute__char(const char *citem, hid_t *hdf5_id,
                             int *n, const char **cattr);
int h5_n_group_members(const char *group_name, hid_t *file_id);
int h5_get_array_size(const char *citem, hid_t *file_id);
bool h5_item_exists(const char *citem_in, hid_t *file_id);
hid_t h5_create_group(hid_t *file_id, const char *cgroup);

#ifdef __cplusplus
}
#endif
#endif /* __H5_C_INTER_H__ */
