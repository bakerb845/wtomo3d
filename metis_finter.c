#include <stdio.h>
#include <stdlib.h>
#include <metis.h>

extern void sort_int_finter(int *values, int *nn);

/*!
 * Converts a mesh to George and Liu style adjacency graph structure
 * with metis 5.  For more information on data structures see metis 5 
 * user manual or George and Liu's classic book.
 *
 * @param[in] job        =1 convert mesh and return size of adjncy (nadj)
 *                       =2 copy xadj/adjncy
 *                       =3 free space
 * @param[in] ne         number of elements in mesh
 * @param[in] nn         number of nodes in mesh
 * @param[inout] nadj    if job is 1, then on output this is the size of 
 *                       the adjacency array
 *                       if job is not 1, then this is the size of the 
 *                       adjacency array but is not used
 * @param[in] iengv_ptr  points from ielem'th element to start index of 
 *                       iengv [ne+1]
 * @param[in] iengv      points from ia'th node of ielem'th to global 
 *                       node number [iengv_ptr[ne]]
 * @param[in] numflag    if 1 then use fortran numbering
 *                       otherwise use c numbering
 *
 * @param[out] xadj      if job is 2 then this is points from the 
 *                       global anchor node to the start of the adjncy
 *                       [nn]
 * @param[out] adjncy    if job is 2 then this is the corresponding 
 *                       column for the i'th node (row) in xadj;
 *                       note this does not include the main diagonal
 *                       [nadj]
 * @param[out] ierr      0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */
void METIS_MeshToNodal_finter(int *job, int *ne, int *nn, int *nadj, 
                              int *iengv_ptr, int *iengv, int *numflag,
                              int *xadj, int *adjncy, int *ierr)
{
    const char *fcnm = "METIS_MeshToNodal_finter\0";
    static int *xadj_ptr = NULL, *adjncy_ptr = NULL; 
    int i, i1, ishift, j, j1, j2, nsort;
    *ierr = 0;
    if (*job == 1)
    {
        *nadj = 0;
        *ierr = METIS_MeshToNodal(ne, nn, iengv_ptr, iengv, numflag,
                                  &xadj_ptr, &adjncy_ptr);
        if (*ierr != METIS_OK)
        {
            printf("%s: Error partitioning mesh\n", fcnm);
            *ierr = 1;
            return;
        }
        *nadj = xadj_ptr[*nn] - 1;
        // Sort
        ishift = 0;
        if (*numflag == 1)
        {
            ishift = 1;
        }
        for (i=0; i<*nn; i++)
        {
            i1 = xadj_ptr[i] - ishift;
            nsort = xadj_ptr[i+1] - xadj_ptr[i];
            if (nsort > 1)
            {
                sort_int_finter(&adjncy_ptr[i1], &nsort);
            }
        }
        *ierr = 0;
    }
    else if (*job == 2)
    {
        for (i=0; i<*nn; i++)
        {
            xadj[i] = xadj_ptr[i]; 
            j1 = xadj_ptr[i] - 1;
            j2 = xadj_ptr[i+1] - 1;
            for (j=j1; j<j2; j++)
            {
                adjncy[j] = adjncy_ptr[j];
            }
        }
        xadj[*nn] = xadj_ptr[*nn];
    }
    else if (*job == 3)
    {
        if (xadj_ptr){free(xadj_ptr);}
        if (adjncy_ptr){free(adjncy_ptr);}
        xadj_ptr = NULL;
        adjncy_ptr = NULL;
    }
    else
    {
        printf("%s: Invalid job\n", fcnm);
        *ierr = 1;
    }
    return;
}

//============================================================================//
/*!
 * Total communication volume minimizing mesh partitioning
 *
 * @param[in] n         number of nodes in mesh 
 * @param[in] nparts    number of parts in which to partition mesh
 * @param[in] xadj      adjacency pointer [n+1]
 * @param[in] adjncy    nodes connected to i'th nodal point (columns of 
 *                      matrix) [xadj[n]-1]
 *
 * @param[out] part     given the i'th nodal point returns the corresponding
 *                      partition number [n]
 * @param[out] ierr     0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
void METIS_PartGraphKway_finter(int *n, int *nparts, int *xadj, int *adjncy, 
                                int *part, int *ierr)
{
    const char *fcnm = "METIS_PartGraphKway_finter\0";
    int options[METIS_NOPTIONS], i, ncon, np, nvtxs, objval;
    // Nothing to do
    *ierr = 0;
    if (*nparts == 1)
    {
        for (i=0; i<*n; i++)
        {
            part[i] = 1;
        }
        return;
    }
    // Partition mesh
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 1;
    nvtxs = *n;
    ncon = 1;
    np = *nparts;
    *ierr = METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy,  
                                NULL, NULL, NULL, &np, NULL,
                                NULL, options, &objval, part);
    if (*ierr != METIS_OK)
    {
        printf("%s: Error partitioning graph\n", fcnm);
        *ierr = 1;
        return;
    }
    *ierr = 0;
    return;
}

//============================================================================//

/*!
 * Partitions a mesh with Metis 5
 *
 * @param[in] ne     number of elements in mesh
 * @param[in] nn     number of nodes in mesh
 * @param[in] eptr   points from e'th element to start of eind [ne+1]
 *                   assumed to have fortran numbering
 * @param[in] eind   points from ia'th node on e'th element to global
 *                   node [eptr[ne]]
 *                   assumed to have fortran numbering
 * @param[in] nparts  number of parts in which to partition mesh
 *
 * @param[out] part   returns the ielem'th elements partition [ne]
 * @param[out] ierr   0 indicates success
 * 
 * @author Ben Baker, ISTI
 *
 */ 
void METIS_PartMeshNodal_finter(int *ne, int *nn, int *eptr, int *eind, 
                                int *nparts, int *part, int *ierr)
{
    const char *fcnm = "METIS_PartMeshNodal_finter\0";
    int options[METIS_NOPTIONS], *npart, i, objval;
    //------------------------------------------------------------------------//
    //
    // nothing to do
    *ierr = 0;
    if (*nparts == 1)
    {
        for (i=0; i<*ne; i++)
        {
            part[i] = 1;
        }
        return;
    }
    // use default options and partition 
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 1;
    npart = (int *)calloc(*nn, sizeof(int));
    *ierr = METIS_PartMeshNodal(ne, nn, eptr, eind, NULL, NULL, 
                                nparts, NULL, options, &objval, part, npart);
    free(npart);
    if (*ierr != METIS_OK)
    {
        printf("%s: Error partitioning graph\n", fcnm);
        *ierr = 1;
        for (i=0; i<*ne; i++)
        {
            part[i] = 1;
        }
        *ierr = 1;
        return;
    }
    *ierr = 0;
    return;
}
