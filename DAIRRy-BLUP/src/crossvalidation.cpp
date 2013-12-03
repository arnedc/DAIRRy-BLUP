#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "shared_var.h"
#include <shared_var.h>
#include <hdf5.h>

extern "C" {
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void blacs_barrier_ ( int*, char* );
    void pdsyrk_ ( char*, char*, int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int* );
    void pdgemm_ ( char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc );
    void pdtran_ ( int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc );
    void pdnrm2_ ( int *n, double *norm2, double *x, int *ix, int *jx, int *descx, int *incx );
    void pdpotrs_ ( char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *info );
    void pddot_( int *n, double *dot, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void pdcopy_( int *n, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void pdscal_( int *n, double *a, double *x, int *ix, int *jx, int *descx, int *incx );
    void dgsum2d_(int *ConTxt, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rdest, int *cdest);
    void dgebs2d_(int *ConTxt, char *scope, char *top, int *m, int *n, double *A, int *lda);
    void dgebr2d_(int *ConTxt, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc);
    void dgesd2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rdest, int *cdest );
    void dgerv2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc );
}

int cv(double * ets, int *dets) {
    FILE *fT;
    int ni, i,j, info;
    int *DT, *DE;
    double *Tb, v, *Ebr;
    int nTb, ns, pTb, sc, lld_T,pcol,ct, lld_E;

    DT= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DT==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }
    DE= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DE==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }

    pcol= * ( pst+1 );
    ns= nt % ( b * * ( ds+1 ) ) ==0 ?  nt / ( b * * ( ds+1 ) ) : ( nt / ( b * * ( ds+1 ) ) ) +1;
    sc= b * * ( ds+1 ); 
    nTb= m%b==0 ? m/b : m/b +1;
    pTb= ( nTb - *pst ) % *ds == 0 ? ( nTb- *pst ) / *ds : ( nTb- *pst ) / *ds +1;
    pTb= pTb <1? 1:pTb;
    lld_T=pTb*b;
    lld_E=ns*b* *(ds+1);

    descinit_ ( DT, &m, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_T, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }

    descinit_ ( DE, &lld_E, &i_one, &lld_E, &i_one, &i_zero, &i_zero, &ICTXT2D, &lld_E, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of EBV returns info: %d\n",info );
        return info;
    }

    Tb= ( double* ) calloc ( pTb*b*b, sizeof ( double ) );
    if ( Tb==NULL ) {
        printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    
    Ebr = ( double* ) calloc ( lld_E,sizeof ( double ) );
    if ( Ebr==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }

    fT=fopen ( Ts,"rb" );
    if ( fT==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }

    for ( ni=0; ni<ns; ++ni ) {
        if ( ni==ns-1 ) {

            free ( Tb );
            Tb= ( double* ) calloc ( pTb*b*b, sizeof ( double ) );
            if ( Tb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
                return -1;
            }
        }
        if ( ( nTb-1 ) % *ds == *pst && m%b !=0 ) {
            for ( i=0; i<b; ++i ) {
                info=fseek ( fT, ( long ) ( ( ( ni * * ( ds+1 ) * b + pcol * b + i ) * ( m+1 ) + b * *pst ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading file" );
                    return -1;
                }
                if ( *pst==0 )
                    fread ( &v,sizeof ( double ),1,fT );
                else
                    info=fseek ( fT,1L * sizeof ( double ), SEEK_CUR );
                for ( j=0; j < pTb-1; ++j ) {
                    fread ( Tb + i*pTb*b + j*b,sizeof ( double ),b,fT );
                    info=fseek ( fT, ( long ) ( ( ( *ds ) -1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading file" );
                        return -1;
                    }
                }
                fread ( Tb + i*pTb*b + j*b,sizeof ( double ),m%b,fT );
            }
        } else {
            for ( i=0; i<b; ++i ) {
                info=fseek ( fT, ( long ) ( ( ( ni * * ( ds+1 ) * b + pcol * b + i ) * ( m+1 ) + b * *pst ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading file" );
                    return -1;
                }
                if ( *pst==0 )
                    fread ( &v,sizeof ( double ),1,fT );
                else
                    info=fseek ( fT,1L * sizeof ( double ), SEEK_CUR );
                for ( j=0; j < pTb; ++j ) {
                    fread ( Tb + i*pTb*b + j*b,sizeof ( double ),b,fT );
                    info=fseek ( fT, ( long ) ( ( * ( ds )-1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading file" );
                        return -1;
                    }
                }
            }
        }
        blacs_barrier_ ( &ICTXT2D,"A" );

        ct=1 + ni * sc;

        pdgemm_ ( "T","N",&sc,&i_one,&m,&d_one,Tb,&i_one,&i_one,DT,ets,&t_plus,&i_one,dets,&d_one,Ebr,&ct,&i_one,DE);
       
    }

    fclose(fT);
    if(*pst==0 && *(pst+1)==0) {
        printdense(nt,1,Ebr,"EBV.txt" );
    }
    blacs_barrier_(&ICTXT2D, "A" );
    free(Ebr);
    free(Tb);
    free(DE);
    free(DT);

    return info;

}

int cv5(double * ets, int *dets) {
    int ni, i,j, info;
    int *DT, *DE;
    double *tb, v, *ebr;
    int ntb, ns, pTb, sc, lld_T,pcol,curtest, lld_E;

    hid_t       fid, dgi, sgi;
    hid_t	pid, mg;
    herr_t	status;
    hsize_t	dm[2], os[2],co[2], st[2],bl[2];

    int mpinfo  = MPI_INFO_NULL;

    pid = H5Pcreate ( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio ( pid, MPI_COMM_WORLD, mpinfo );

    fid = H5Fopen ( fDn, H5F_ACC_RDWR, pid );
    dgi = H5Dopen ( fid, Ts, H5P_DEFAULT );
    sgi=H5Dget_space(dgi);

    DT= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DT==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }

    DE= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DE==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }

    pcol= * ( pst+1 );
    ns= nt % ( b * * ( ds+1 ) ) ==0 ?  nt / ( b * * ( ds+1 ) ) : ( nt / ( b * * ( ds+1 ) ) ) +1;
    sc= b * * ( ds+1 ); 
    ntb= m%b==0 ? m/b : m/b +1;
    pTb= ( ntb - *pst ) % *ds == 0 ? ( ntb- *pst ) / *ds : ( ntb- *pst ) / *ds +1;
    pTb= pTb <1? 1:pTb;
    lld_T=pTb*b;
    lld_E=ns*b* *(ds+1);

    descinit_ ( DT, &m, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_T, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }

    descinit_ ( DE, &lld_E, &i_one, &lld_E, &i_one, &i_zero, &i_zero, &ICTXT2D, &lld_E, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of EBV returns info: %d\n",info );
        return info;
    }

    tb= ( double* ) calloc ( pTb*b*b, sizeof ( double ) );
    if ( tb==NULL ) {
        printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    dm[0]=b;
    dm[1]=pTb*b;
    mg = H5Screate_simple(2,dm,NULL);

    ebr = ( double* ) calloc ( lld_E,sizeof ( double ) );
    if ( ebr==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }

    pid = H5Pcreate ( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio ( pid, H5FD_MPIO_INDEPENDENT );

    for ( ni=0; ni<ns; ++ni ) {
        if(*pst >= ntb)
            goto CALC;
        if ( ni==ns-1 ) {

            free ( tb );
            tb= ( double* ) calloc ( pTb*b*b, sizeof ( double ) );
            if ( tb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
                return -1;
            }
            if((pcol + 1 + (ns-1) * *(ds+1))*b <= nt)
                bl[0]=b;
            else if ((pcol + (ns-1) * *(ds+1))*b >= nt)
                bl[0]=0;
            else
                bl[0]=nt%b;
        }
        else {
            bl[0]=b;
        }
        if ( ( ntb-1 ) % *ds == *pst && m%b !=0 ) {
            os[0] = ni * *(ds+1) * b + pcol * b;
            os[1] = *pst * b;
            co[0] = 1;
            co[1] = pTb-1;
            st[0] = b * *(ds+1);
            st[1] = b * *ds;
            bl[1] = b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_SET, os, st, co, bl );
            if (status<0) {
                printf("selection of geno hyperslab in file was unsuccesful, strip: %d\n",ni);
                return status;
            }
            os[0] = 0;
            os[1] = 0;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( mg, H5S_SELECT_SET, os, st, co, bl );
            if (status<0) {
                printf("selection of hyperslab in memory was unsuccesful, strip: %d\n",ni);
                return status;
            }

            os[0] = ni * *(ds+1) * b + pcol * b;
            os[1] = (ntb-1) * b;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * *(ds+1);
            st[1] = b * *ds;
            bl[1] = m%b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_OR, os, st, co, bl );
            if (status<0) {
                printf("selection of geno extended hyperslab in file was unsuccesful, strip: %d\n",ni);
                return status;
            }

            os[0] = 0;
            os[1] = (pTb-1) * b;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( mg, H5S_SELECT_OR, os, st, co, bl );
            if (status<0) {
                printf("selection of hyperslab in memory was unsuccesful, strip: %d\n",ni);
                return status;
            }
        }
        else {
            os[0] = ni * *(ds+1) * b + pcol * b;
            os[1] = *pst * b;
            co[0] = 1;
            co[1] = pTb;
            st[0] = b * *(ds+1);
            st[1] = b * *ds;
            bl[1] = b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_SET, os, st, co, bl );
            if (status<0) {
                printf("selection of geno hyperslab in file was unsuccesful\n");
                return status;
            }

            os[0] = 0;
            os[1] = 0;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( mg, H5S_SELECT_SET, os, st, co, bl );
            if (status<0) {
                printf("selection of hyperslab in memory was unsuccesful\n");
                return status;
            }
        }
        status= H5Dread ( dgi,H5T_NATIVE_DOUBLE_g,mg,sgi,pid,tb );
        if (status<0) {
            printf("reading of geno hyperslab was unsuccesful\n");
            return status;
        }
CALC:
        blacs_barrier_ ( &ICTXT2D,"A" );

        curtest=1 + ni * sc;

        pdgemm_ ( "T","N",&sc,&i_one,&m,&d_one,tb,&i_one,&i_one,DT,ets,&t_plus,&i_one,dets,&d_one,ebr,&curtest,&i_one,DE);
        
    }

    if(*pst==0 && *(pst+1)==0) {
        printdense(nt,1,ebr,"EBV.txt" );
    }
    blacs_barrier_(&ICTXT2D, "A" );
    free(ebr);
    free(tb);
    free(DE);
    free(DT);
    
    H5Dclose ( dgi );
    H5Sclose ( mg );
    H5Sclose ( sgi );

    H5Pclose ( pid );

    H5Fclose ( fid );

    return info;

}
