#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <string.h>
#include <math.h>
#include <shared_var.h>

extern "C" {
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void blacs_barrier_ ( int*, char* );
    void pdsyrk_ ( char*, char*, int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int* );
    void pdgemm_ ( char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, double *beta, double *c, int *ic, int *jc, int *descc );
    void pdtran_ ( int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *beta, double *c, int *ic, int *jc, int *descc );
    void pdnrm2_ ( int *n, double *norm2, double *x, int *ix, int *jx, int *descx, int *incx );
    void pdpotrs_ ( char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *info );
    void pddot_ ( int *n, double *dot, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void pdcopy_ ( int *n, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void pdscal_ ( int *n, double *a, double *x, int *ix, int *jx, int *descx, int *incx );
    H5_DLL hid_t H5Pcreate ( hid_t cls_id );
}

#define MPI_INFO_NULL         ((MPI_Info)0x1c000000)

int Csu5 ( int * dco, double * mco, int * dy, double * ty, double *rn) {
    FILE *fX;
    int ni, i,j, info;
    int *dz, *dey, *dx;
    double *zb, *xb, *yb, *nb, *temp;
    int nzb, nxb, nst, pzb, pxb, sc, lld_Z, lld_X, pcol, ccu,rcu;

    hid_t       fid, dgi, dpi, sgi;
    hid_t	pid, msg, spi, msp;
    herr_t	status;
    hsize_t	dm[2], os[2],co[2], st[2],bl[2];

    int mpinfo  = MPI_INFO_NULL;

    pid = H5Pcreate ( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio ( pid, MPI_COMM_WORLD, mpinfo );
    if (mpinfo<0) {
        printf("Something went wrong with setting IO options for HDF5-file, error: %d \n",mpinfo);
        return mpinfo;
    }

    fid = H5Fopen ( fDn, H5F_ACC_RDONLY, pid );
    if (fid <0) {
        printf("Something went wrong with opening HDF5-file, error: %d \n",fid);
        return fid;
    }
    dgi = H5Dopen ( fid, Sd, H5P_DEFAULT );
    if (dgi <0) {
        printf("Something went wrong with opening dataset in HDF5-file, error: %d \n",dgi);
        return dgi;
    }
    dpi = H5Dopen ( fid, pd, H5P_DEFAULT );
    if (dpi <0) {
        printf("Something went wrong with opening dataset in HDF5-file, error: %d \n",dpi);
        return dpi;
    }
    sgi=H5Dget_space ( dgi );
    if (sgi <0) {
        printf("Something went wrong with opening dataset in HDF5-file, error: %d \n",sgi);
        return sgi;
    }
    spi=H5Dget_space ( dpi );
    if (spi <0) {
        printf("Something went wrong with opening dataset in HDF5-file, error: %d \n",spi);
        return spi;
    }

    dz= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dz==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }
    dey= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dey==NULL ) {
        printf ( "unable to allocate memory for descriptor for Y\n" );
        return -1;
    }
    dx= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dx==NULL ) {
        printf ( "unable to allocate memory for descriptor for Y\n" );
        return -1;
    }

    pcol= * ( pst+1 );
    nst= n % ( b * * ( ds+1 ) ) ==0 ?  n / ( b * * ( ds+1 ) ) : ( n / ( b * * ( ds+1 ) ) ) +1; 
    sc= b * * ( ds+1 ); 
    nzb= m%b==0 ? m/b : m/b +1;
    pzb= ( nzb - *pst ) % *ds == 0 ? ( nzb- *pst ) / *ds : ( nzb- *pst ) / *ds +1;
    pzb= pzb <1? 1:pzb;
    lld_Z=pzb*b;		
    nxb= t%b==0 ? t/b : t/b +1;	
    pxb= ( nxb - *pst ) % *ds == 0 ? ( nxb- *pst ) / *ds : ( nxb- *pst ) / *ds +1;
    pxb= pxb <1? 1:pxb;
    lld_X=pxb*b;		

    descinit_ ( dz, &m, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_Z, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }

    descinit_ ( dey, &i_one, &sc, &i_one, &b, &i_zero, &i_zero, &ICTXT2D, &i_one, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Y returns info: %d\n",info );
        return info;
    }
    descinit_ ( dx, &t, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_X, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix X returns info: %d\n",info );
        return info;
    }

    zb= ( double* ) calloc ( pzb*b*b, sizeof ( double ) );
    if ( zb==NULL ) {
        printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    dm[0]=b;
    dm[1]=pzb*b;
    msg = H5Screate_simple ( 2,dm,NULL );

    yb = ( double* ) calloc ( b,sizeof ( double ) );
    if ( yb==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }
    dm[0]=b;
    dm[1]=1;
    msp = H5Screate_simple ( 1,dm,NULL );

    xb= ( double* ) calloc ( pxb*b*b, sizeof ( double ) );
    if ( xb==NULL ) {
        printf ( "Error in allocating memory for a strip of X in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    nb = ( double* ) calloc ( 1,sizeof ( double ) );
    if ( nb==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }

    temp=mco;
    for ( i=0,rcu=0,ccu=0; i<Cb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds ) {
            rcu=0;
            temp += b;
        }
        if ( ccu==* ( ds+1 ) ) {
            ccu=0;
            temp += b*lld_C;
        }
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            for ( j=0; j<b; ++j ) {
                * ( temp + j  * lld_C +j ) =l;
            }
            if ( i==Cb-1 && Cd % b != 0 ) {
                for ( j=b-1; j>= Cd % b; --j ) {
                    * ( temp + j * lld_C + j ) =0.0;
                }
            }
        }

    }
    temp=mco;
    for ( i=0,rcu=0,ccu=0; i<nxb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds ) {
            rcu=0;
            temp += b;
        }
        if ( ccu==* ( ds+1 ) ) {
            ccu=0;
            temp += b*lld_C;
        }
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            if ( i<nxb-1 ) {
                for ( j=0; j<b; ++j ) {
                    * ( temp + j * lld_C + j ) =0.0;
                }
            } else {
                for ( j=0; j<= ( t-1 ) %b; ++j ) {
                    * ( temp + j * lld_C + j ) =0.0;
                }
            }
        }
    }

    fX=fopen ( fXn,"rb" );
    if ( fX==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }
    *rn=0.0;
    *nb=0.0;

    pid = H5Pcreate ( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio ( pid, H5FD_MPIO_INDEPENDENT );  

    for ( ni=0; ni<nst; ++ni ) {
        if ( *pst >= nzb )
            goto CALC;
        if ( ni==nst-1 ) {

            free ( zb );
            free ( yb );
            free ( xb );
            zb= ( double* ) calloc ( pzb*b*b, sizeof ( double ) );
            if ( zb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)\n",*pst,* ( pst+1 ) );
                return -1;
            }
            yb = ( double* ) calloc ( b,sizeof ( double ) );
            if ( yb==NULL ) {
                printf ( "unable to allocate memory for Matrix Y\n" );
                return EXIT_FAILURE;
            }
            xb= ( double* ) calloc ( pxb*b*b, sizeof ( double ) );
            if ( xb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)\n",*pst,* ( pst+1 ) );
                return -1;
            }

            if ( ( pcol + 1 + ( nst-1 ) * * ( ds+1 ) ) *b <= n )
                bl[0]=b;
            else if ( ( pcol + ( nst-1 ) * * ( ds+1 ) ) *b >= n )
                bl[0]=0;
            else
                bl[0]=n%b;
        } else {
            bl[0]=b;
        }
        if ( ( nzb-1 ) % *ds == *pst && m%b !=0 ) {
            os[0] = ni * * ( ds+1 ) * b + pcol * b;
            os[1] = *pst * b;
            co[0] = 1;
            co[1] = pzb-1;
            st[0] = b * * ( ds+1 );
            st[1] = b * *ds;
            bl[1] = b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of geno hyperslab in file was unsuccesful, strip: %d\n",ni );
                return status;
            }
            os[0] = 0;
            os[1] = 0;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( msg, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of hyperslab in memory was unsuccesful, strip: %d\n",ni );
                return status;
            }

            os[0] = ni * * ( ds+1 ) * b + pcol * b;
            os[1] = ( nzb-1 ) * b;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * * ( ds+1 );
            st[1] = b * *ds;
            bl[1] = m%b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_OR, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of geno extended hyperslab in file was unsuccesful, strip: %d\n",ni );
                return status;
            }

            os[0] = 0;
            os[1] = ( pzb-1 ) * b;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( msg, H5S_SELECT_OR, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of hyperslab in memory was unsuccesful, strip: %d\n",ni );
                return status;
            }
        } else {
            os[0] = ni * * ( ds+1 ) * b + pcol * b;
            os[1] = *pst * b;
            co[0] = 1;
            co[1] = pzb;
            st[0] = b * * ( ds+1 );
            st[1] = b * *ds;
            bl[1] = b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of geno hyperslab in file was unsuccesful\n" );
                return status;
            }

            os[0] = 0;
            os[1] = 0;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( msg, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of hyperslab in memory was unsuccesful\n" );
                return status;
            }
        }
        status= H5Dread ( dgi,H5T_NATIVE_DOUBLE_g,msg,sgi,pid,zb );
        if ( status<0 ) {
            printf ( "reading of geno hyperslab was unsuccesful\n" );
            return status;
        }
        if ( *pst==0 ) {

            os[0] = ni * b * * ( ds+1 ) + pcol * b;
            os[1] = 0;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * *ds;
            st[1] = 1;
            bl[1] = 1;

            status = H5Sselect_hyperslab ( spi, H5S_SELECT_SET, os, st, co,bl );
            if ( status<0 ) {
                printf ( "selection of pheno hyperslab in file was unsuccesful\n" );
                return -1;
            }
            os[0] = 0;
            os[1] = 0;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * *ds;
            st[1] = 1;
            bl[1] = 1;

            status = H5Sselect_hyperslab ( msp, H5S_SELECT_SET, os, st, co,bl );
            if ( status<0 ) {
                printf ( "selection of pheno hyperslab in file was unsuccesful\n" );
                return -1;
            }

            status=H5Dread ( dpi,H5T_NATIVE_DOUBLE_g,msp,spi,pid,yb );
            if ( status<0 ) {
                printf ( "reading of pheno hyperslab was unsuccesful\n" );
                return -1;
            }

        }

        if ( ( nxb-1 ) % *ds == *pst && t%b !=0 ) {
            if ( ni==0 ) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            } else {
                info=fseek ( fX, ( long ) ( b * ( * ( ds+1 )-1 ) * t * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fX, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                for ( j=0; j < pxb-1; ++j ) {
                    fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                    info=fseek ( fX, ( long ) ( ( ( *ds ) -1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( xb + i*pxb*b + j*b,sizeof ( double ),t%b,fX );
            }
        } else {
            if ( ni==0 ) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            } else {
                info=fseek ( fX, ( long ) ( b * ( * ( ds+1 )-1 ) * t * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fX, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                for ( j=0; j < pxb-1; ++j ) {
                    fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                    info=fseek ( fX, ( long ) ( ( * ( ds )-1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                if ( t>*pst * b ) {
                    info=fseek ( fX, ( long ) ( ( t - b * ( ( pxb-1 ) * *ds + *pst +1 ) ) * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
            }
        }

CALC:
        blacs_barrier_ ( &ICTXT2D,"A" );

        pdsyrk_ ( "U","N",&m,&sc,&d_one, zb,&i_one, &i_one,dz, &d_one, mco, &t_plus, &t_plus, dco ); 

        pdgemm_ ( "N","T",&m,&i_one,&sc,&d_one,zb,&i_one, &i_one, dz,yb,&i_one,&i_one,dey,&d_one,ty,&t_plus,&i_one,dey ); 

        pdsyrk_ ( "U","N",&t,&sc,&d_one, xb,&i_one, &i_one,dx, &d_one, mco, &i_one, &i_one, dco ); 

        pdgemm_ ( "N","T",&t,&i_one,&sc,&d_one,xb,&i_one, &i_one, dx,yb,&i_one,&i_one,dey,&d_one,ty,&i_one,&i_one,dey ); 

        pdgemm_ ( "N","T",&t,&m,&sc,&d_one,xb,&i_one, &i_one, dx,zb,&i_one,&i_one,dz,&d_one,mco,&i_one,&t_plus,dco ); 

        pdnrm2_ ( &sc,nb,yb,&i_one,&i_one,dey,&i_one );
        *rn += *nb * *nb;

        blacs_barrier_ ( &ICTXT2D,"A" );
    }

    info=fclose ( fX );
    if ( info!=0 ) {
        printf ( "Error in closing open streams" );
        return -1;
    }

    H5Dclose ( dgi );
    H5Dclose ( dpi );
    H5Sclose ( msg );
    H5Sclose ( msp );
    H5Sclose ( sgi );
    H5Sclose ( spi );

    H5Pclose ( pid );

    H5Fclose ( fid );
    free ( dx );
    free ( dey );
    free ( dz );
    free ( zb );
    free ( xb );
    free ( yb );
    free ( nb );
    return 0;

}

int Asu5 ( double * mai, int * dai,int * dyt, double * yt, int * dc, double * mc, double sa) {

    FILE *fX;
    int ni, i,j, info;
    int *dz, *dy, *dx, *dzu, *dqr, *dqs;
    double *zb, *xb, *yb, *zub, *qr, *qs,*nb, sr;
    int nzb, nxb, nst, pzb, pxb, sc, lld_Z, lld_X, pcol, ccu,rcu;

    hid_t       fid, dgi, dpi, sgi;    
    hid_t	pid, msg, spi, msp;    
    herr_t	status;
    hsize_t	dm[2], os[2],co[2], st[2],bl[2];

    MPI_Info mpinfo  = MPI_INFO_NULL;

    pid = H5Pcreate ( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio ( pid, MPI_COMM_WORLD, mpinfo );

    fid = H5Fopen ( fDn, H5F_ACC_RDWR, pid );
    dgi = H5Dopen ( fid, Sd, H5P_DEFAULT );
    dpi = H5Dopen ( fid, pd, H5P_DEFAULT );
    sgi=H5Dget_space ( dgi );
    spi=H5Dget_space ( dpi );

    dz= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dz==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }
    dy= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dy==NULL ) {
        printf ( "unable to allocate memory for descriptor for Y\n" );
        return -1;
    }
    dx= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dx==NULL ) {
        printf ( "unable to allocate memory for descriptor for X\n" );
        return -1;
    }
    dzu= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dzu==NULL ) {
        printf ( "unable to allocate memory for descriptor for Zu\n" );
        return -1;
    }
    dqr= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dqr==NULL ) {
        printf ( "unable to allocate memory for descriptor for QRHS\n" );
        return -1;
    }
    dqs= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( dqs==NULL ) {
        printf ( "unable to allocate memory for descriptor for QSOL\n" );
        return -1;
    }

    pcol= * ( pst+1 );
    nst= n % ( b * * ( ds+1 ) ) ==0 ?  n / ( b * * ( ds+1 ) ) : ( n / ( b * * ( ds+1 ) ) ) +1; 
    sc= b * * ( ds+1 ); 
    nzb= m%b==0 ? m/b : m/b +1;	
    pzb= ( nzb - *pst ) % *ds == 0 ? ( nzb- *pst ) / *ds : ( nzb- *pst ) / *ds +1;
    pzb= pzb <1? 1:pzb;
    lld_Z=pzb*b;	
    nxb= t%b==0 ? t/b : t/b +1;	
    pxb= ( nxb - *pst ) % *ds == 0 ? ( nxb- *pst ) / *ds : ( nxb- *pst ) / *ds +1;
    pxb= pxb <1? 1:pxb;
    lld_X=pxb*b;		
    sr=1/sa;

    descinit_ ( dz, &m, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_Z, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }
    descinit_ ( dy, &i_one, &sc, &i_one, &b, &i_zero, &i_zero, &ICTXT2D, &i_one, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Y returns info: %d\n",info );
        return info;
    }
    descinit_ ( dx, &t, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_X, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix X returns info: %d\n",info );
        return info;
    }
    descinit_ ( dzu, &i_one, &sc, &i_one, &b, &i_zero, &i_zero, &ICTXT2D, &i_one, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }
    descinit_ ( dqr, &Cd, &i_two, &b, &i_two, &i_zero, &i_zero, &ICTXT2D, &lld_C, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }
    descinit_ ( dqs, &Cd, &i_two, &b, &i_two, &i_zero, &i_zero, &ICTXT2D, &lld_C, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }

    zb= ( double* ) calloc ( pzb*b*b, sizeof ( double ) );
    if ( zb==NULL ) {
        printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    dm[0]=b;
    dm[1]=pzb*b;
    msg = H5Screate_simple ( 2,dm,NULL );

    yb = ( double* ) calloc ( b,sizeof ( double ) );
    if ( yb==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }
    dm[0]=b;
    dm[1]=1;
    msp = H5Screate_simple ( 1,dm,NULL );

    zub = ( double* ) calloc ( b,sizeof ( double ) );
    if ( zub==NULL ) {
        printf ( "unable to allocate memory for Matrix Zu\n" );
        return EXIT_FAILURE;
    }
    xb= ( double* ) calloc ( pxb*b*b, sizeof ( double ) );
    if ( xb==NULL ) {
        printf ( "Error in allocating memory for a strip of X in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    qr= ( double * ) calloc ( Cr * b * 2,sizeof ( double ) );
    if ( qr==NULL ) {
        printf ( "Error in allocating memory for QRHS in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    qs= ( double * ) calloc ( Cr * b * 2,sizeof ( double ) );
    if ( qs==NULL ) {
        printf ( "Error in allocating memory for QRHS in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    nb = ( double* ) calloc ( 1,sizeof ( double ) );
    if ( nb==NULL ) {
        printf ( "unable to allocate memory for norm\n" );
        return EXIT_FAILURE;
    }

    fX=fopen ( fXn,"rb" );
    if ( fX==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }
    *nb=0.0;

    pid = H5Pcreate ( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio ( pid, H5FD_MPIO_INDEPENDENT );

    for ( ni=0; ni<nst; ++ni ) {

        if ( *pst >= nzb )
            goto CALC;
        if ( ni==nst-1 ) {

            free ( zb );
            free ( yb );
            free ( xb );
            zb= ( double* ) calloc ( pzb*b*b, sizeof ( double ) );
            if ( zb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)\n",*pst,* ( pst+1 ) );
                return -1;
            }
            yb = ( double* ) calloc ( b,sizeof ( double ) );
            if ( yb==NULL ) {
                printf ( "unable to allocate memory for Matrix Y\n" );
                return EXIT_FAILURE;
            }
            xb= ( double* ) calloc ( pxb*b*b, sizeof ( double ) );
            if ( xb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)\n",*pst,* ( pst+1 ) );
                return -1;
            }

            if ( ( pcol + 1 + ( nst-1 ) * * ( ds+1 ) ) *b <= n )
                bl[0]=b;
            else if ( ( pcol + ( nst-1 ) * * ( ds+1 ) ) *b >= n )
                bl[0]=0;
            else
                bl[0]=n%b;
        } else {
            bl[0]=b;
        }
        if ( ( nzb-1 ) % *ds == *pst && m%b !=0 ) {
            os[0] = ni * * ( ds+1 ) * b + pcol * b;
            os[1] = *pst * b;
            co[0] = 1;
            co[1] = pzb-1;
            st[0] = b * * ( ds+1 );
            st[1] = b * *ds;
            bl[1] = b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of geno hyperslab in file was unsuccesful, strip: %d\n",ni );
                return status;
            }
            os[0] = 0;
            os[1] = 0;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( msg, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of hyperslab in memory was unsuccesful, strip: %d\n",ni );
                return status;
            }

            os[0] = ni * * ( ds+1 ) * b + pcol * b;
            os[1] = ( nzb-1 ) * b;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * * ( ds+1 );
            st[1] = b * *ds;
            bl[1] = m%b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_OR, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of geno extended hyperslab in file was unsuccesful, strip: %d\n",ni );
                return status;
            }

            os[0] = 0;
            os[1] = ( pzb-1 ) * b;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( msg, H5S_SELECT_OR, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of hyperslab in memory was unsuccesful, strip: %d\n",ni );
                return status;
            }
        } else {
            os[0] = ni * * ( ds+1 ) * b + pcol * b;
            os[1] = *pst * b;
            co[0] = 1;
            co[1] = pzb;
            st[0] = b * * ( ds+1 );
            st[1] = b * *ds;
            bl[1] = b;

            status = H5Sselect_hyperslab ( sgi, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of geno hyperslab in file was unsuccesful\n" );
                return status;
            }

            os[0] = 0;
            os[1] = 0;
            st[0] = b;
            st[1] = b;

            status = H5Sselect_hyperslab ( msg, H5S_SELECT_SET, os, st, co, bl );
            if ( status<0 ) {
                printf ( "selection of hyperslab in memory was unsuccesful\n" );
                return status;
            }
        }
        status= H5Dread ( dgi,H5T_NATIVE_DOUBLE_g,msg,sgi,pid,zb );
        if ( status<0 ) {
            printf ( "reading of geno hyperslab was unsuccesful\n" );
            return status;
        }
        if ( *pst==0 ) {

            os[0] = ni * b * * ( ds+1 ) + pcol * b;
            os[1] = 0;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * *ds;
            st[1] = 1;
            bl[1] = 1;

            status = H5Sselect_hyperslab ( spi, H5S_SELECT_SET, os, st, co,bl );
            if ( status<0 ) {
                printf ( "selection of pheno hyperslab in file was unsuccesful\n" );
                return -1;
            }
            os[0] = 0;
            os[1] = 0;
            co[0] = 1;
            co[1] = 1;
            st[0] = b * *ds;
            st[1] = 1;
            bl[1] = 1;

            status = H5Sselect_hyperslab ( msp, H5S_SELECT_SET, os, st, co,bl );
            if ( status<0 ) {
                printf ( "selection of pheno hyperslab in file was unsuccesful\n" );
                return -1;
            }

            status=H5Dread ( dpi,H5T_NATIVE_DOUBLE_g,msp,spi,pid,yb );
            if ( status<0 ) {
                printf ( "reading of pheno hyperslab was unsuccesful\n" );
                return -1;
            }

        }

        if ( ( nxb-1 ) % *ds == *pst && t%b !=0 ) {
            if ( ni==0 ) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            } else {
                info=fseek ( fX, ( long ) ( b * ( * ( ds+1 )-1 ) * t * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fX, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                for ( j=0; j < pxb-1; ++j ) {
                    fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                    info=fseek ( fX, ( long ) ( ( ( *ds ) -1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( xb + i*pxb*b + j*b,sizeof ( double ),t%b,fX );
            }
        } else {
            if ( ni==0 ) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            } else {
                info=fseek ( fX, ( long ) ( b * ( * ( ds+1 )-1 ) * t * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fX, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                for ( j=0; j < pxb-1; ++j ) {
                    fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                    info=fseek ( fX, ( long ) ( ( * ( ds )-1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                if ( t>*pst * b ) {
                    info=fseek ( fX, ( long ) ( ( t - b * ( ( pxb-1 ) * *ds + *pst +1 ) ) * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
            }
        }
CALC:
        blacs_barrier_ ( &ICTXT2D,"A" );

        pdgemm_ ( "T","N", &i_one, &sc,&m,&l, yt, &t_plus,&i_one,dyt,zb,&i_one,&i_one,dz,&d_zero,zub,&i_one,&i_one,dzu ); 

        pdgemm_ ( "N","T",&m,&i_one,&sc,&sr,zb,&i_one, &i_one, dz,yb,&i_one,&i_one,dy,&d_one,qr,&t_plus,&i_one,dqr ); 

        pdgemm_ ( "N","T",&t,&i_one,&sc,&sr,xb,&i_one, &i_one, dx,yb,&i_one,&i_one,dy,&d_one,qr,&i_one,&i_one,dqr ); 

        pdgemm_ ( "N","T",&m,&i_one,&sc,&d_one,zb,&i_one, &i_one, dz,zub,&i_one,&i_one,dzu,&d_one,qr,&t_plus,&i_two,dqr ); 

        pdgemm_ ( "N","T",&t,&i_one,&sc,&d_one,xb,&i_one, &i_one, dx,zub,&i_one,&i_one,dzu,&d_one,qr,&i_one,&i_two,dqr ); 

        pdnrm2_ ( &sc,nb,yb,&i_one,&i_one,dy,&i_one );
        *mai += *nb * *nb/sa/sa; 
        pdnrm2_ ( &sc,nb,zub,&i_one,&i_one,dzu,&i_one );
        * ( mai + 3 ) += *nb * *nb;
        pddot_ ( &sc,nb,zub,&i_one,&i_one,dzu,&i_one,yb,&i_one,&i_one,dy,&i_one );
        * ( mai + 1 ) += *nb /sa;
        * ( mai + 2 ) += *nb /sa;
        blacs_barrier_ ( &ICTXT2D,"A" );
    }

    pdcopy_ ( &Cd,qr,&i_one,&i_two,dqr,&i_one,qs,&i_one,&i_two,dqs,&i_one );
    pdcopy_ ( &Cd,yt,&i_one,&i_one,dyt,&i_one,qs,&i_one,&i_one,dqs,&i_one );
    pdscal_ ( &Cd,&sr,qs,&i_one,&i_one,dqs,&i_one );
    pdpotrs_ ( "U",&Cd,&i_one,mc,&i_one,&i_one,dc,qs,&i_one,&i_two,dqs,&info );
    if ( info!=0 )
        printf ( "Parallel Cholesky solution for Q was unsuccesful, error returned: %d\n",info );

    pdgemm_ ( "T","N",&i_two,&i_two,&Cd,&d_negone,qr,&i_one,&i_one,dqr,qs,&i_one,&i_one,dqs,&d_one, mai,&i_one,&i_one,dai );

    for ( i=0; i<4; ++i )
        * ( mai + i ) = * ( mai + i ) / 2 / sa;

    info=fclose ( fX );
    if ( info!=0 ) {
        printf ( "Error in closing open streams" );
        return -1;
    }

    H5Dclose ( dgi );
    H5Dclose ( dpi );
    H5Sclose ( msg );
    H5Sclose ( msp );
    H5Sclose ( sgi );
    H5Sclose ( spi );
    H5Pclose ( pid );
    H5Fclose ( fid );

    free ( dqr );
    free ( dqs );
    free ( dx );
    free ( dy );
    free ( dz );
    free ( dzu );
    free ( zb );
    free ( xb );
    free ( yb );
    free ( nb );
    free ( qr );
    free ( qs );
    free ( zub );

    return 0;
}
