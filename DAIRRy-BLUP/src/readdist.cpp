#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "shared_var.h"
#include <shared_var.h>
#include <limits.h>
#define MPI_SUCCESS          0
#define DLEN_ 		     9
typedef int MPI_Comm;


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
}

int Csu ( int * DC, double * mc, int * DYT, double * yt, double *rn ) {

    FILE *fZ, *fX;
    int ni, i,j, info;
    int *DZ, *DY, *DX;
    double *zb, *Xb, *yb, *nb, *temp;
    int nzb, nxb, ns, pzb, pxb, sc, lld_Z, lld_X, pcol, ccu,rcu;

    DZ= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DZ==NULL ) {
        printf ( "unable to allocate memory for descriptor for Z\n" );
        return -1;
    }
    DY= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DY==NULL ) {
        printf ( "unable to allocate memory for descriptor for Y\n" );
        return -1;
    }
    DX= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DX==NULL ) {
        printf ( "unable to allocate memory for descriptor for Y\n" );
        return -1;
    }

    pcol= * ( pst+1 );
    ns= n % ( b * * ( ds+1 ) ) ==0 ?  n / ( b * * ( ds+1 ) ) : ( n / ( b * * ( ds+1 ) ) ) +1;
    sc= b * * ( ds+1 );
    nzb= m%b==0 ? m/b : m/b +1;
    pzb= ( nzb - *pst ) % *ds == 0 ? ( nzb- *pst ) / *ds : ( nzb- *pst ) / *ds +1;
    pzb= pzb <1? 1:pzb;
    lld_Z=pzb*b;
    nxb= t%b==0 ? t/b : t/b +1;
    pxb= ( nxb - *pst ) % *ds == 0 ? ( nxb- *pst ) / *ds : ( nxb- *pst ) / *ds +1;
    pxb= pxb <1? 1:pxb;
    lld_X=pxb*b;

    descinit_ ( DZ, &m, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_Z, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Z returns info: %d\n",info );
        return info;
    }

    descinit_ ( DY, &i_one, &sc, &i_one, &b, &i_zero, &i_zero, &ICTXT2D, &i_one, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix Y returns info: %d\n",info );
        return info;
    }
    descinit_ ( DX, &t, &sc, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_X, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix X returns info: %d\n",info );
        return info;
    }

    zb= ( double* ) calloc ( pzb*b*b, sizeof ( double ) );
    if ( zb==NULL ) {
        printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }

    yb = ( double* ) calloc ( b,sizeof ( double ) );
    if ( yb==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }
    Xb= ( double* ) calloc ( pxb*b*b, sizeof ( double ) );
    if ( Xb==NULL ) {
        printf ( "Error in allocating memory for a strip of X in processor (%d,%d)",*pst,* ( pst+1 ) );
        return -1;
    }
    nb = ( double* ) calloc ( 1,sizeof ( double ) );
    if ( nb==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }

    temp=mc;
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
    temp=mc;
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


    fZ=fopen ( fDn,"rb" );
    if ( fZ==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }

    fX=fopen ( fXn,"rb" );
    if ( fX==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }
    *rn=0.0;
    *nb=0.0;

    for ( ni=0; ni<ns; ++ni ) {
        if ( ni==ns-1 ) {

            free ( zb );
            free ( yb );
            free ( Xb );
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
            Xb= ( double* ) calloc ( pxb*b*b, sizeof ( double ) );
            if ( Xb==NULL ) {
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)\n",*pst,* ( pst+1 ) );
                return -1;
            }
        }

        if ( ( nzb-1 ) % *ds == *pst && m%b !=0 ) {
            if (ni==0) {
                info=fseek ( fZ, ( long ) ( pcol * b * ( m+1 ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fZ, ( long ) ( b * (*( ds+1 )-1) * ( m+1 ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fZ, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                if ( *pst==0 )
                    fread ( yb + i,sizeof ( double ),1,fZ );
                else
                    info=fseek ( fZ,1L * sizeof ( double ), SEEK_CUR );
                for ( j=0; j < pzb-1; ++j ) {
                    fread ( zb + i*pzb*b + j*b,sizeof ( double ),b,fZ );
                    info=fseek ( fZ, ( long ) ( ( ( *ds ) -1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( zb + i*pzb*b + j*b,sizeof ( double ),m%b,fZ );
            }
        }
        else {
            if (ni==0) {
                info=fseek ( fZ, ( long ) ( pcol * b * ( m+1 ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fZ, ( long ) ( b * (*( ds+1 )-1) * ( m+1 ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fZ, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                if ( *pst==0 )
                    fread ( yb + i,sizeof ( double ),1,fZ );
                else
                    info=fseek ( fZ,1L * sizeof ( double ), SEEK_CUR );
                for ( j=0; j < pzb-1; ++j ) {
                    fread ( zb + i*pzb*b + j*b,sizeof ( double ),b,fZ );
                    info=fseek ( fZ, ( long ) ( ( * ( ds )-1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( zb + i*pzb*b + j*b,sizeof ( double ),b,fZ );
                info=fseek ( fZ, ( long ) ( (m - b * ((pzb-1) * *ds + *pst +1 )) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
        }

        if ( ( nxb-1 ) % *ds == *pst && t%b !=0 ) {
            if (ni==0) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fX, ( long ) ( b * (*( ds+1 )-1) * t * sizeof ( double ) ),SEEK_CUR );
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
                    fread ( Xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                    info=fseek ( fX, ( long ) ( ( ( *ds ) -1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( Xb + i*pxb*b + j*b,sizeof ( double ),t%b,fX );
            }
        } else {
            if (ni==0) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fX, ( long ) ( b * (*( ds+1 )-1) * t * sizeof ( double ) ),SEEK_CUR );
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
                    fread ( Xb + i*pxb*b + j*b,sizeof ( double ),b,fX );
                    info=fseek ( fX, ( long ) ( ( * ( ds )-1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( Xb + i*pxb*b + j*b,sizeof ( double ),b,fX  );
                info=fseek ( fX, ( long ) ( (t - b * ((pxb-1) * *ds + *pst +1 )) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
        }

        blacs_barrier_ ( &ICTXT2D,"A" );

        pdsyrk_ ( "U","N",&m,&sc,&d_one, zb,&i_one, &i_one,DZ, &d_one, mc, &t_plus, &t_plus, DC );

        pdgemm_ ( "N","T",&m,&i_one,&sc,&d_one,zb,&i_one, &i_one, DZ,yb,&i_one,&i_one,DY,&d_one,yt,&t_plus,&i_one,DYT );

        pdsyrk_ ( "U","N",&t,&sc,&d_one, Xb,&i_one, &i_one,DX, &d_one, mc, &i_one, &i_one, DC );

        pdgemm_ ( "N","T",&t,&i_one,&sc,&d_one,Xb,&i_one, &i_one, DX,yb,&i_one,&i_one,DY,&d_one,yt,&i_one,&i_one,DYT );

        pdgemm_ ( "N","T",&t,&m,&sc,&d_one,Xb,&i_one, &i_one, DX,zb,&i_one,&i_one,DZ,&d_one,mc,&i_one,&t_plus,DC );

        pdnrm2_ ( &sc,nb,yb,&i_one,&i_one,DY,&i_one );
        *rn += *nb * *nb;

        blacs_barrier_ ( &ICTXT2D,"A" );
    }

    info=fclose ( fX );
    if ( info!=0 ) {
        printf ( "Error in closing open streams" );
        return -1;
    }
    info=fclose ( fZ );
    if ( info!=0 ) {
        printf ( "Error in closing open streams" );
        return -1;
    }

    free ( DX );
    free ( DY );
    free ( DZ );
    free ( zb );
    free ( Xb );
    free ( yb );
    free ( nb );
    return 0;
}

int Cu ( int * DC, double * mc, double du) {

    int i,j, rcu,ccu,nxb;

    nxb= t%b==0 ? t/b : t/b +1;

    for ( i=0,rcu=0,ccu=0; i<Cb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds )
            rcu=0;
        if ( ccu==* ( ds+1 ) )
            ccu=0;
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            if ( i< ( Cb -1 ) ) {
                for ( j=0; j<b; ++j ) {
                    * ( mc+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j ) +=du;
                }
            } else {
                for ( j=0; j< Cd % b; ++j ) {
                    * ( mc+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j ) +=du;
                }
            }
        }
    }

    for ( i=0,rcu=0,ccu=0; i<nxb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds )
            rcu=0;
        if ( ccu==* ( ds+1 ) )
            ccu=0;
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            if ( i<nxb-1 ) {
                for ( j=0; j<b; ++j ) {
                    * ( mc+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j ) -=du;
                }
            } else {
                for ( j=0; j< t%b; ++j ) {
                    * ( mc+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j ) -=du;
                }
            }
        }
    }
}


int Asu ( double * mai, int * dai,int * dyt, double * yt, int * dc, double * mc, double sa ) {

    FILE *fZ, *fX;
    int ni, i,j, info;
    int *dz, *dy, *dx, *dzu, *dqr, *dqs;
    double *zb, *xb, *yb, *zub, *qr, *qs,*nb, sr;
    int nzb, nxb, nst, pzb, pxb, sc, lld_Z, lld_X, pcol, ccu,rcu;

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
    yb = ( double* ) calloc ( b,sizeof ( double ) );
    if ( yb==NULL ) {
        printf ( "unable to allocate memory for Matrix Y\n" );
        return EXIT_FAILURE;
    }
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

    fZ=fopen ( fDn,"rb" );
    if ( fZ==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }

    fX=fopen ( fXn,"rb" );
    if ( fX==NULL ) {
        printf ( "Error opening file\n" );
        return -1;
    }
    *nb=0.0;

    for ( ni=0; ni<nst; ++ni ) {
        if ( ni==nst-1 ) {
            free ( zb );
            free ( zub );
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
                printf ( "Error in allocating memory for a strip of Z in processor (%d,%d)",*pst,* ( pst+1 ) );
                return -1;
            }
            zub = ( double* ) calloc ( b,sizeof ( double ) );
            if ( zub==NULL ) {
                printf ( "unable to allocate memory for Matrix Y\n" );
                return EXIT_FAILURE;
            }
        }


        if ( ( nzb-1 ) % *ds == *pst && m%b !=0 ) {
            if (ni==0) {
                info=fseek ( fZ, ( long ) ( pcol * b * ( m+1 ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fZ, ( long ) ( b * (*( ds+1 )-1) * ( m+1 ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fZ, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                if ( *pst==0 )
                    fread ( yb + i,sizeof ( double ),1,fZ );
                else
                    info=fseek ( fZ,1L * sizeof ( double ), SEEK_CUR );
                for ( j=0; j < pzb-1; ++j ) {
                    fread ( zb + i*pzb*b + j*b,sizeof ( double ),b,fZ );
                    info=fseek ( fZ, ( long ) ( ( ( *ds ) -1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( zb + i*pzb*b + j*b,sizeof ( double ),m%b,fZ );
            }
        }
        else {
            if (ni==0) {
                info=fseek ( fZ, ( long ) ( pcol * b * ( m+1 ) * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fZ, ( long ) ( b * (*( ds+1 )-1) * ( m+1 ) * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            for ( i=0; i<b; ++i ) {
                info=fseek ( fZ, ( long ) ( b * *pst * sizeof ( double ) ),SEEK_CUR );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
                if ( *pst==0 )
                    fread ( yb + i,sizeof ( double ),1,fZ );
                else
                    info=fseek ( fZ,1L * sizeof ( double ), SEEK_CUR );
                for ( j=0; j < pzb-1; ++j ) {
                    fread ( zb + i*pzb*b + j*b,sizeof ( double ),b,fZ );
                    info=fseek ( fZ, ( long ) ( ( * ( ds )-1 ) * b * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
                fread ( zb + i*pzb*b + j*b,sizeof ( double ),b,fZ );
                if (m>*pst * b) {
                    info=fseek ( fZ, ( long ) ( (m - b * ((pzb-1) * *ds + *pst +1 )) * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading Z file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
            }
        }

        if ( ( nxb-1 ) % *ds == *pst && t%b !=0 ) {
            if (ni==0) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fX, ( long ) ( b * (*( ds+1 )-1) * t * sizeof ( double ) ),SEEK_CUR );
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
            if (ni==0) {
                info=fseek ( fX, ( long ) ( pcol * b *  t * sizeof ( double ) ),SEEK_SET );
                if ( info!=0 ) {
                    printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                    return -1;
                }
            }
            else {
                info=fseek ( fX, ( long ) ( b * (*( ds+1 )-1) * t * sizeof ( double ) ),SEEK_CUR );
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
                fread ( xb + i*pxb*b + j*b,sizeof ( double ),b,fX  );
                if (t>*pst * b) {
                    info=fseek ( fX, ( long ) ( (t - b * ((pxb-1) * *ds + *pst +1 )) * sizeof ( double ) ),SEEK_CUR );
                    if ( info!=0 ) {
                        printf ( "Error in setting correct begin position for reading X file\nprocessor (%d,%d), error: %d \n", *pst,pcol,info );
                        return -1;
                    }
                }
            }
        }
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
    info=fclose ( fZ );
    if ( info!=0 ) {
        printf ( "Error in closing open streams" );
        return -1;
    }
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

double CZt ( double *cm, int * dcm ) {

    double tp;
    int i, j, rcu,ccu, nxb;

    tp=0.0;

    nxb= t%b==0 ? t/b : t/b +1;

    for ( i=0,rcu=0,ccu=0; i<Cb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds )
            rcu=0;
        if ( ccu==* ( ds+1 ) )
            ccu=0;
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            if ( i< ( Cb -1 ) ) {
                for ( j=0; j<b; ++j ) {
                    tp += * ( cm+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j );
                }
            } else {
                for ( j=0; j< Cd % b; ++j ) {
                    tp += * ( cm+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j );
                }
            }
        }
    }

    for ( i=0,rcu=0,ccu=0; i<nxb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds )
            rcu=0;
        if ( ccu==* ( ds+1 ) )
            ccu=0;
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            if ( i<nxb-1 ) {
                for ( j=0; j<b; ++j ) {
                    tp -= * ( cm+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j );
                }
            } else {
                for ( j=0; j<= (t-1)%b; ++j ) {
                    tp -= * ( cm+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j );
                }
            }
        }
    }
    return tp;
}

double Cld ( double *cm, int * dcm ) {

    double ldp;
    int i, j, rcu,ccu, nxb;

    ldp=0.0;

    for ( i=0,rcu=0,ccu=0; i<Cb; ++i, ++ccu, ++rcu ) {
        if ( rcu==*ds )
            rcu=0;
        if ( ccu==* ( ds+1 ) )
            ccu=0;
        if ( *pst==rcu && * ( pst+1 ) == ccu ) {
            if ( i== ( Cb -1 ) && Cd %b !=0 ) {
                for ( j=0; j< (Cd+1) % b; ++j ) {
                    ldp += log(* ( cm+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j ));
                }

            } else {
                for ( j=0; j<b; ++j ) {
                    ldp += log( * ( cm+ ( j + i / * ( ds+1 ) * b ) * lld_C + i / *ds *b +j ));
                }
            }
        }
    }
    return ldp;
}
