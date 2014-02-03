#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "src/shared_var.h"
#include <shared_var.h>
#define MPI_SUCCESS          0      
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)


typedef int MPI_Comm;
typedef int MPI_Datatype;

double d_one=1.0, d_zero=0.0, d_negone=-1.0;
int DLEN_=9, i_negone=-1, i_zero=0, i_one=1, i_two=2, i_four=4; 
int m,n,t, b, m_plus, t_plus; 
int lld_C, Cb, Cd;
long Cr, Cc;
int size, *ds, *pst, ICTXT2D, iam;
int nt, mi,dh5, Ccop;
char *Sd, *pd;
char *fDn, *fXn, *Ts;
double l, e;

extern "C" {
    void blacs_pinfo_ ( int *mypnum, int *nprocs );
    void blacs_setup_ ( int *mypnum, int *nprocs );
    void blacs_get_ ( int *ConTxt, int *what, int *val );
    void blacs_gridinit_ ( int *ConTxt, char *order, int *nprow, int *npcol );
    void blacs_gridexit_ ( int *ConTxt );
    void blacs_pcoord_ ( int *ConTxt, int *nodenum, int *prow, int *pcol );
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void blacs_barrier_ ( int*, char* );
    void igebs2d_(int *ConTxt, char *scope, char *top, int *m, int *n, int *A, int *lda);
    void igebr2d_(int *ConTxt, char *scope, char *top, int *m, int *n, int *A, int *lda, int *rsrc, int *csrc);
    void dgebs2d_(int *ConTxt, char *scope, char *top, int *m, int *n, double *A, int *lda);
    void dgebr2d_(int *ConTxt, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc);
    void dgsum2d_(int *ConTxt, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rdest, int *cdest);
    void dgesd2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rdest, int *cdest );
    void dgerv2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc );
    void pdcopy_ ( int *n, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void pddot_ ( int *n, double *dot, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    double pdlansy_ ( char *norm, char *uplo, int *n, double *a, int *ia, int *ja, int *desca, double *work );
    void pdlacpy_ (char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
    void pdpotrf_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdpotrs_ ( char *uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *info );
    void pdpotri_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info);
    void pdnrm2_( int *n, double *norm2, double *x, int *ix, int *jx, int *descx, int *incx );
    void dsyrk(const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *beta, double *c, const int *ldc);
    void dgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
    double  dnrm2(const int *n, const double *x, const int *incx);
    void    dcopy(const int *n, const double *x, const int *incx, double *y, const int *incy);
    double  ddot (const int *n, const double *x, const int *incx, const double *y, const int *incy);
    void dpotrf_( const char* uplo, const int* n, double* a, const int* lda, int* info );
    void dpotrs_( const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda, double* b, const int* ldb, int* info );
    void dpotri_( const char* uplo, const int* n, double* a, const int* lda, int* info );
    int MPI_Init ( int *, char *** );
    int MPI_Dims_create ( int, int, int * );
    int MPI_Finalize ( void );
    int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm );
}

int main ( int argc, char **argv ) {
    int info;
    info = MPI_Init ( &argc, &argv );
    if ( info != MPI_SUCCESS ) {
        printf ( "Error in MPI initialisation: %d",info );
        return info;
    }
    int i,j,pcol, co, bv;
    double *mc, *ytot, *RHS, *rn, *ran, *mai, *ets, *mcc;
    double si, dot, tp, ccrit, llh,pll,ull;
    int *DESCC, *DESCYTOT, *DESCRHS, *DESCAI, *DESCCCOPY;
    double c0, c1, c2, c3, c4;
    struct timeval tz0,tz1, tz2,tz3;
    double vm_usage, resident_set, cpu_sys, cpu_user;
    double *work, normC,norm1C, norminv, norm1inv, Cmax, colmax;

    DESCC= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCC==NULL ) {
        printf ( "unable to allocate memory for descriptor for C\n" );
        return -1;
    }
    DESCYTOT= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCYTOT==NULL ) {
        printf ( "unable to allocate memory for descriptor for Ytot\n" );
        return -1;
    }
    DESCRHS= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCRHS==NULL ) {
        printf ( "unable to allocate memory for descriptor for RHS\n" );
        return -1;
    }
    DESCAI= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
    if ( DESCAI==NULL ) {
        printf ( "unable to allocate memory for descriptor for AI\n" );
        return -1;
    }

    pst= ( int* ) calloc ( 2,sizeof ( int ) );
    if ( pst==NULL ) {
        printf ( "unable to allocate memory for processor position coordinate\n" );
        return EXIT_FAILURE;
    }
    ds= ( int* ) calloc ( 2,sizeof ( int ) );
    if ( ds==NULL ) {
        printf ( "unable to allocate memory for grid dimensions coordinate\n" );
        return EXIT_FAILURE;
    }

    blacs_pinfo_ ( &iam,&size );	
    blacs_setup_ ( &iam,&size );
    if ( iam ==-1 ) {
        printf ( "Error in initialisation of proces grid" );
        return -1;
    }
    info=MPI_Dims_create ( size, 2, ds );
    if ( info != MPI_SUCCESS ) {
        printf ( "Error in MPI creation of dimensions: %d",info );
        return info;
    }
    blacs_get_ ( &i_negone,&i_zero,&ICTXT2D );

    blacs_gridinit_ ( &ICTXT2D,"R",ds, ds+1 );

    blacs_pcoord_ ( &ICTXT2D,&iam,pst, pst+1 );
    if ( *pst ==-1 ) {
        printf ( "Error in proces grid" );
        return -1;
    }

    if ( argc !=2 ) {
        if ( * ( pst+1 ) ==0 && *pst==0 ) {

            printf ( "The correct use of DAIRRy-BLUP is:\n ./DAIRRy-BLUP <input_file>\n" );
            return -1;
        }
        else
            return -1;
    }
    
    info=ri(*++argv);

    if(info!=0) {
        printf("Something went wrong when reading input file for processor %d\n",iam);
        return -1;
    }
    blacs_barrier_(&ICTXT2D,"ALL");
    if ( * ( pst+1 ) ==0 && *pst==0 )
        printf("Reading of input-file succesful\n");

    m_plus=m+1;
    t_plus=t+1;

    Cd=m+t;
    pcol= * ( pst+1 );
    Cb= Cd%b==0 ? Cd/b : Cd/b +1;
    Cr= ( Cb - *pst ) % *ds == 0 ? ( Cb- *pst ) / *ds : ( Cb- *pst ) / *ds +1;
    Cr= Cr<1? 1 : Cr;
    Cc= ( Cb - pcol ) % * ( ds+1 ) == 0 ? ( Cb- pcol ) / * ( ds+1 ) : ( Cb- pcol ) / * ( ds+1 ) +1;
    Cc=Cc<1? 1 : Cc;
    lld_C=Cr*b;

    if(Ccop) {
        DESCCCOPY= ( int* ) malloc ( DLEN_ * sizeof ( int ) );
        if ( DESCCCOPY==NULL ) {
            printf ( "unable to allocate memory for descriptor for copy of C\n" );
            return -1;
        }
        descinit_ ( DESCCCOPY, &Cd, &Cd, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_C, &info );
        if ( info!=0 ) {
            printf ( "Descriptor of copy of matrix C returns info: %d\n",info );
            return info;
        }
        mcc= ( double* ) calloc ( Cr * (long) b * Cc * (long) b,sizeof ( double ) );
        if ( mcc==NULL ) {
            printf ( "unable to allocate memory for copy of Matrix C (required: %dl bytes)\n", Cr * (long) b * Cc * (long) b );
            return EXIT_FAILURE;
        }
    }

    descinit_ ( DESCC, &Cd, &Cd, &b, &b, &i_zero, &i_zero, &ICTXT2D, &lld_C, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of matrix C returns info: %d\n",info );
        return info;
    }
    descinit_ ( DESCYTOT, &Cd, &i_one, &b, &i_one, &i_zero, &i_zero, &ICTXT2D, &lld_C, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of response matrix returns info: %d\n",info );
        return info;
    }
    descinit_ ( DESCRHS, &Cd, &i_one, &b, &i_one, &i_zero, &i_zero, &ICTXT2D, &lld_C, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of RHS matrix returns info: %d\n",info );
        return info;
    }
    descinit_ ( DESCAI, &i_two, &i_two, &i_two, &i_two, &i_zero, &i_zero, &ICTXT2D, &i_two, &info );
    if ( info!=0 ) {
        printf ( "Descriptor of AI matrix returns info: %d\n",info );
        return info;
    }

    ccrit=0;
    co=0;

    if ( * ( pst+1 ) ==0 && *pst==0 ) {
        printf("\nA linear mixed model with %d observations, %d random effects and %d fixed effects\n", n,m,t);
        printf("was analyzed using %d (%d x %d) processors\n",size,*ds,*(ds+1));
        gettimeofday ( &tz2,NULL );
        c2= tz2.tv_sec*1000000 + ( tz2.tv_usec );
        ets=(double *) calloc(Cb * b, sizeof(double));
        if(ets==NULL) {
            printf("unable to allocate memory for solution matrix\n");
            return EXIT_FAILURE;
        }
    }

    mc= ( double* ) calloc ( Cr * (long) b * Cc * (long) b,sizeof ( double ) );
    if ( mc==NULL ) {
        printf ( "unable to allocate memory for Matrix C  (required: %dl bytes)\n", Cr * (long) b * Cc * (long) b * sizeof ( double ) );
        return EXIT_FAILURE;
    }
    ytot = ( double* ) calloc ( Cr * b,sizeof ( double ) );
    if ( ytot==NULL ) {
        printf ( "unable to allocate memory for Matrix Y (required: %d bytes)\n", Cr * b*sizeof ( double )  );
        return EXIT_FAILURE;
    }
    RHS = ( double* ) calloc ( Cr * b,sizeof ( double ) );
    if ( RHS==NULL ) {
        printf ( "unable to allocate memory for RHS (required: %d bytes)\n", Cr * b * sizeof ( double )  );
        return EXIT_FAILURE;
    }
    mai = ( double* ) calloc ( 2*2,sizeof ( double ) );
    if ( RHS==NULL ) {
        printf ( "unable to allocate memory for AI matrix (required: %d bytes)\n",2*2*sizeof ( double ) );
        return EXIT_FAILURE;
    }
    rn= ( double * ) calloc ( 1,sizeof ( double ) );
    if ( rn==NULL ) {
        printf ( "unable to allocate memory for norm\n" );
        return EXIT_FAILURE;
    }
    ran= ( double * ) calloc ( 1,sizeof ( double ) );
    if ( ran==NULL ) {
        printf ( "unable to allocate memory for norm\n" );
        return EXIT_FAILURE;
    }

    while ( fabs ( ccrit ) >e || fabs(ull/llh) > e || co<2 ) {
        ++co;

        if (co > mi) {
            if ( * ( pst+1 ) ==0 && *pst==0 ) {
                printf("maximum number of iterations reached, AI-REML has not converged\n");
            }
            break;
        }

        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            printf ( "\nParallel results: loop %d\n",co );
            printf ( "=========================\n" );
            gettimeofday ( &tz3,NULL );
            c3= tz3.tv_sec*1000000 + ( tz3.tv_usec );
        }

        if (co > 1 ) {
            free ( ytot );
            //free ( Cmat );
            free ( mai );

            mc= ( double* ) calloc ( Cr*b*Cc*b,sizeof ( double ) );
            if ( mc==NULL ) {
                printf ( "unable to allocate memory for Matrix C (required: %d bytes)\n", Cr*b*Cc*b*sizeof ( double ) );
                return EXIT_FAILURE;
            }
            ytot = ( double* ) calloc ( Cr * b,sizeof ( double ) );
            if ( ytot==NULL ) {
                printf ( "unable to allocate memory for Matrix Y\n" );
                return EXIT_FAILURE;
            }
            mai = ( double* ) calloc ( 2*2,sizeof ( double ) );
            if ( RHS==NULL ) {
                printf ( "unable to allocate memory for AI matrix\n" );
                return EXIT_FAILURE;
            }
            blacs_barrier_ ( &ICTXT2D,"A" );

            if ( * ( pst+1 ) ==0 && *pst==0 ) {
                gettimeofday ( &tz1,NULL );
                c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
                printf ( "\t elapsed wall time allocation of memory:		%10.3f s\n", ( c1 - c3 ) /1000000.0 );
            }

            if (Ccop) {
		Cu(DESCCCOPY,mcc,l/ ( 1+ccrit ) - l);
                pdlacpy_ ( "U", &Cd, &Cd, mcc, &i_one, &i_one, DESCCCOPY, mc, &i_one, &i_one, DESCC);
                pdcopy_ ( &Cd, RHS,&i_one,&i_one,DESCRHS,&i_one,ytot,&i_one,&i_one,DESCYTOT,&i_one );
                if ( * ( pst+1 ) ==0 && *pst==0 ) {
                    gettimeofday ( &tz1,NULL );
                    c0= tz1.tv_sec*1000000 + ( tz1.tv_usec );
                    printf ( "\t elapsed wall time copy of Y and C:			%10.3f s\n", ( c0 - c1 ) /1000000.0 );
                }
                l=l/ ( 1+ccrit ); 
            }
            else {
	      free ( rn );
	      rn= ( double * ) calloc ( 1,sizeof ( double ) );
	      if ( rn==NULL ) {
                printf ( "unable to allocate memory for norm\n" );
                return EXIT_FAILURE;
            }
	      l=l/ ( 1+ccrit );
                if(dh5)
                    info = Csu5 ( DESCC, mc, DESCYTOT, ytot, rn);
                else
                    info = Csu ( DESCC, mc, DESCYTOT, ytot, rn);
                if ( info!=0 ) {
                    printf ( "Something went wrong with set-up of matrix C, error nr: %d\n",info );
                    return info;
                }
                if ( * ( pst+1 ) ==0 && *pst==0 ) {
                    gettimeofday ( &tz0,NULL );
                    c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
                    printf ( "\t elapsed wall time set-up of C and Y:			%10.3f s\n", ( c0 - c1 ) /1000000.0 );
                }

            }
        }
        else {
            if ( * ( pst+1 ) ==0 && *pst==0 ) {
                gettimeofday ( &tz1,NULL );
                c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
                printf ( "\t elapsed wall time allocation of memory:		%10.3f s\n", ( c1 - c3 ) /1000000.0 );
            }
            if(dh5)
                info = Csu5 ( DESCC, mc, DESCYTOT, ytot, rn);
            else
                info = Csu ( DESCC, mc, DESCYTOT, ytot, rn);
            if ( info!=0 ) {
                printf ( "Something went wrong with set-up of matrix C, error nr: %d\n",info );
                return info;
            }
            if ( * ( pst+1 ) ==0 && *pst==0 ) {
                gettimeofday ( &tz0,NULL );
                c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
                printf ( "\t elapsed wall time set-up of C and Y:			%10.3f s\n", ( c0 - c1 ) /1000000.0 );
            }

            pdcopy_ ( &Cd, ytot,&i_one,&i_one,DESCYTOT,&i_one,RHS,&i_one,&i_one,DESCRHS,&i_one );
            if(Ccop)
                pdlacpy_ ( "U", &Cd, &Cd, mc, &i_one, &i_one, DESCC, mcc, &i_one, &i_one, DESCCCOPY);
	    
            if ( * ( pst+1 ) ==0 && *pst==0 ) {
                gettimeofday ( &tz1,NULL );
                c4= tz1.tv_sec*1000000 + ( tz1.tv_usec );
                printf ( "\t elapsed wall time copy of Y (and C):			%10.3f s\n", ( c4 - c0 ) /1000000.0 );
            }
        }

        work= ( double * ) calloc ( 2*b*(Cc+Cr),sizeof ( double ) );
	if ( work==NULL ) {
                printf ( "unable to allocate memory for work (norm)\n" );
                return EXIT_FAILURE;
            }
        normC=pdlansy_ ( "F","U",&Cd,mc,&i_one,&i_one,DESCC,work );
        norm1C=pdlansy_ ( "1","U",&Cd,mc,&i_one,&i_one,DESCC,work );
        Cmax=pdlansy_ ( "M","U",&Cd,mc,&i_one, &i_one, DESCC,work );
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz0,NULL );
            c1= tz0.tv_sec*1000000 + ( tz0.tv_usec );
            printf ( "\t elapsed wall time norm and max of C:			%10.3f s\n", ( c1 - c0 ) /1000000.0 );
            printf ( "The new parallel lambda is: %15.10g\n",l );
            printf ( "norm of y-vector is: %g\n",*rn );
        }

        pdpotrf_ ( "U",&Cd,mc,&i_one, &i_one,DESCC,&info );
        if ( info!=0 ) {
            printf ( "Parallel Cholesky decomposition of C was unsuccesful, error returned: %d\n",info );
            return -1;
        }
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz1,NULL );
            c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
            printf ( "\t elapsed wall time Cholesky decomposition of C:		%10.3f s\n", ( c1 - c0 ) /1000000.0 );
        }

        pdpotrs_ ( "U",&Cd,&i_one,mc,&i_one,&i_one,DESCC,ytot,&i_one,&i_one,DESCYTOT,&info );
        if ( info!=0 )
            printf ( "Parallel Cholesky solution was unsuccesful, error returned: %d\n",info );
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz0,NULL );
            c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
            printf ( "\t elapsed wall time estimation of effects:		%10.3f s\n", ( c0 - c1 ) /1000000.0 );
        }

        pddot_ ( &Cd,&dot,RHS,&i_one,&i_one,DESCRHS,&i_one,ytot,&i_one,&i_one,DESCYTOT,&i_one );

        si= ( *rn - dot ) / ( n-t );
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            dgebs2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one,&si,&i_one );
        } else
            dgebr2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one, &si,&i_one,&i_zero,&i_zero );

        llh=Cld(mc,DESCC);
        dgsum2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one,&llh,&i_one,&i_negone,&i_negone );

        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz1,NULL );
            c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
            printf ( "\t elapsed wall time calculation and sending of sigma and log(det(C)):	%10.3f s\n", ( c1 - c0 ) /1000000.0 );
        }

        if(dh5)
            info = Asu5 ( mai, DESCAI,DESCYTOT, ytot, DESCC, mc, si) ;
        else
            info = Asu ( mai, DESCAI,DESCYTOT, ytot, DESCC, mc, si) ;

        if ( info!=0 ) {
            printf ( "Something went wrong with set-up of AI-matrix, error nr: %d\n",info );
            return EXIT_FAILURE;
        }

        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz0,NULL );
            c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
            printf ( "\t elapsed wall time set up of AI matrix:			%10.3f s\n", ( c0 - c1 ) /1000000.0 );
        }

        pdpotri_ ( "U",&Cd,mc,&i_one,&i_one,DESCC,&info );
        if ( info!=0 )
            printf ( "Parallel Cholesky inverse was unsuccesful, error returned: %d\n",info );
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz1,NULL );
            c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
            printf ( "\t elapsed wall time inverse of C:			%10.3f s\n", ( c1 - c0 ) /1000000.0 );
        }

        norminv=pdlansy_ ( "F","U",&Cd,mc,&i_one,&i_one,DESCC,work );
        norm1inv=pdlansy_ ( "1","U",&Cd,mc,&i_one,&i_one,DESCC,work );
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz0,NULL );
            c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
            printf ( "\t elapsed wall time set norm of inverse of C:		%10.3f s\n", ( c0 - c1 ) /1000000.0 );
	    process_mem_usage ( vm_usage, resident_set, cpu_user, cpu_sys );
        }

        free ( work );

        tp=CZt ( mc,DESCC );
	
	free ( mc );

        dgsum2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one,&tp,&i_one,&i_negone,&i_negone );

        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz1,NULL );
            c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
            printf ( "\t elapsed wall time trace of inverse of C:		%10.3f s\n", ( c1 - c0 ) /1000000.0 );
        }

        pdnrm2_ ( &m,ran,ytot,&t_plus,&i_one,DESCYTOT,&i_one );

        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz0,NULL );
            c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
            printf ( "\t elapsed wall time set norm of estimation of u:		%10.3f s\n", ( c0 - c1 ) /1000000.0 );
            double *score;
            printf ( "dot product = %15.10g \n",dot );
            printf ( "parallel sigma = %15.10g\n",si );
            printf ( "The trace of CZZ is: %15.10g \n",tp );
            printf ( "The norm of the estimation of u is: %g \n",*ran );
            llh += (m * log(1/l) + (n-t) * log(si) + n-t)/2;
            llh *= -1.0;


            score= ( double * ) calloc ( 2,sizeof ( double ) );
	    if ( score==NULL ) {
                printf ( "unable to allocate memory for score function\n" );
                return EXIT_FAILURE;
            }
            * ( score+1 ) = - ( m-tp*l- *ran * *ran * l / si ) * l / 2;
            printf ( "The score function is: %g\n",* ( score+1 ) );
            printdense ( 2,2, mai, "AI_par.txt" );
            bv=0;
            if (fabs(*(score+1))< e * e) {
                printf("Score function too close to zero to go further, solution may not have converged\n ");
                bv=1;
                igebs2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one,&bv,&i_one );
                break;
            }
            igebs2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one,&bv,&i_one );

            if (co==1) {
                tp= *mai + *(mai+3);
                pll=llh;
                printf ( "The loglikelihood is: %g\n",llh );

            }
            else {
                ull=llh-pll;
                pll=llh;
                printf ( "The update for the loglikelihood is: %g\n",ull );
                printf ( "The new loglikelihood is: %g\n",llh );
            }
            dgebs2d_ (&ICTXT2D,"ALL","1-tree",&i_one,&i_one,&ull,&i_one);
            dpotrf_ ( "U", &i_two, mai, &i_two, &info );
            if ( info!=0 ) {
                printf ( "Cholesky decomposition of AI matrix was unsuccesful, error returned: %d\n",info );
                return -1;
            }
            dpotrs_ ( "U",&i_two,&i_one,mai,&i_two,score,&i_two,&info );
            if ( info!=0 ) {
                printf ( "Parallel solution for AI matrix was unsuccesful, error returned: %d\n",info );
                return -1;
            }
            gettimeofday ( &tz1,NULL );
            c1= tz1.tv_sec*1000000 + ( tz1.tv_usec );
            printf ( "\t elapsed wall time update for lambda:			%10.3f s\n", ( c1 - c0 ) /1000000.0 );
            printf ( "The update for gamma is: %g \n", *(score+1) );
            while (*(score+1)+1/l <0) {
                *(score+1)=*(score+1)/2;
                printf("Half a step is used to avoid negative gamma\n");
            }
            ccrit=l * * ( score+1 );
            dgebs2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one,&ccrit,&i_one );
            free ( score );
            printf ( "The eventual relative update for gamma is: %g \n", ccrit );

        }
        else {
            igebr2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one, &bv,&i_one,&i_zero,&i_zero );
            if (bv >0) {
                break;
            }
            dgebr2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one, &ull,&i_one,&i_zero,&i_zero );
            dgebr2d_ ( &ICTXT2D,"ALL","1-tree",&i_one,&i_one, &ccrit,&i_one,&i_zero,&i_zero );
        }
        if ( * ( pst+1 ) ==0 && *pst==0 ) {
            gettimeofday ( &tz0,NULL );
            c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
            printf ( "\t elapsed wall time sending and receiving update lambda:	%10.3f s\n", ( c0 - c1 ) /1000000.0 );
            printf ( "\t elapsed wall time iteration loop %d:			%10.3f s\n", co, ( c0 - c3 ) /1000000.0 );
        }

    }
    blacs_barrier_ ( &ICTXT2D,"A" );
    
    
    free ( mai );
    free ( RHS );
    free ( rn );
    free ( ran );
    free ( DESCC );
    free ( DESCAI );
    free ( DESCRHS );
    
    if (Ccop){
      free ( mcc );
      free ( DESCCCOPY);
    }

    if (nt>0) {
        if(dh5)
            info=cv5(ytot,DESCYTOT);
        else
            info=cv( ytot, DESCYTOT);
        if ( info!=0 ) {
            printf ( "Cross-validation was unsuccesful, error returned: %d\n",info );
            return -1;
        }
    }

    if ( * ( pst+1 ) ==0 ) {
        for ( i=0,j=0; i<Cb; ++i,++j ) {
            if ( j==*ds )
                j=0;
            if ( *pst==j ) {
                dgesd2d_ ( &ICTXT2D,&b,&i_one,ytot+ i / *ds *b,&b,&i_zero,&i_zero );
            }
            if ( *pst==0 ) {
                dgerv2d_ ( &ICTXT2D,&b,&i_one,ets+b*i,&b,&j,&i_zero );
            }
        }
    }

    blacs_barrier_ ( &ICTXT2D, "A" );

    if ( * ( pst+1 ) ==0 && *pst==0 ) {
        gettimeofday ( &tz0,NULL );
        c0= tz0.tv_sec*1000000 + ( tz0.tv_usec );
        printf ( "\n\tOverall results:\n");
        printf ( "\t================\n");
        printf ( "\tThe maximum element in C is:          %10.5f\n", Cmax );
        printf ( "\tThe Frobenius norm of C is:           %15.10e\n", normC );
        printf ( "\tThe 1-norm of C is:                   %15.10e\n", norm1C );
        printf ( "\tThe Frobenius norm of Cinv is:        %15.10e\n", norminv );
        printf ( "\tThe 1-norm of Cinv is:                %15.10e\n", norm1inv );
        printf ( "\tThe Frobenius condition number is:    %15.10e\n", norminv*normC );
        printf ( "\tThe condition number (1-norm) is:     %15.10e\n", norm1inv*norm1C );
        printf ( "\tThe accuracy is:                      %15.10e\n", norminv*normC*Cmax/pow ( 2,53 ) );
        printf ( "\tThe ultimate lambda is:               %15.10g\n",l );
        printf ( "\tThe ultimate sigma is:                %15.10g\n", si );

        printf ( "\telapsed total wall time:              %10.3f s\n", ( c0 - c2 ) /1000000.0 );

        printf ( "\tProcessor: %d \n\t ========================\n", iam );
        printf ( "\tVirtual memory used:                  %10.0f kb\n", vm_usage );
        printf ( "\tResident set size:                    %10.0f kb\n", resident_set );
        printf ( "\tCPU time (user):                      %10.3f s\n", cpu_user );
        printf ( "\tCPU time (system):                    %10.3f s\n", cpu_sys );
        printdense(Cd,1,ets,"estimates.txt" );
	free(ets);
    }
    
    free ( DESCYTOT );
    free ( pst ),free ( ds );
    free ( fDn );
    free ( fXn );
    free ( Ts );
    free ( ytot );
    free (Sd);
    free (pd);

    blacs_barrier_ ( &ICTXT2D, "A" );
    blacs_gridexit_ ( &ICTXT2D );
    
    return 0;
    
    MPI_Finalize();
}
