#define MPI_COMM_WORLD ((MPI_Comm)0x44000000)

void process_mem_usage(double& , double& , double& , double& );
void printdense ( int , int , double *, char * );
int Csu (int* , double* , int* , double* , double* ) ;
int Csu5 ( int* , double* , int*, double* , double*  ) ;
int Cu ( int * , double * , double ) ;
int Asu (double* , int* , int* , double* , int* , double* , double );
int Asu5 ( double* , int* , int* , double* , int* , double* , double  );
double CZt(double *, int * );
double Cld ( double *, int *  ) ;
int ri(char* ) ;
int cv(double* , int* );
int cv5(double * , int *) ;

extern double d_one, d_zero, d_negone;
extern int DLEN_, i_negone, i_zero, i_one, i_two, i_four;
extern int m,n,t, b, m_plus;
extern int lld_C, Cb, Cd, t_plus;
extern long Cr, Cc;
extern int size, *ds, * pst, ICTXT2D, iam;
extern int nt, mi,dh5, Ccop;
extern char *Sd, *pd;
extern char *fDn, *fXn, *Ts;
extern double l, e;
