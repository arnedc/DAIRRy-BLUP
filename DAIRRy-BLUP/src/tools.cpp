#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printdense ( int m, int n, double *mat, char *filename ) {
    FILE *fd;
    fd = fopen ( filename,"w" );
    if ( fd==NULL )
        printf ( "error creating file" );
    int i,j;
    for ( i=0; i<m; ++i ) {
        fprintf ( fd,"[\t" );
        for ( j=0; j<n; ++j ) {
            fprintf ( fd,"%12.8g\t",*(mat+i*n +j));
        }
        fprintf ( fd,"]\n" );
    }
    fclose ( fd );
}


//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set, double& cpu_user, double& cpu_sys)
{
    using std::ios_base;
    using std::ifstream;
    using std::string;

    vm_usage     = 0.0;
    resident_set = 0.0;
    cpu_user     = 0.0;
    cpu_sys      = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize, utime, stime;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage     = vsize / 1024.0;
    resident_set = rss * page_size_kb;
    cpu_user     = utime / (float) sysconf(_SC_CLK_TCK);
    cpu_sys      = stime / (float) sysconf(_SC_CLK_TCK);
}

