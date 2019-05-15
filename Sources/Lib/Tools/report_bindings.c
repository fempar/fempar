#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#include <mpi.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#ifdef __APPLE__

#include <stdint.h>

#ifndef CPU_SETSIZE
#define CPU_SETSIZE 1024
#endif

#ifndef SYSCTL_CORE_COUNT 
#define SYSCTL_CORE_COUNT   "machdep.cpu.core_count"
#endif

typedef struct cpu_set {
  uint32_t    count;
} cpu_set_t;

static inline int CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set) {
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int     i;
  int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }

  return 0;
}

#endif


/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
    char *ptr = str;
    int i, j, entry_made = 0;
    for (i = 0; i < CPU_SETSIZE; i++) {
        if (CPU_ISSET(i, mask)) {
            int run = 0;
            entry_made = 1;
            for (j = i + 1; j < CPU_SETSIZE; j++) {
                if (CPU_ISSET(j, mask)) run++;
                else break;
            }
            if (!run)
                sprintf(ptr, "%d,", i);
            else if (run == 1) {
                sprintf(ptr, "%d,%d,", i, i + 1);
                i++;
            } else {
                sprintf(ptr, "%d-%d,", i, i + run);
                i += run;
            }
            while (*ptr != 0) ptr++;
        }
    }
    ptr -= entry_made;
    *ptr = 0;
    return(str);
}

void report_bindings(const MPI_Fint Fcomm)
{
    int rank, thread;
    cpu_set_t coremask;
    char clbuf[7 * CPU_SETSIZE], hnbuf[64];

    MPI_Comm Ccomm;
    Ccomm = MPI_Comm_f2c(Fcomm); // Convert Fortran->C communicator
    
    MPI_Comm_rank(Ccomm, &rank);
    memset(clbuf, 0, sizeof(clbuf));
    memset(hnbuf, 0, sizeof(hnbuf));
    (void)gethostname(hnbuf, sizeof(hnbuf));
#ifndef _OPENMP
    thread = 0;
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);
    printf("Hello from rank %d, thread %d, on %s. (core affinity = %s)\n",
           rank, thread, hnbuf, clbuf);
#else
    thread = omp_get_thread_num();
    (void)sched_getaffinity(0, sizeof(coremask), &coremask);
    cpuset_to_cstr(&coremask, clbuf);
    #pragma omp barrier
    printf("Hello from rank %d, thread %d, on %s. (core affinity = %s)\n",
            rank, thread, hnbuf, clbuf);
#endif    
    
    
}
