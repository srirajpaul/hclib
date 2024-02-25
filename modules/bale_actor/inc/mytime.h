#include <time.h>
#include <sys/time.h>

double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

