#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include "gene_core.h"
#include "exsort.h"

#undef DEBUG

#define SHELL 24
#define SMAX  20
#define MOVE(A,B) memcpy(A,B,rsize)

static int64 Shell;

#ifdef DEBUG

static inline void sorted(uint8 *array, int asize, int rsize, int ksize)
{ int p;
  for (p = rsize; p < asize; p += rsize)
    if (memcmp(array+(p-rsize),array+p,ksize) > 0)
      printf("Not sorted\n");
}

#endif

static inline void gap_sort(uint8 *array, int asize, int rsize, int ksize, int gap)
{ int    i, j, step;
  uint8  temp[rsize];
  uint8 *garray;

  step   = gap*rsize;
  garray = array + step;
  for (i = step; i < asize; i += rsize)
    { j = i-step;
      if (memcmp(array+j,array+i,ksize) <= 0)
        continue;
      MOVE(temp,array+i);
      MOVE(array+i,array+j);
      for(j -= step; j >= 0; j -= step)
        { if (memcmp(array+j,temp,ksize) <= 0)
            break;
          MOVE(garray+j,array+j);
        }
      MOVE(garray+j,temp);
    }
}

static inline void shell_sort(uint8 *array, int asize, int rsize, int ksize)
{ if (memcmp(array,array+(asize-rsize),ksize) == 0)
    return;
  gap_sort(array,asize,rsize,ksize,10);
  gap_sort(array,asize,rsize,ksize,4);
  gap_sort(array,asize,rsize,ksize,1);
}

typedef struct
  { char  *path;
    int64  asize;
    int    rsize;
    int    ksize;
    int64  offset;
    int    npart;
    int64 *parts;
  } Arg;

static Arg   *Parms;
static int    Nthreads;
static int64  Parts[256];

static void radix_sort(uint8 *array, int64 asize, int rsize, int ksize, int digit)
{ int64  len[256];
  int64  beg[256];
  int    x;

  { int64  end[256];
    uint8 *darray = array + digit;

    int64  off[256];
    uint8  temp[rsize];
    int64  stack[SMAX];

    // fprintf(stderr,"digit: %d asize=%lld, digit=%d, rsize=%d, ksize=%d)\n",
                 // digit,asize,digit,rsize,ksize);

    { int64 o;

      for (x = 0; x < 256; x++)
        len[x] = 0;
      for (o = 0; o < asize; o += rsize)
        len[darray[o]] += rsize;

      o = 0;
      for (x = 0; x < 256; x++)
        { beg[x] = off[x] = o;
          end[x] = o += len[x];
// printf("  %3d: %12lld .. %12lld = %12lld\n",x,beg[x],end[x],len[x]);
        }
    }

    for (x = 0; x < 256; x++)
      { uint8 *o, *p;
        int    t, s;
        int64  u;

        while (off[x] < end[x])
          { t = darray[off[x]];
            if (t == x)
              off[x] += rsize;
            else
              { s = 0;
                stack[s++] = off[x];
                while (s < SMAX)
                  { u = off[t];
                    off[t] = u + rsize;
                    if (t == x)
                      break;
                    stack[s++] = u;
                    t = darray[u];
                  }
  
                o = array + stack[--s];
                MOVE(temp,o);
	        while (s > 0)
                  { p = array + stack[--s];
                    MOVE(o,p);
                    o = p;
                  }
                MOVE(o,temp);
              }
          }
      }

// printf("Pass finished\n");
  }
  
  digit += 1;
  if (digit >= ksize)
    return;

  if (digit == 1)
    { int   n, beg;
      int64 sum, thr, off;

      n   = 0;
      thr = asize / Nthreads;
#ifdef DEBUG
      printf("thr = %lld\n",thr);
#endif
      off = 0;
      sum = 0;
      beg = 0;
      for (x = 0; x < 256; x++)
        { Parts[x] = len[x];
          sum += len[x];
#ifdef DEBUG
          printf("%d: sum %lld\n",x,sum);
#endif
          if (sum >= thr)
            { Parms[n].offset = off;
              Parms[n].npart  = (x+1) - beg;
              Parms[n].parts  = Parts + beg;
#ifdef DEBUG
              printf("  %d - %d\n",beg,x);
#endif
              n  += 1;
              thr = (asize * (n+1))/Nthreads;
              beg = x+1;
              off = sum;
#ifdef DEBUG
              printf("thr = %lld\n",thr);
#endif
            }
	}
    }
  else
    { for (x = 0; x < 256; x++)
        if (len[x] > Shell)
          radix_sort(array + beg[x], len[x], rsize, ksize, digit);
        else if (len[x] > rsize)
          shell_sort(array + beg[x], len[x], rsize, ksize);
    }
}

static void *sort_thread(void *arg) 
{ Arg *param = (Arg *) arg;

  char  *path  = param->path;
  int64  asize = param->asize;
  int    rsize = param->rsize;
  int    ksize = param->ksize;
  int    npart = param->npart;
  int64 *parts = param->parts;

  uint8 *array;
  int    fd;
  int64  off;
  int    x;

  fd = open(path, O_RDWR);
  if (fd == -1)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,path);
      exit (1);
    }

  array = (uint8 *) mmap(NULL,asize,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);
  if (array == NULL)
    { close(fd);
      fprintf(stderr,"%s: Cannot memory map %s\n",Prog_Name,path);
      exit (1);
    }

  madvise(array, asize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

#ifdef DEBUG

   printf("Thread %d: %lld",npart,param->offset);
   for (x = 0; x < npart; x++)
     printf(" %lld",parts[x]);
   printf("\n");

#endif

  off = param->offset;
  for (x = 0; x < npart; x++)
    { if (parts[x] > Shell)
        radix_sort(array + off, parts[x], rsize, ksize, 1);
      else if (parts[x] > rsize)
        shell_sort(array + off, parts[x], rsize, ksize);
      off += parts[x];
    }

  munmap(array,asize);
  close(fd);

  return (NULL);
}


void Ex_sort(char *path, int rsize, int ksize, int nthreads)
{ pthread_t   threads[nthreads];
  Arg         parms[nthreads];
  int64       asize;
  int         x;

  { uint8      *array;
    struct stat stats;
    int         fd;

    fd = open(path, O_RDWR);
    if (fd == -1)
      { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,path);
        exit (1);
      }

    if (fstat(fd, &stats) == -1)
      { close(fd);
        fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,path);
        exit (1);
      }
    asize = stats.st_size;

    array = (uint8 *) mmap(NULL,asize,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);
    if (array == NULL)
      { close(fd);
        fprintf(stderr,"%s: Cannot memory map %s\n",Prog_Name,path);
        exit (1);
      }

    madvise(array, asize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

    Parms    = parms;
    Shell    = SHELL * rsize;
    Nthreads = nthreads;

    radix_sort(array, asize, rsize, ksize, 0);

    munmap(array,asize);
    close(fd);
  }

  if (ksize == 1)
    return;

  for (x = 0; x < nthreads; x++)
    { parms[x].path  = path;
      parms[x].asize = asize;
      parms[x].rsize = rsize;
      parms[x].ksize = ksize;
    }

  for (x = 0; x < nthreads; x++)
    pthread_create(threads+x,NULL,sort_thread,parms+x);

  for (x = 0; x < nthreads; x++)
    pthread_join(threads[x],NULL);

#ifdef DEBUG
  { uint8      *array;
    int         fd;

    fd = open(path, O_RDWR);
    if (fd == -1)
      { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,path);
        exit (1);
      }

    array = (uint8 *) mmap(NULL,asize,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);
    if (array == NULL)
      { close(fd);
        fprintf(stderr,"%s: Cannot memory map %s\n",Prog_Name,path);
        exit (1);
      }

    madvise(array, asize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

    sorted(array,asize,rsize,ksize);

    munmap(array,asize);
    close(fd);
  }
#endif
  
}

int main(int argc, char** argv) {
    int opt;
    int char_start = 0;
    int char_stop = 255;
    int record_size=100;
    int key_size=10;
    int stack_size=5;
    int cut_off = 4;
    int switch_to_shell = 20;
    int verbosity = 0;
    int nthreads = 1;
    while ((opt = getopt(argc, argv, "var:k:s:c:t:w:")) != -1) {
        switch (opt) {
        case 'v':
            verbosity += 1;
            break;
        case 'a':
            char_start = 32;
            char_stop = 128;
            break;
        case 'r':
            record_size = atoi(optarg);
            break;
        case 'k':
            key_size = atoi(optarg);
            break;
        case 's':
            stack_size = atoi(optarg);
            break;
        case 'w':
            switch_to_shell = atoi(optarg);
            break;
        case 't':
            nthreads = atoi(optarg);
            break;
        case 'c':
            cut_off = atoi(optarg);
            break;
        default:
            fprintf(stderr, "Invalid parameter: -%c\n", opt);
            return 1;
            break;
        }
    }

    while(optind < argc) {
        Ex_sort(argv[optind++], record_size, key_size, nthreads);
    }

    return 0;
}
