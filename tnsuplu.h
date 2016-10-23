/*  activate only one   */
#define SUPERLU_USE 0
#define SUPERLU_MT_USE 0
#define SUPERLU_DIST_USE 0
/* don't change the next lines */
#if SUPERLU_USE+SUPERLU_MT_USE+SUPERLU_DIST_USE > 1
 ERROR___  just pick one 
#endif
#if SUPERLU_DIST_USE
#undef MPI_USE
#define MPI_USE 1
#endif
