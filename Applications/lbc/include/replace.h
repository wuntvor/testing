#ifdef LAYOUT_LIJK
#define NDX(l,i,j,k) l,i,j,k
#endif
#ifdef LAYOUT_IJKL
#define NDX(l,i,j,k) i,j,k,l
#endif
