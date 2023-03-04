
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define max(a,b) ( ( (a) > (b) ) ? (a) : (b) )
#define min(a,b) ( ( (a) < (b) ) ? (a) : (b) )

#define	g  (double)1.4
#define	ga (double)(0.5*(g+1.)/g)
#define	gb (double)(0.5*(g-1.)/g)
#define	gc (double)(2./(g-1.))
#define	gg (double)(1./g)

#define maxit 100
#define ee 0.0000001
#define estr 1.0000001 

