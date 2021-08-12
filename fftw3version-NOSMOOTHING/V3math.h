/* V3math.h */
#ifndef H_V3MATH
#define H_V3MATH
#include <math.h>
#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif
#define Sqr(x)	((x) * (x))
#define Cub(x)	((x) * (x) * (x))
#define VAdd(v1, v2, v3)   \
	(v1).x = (v2).x + (v3).x,   \
	(v1).y = (v2).y + (v3).y,	\
	(v1).z = (v2).z + (v3).z
#define VSub(v1, v2, v3)   \
	(v1).x = (v2).x - (v3).x,   \
	(v1).y = (v2).y - (v3).y,	\
	(v1).z = (v2).z - (v3).z
#define VDot(v1, v2)   \
	((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)
#define VSAdd(v1, v2, s3, v3)   \
	(v1).x = (v2).x + (s3) * (v3).x,	\
	(v1).y = (v2).y + (s3) * (v3).y,	\
	(v1).z = (v2).z + (s3) * (v3).z	
#define VSet(v, sx, sy, sz)   \
	(v).x = sx,   \
	(v).y = sy,	  \
	(v).z = sz
#define VMul(v1, v2, v3)   \
	(v1).x = (v2).x * (v3).x,   \
	(v1).y = (v2).y * (v3).y,	\
	(v1).z = (v2).z * (v3).z
#define VDiv(v1, v2, v3)   \
	(v1).x = (v2).x / (v3).x,   \
	(v1).y = (v2).y / (v3).y,	\
	(v1).z = (v2).z / (v3).z
#define VScale(v, s)  \
	(v).x *= s,		\
	(v).y *= s,		\
	(v).z *= s
#define VSCopy(v2, s1, v1)   \
	(v2).x = (s1)*(v1).x,   \
	(v2).y = (s1)*(v1).y,	\
	(v2).z = (s1)*(v1).z
#define VCSum(v)   \
	((v).x + (v).y + (v).z)
#define VProd(v)	((v).x * (v).y * (v).z)
#define VVAdd(v1, v2) VAdd(v1, v1, v2)
#define VSetAll(v,s)		VSet(v, s, s, s)
#define VZero(v)			VSetAll (v,0.0)
#define VVSAdd(v1, s2, v2)	VSAdd(v1, v1, s2, v2)
#define VLenSq(v)			VDot(v, v)
#define VSecNorm(v)	\
	(sqrt(Sqr((v).x) + Sqr((v).y) + Sqr((v).z)))
#define Max(x1, x2)   \
	(((x1) > (x2)) ? (x1) : (x2))
#define Min(x1, x2)   \
	(((x1) < (x2)) ? (x1) : (x2))
#define VCross(v0, v1, v2)	\
	(v0).x = (v1).y * (v2).z - (v1).z * (v2).y,	\
	(v0).y = (v1).z * (v2).x - (v1).x * (v2).z,	\
	(v0).z = (v1).x * (v2).y - (v1).y * (v2).x

/* complex operators */
#define CmplxSet(v0, x, y)	\
	(v0).re = x, (v0).im = y
#define CmplxCopy(v0, v1)	\
	(v0).re = (v1).re, (v0).im = (v1).im
#define CmplxAdd(c0, c1, c2)	\
	(c0).re = (c1).re + (c2).re,	\
	(c0).im = (c1).im + (c2).im
#define Cmplx3Add(c0, c1, c2, c3)	\
	(c0).re = (c1).re + (c2).re + (c3).re,	\
	(c0).im = (c1).im + (c2).im + (c3).im
#define Cmplx3SAdd(c0, s1, c1, s2, c2, s3, c3)	\
	(c0).re = (s1)*(c1).re + (s2)*(c2).re + (s3)*(c3).re,	\
	(c0).im = (s1)*(c1).im + (s2)*(c2).im + (s3)*(c3).im
#define CmplxSAdd(c0, c1, s, c2)	\
	(c0).re = (c1).re + (s)*(c2).re,	\
	(c0).im = (c1).im + (s)*(c2).im
#define Cmplx2SAdd(c0, s1, c1, s2, c2)	\
	(c0).re = (s1)*(c1).re + (s2)*(c2).re,	\
	(c0).im = (s1)*(c1).im + (s2)*(c2).im
#define CmplxScale(c0, s)	\
	(c0).re *= (s); (c0).im *= (s)

#endif
