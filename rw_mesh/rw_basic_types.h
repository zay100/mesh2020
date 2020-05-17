/*changed 20.06.2012*/
#pragma warning(disable:4996)
#ifndef RW_BASIC_TYPES_H_
#define RW_BASIC_TYPES_H_

#ifndef REAL
#define REAL double
#endif

#ifndef REAL2
typedef REAL REAL2[2];
#define REAL2 REAL2 
#endif

#ifndef REAL3
typedef REAL REAL3[3];
#define REAL3 REAL3 
#endif

#ifndef REAL4
typedef REAL REAL4[4];
#define REAL4 REAL4 
#endif

#ifndef REAL6
typedef REAL REAL6[6];
#define REAL6 REAL6 
#endif

#ifndef REAL7
typedef REAL REAL7[7];
#define REAL7 REAL7 
#endif

#ifndef REAL9
typedef REAL REAL9[9];
#define REAL9 REAL9 
#endif

#ifndef REAL22
typedef REAL REAL22[2][2];
#define REAL22 REAL22 
#endif

#ifndef REAL33
typedef REAL REAL33[3][3];
#define REAL33 REAL33 
#endif

#ifndef REAL23
typedef REAL REAL23[2][3];
#define REAL23 REAL23 
#endif

#ifndef INT2
typedef int INT2[2];
#define INT2 INT2 
#endif

#ifndef INT3
typedef int INT3[3];
#define INT3 INT3 
#endif

#ifndef INT4
typedef int INT4[4];
#define INT4 INT4 
#endif

#ifndef INT5
typedef int INT5[5];
#define INT5 INT5 
#endif

#ifndef INT6
typedef int INT6[6];
#define INT6 INT6 
#endif

#ifndef INT_8
typedef int INT_8[8];
#define INT_8 INT_8 
#endif

#endif /*RW_BASIC_TYPES_H_*/
