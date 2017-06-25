 
#pragma once

#include "targetver.h"


#define MAX_ITER 150
#define N_CONSTRAINTS 2
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN             //  从 Windows 头文件中排除极少使用的信息
// Windows 头文件:


#include <windows.h>
#include <tchar.h>
#endif
// C 运行时头文件
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>

#include <ctime>
#ifdef __unix__
#define WINAPI 

typedef  void * LPVOID;
#ifdef __x86_64__
typedef unsigned int DWORD;
#elif __i386__
typedef unsigned long DWORD;
#endif
#include <pthread.h>
#define  HANDLE pthread_t
#define THREAD_RETURN void*
#define THREAD_PARAM void *
#else
#define THREAD_RETURN DWORD
#define THREAD_PARAM LPVOID
#endif

typedef DWORD cuDWORD;

typedef int(*Device_Stub)(void *, float*, float*, float *, float *);// device function format

enum Stubs { CROSSINGREMOVAL, STRESSCONSTRAINTS, CIRCLE, EUQALANGLE, NOVERLAP, SYMMETRYCONSTRAINTS, ATOMCONSTRAINTS };

enum LensType{ UNKNOWTYPE, PATH, REMOVECROSSING, CIRCLES, ANGLE, ORTHO, AREA_ANGLE, 
	AREA_HORI, AREA_VERTI, AREA_ORTHO, SELECTNODE, SYMMETRY, SHAREDNODESLENS, AREA_ASSIGNEDDIRE
};
