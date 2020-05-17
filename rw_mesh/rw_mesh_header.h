/*
 * rw_mesh_header.h
 * Общий заголовочный файл для всех форматов
 *
 *  Created on: 26.06.2013
 *      Author: moric
 */

#ifndef RW_MESH_HEADER_H_
#define RW_MESH_HEADER_H_

#include <stdio.h>
#include <locale.h>
#include "rw_basic_types.h"

#ifndef NULL
#define NULL 0
#endif

typedef int INT;
typedef INT INT7[7];
//typedef INT INT8[8];


#ifndef log_write
#define log_write
#endif

#define __save_locale char*oldlocale = setlocale(LC_NUMERIC, "C")
#define __load_locale setlocale(LC_NUMERIC, oldlocale)

#ifndef ffree
#define ffree(pointer) if(pointer){free(pointer);pointer=NULL;}
#endif

void rw_mesh_set_filename(const char*filename);
void rw_mesh_set_error(int lineNumber,const char*error);
void rw_mesh_get_error(int*lineNumber,char*error,char*filename);
void rw_mesh_print_error();

#endif /* RW_MESH_HEADER_H_ */
