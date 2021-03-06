/*
 * rw_mesh_header.c
 *
 *  Created on: 01.06.2014
 *      Author: moric
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rw_mesh_header.h"

int __error_lineNumber;
char __error_filename[256];
char __error_comment[256];

void rw_mesh_set_filename(const char*filename){
	strncpy(__error_filename,filename,256);
}

void rw_mesh_set_error(int lineNumber,const char *error){
	__error_lineNumber = lineNumber;
	strncpy(__error_comment,error,256);
}

void rw_mesh_get_error(int*lineNumber,char*error,char*filename){
	if(lineNumber)*lineNumber = __error_lineNumber;
	if(error)strncpy(error,__error_comment,256);
	if(filename)strncpy(filename,__error_filename,256);
}

void rw_mesh_print_error(){
	printf("Error parsing file '%s' on line %d: %s;\n",__error_filename,__error_lineNumber,__error_comment);
}


