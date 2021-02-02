#ifndef __MAKEPATH_H__
#define __MAKEPATH_H__

#include <string.h>
#include <stdio.h>

void _splitpath(char* str, char* drive, char* dir, char* fname, char* ext);
void _makepath(char* str, char* drive, char* dir, char* fname, char* ext);

/*int main(int argc, char* argv[]){
	char str[3000];
	char drive[20];
	char dir[2000];
	char fname[100];
	char ext[30];
	
	puts(argv[1]);
	_splitpath(argv[1], drive, dir, fname, ext);
	puts(ext);
	puts(fname);
	puts(dir);
	
	_makepath(str, drive, dir, fname, ext);
	puts(str);
	return 0;
}*/

#endif
