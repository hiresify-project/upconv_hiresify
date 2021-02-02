#include <string.h>
#include <stdio.h>
#include "makepath.h"
//void _splitpath(char* str, char* drive, char* dir, char* fname, char* ext);
//void _makepath(char* str, char* drive, char* dir, char* fname, char* ext);

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
void _splitpath(char* str, char* drive, char* dir, char* fname, char* ext){
	int dotcount = 0;
	int slashcount = 0;
	char tmp[4096];
	tmp[0] = '\0';
	int i;
	for(i=strlen(str);i>=0;i--){
		if( str[i] == '.' && dotcount == 0){
			strncpy( ext, &(str[i+1]), strlen(str) - i -1);
			ext[strlen(str) - i] = '\0';
			dotcount = i;
		}
		if( str[i] == '/' && slashcount == 0){
			if(dotcount != 0){
				strncpy( fname, &(str[i+1]), dotcount - i -1);
				fname[dotcount - i -1 ] = '\0';
				slashcount = i;
				
				strncpy(dir, str, slashcount);
				dir[slashcount] = '\0';
				drive[0] = '\0';
			}
		}
	}

//	printf("%s,%s,%s\n", dir, fname, ext);
}

void _makepath(char* str, char* drive, char* dir, char* fname, char* ext){
	str[0] = '\0';
	strcat(str, dir);
	strcat(str, "/");
	strcat(str, fname);
	strcat(str, ".");
	strcat(str, ext);
}
