/****************************************************************************/
/* raw2wav (C) 2011-2012 By 59414d41										*/
/*																			*/
/*																			*/
/****************************************************************************/

/*--- Log ------------------------------------------------------------------
 * Ver 0.10 <09/07/15> - upconvから分離
 * Ver 0.20 <09/11/01> - パラメータファイルの採用
 * Ver 0.21 <09/11/16> - エラー情報をファイルへ出力するようにした
 * Ver 0.22 <09/12/08> - 32k のサンプリングレートに対応
 * Ver 0.50 <10/11/02> - 処理修正
 * Ver 0.80 <12/02/11> - コンパイラをmingwに変更
 *						 大きなファイルに対応
 *						 split処理をやめたことによる修正
 *						 マルチチャンネルに対応
 * Ver 0.90 <12/07/02> - multconv対応、multconvのマルチチャンネル出力、
 *						 個別にモノラルファイルとして出力する機能対応
 */

#define STR_COPYRIGHT	"raw2wav (c) 2011-2012 Ver 0.90 By 59414d41\n\n"
#define STR_USAGE		"raw2wav fromfile rawfile\nraw2wav fromfile rawfile flacfile"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./../upconv/fileio.h"
#include "PLG_AudioIO.h"

#define TRUE 1
#define FALSE 0

// サンプルを処理するデータ型
#define SSIZE	signed long long
#define USSIZE	unsigned long long

typedef struct {
	long sampling;
	long bitwidth;
	long mode;
	long norm;
	long split;
	long downmix;
	int ch;
	int chC;
	int chS;
	int chLFE;
	int bwf;
	int raw;
	int rf64;
	int abe;
	int hfa;
	int flac;
	int thread;
	int fio;
	int	errLine;
	double	nx1;
	double	nx2;
	double	nx3;
	double	nx4;
	double	nx5;
	double	nx6;
	SSIZE	sum1;
	SSIZE	sum2;
	SSIZE	sum3;
	SSIZE	sum4;
	SSIZE	sum5;
	SSIZE	sum6;
	SSIZE	sum1_2nd;
	SSIZE	sum2_2nd;
	SSIZE	sum3_2nd;
	SSIZE	sum4_2nd;
	SSIZE	sum5_2nd;
	SSIZE	sum6_2nd;
	SSIZE	sum1_3rd;
	SSIZE	sum2_3rd;
	SSIZE	sum3_3rd;
	SSIZE	sum4_3rd;
	SSIZE	sum5_3rd;
	SSIZE	sum6_3rd;
	SSIZE	old_s1;
	SSIZE	old_s2;
	SSIZE	old_s3;
	SSIZE	old_s4;
	SSIZE	old_s5;
	SSIZE	old_s6;
	SSIZE	s_b1,s_a1;
	SSIZE	s_b2,s_a2;
	SSIZE	s_b3,s_a3;
	SSIZE	s_b4,s_a4;
	SSIZE	s_b5,s_a5;
	SSIZE	s_b6,s_a6;
	int ditherLv;
	char fromfile[_MAX_PATH];
	char flacfile[_MAX_PATH];
	char tofile[_MAX_PATH];
	FIO *fp_w[6];
} PARAM_INFO;

// BWF の link チャンク対応
char link_start[]="<LINK>\r\n";
char link_file[]=						\
	"\t<FILE type=\"%s\">\r\n"				\
	"\t\t<FILENUMBER>%d</FILENUMBER>\r\n"	\
	"\t\t<FILENAME>%s</FILENAME>\r\n"		\
	"\t</FILE>\r\n";
char link_end[]="</LINK>";

#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("/tmp/raw2wav.log","a");									\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n",__FUNCTION__,__LINE__,s);			\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif

int Normalize(int *nCount,SOUNDFMT *inFmt,SOUNDFMT *outFmt,FIO *fp_r1,FIO *fp_r2,FIO *fp_r3,FIO *fp_r4,FIO *fp_r5,FIO *fp_r6,DWORD *start,DWORD nSample,PARAM_INFO *param);
int Normalize_Mx(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M0(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M1(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M2(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int Normalize_M3(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param);
int UpdateBext(BROADCAST_EXT *bext,SOUNDFMT *inFmt,SOUNDFMT *outFmt,PARAM_INFO *param,long bwf_size);

int WriteData(FIO *fp_w,int bit,SSIZE data);
double normalNoise(void);

#if 1
#define CLIP_NX(v,nx)	(v) == 0 ? ((SSIZE)0) : \
							((v) > 0 ?	\
								(((SSIZE)(v) * nx) >= (SSIZE)(0x007FFFFFFFFFFFFF) ?		(SSIZE)(0x7FFFFFFFFFFFFF)		: ((SSIZE)((v) * nx))) : \
								(((SSIZE)(v) * nx) <= ((SSIZE)0x007FFFFFFFFFFFFF * -1) ? ((SSIZE)0x7FFFFFFFFFFFFF * -1) : ((SSIZE)((v) * nx))))
#else
#define CLIP_NX(v,nx)	((SSIZE)((v) * (nx)) << 6)
#endif

#define ROUND_NBIT(b)	((SSIZE)1 << (b) - 1) 

#define CLIP_ADD(v,a)	(SSIZE)(v) + (a)

#define CLIP_MAX(v)	(v) == 0 ? ((SSIZE)0) : \
							((v) > 0 ?	\
								((SSIZE)(v) >= (SSIZE)(0x007FFFFFFFFFFFFF) ?		(SSIZE)(0x7FFFFFFFFFFFFFFF)		: ((SSIZE)((v) << 8))) : \
								(((SSIZE)(v)) <= ((SSIZE)0x007FFFFFFFFFFFFF * -1) ? ((SSIZE)0x7FFFFFFFFFFFFFFF * -1) : ((SSIZE)((v) << 8))))

#define CLIP_MAX_N(v,b)	(v) == 0 ? ((SSIZE)0) : \
							((v) > 0 ?	\
								((SSIZE)(v) >= ((((SSIZE)1 << (b)) - 1) >> 1) ?		((((SSIZE)1 << (b)) - 1) >> 1)		: (SSIZE)(v)) : \
								(((SSIZE)(v)) <= (((((SSIZE)1 << (b)) - 1) * -1) >> 1) ? (((((SSIZE)1 << (b)) - 1) * -1) >> 1) : ((SSIZE)(v))))

PARAM_INFO paramInfo;

//---------------------------------------------------------------------------
// Function   : main
// Description: 引数を処理し変換関数を呼び出す
//
//
int main(int argc, char *argv[])
{
	char tmppath[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	char param[512];
	char fn[6][10]={"_1","_2","_3","_4","_5","_6"};
	char s[300];
	int paramFlag;
	int ccc;
	FILE *fp;
	FIO fp_r1,fp_r2,fp_r3,fp_r4,fp_r5,fp_r6,fp_w;
	FIO fp_ws[6];
	FIO *p_fp1,*p_fp2,*p_fp3,*p_fp4,*p_fp5,*p_fp6;
	SOUNDFMT inFmt;
	SOUNDFMT outFmt;
	SOUNDFMT flacFmt;
	char *p1,*p2;
	DWORD inSample,outSample,flacSample,startSample;
	SSIZE size;
	double nx;
	long rd;
	int i,count;
	int retCode;
	int thread;
	int fio;
	int genCh;
	double *p_nx[6];
	SSIZE  *p_max[6];
	SSIZE  *p_avg[6];
	FILEINFO fileInfo;
	long temp;
	double perR1,perR2,perR3,perR4,perR5,perR6;
	SSIZE maxR1,maxR2,maxR3,maxR4,maxR5,maxR6;
	SSIZE avgR1,avgR2,avgR3,avgR4,avgR5,avgR6;
	double per;
	SSIZE max,avg;
	thread = 1;
	genCh = 0;
	do {
		memset(&paramInfo,0,sizeof (PARAM_INFO));
		memset(&fileInfo,0,sizeof (FILEINFO));
		memset(&fp_ws,0,sizeof (FIO) * 6);
		paramFlag = 0;
		paramInfo.chC = 0;
		paramInfo.chS = 0;
		paramInfo.chLFE = 0;
		paramInfo.fio = -1;
		// Wave ファイルのは違う順番
		p_nx[0]  = &paramInfo.nx1;	// Left
		p_nx[1]  = &paramInfo.nx2;	// Right
		p_nx[2]  = &paramInfo.nx3;	// Center
		p_nx[3]  = &paramInfo.nx4;	// Surround Left
		p_nx[4]  = &paramInfo.nx5;	// Surround Right
		p_nx[5]  = &paramInfo.nx6;	// LFE
		p_max[0] = &maxR1;
		p_max[1] = &maxR2;
		p_max[2] = &maxR3;
		p_max[3] = &maxR4;
		p_max[4] = &maxR5;
		p_max[5] = &maxR6;
		p_avg[0] = &avgR1;
		p_avg[1] = &avgR2;
		p_avg[2] = &avgR3;
		p_avg[3] = &avgR4;
		p_avg[4] = &avgR5;
		p_avg[5] = &avgR6;

		p_fp1 = p_fp2 = p_fp3 = p_fp4 = p_fp5 = p_fp6 = NULL;
		
		retCode = STATUS_SUCCESS;
		paramInfo.errLine = __LINE__;
		if (argc == 3 || argc == 4) {
			strcpy(paramInfo.fromfile,argv[1]);
			// param ファイル
			_splitpath(argv[2],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"param");
			// ファイルオープン
			fp = fopen(tmppath,"r");
			if (fp == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			
			if (fgets(param,512,fp) == NULL) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			perR1 = perR2 = perR3 = perR4 = perR5 = perR6 = 0;
			if (fscanf(fp,"r1=%lf,%llx\n",&perR1,&avgR1) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r2=%lf,%llx\n",&perR2,&avgR2) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r3=%lf,%llx\n",&perR3,&avgR3) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r4=%lf,%llx\n",&perR4,&avgR4) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r5=%lf,%llx\n",&perR5,&avgR5) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			if (fscanf(fp,"r6=%lf,%llx\n",&perR6,&avgR6) != 2) {
				retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
				break;
			}
			maxR1 = (SSIZE)(perR1 * 0x7FFFFFFFFFFFFF);
			maxR2 = (SSIZE)(perR2 * 0x7FFFFFFFFFFFFF);
			maxR3 = (SSIZE)(perR3 * 0x7FFFFFFFFFFFFF);
			maxR4 = (SSIZE)(perR4 * 0x7FFFFFFFFFFFFF);
			maxR5 = (SSIZE)(perR5 * 0x7FFFFFFFFFFFFF);
			maxR6 = (SSIZE)(perR6 * 0x7FFFFFFFFFFFFF);

			fclose(fp);
			p1 = param;
			p2 = strchr(p1,(int)' ');

			for (;p1 != NULL;) {
				if (p2 != NULL) {
					*p2 = '\0';
				}

				if (sscanf(p1,"-s:%ld",&temp) == 1) {
					switch (temp) {
						case 32000:
						case 44100:
						case 48000:
						case 88200:
						case 96000:
						case 176400:
						case 192000:
						case 352800:
						case 384000:
							paramInfo.sampling = temp;
							paramFlag |= 0x01;
							break;
					}
				}
				if (sscanf(p1,"-w:%ld",&temp) == 1) {
					switch (temp) {
						case 16:
						case 24:
						case 32:
						case 64:
							paramInfo.bitwidth = temp;
							paramFlag |= 0x02;
							break;
					}
				}

				if (sscanf(p1,"-m:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
						case 2:
						case 3:
							paramInfo.mode = temp;
							paramFlag |= 0x04;
							break;
					}
				}
				if (sscanf(p1,"-n:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
							paramInfo.norm = temp;
							paramFlag |= 0x08;
							break;
					}
				}
				if (sscanf(p1,"-ditherLv:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= 16) {
						paramInfo.ditherLv = temp;
						paramFlag |= 0x10;
					}
				}
				if (sscanf(p1,"-ch:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 6) {
						paramInfo.ch = temp;
						paramFlag |= 0x20;
					}
				}
				if (strcmp(p1,"-bwf") == 0) {
					paramInfo.bwf = 1;
				}
				if (strcmp(p1,"-raw") == 0) {
					paramInfo.raw = 1;
				}
				if (strcmp(p1,"-rf64") == 0) {
					paramInfo.rf64 = 1;
				}
				if (strcmp(p1,"-C") == 0) {
					paramInfo.chC = 1;
					genCh++;
				}
				if (strcmp(p1,"-mcdownmix") == 0) {
					paramInfo.downmix = 1;
				}
				if (strcmp(p1,"-SLR") == 0) {
					paramInfo.chS = 1;
					genCh += 2;
				}
				if (strcmp(p1,"-LFE") == 0) {
					paramInfo.chLFE = 1;
					genCh++;
				}
				if (strcmp(p1,"-split") == 0) {
					paramInfo.split = 1;
				}
				if (sscanf(p1,"-hfa:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
						case 2:
						case 3:
							paramInfo.hfa = temp;
							break;
					}
				}
				if (strcmp(p1,"-adjBE") == 0) {
					paramInfo.abe = 1;
				}
				if (strcmp(p1,"-flac") == 0 && argc == 4) {
					paramInfo.flac = 1;
				}
				if (sscanf(p1,"-thread:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 24) {
						thread = (int)temp;
					}
				}
				if (sscanf(p1,"-fio:%ld",&temp) == 1) {
					if (temp >= 30 && temp <= 16000) {
						paramInfo.fio = temp / 10;
					}
				}

				if (p2 == NULL) {
					break;
				}
				p1 = p2 + 1;
				p2 = strchr(p1,(int)' ');
			}
			if (paramFlag == 0x3F) {
#ifdef _OPENMP
	omp_set_num_threads(thread);
#endif
				if (paramInfo.downmix) {
					genCh = 0;
				}
				paramInfo.ch += genCh;

				retCode = PLG_InfoAudioData(argv[1],&inFmt,&inSample,&fileInfo);
				if (retCode != STATUS_SUCCESS) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				if (paramInfo.flac == 1) {
					if (PLG_InfoAudioData(argv[3],&flacFmt,&flacSample,NULL) || strcmp(flacFmt.fmt,"flac")) {
						paramInfo.flac = 0;
					} else {
						strcpy(paramInfo.flacfile,argv[3]);
					}
				}

				outFmt.sample  = paramInfo.sampling;
				if (paramInfo.split == 0) {
					outFmt.channel = paramInfo.ch;
				} else {
					outFmt.channel = 1;
				}
				outFmt.bitsPerSample = (unsigned char)paramInfo.bitwidth;

				if (paramInfo.ch == 3) {
					if (paramInfo.chC == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_max[2] = &maxR3;
						p_avg[2] = &avgR3;
						strcpy(fn[2],"_C");
					} else if (paramInfo.chLFE == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_max[2] = &maxR6;
						p_avg[2] = &avgR6;
						strcpy(fn[2],"_LFE");
					}
				}
				if (paramInfo.ch == 4) {
					if (paramInfo.chS == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_max[2] = &maxR4;
						p_max[3] = &maxR5;
						p_avg[2] = &avgR4;
						p_avg[3] = &avgR5;
						strcpy(fn[2],"_SL");
						strcpy(fn[3],"_SR");
					}
					if (paramInfo.chC == 1 && paramInfo.chLFE == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_max[2] = &maxR3;
						p_max[3] = &maxR6;
						p_avg[2] = &avgR3;
						p_avg[3] = &avgR6;
						strcpy(fn[2],"_C");
						strcpy(fn[3],"_LFE");
					}
				}
				if (paramInfo.ch == 5) {
					if (paramInfo.chC == 1 && paramInfo.chS == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_nx[4]  = &paramInfo.nx5;
						p_max[2] = &maxR3;
						p_max[3] = &maxR4;
						p_max[4] = &maxR5;
						p_avg[2] = &avgR3;
						p_avg[3] = &avgR4;
						p_avg[4] = &avgR5;
						strcpy(fn[2],"_C");
						strcpy(fn[3],"_SL");
						strcpy(fn[4],"_SR");
					}
					if (paramInfo.chLFE == 1 && paramInfo.chS == 1) {
						p_nx[2]  = &paramInfo.nx3;
						p_nx[3]  = &paramInfo.nx4;
						p_nx[4]  = &paramInfo.nx5;
						p_max[2] = &maxR6;
						p_max[3] = &maxR4;
						p_max[4] = &maxR5;
						p_avg[2] = &avgR6;
						p_avg[3] = &avgR4;
						p_avg[4] = &avgR5;
						strcpy(fn[2],"_LFE");
						strcpy(fn[3],"_SL");
						strcpy(fn[4],"_SR");
					}
				}
				if (paramInfo.ch == 6) {
					p_nx[2]  = &paramInfo.nx3;
					p_nx[5]  = &paramInfo.nx6;
					p_nx[3]  = &paramInfo.nx4;
					p_nx[4]  = &paramInfo.nx5;
					strcpy(fn[2],"_C");
					strcpy(fn[3],"_LFE");
					strcpy(fn[4],"_SL");
					strcpy(fn[5],"_SR");
				}

				_splitpath(argv[2],drive,dir,fname,ext);
				_makepath(tmppath,drive,dir,fname,"r1.param");
				fp = fopen(tmppath,"rb");
				if (fp == NULL) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				rd = fscanf(fp,"%lf,%llx",&per,&avg);
				if (rd != 2) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				fclose(fp);
				nx = 0;
				if (avg > 10000000000) {
					nx = (double)*p_avg[0] / avg;
					if (nx > 0) {
						max = (SSIZE)((double)*p_max[0] / nx);
					}
				} else {
					max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
					max = (SSIZE)((double)max * 1.02);
				}
				nx = 0;
				if (max > 0) {
					if (paramInfo.norm == 1) {
						nx = (double)0x7FFFFFFFFFFFFF / max;
					} else {
						nx = (double)*p_max[0] / max;
					}
					nx *= 0.98;
				}
				*p_nx[0] = nx;

				if (paramInfo.ch >= 2) {
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r2.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[1] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[1] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[1] / max;
						}
						nx *= 0.98;
					}
					*p_nx[1] = nx;
				}
				if (paramInfo.ch == 3) {
					if (paramInfo.chC == 1) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r3.param");
					} else if (paramInfo.chLFE) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r6.param");
					} else {	// not genCh
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r3.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;
				}
				if (paramInfo.ch == 4 && paramInfo.chS == 1) {
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r4.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r5.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
				} else if (paramInfo.ch == 4 && paramInfo.chC == 1 && paramInfo.chLFE == 1) {
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r3.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r6.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
				} else if (paramInfo.ch == 4) {
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r3.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r4.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
				}

				if (paramInfo.ch == 5) {
					if (paramInfo.chC == 1) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r3.param");
					} else if (paramInfo.chLFE) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r6.param");
					} else {	// not genCh
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r3.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r4.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r5.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[4] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[4] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[4] / max;
						}
						nx *= 0.98;
					}
					*p_nx[4] = nx;
				}

				if (paramInfo.ch == 6) {
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r3.param");
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[2] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[2] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[2] / max;
						}
						nx *= 0.98;
					}
					*p_nx[2] = nx;

					if (genCh > 0) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r6.param");
					} else {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r4.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[3] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[3] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[3] / max;
						}
						nx *= 0.98;
					}
					*p_nx[3] = nx;

					if (genCh > 0) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r4.param");
					} else {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r5.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[4] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[4] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[4] / max;
						}
						nx *= 0.98;
					}
					*p_nx[4] = nx;

					if (genCh > 0) {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r5.param");
					} else {
						_splitpath(argv[2],drive,dir,fname,ext);
						_makepath(tmppath,drive,dir,fname,"r6.param");
					}
					fp = fopen(tmppath,"rb");
					if (fp == NULL) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					rd = fscanf(fp,"%lf,%llx",&per,&avg);
					if (rd != 2) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fclose(fp);
					nx = 0;
					if (avg > 10000000000) {
						nx = (double)*p_avg[5] / avg;
						if (nx > 0) {
							max = (SSIZE)((double)*p_max[5] / nx);
						}
					} else {
						max = (SSIZE)(per * 0x7FFFFFFFFFFFFF);
						max = (SSIZE)((double)max * 1.02);
					}
					nx = 0;
					if (max > 0) {
						if (paramInfo.norm == 1) {
							nx = (double)0x7FFFFFFFFFFFFF / max;
						} else {
							nx = (double)*p_max[5] / max;
						}
						nx *= 0.98;
					}
					*p_nx[5] = nx;
				}

sprintf(s,"nx0:%lf",*p_nx[0]);
PRINT_LOG(s);
sprintf(s,"nx1:%lf",*p_nx[1]);
PRINT_LOG(s);
sprintf(s,"nx2:%lf",*p_nx[2]);
PRINT_LOG(s);
sprintf(s,"nx3:%lf",*p_nx[3]);
PRINT_LOG(s);
sprintf(s,"nx4:%lf",*p_nx[4]);
PRINT_LOG(s);
sprintf(s,"nx5:%lf",*p_nx[5]);
PRINT_LOG(s);

PRINT_LOG("");
				// R1 ファイル
				_splitpath(argv[2],drive,dir,fname,ext);
				_makepath(tmppath,drive,dir,fname,"r1");

				// ファイルオープン
				fio_open(&fp_r1,tmppath,FIO_MODE_R);
				if (fp_r1.error) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				p_fp1 = &fp_r1;
				if (paramInfo.ch >= 2) {
					// R2 ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r2");

					// ファイルオープン
					fio_open(&fp_r2,tmppath,FIO_MODE_R);
					if (fp_r2.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp2 = &fp_r2;
				}

				if (paramInfo.ch >= 3 && (paramInfo.chC == 1 || genCh == 0)) {
					// R3 ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r3");

					// ファイルオープン
					fio_open(&fp_r3,tmppath,FIO_MODE_R);
					if (fp_r3.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp3 = &fp_r3;
				}

				if (paramInfo.ch >= 4 && (paramInfo.chS == 1 || genCh == 0)) {
					// R4 ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r4");

					// ファイルオープン
					fio_open(&fp_r4,tmppath,FIO_MODE_R);
					if (fp_r4.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp4 = &fp_r4;

					// R5 ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r5");

					// ファイルオープン
					fio_open(&fp_r5,tmppath,FIO_MODE_R);
					if (fp_r5.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp5 = &fp_r5;

				}

				if ((paramInfo.ch >= 3 && paramInfo.chLFE == 1) || (paramInfo.ch == 6 && genCh == 0)) {
					// R6 ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					_makepath(tmppath,drive,dir,fname,"r6");

					// ファイルオープン
					fio_open(&fp_r6,tmppath,FIO_MODE_R);
					if (fp_r6.error) {
						retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					p_fp6 = &fp_r6;
				}

				outSample = 0;
				fio_get_filesize(&fp_r1,&size);
				if (fp_r1.error || size < sizeof (SSIZE)) {
					retCode = STATUS_FILE_READ_ERR;paramInfo.errLine = __LINE__;
					break;
				}
				outSample = size / sizeof (SSIZE);

				fprintf(stdout,"NORM\n");
				fflush(stdout);
				if (paramInfo.split == 0) {
					// 出力ファイル
					_splitpath(argv[2],drive,dir,fname,ext);
					if (paramInfo.raw) {
						_makepath(tmppath,drive,dir,fname,"raw");
					} else {
						_makepath(tmppath,drive,dir,fname,"wav");
					}
					strcpy(paramInfo.tofile,fname);
					strcat(paramInfo.tofile,".wav");
					// ファイルオープン
					fio_open(&fp_w,tmppath,FIO_MODE_W);
					if (fp_w.error) {
						retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
						break;
					}
					fio_set_memory_limit(&fp_w,paramInfo.fio);
					paramInfo.fp_w[0] = &fp_w;
					for (count = 1,startSample = 0;count > 0;) {
						retCode = Normalize(&count,&inFmt,&outFmt,p_fp1,p_fp2,p_fp3,p_fp4,p_fp5,p_fp6,&startSample,outSample,&paramInfo);
						fio_close(paramInfo.fp_w[0]);
						if (paramInfo.fp_w[0]->error) {
							break;
						}
						if (count == 0) {
							break;
						}
						count++;
						// 出力ファイル
						_splitpath(argv[2],drive,dir,fname,ext);
						sprintf(fname,"%s_%d.wav",fname,count);
						_makepath(tmppath,drive,dir,fname,NULL);
						// ファイルオープン
						memset(&fp_w,0,sizeof (FIO));
						fio_open(&fp_w,tmppath,FIO_MODE_W);
						if (fp_w.error) {
							retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
							break;
						}
						fio_set_memory_limit(&fp_w,paramInfo.fio);
						strcpy(paramInfo.tofile,fname);
						paramInfo.fp_w[0] = &fp_w;
					}
				} else {
					// 出力ファイル
					for (ccc = 0;ccc < paramInfo.ch;ccc++) {
						_splitpath(argv[2],drive,dir,fname,ext);
						strcat(fname,fn[ccc]);
						if (paramInfo.raw) {
							_makepath(tmppath,drive,dir,fname,"raw");
						} else {
							_makepath(tmppath,drive,dir,fname,"wav");
						}
						if (ccc == 0) {
							strcpy(paramInfo.tofile,fname);
							strcat(paramInfo.tofile,"wav");
						}
						// ファイルオープン
						fio_open(&fp_ws[ccc],tmppath,FIO_MODE_W);
						if (fp_ws[ccc].error) {
							retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
							break;
						}
						fio_set_memory_limit(&fp_ws[ccc],paramInfo.fio);
						paramInfo.fp_w[ccc] = &fp_ws[ccc];
					}
					for (count = 1,startSample = 0;count > 0;) {
						retCode = Normalize(&count,&inFmt,&outFmt,p_fp1,p_fp2,p_fp3,p_fp4,p_fp5,p_fp6,&startSample,outSample,&paramInfo);
						if (retCode != STATUS_SUCCESS) {
							break;
						}
						for (ccc = 0;ccc < paramInfo.ch;ccc++) {
							fio_close(paramInfo.fp_w[ccc]);
							if (paramInfo.fp_w[ccc]->error) {
								retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
								count = 0;
								break;
							}
						}
						if (count == 0) {
							break;
						}
						count++;
						for (ccc = 0;ccc < paramInfo.ch;ccc++) {
							// 出力ファイル
							_splitpath(argv[2],drive,dir,fname,ext);
							sprintf(fname,"%s_%d%s",fname,count,fn[ccc]);
							_makepath(tmppath,drive,dir,fname,"wav");
							// ファイルオープン
							fio_open(&fp_w,tmppath,FIO_MODE_W);
							if (fp_w.error) {
								retCode = STATUS_FILE_WRITE_ERR;paramInfo.errLine = __LINE__;
								break;
							}
							fio_set_memory_limit(&fp_w,paramInfo.fio);
							strcpy(paramInfo.tofile,fname);
							paramInfo.fp_w[ccc] = &fp_w;
							memset(&fp_w,0,sizeof (FIO));
						}
					}
				}
				fio_close(&fp_r1);
				if (p_fp2 != NULL) {
					fio_close(&fp_r2);
				}
				if (p_fp3 != NULL) {
					fio_close(&fp_r3);
				}
				if (p_fp4 != NULL) {
					fio_close(&fp_r4);
				}
				if (p_fp5 != NULL) {
					fio_close(&fp_r5);
				}
				if (p_fp6 != NULL) {
					fio_close(&fp_r6);
				}
			}
		}
		if (!(argc == 3 || argc == 4) || paramFlag != 0x3F) {
			printf(STR_COPYRIGHT);
			printf(STR_USAGE);
			exit(0);
		}
	} while (0);

	if (retCode != STATUS_SUCCESS) {
		_splitpath(argv[2],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"err");
		fp = fopen(tmppath,"a");
		if (fp) {
			switch (retCode) {
				case STATUS_FILE_READ_ERR:
					fprintf(fp,"raw2wav:[%04d] File read error.\n",paramInfo.errLine);
					break;
				case STATUS_FILE_WRITE_ERR:
					fprintf(fp,"raw2wav:[%04d] File write error.\n",paramInfo.errLine);
					break;
				case STATUS_MEM_ALLOC_ERR:
					fprintf(fp,"raw2wav:[%04d] Memory Allocation error.\n",paramInfo.errLine);
					break;
				default:
					fprintf(fp,"raw2wav:[%04d] Other error.\n",paramInfo.errLine);
			}
			fclose(fp);
		}
	}

	exit(retCode);
}

//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理
// ---
//	pcount	:ファイル番号のアドレス
//	inFmt	:入力ファイル音声形式情報
//	outFmt	:出力ファイル音声形式情報
//	fp_r1	:音声データのFIO構造体
//	fp_r2	:音声データのFIO構造体
//	fp_r3	:音声データのFIO構造体
//	fp_r4	:音声データのFIO構造体
//	fp_r5	:音声データのFIO構造体
//	fp_r6	:音声データのFIO構造体
//	startSample : 開始サンプル
//	nSample :処理をするサンプル数のアドレス
//	param	:パラメーター構造体
//
int Normalize(int *pCount,SOUNDFMT *inFmt,SOUNDFMT *outFmt,FIO *fp_r1,FIO *fp_r2,FIO *fp_r3,FIO *fp_r4,FIO *fp_r5,FIO *fp_r6,DWORD *startSample,DWORD nSample,PARAM_INFO *param)
{
	DWORD i,j,outSample;
	DWORD chMask;
	BYTE *header;
	BYTE *bwf_chunk;
	BYTE *buf1,*buf2,*buf3,*buf4,*buf5,*buf6;
	FIO *p_fp1,*p_fp2,*p_fp3,*p_fp4,*p_fp5,*p_fp6;
	long header_size;
	fio_size write_size,ws,data_size,file_size;
	fio_size data_start,data_end;
	BROADCAST_EXT *bext;
	BYTE *chunk_data;
	long  chunk_size;
	int retCode,retCode1,retCode2,retCode3,retCode4,retCode5,retCode6;
	int nChunk;
	int bwf_enable;
	long bwf_size;
	int div_flag;
	double persent,per;
	char ppp[50];
	int ccc;
	
	do {
		div_flag = 0;
		header = NULL;
		buf1 = NULL;
		buf2 = NULL;
		buf3 = NULL;
		buf4 = NULL;
		buf5 = NULL;
		buf6 = NULL;
		p_fp1 = p_fp2 = p_fp3 = p_fp4 = p_fp5 = p_fp6 = NULL;

		retCode1 = retCode2 = retCode3 = retCode4 = retCode5 = retCode6 = 0;
		retCode = STATUS_SUCCESS;
		// ヘッダ処理
		if (param->raw == 0) {
			header = malloc(1 * 1024 * 1024);
			if (header == NULL) {
				retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
				break;
			}
			if (param->split == 0) {
				data_size = (param->sampling * (param->bitwidth / 8)) * param->ch;
				data_size *= 10;	// 10秒
			} else {
				data_size = (param->sampling * (param->bitwidth / 8)) * 1;
				data_size *= 10;	// 10秒
			}
			buf1 = malloc(data_size);
			if (buf1 == NULL) {
				retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
				break;
			}
			if (param->ch >= 2) {
				buf2 = malloc(data_size);
				if (buf2 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 3) {
				buf3 = malloc(data_size);
				if (buf3 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 4) {
				buf4 = malloc(data_size);
				if (buf4 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 5) {
				buf5 = malloc(data_size);
				if (buf5 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->ch >= 6) {
				buf6 = malloc(data_size);
				if (buf6 == NULL) {
					retCode = STATUS_MEM_ALLOC_ERR;param->errLine = __LINE__;
					break;
				}
			}
			if (param->rf64 == 0) {
				retCode = PLG_MakeHeaderWAV(inFmt,outFmt,header,1 * 1024 * 1024,&header_size);
				if (retCode != STATUS_SUCCESS) {
					break;
				}
				if (param->ch >= 3 && param->split == 0) {
					chMask	= 0;
					chMask |= SPEAKER_FRONT_LEFT;
					chMask |= SPEAKER_FRONT_RIGHT;
					if (param->chC == 1) {
						chMask |= SPEAKER_FRONT_CENTER;
					}
					if (param->chS == 1) {
						chMask |= SPEAKER_BACK_LEFT;
						chMask |= SPEAKER_BACK_RIGHT;
					}
					if (param->chLFE == 1) {
						chMask |= SPEAKER_LOW_FREQUENCY;
					}
					(*(DWORD *)(&header[40])) = chMask;
				}
			} else {
				retCode = PLG_MakeHeaderRF64(inFmt,outFmt,header,1 * 1024 * 1024,&header_size);
				if (retCode != STATUS_SUCCESS) {
					break;
				}
			}
			if (param->split == 0) {
				ws = fio_write(header,1,header_size,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != header_size) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
				fio_get_datasize(param->fp_w[0],&data_start);
			} else {
				for (ccc = 0;ccc < param->ch;ccc++) {
					ws = fio_write(header,1,header_size,param->fp_w[ccc]);
					if (param->fp_w[ccc]->error || ws != header_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
					fio_get_datasize(param->fp_w[ccc],&data_start);
				}
				if (retCode != STATUS_SUCCESS) {
					break;
				}
			}
		}
		p_fp1 = fp_r1;
		p_fp2 = fp_r2;
		if (param->ch == 3) {
			if (param->chC == 1) {
				p_fp3 = fp_r3;
			} else if (param->chLFE == 1) {
				p_fp3 = fp_r6;
			}
		}
		if (param->ch == 4) {
			if (param->chS == 1) {
				p_fp3 = fp_r4;
				p_fp4 = fp_r5;
			}
			if (param->chC == 1 && param->chLFE == 1) {
				p_fp3 = fp_r3;
				p_fp4 = fp_r6;
			}
		}
		if (param->ch == 5) {
			if (param->chC == 1 && param->chS == 1) {
				p_fp3 = fp_r3;
				p_fp4 = fp_r4;
				p_fp5 = fp_r5;
			}
			if (param->chLFE == 1 && param->chS == 1) {
				p_fp3 = fp_r6;
				p_fp4 = fp_r4;
				p_fp5 = fp_r5;
			}
		}
		if (param->ch == 6) {
			p_fp3 = fp_r3;
			p_fp4 = fp_r6;
			p_fp5 = fp_r4;
			p_fp6 = fp_r5;
		}
		per = -1;
		for (i = *startSample;i < nSample;i += param->sampling * 10) {
			persent = ((double)i / nSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;

			// 10秒ずつ処理する
			outSample = param->sampling * 10;
			if (i + outSample > nSample) {
				outSample = nSample - i;
			}
			if (buf1) memset(buf1,0,data_size);
			if (buf2) memset(buf2,0,data_size);
			if (buf3) memset(buf3,0,data_size);
			if (buf4) memset(buf4,0,data_size);
			if (buf5) memset(buf5,0,data_size);
			if (buf6) memset(buf6,0,data_size);

			// データ順
			// Left,Right,Center,LFE,Back Left,Back Right

			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					{
						retCode1 = Normalize_Mx(param->ch,0,param->bitwidth,p_fp1,outSample,buf1,param);
					}
					#pragma omp section
					{
						if (param->ch >= 2) {
							retCode2 = Normalize_Mx(param->ch,1,param->bitwidth,p_fp2,outSample,buf2,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 3) {
							retCode3 = Normalize_Mx(param->ch,2,param->bitwidth,p_fp3,outSample,buf3,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 4) {
							retCode4 = Normalize_Mx(param->ch,3,param->bitwidth,p_fp4,outSample,buf4,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 5) {
							retCode5 = Normalize_Mx(param->ch,4,param->bitwidth,p_fp5,outSample,buf5,param);
						}
					}
					#pragma omp section
					{
						if (param->ch >= 6) {
							retCode6 = Normalize_Mx(param->ch,5,param->bitwidth,p_fp6,outSample,buf6,param);
						}
					}
				}
			}
			retCode = retCode1 | retCode2 | retCode3 | retCode4 | retCode5 | retCode6;
			if (retCode != STATUS_SUCCESS) {
				break;
			}
			if (param->split == 0) {
				if (param->ch == 2) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j];
					}
				} else if (param->ch == 3) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j];
					}
				} else if (param->ch == 4) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j] | buf4[j];
					}
				} else if (param->ch == 5) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j] | buf4[j] | buf5[j];
					}
				} else if (param->ch == 6) {
					#pragma omp parallel for
					for (j = 0;j < data_size;j++) {
						buf1[j] = buf1[j] | buf2[j] | buf3[j] | buf4[j] | buf5[j] | buf6[j];
					}
				}
				if (outSample == param->sampling * 10) {
					ws = fio_write(buf1,1,data_size,param->fp_w[0]);
					if (param->fp_w[0]->error || ws != data_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
				} else {
					write_size = (outSample * (param->bitwidth / 8)) * param->ch;
					ws = fio_write(buf1,1,write_size,param->fp_w[0]);
					if (param->fp_w[0]->error || ws != write_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
				}
			} else {
				if (outSample == param->sampling * 10) {
					BYTE *pb[6];
					pb[0] = buf1;
					pb[1] = buf2;
					pb[2] = buf3;
					pb[3] = buf4;
					pb[4] = buf5;
					pb[5] = buf6;
					for (ccc = 0;ccc < param->ch;ccc++) {
						ws = fio_write(pb[ccc],1,data_size,param->fp_w[ccc]);
						if (param->fp_w[ccc]->error || ws != data_size) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					}
					if (retCode != STATUS_SUCCESS) {
						break;
					}
				} else {
					BYTE *pb[6];
					pb[0] = buf1;
					pb[1] = buf2;
					pb[2] = buf3;
					pb[3] = buf4;
					pb[4] = buf5;
					pb[5] = buf6;
					for (ccc = 0;ccc < param->ch;ccc++) {
						write_size = (outSample * (param->bitwidth / 8)) * 1;
						ws = fio_write(pb[ccc],1,write_size,param->fp_w[ccc]);
						if (param->fp_w[ccc]->error || ws != write_size) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					}
					if (retCode != STATUS_SUCCESS) {
						break;
					}
				}
			}
			/*if (param->rf64 == 0) {
				fio_get_datasize(param->fp_w[0],&file_size);
				if (file_size > (fio_size)1 * 1024 * 1024 * 1024 && outSample == param->sampling * 10 && (nSample - i) >= param->sampling * 20) {
					// データサイズが大きいので分割する(rf64以外)
					*startSample = i + outSample;
					div_flag = 1;
					break;
				}
			}
			*/
		}
		if (retCode != STATUS_SUCCESS) {
			break;
		}
		fio_get_datasize(param->fp_w[0],&data_end);
		bwf_enable = 0;
		if (param->bwf) {
			bwf_enable = 1;
		}
		if (param->raw == 0) {
			chunk_data = NULL;
			for (nChunk = 1;;) {
				if (param->flac == 0) {
					PLG_GetExtraChunk(param->fromfile,nChunk,&chunk_data,&chunk_size);
				} else {
					PLG_GetExtraChunk(param->flacfile,nChunk,&chunk_data,&chunk_size);
				}
				if (chunk_size == 0) {
					break;
				}
				if (chunk_data[0] == 'b' && chunk_data[1] == 'e' && chunk_data[2] == 'x' && chunk_data[3] == 't') {
					bwf_enable = 0;
					if (param->bwf) {
						int alloc_flag = 0;
						bwf_size = chunk_data[7];bwf_size <<= 8;
						bwf_size |= chunk_data[6];bwf_size <<= 8;
						bwf_size |= chunk_data[5];bwf_size <<= 8;
						bwf_size |= chunk_data[4];
						if (bwf_size <= sizeof (BROADCAST_EXT) + 128) {
							alloc_flag = 1;
						} else {
							bext = (BROADCAST_EXT *)&chunk_data[8];
							if (chunk_size < sizeof (BROADCAST_EXT) + strlen(bext->codingHistory) + 128) {
								alloc_flag = 1;
							}
						}
						if (alloc_flag == 1) {
							bwf_size += 256;
							bext = (BROADCAST_EXT *)malloc(bwf_size);
						}
						if (bext) {
							if (alloc_flag == 1) {
								memset(bext,0,bwf_size);
								memcpy(bext,chunk_data + 8,chunk_size);
							}
							UpdateBext(bext,inFmt,outFmt,param,bwf_size);
							if (param->split == 0) {
								ws = fio_write("bext",1,4,param->fp_w[0]);
								if (param->fp_w[0]->error || ws != 4) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
								ws = fio_write(&bwf_size,1,4,param->fp_w[0]);
								if (param->fp_w[0]->error || ws != 4) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
								ws = fio_write(bext,1,bwf_size,param->fp_w[0]);
								if (param->fp_w[0]->error || ws != bwf_size) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
							} else {
								for (ccc = 0;ccc < param->ch;ccc++) {
									ws = fio_write("bext",1,4,param->fp_w[ccc]);
									if (param->fp_w[ccc]->error || ws != 4) {
										retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
										break;
									}
									ws = fio_write(&bwf_size,1,4,param->fp_w[ccc]);
									if (param->fp_w[ccc]->error || ws != 4) {
										retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
										break;
									}
									ws = fio_write(bext,1,bwf_size,param->fp_w[ccc]);
									if (param->fp_w[ccc]->error || ws != bwf_size) {
										retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
										break;
									}
								}
								if (retCode != STATUS_SUCCESS) {
									break;
								}
							}
							if (alloc_flag) {
								free(bext);
							}
						}
					}
				} else {
					if (param->split == 0) {
						ws = fio_write(chunk_data,1,chunk_size,param->fp_w[0]);
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(chunk_data,1,chunk_size,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != chunk_size) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
				}
				nChunk++;
			}
			if (bwf_enable) {
				bwf_chunk = malloc(512 + 256);
				if (bwf_chunk) {
					memset(bwf_chunk,0,512 + 256);
					bwf_chunk[0] = 'b';
					bwf_chunk[1] = 'e';
					bwf_chunk[2] = 'x';
					bwf_chunk[3] = 't';
					bwf_chunk[4] = (BYTE)((512 + 256) - 8);
					bwf_chunk[5] = (BYTE)(((512 + 256) - 8) >> 8);
					bwf_chunk[6] = (BYTE)(((512 + 256) - 8) >> 16);
					bwf_chunk[7] = (BYTE)(((512 + 256) - 8) >> 24);
					bext = (BROADCAST_EXT *)&bwf_chunk[8];
					UpdateBext(bext,inFmt,outFmt,param,512 + 256);
					if (param->split == 0) {
						ws = fio_write(bwf_chunk,1,512 + 256,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != 512 + 256) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(bwf_chunk,1,512 + 256,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != 512 + 256) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
					free(bwf_chunk);
					bwf_chunk = NULL;
				}
			}
			if (param->bwf && (div_flag == 1 || *pCount > 1)) {
				// bwf の link チャンク
				char *link,*wk_str;
				long link_size;
				link = (char *)malloc(strlen(param->tofile) + 128);
				wk_str = malloc(_MAX_PATH + 128);
				if (link != NULL && wk_str != NULL) {
					link[0] = '\0';
					strcat(link,link_start);
					if (*pCount == 1) {
						sprintf(wk_str,link_file,"actual",(int)*pCount,param->tofile);
					} else {
						sprintf(wk_str,link_file,"other",(int)*pCount,param->tofile);
					}
					strcat(link,wk_str);
					strcat(link,link_end);
					if (param->split == 0) {
						ws = fio_write("link",1,4,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != 4) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write("link",1,4,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != 4) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
					link_size = strlen(link);
					if (param->split == 0) {
						ws = fio_write(&link_size,1,4,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != 4) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(&link_size,1,4,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != 4) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
					}
					if (param->split == 0) {
						ws = fio_write(link,1,link_size,param->fp_w[0]);
						if (param->fp_w[0]->error || ws != link_size) {
							retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
							break;
						}
					} else {
						for (ccc = 0;ccc < param->ch;ccc++) {
							ws = fio_write(link,1,link_size,param->fp_w[ccc]);
							if (param->fp_w[ccc]->error || ws != link_size) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						}
						if (retCode != STATUS_SUCCESS) {
							break;
						}
					}
					free(link);
					free(wk_str);
					if (link_size & 0x01) {
						// padding
						if (param->split == 0) {
							ws = fio_write("",1,1,param->fp_w[0]);
							if (param->fp_w[0]->error || ws != 1) {
								retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
								break;
							}
						} else {
							for (ccc = 0;ccc < param->ch;ccc++) {
								ws = fio_write("",1,1,param->fp_w[ccc]);
								if (param->fp_w[ccc]->error || ws != 1) {
									retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
									break;
								}
							}
						}
					}
				}
			}
			if (param->split == 0) {
				fio_get_datasize(param->fp_w[0],&data_size);
				fio_flush(param->fp_w[0]);
				if (param->fp_w[0]->error) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
				fio_rewind(param->fp_w[0]);
				if (param->rf64 == 0) {
					retCode = PLG_UpdateHeaderWAV(outFmt,data_size,data_end - data_start,header,header_size);
				} else {
					retCode = PLG_UpdateHeaderRF64(outFmt,data_size,data_end - data_start,header,header_size);
				}
				ws = fio_write(header,1,header_size,param->fp_w[0]);
				if (param->fp_w[0]->error || ws != header_size) {
					retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
					break;
				}
			} else {
				fio_get_datasize(param->fp_w[0],&data_size);
				for (ccc = 0;ccc < param->ch;ccc++) {
					fio_flush(param->fp_w[ccc]);
					if (param->fp_w[ccc]->error) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
					fio_rewind(param->fp_w[ccc]);
					if (param->rf64 == 0) {
						retCode = PLG_UpdateHeaderWAV(outFmt,data_size,data_end - data_start,header,header_size);
					} else {
						retCode = PLG_UpdateHeaderRF64(outFmt,data_size,data_end - data_start,header,header_size);
					}
					ws = fio_write(header,1,header_size,param->fp_w[ccc]);
					if (param->fp_w[ccc]->error || ws != header_size) {
						retCode = STATUS_FILE_WRITE_ERR;param->errLine = __LINE__;
						break;
					}
				}
			}
		}
	} while (0);
	if (div_flag == 0) {
		*pCount = 0;
	}
	return retCode;
}
//---------------------------------------------------------------------------
// Function   : Normalize_Mx
// Description: ノーマライズ処理
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_Mx(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	int retCode;

	// 音声データ出力
	switch (param->mode) {
		case 0:
			retCode = Normalize_M0(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		case 1:
			retCode = Normalize_M1(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		case 2:
			retCode = Normalize_M2(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		case 3:
			retCode = Normalize_M3(nCh,ch,bit,fp_r,nSample,buffer,param);
			break;
		default:
			retCode = -1;
			break;
	}
	return retCode;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(カット)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M0(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	DWORD i,w_off;
	int next;
	double nx;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	char p[50];
	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
	} else if (ch == 1) {
		nx = param->nx2;
	} else if (ch == 2) {
		nx = param->nx3;
	} else if (ch == 3) {
		nx = param->nx4;
	} else if (ch == 4) {
		nx = param->nx5;
	} else if (ch == 5) {
		nx = param->nx6;
	}
	for (i = 0;i < nSample;i++) {
		rs = fio_read(&s,sizeof (SSIZE),1,fp_r);
		if (fp_r->error || rs != 1) {
			param->errLine = __LINE__;
			return STATUS_FILE_READ_ERR;
		}
		s = CLIP_NX(s,nx);
		
		if (bit == 16) {
			s = CLIP_MAX(s);
			s >>= (64 - 16);
			s_os = (short *)&buffer[w_off];
			*s_os = (short)s;
			w_off += next;
		} else if (bit == 24) {
			s = CLIP_MAX(s);
			s >>= (64 - 24);
			buffer[w_off + 0] = (unsigned char)s;
			buffer[w_off + 1] = (unsigned char)(s >> 8);
			buffer[w_off + 2] = (unsigned char)(s >> 16);
			w_off += next;
		} else if (bit == 32) {
			f_os = (float *)&buffer[w_off];
			s = CLIP_MAX(s);
			*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		} else if (bit == 64) {
			d_os = (double *)&buffer[w_off];
			s = CLIP_MAX(s);
			*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		}
	}
	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(ディザ)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M1(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	DWORD i;
	int ignore_s;
	double nx;
	double noise;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	int next;
	int w_off;
	
	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
	} else if (ch == 1) {
		nx = param->nx2;
	} else if (ch == 2) {
		nx = param->nx3;
	} else if (ch == 3) {
		nx = param->nx4;
	} else if (ch == 4) {
		nx = param->nx5;
	} else if (ch == 5) {
		nx = param->nx6;
	}

	for (i = 0;i < nSample;i++) {
		rs = fio_read(&s,sizeof (SSIZE),1,fp_r);
		if (fp_r->error || rs != 1) {
			param->errLine = __LINE__;
			return STATUS_FILE_READ_ERR;
		}

		s = CLIP_NX(s,nx);

		noise = normalNoise();
		if (bit == 16) {
			if (param->ditherLv > 0) {
				noise *= (0x10000000000 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-16));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			s >>= (64 - 16);
			s_os = (short *)&buffer[w_off];
			*s_os = (short)s;
			w_off += next;

		} else if (bit == 24) {
			if (param->ditherLv > 0) {
				noise *= (0x100000000 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-24));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			s >>= (64 - 24);
			buffer[w_off + 0] = (unsigned char)s;
			buffer[w_off + 1] = (unsigned char)(s >> 8);
			buffer[w_off + 2] = (unsigned char)(s >> 16);
			w_off += next;
		} else if (bit == 32) {
			if (param->ditherLv > 0) {
				noise *= (0x1000000 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-32));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			f_os = (float *)&buffer[w_off];
			*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		} else if (bit == 64) {
			if (param->ditherLv > 0) {
				noise *= (0x100 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-48));
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX(s);

			d_os = (double *)&buffer[w_off];
			*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		}
	}

	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(誤差累積)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M2(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	SSIZE old_s;
	SSIZE sum,sum_2nd,sum_3rd,val;
	DWORD i;
	int ignore_s;
	double nx;
	double noise;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	int next;
	int w_off;
	char p[100];


	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
		sum = param->sum1;
		sum_2nd = param->sum1_2nd;
		sum_3rd = param->sum1_3rd;
		old_s	= param->old_s1;
	} else if (ch == 1) {
		nx = param->nx2;
		sum = param->sum2;
		sum_2nd = param->sum2_2nd;
		sum_3rd = param->sum2_3rd;
		old_s	= param->old_s2;
	} else if (ch == 2) {
		nx = param->nx3;
		sum = param->sum3;
		sum_2nd = param->sum3_2nd;
		sum_3rd = param->sum3_3rd;
		old_s	= param->old_s3;
	} else if (ch == 3) {
		nx = param->nx4;
		sum = param->sum4;
		sum_2nd = param->sum4_2nd;
		sum_3rd = param->sum4_3rd;
		old_s	= param->old_s4;
	} else if (ch == 4) {
		nx = param->nx5;
		sum = param->sum5;
		sum_2nd = param->sum5_2nd;
		sum_3rd = param->sum5_3rd;
		old_s	= param->old_s5;
	} else if (ch == 5) {
		nx = param->nx6;
		sum = param->sum6;
		sum_2nd = param->sum6_2nd;
		sum_3rd = param->sum6_3rd;
		old_s	= param->old_s6;
	}

	for (i = 0;i < nSample;i++) {
		rs = fio_read(&s,sizeof (SSIZE),1,fp_r);
		if (fp_r->error || rs != 1) {
			param->errLine = __LINE__;
			return STATUS_FILE_READ_ERR;
		}

		s = CLIP_NX(s,nx);

		noise = normalNoise();
		if (bit == 16) {
			// 四捨五入

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-16));
			
			s -= old_s;
			sum += s;

			s = sum >> ((64-8) - 16);
			s <<= ((64-8) - 16);
			s -= old_s;
			sum_2nd += s;
			s = sum_2nd >> ((64-8) - 16);

			s <<= ((64-8) - 16);
			s -= old_s;
			sum_3rd += s;
			s = sum_3rd >> ((64-8) - 16);

			old_s = s << ((64-8) - 16);

			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,16);
			s_os = (short *)&buffer[w_off];
			*s_os = (short)s;
			w_off += next;
		} else if (bit == 24) {

			s = CLIP_ADD(s,ROUND_NBIT((64-8)-24));
			s -= old_s;
			sum += s;

			s = sum >> ((64-8) - 24);
			s <<= ((64-8) - 24);
			s -= old_s;
			sum_2nd += s;
			s = sum_2nd >> ((64-8) - 24);

			s <<= ((64-8) - 24);
			s -= old_s;
			sum_3rd += s;
			s = sum_3rd >> ((64-8) - 24);

			old_s = s << ((64-8) - 24);

			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,24);
			buffer[w_off + 0] = (unsigned char)s;
			buffer[w_off + 1] = (unsigned char)(s >> 8);
			buffer[w_off + 2] = (unsigned char)(s >> 16);
			w_off += next;
		} else if (bit == 32) {
			s = CLIP_ADD(s,ROUND_NBIT((64-8)-32));

			// 誤差の累積
			sum += s % 0x1000000;

			// 誤差の判定
			val = sum / 0x1000000;
			if (ignore_s == 0 && val > 0) {
				s += 0x1000000;
				sum %= 0x1000000;
			} else if (ignore_s == 0 && val < 0) {
				s -= 0x1000000;
				sum %= 0x1000000;
			}

			// 2次 誤差の累積
			sum_2nd += s % 0x1000000;
			// 誤差の判定
			val = sum_2nd / 0x1000000;
			if (val > 0) {
				s += 0x1000000;
				sum_2nd %= 0x1000000;
			} else if (val < 0) {
				s -= 0x1000000;
				sum_2nd %= 0x1000000;
			}

			// 32bit 化する
			s >>= ((64-8) - 32);
			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}

			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,32);
			s <<= (64 - 32);
			f_os = (float *)&buffer[w_off];
			*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
			w_off += next;
		} else if (bit == 64) {
			// 48bit
			s = CLIP_ADD(s,ROUND_NBIT((64-8)-48));

			// 誤差の累積
			sum += s % 0x10000;

			// 誤差の判定
			val = sum / 0x10000;
			if (val > 0) {
				s += 0x10000;
				sum %= 0x10000;
			} else if (val < 0) {
				s -= 0x10000;
				sum %= 0x10000;
			}

			// 2次 誤差の累積
			sum_2nd += s % 0x10000;
			// 誤差の判定
			val = sum_2nd / 0x10000;
			if (val > 0) {
				s += 0x10000;
				sum_2nd %= 0x10000;
			} else if (val < 0) {
				s -= 0x10000;
				sum_2nd %= 0x10000;
			}

			// 64bit 化する
			s >>= ((64-8) - 48);
			if (param->ditherLv > 0) {
				noise *= (1 << (param->ditherLv - 1));
			} else {
				noise = 0;
			}
			s = CLIP_ADD(s,(SSIZE)noise);
			s = CLIP_MAX_N(s,48);
			s <<= (64 - 48);
			d_os = (double *)&buffer[w_off];
			
			if(s < 0x8000000000000000){
				*d_os = (double)((__float128)s / 9223372036854775808.0q);
			} else {
				unsigned long long temp_value = ( ( (~s) + 1 ) & 0xFFFFFFFFFFFFFFFF );
				*d_os = (double)( -1.0q * ( (__float128)temp_value / 9223372036854775808.0q ) );
			}
			
			//*d_os = (double)((__float128)s / (__float128)0x7FFFFFFFFFFFFFFF);

			//printf("%f\n", *d_os);
			
			w_off += next;
		}
	}
	if (ch == 0) {
		param->sum1 = sum;
		param->sum1_2nd = sum_2nd;
		param->sum1_3rd = sum_3rd;
		param->old_s1 = old_s;
	} else if (ch == 1) {
		param->sum2 = sum;
		param->sum2_2nd = sum_2nd;
		param->sum2_3rd = sum_3rd;
		param->old_s2 = old_s;
	} else if (ch == 2) {
		param->sum3 = sum;
		param->sum3_2nd = sum_2nd;
		param->sum3_3rd = sum_3rd;
		param->old_s3 = old_s;
	} else if (ch == 3) {
		param->sum4 = sum;
		param->sum4_2nd = sum_2nd;
		param->sum4_3rd = sum_3rd;
		param->old_s4 = old_s;
	} else if (ch == 4) {
		param->sum5 = sum;
		param->sum5_2nd = sum_2nd;
		param->sum5_3rd = sum_3rd;
		param->old_s5 = old_s;
	} else if (ch == 5) {
		param->sum6 = sum;
		param->sum6_2nd = sum_2nd;
		param->sum6_3rd = sum_3rd;
		param->old_s6 = old_s;
	}

	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : Normalize
// Description: ノーマライズ処理(誤差拡散)
// ---
//	nCh		:チャンネル数
//	ch		:チャンネルの何番目を処理するかの指定
//	bit 	:入力データのビット数
//	fp_r	:音声データのFIO構造体
//	nSample :処理をするサンプル数
//	buffer	:出力先のバッファ
//
int Normalize_M3(int nCh,int ch,int bit,FIO *fp_r,DWORD nSample,BYTE *buffer,PARAM_INFO *param)
{
	SSIZE s;
	SSIZE ss[3];
	SSIZE sd[3];
	SSIZE val;
	SSIZE *buf;
	DWORD pastSample;
	long i,remain;
	int ignore_s;
	double nx;
	double noise;
	short *s_os;
	float *f_os;
	double *d_os;
	fio_size rs;
	int next;
	int w_off;

	if (bit == 16) {
		next = 2 * nCh;
		w_off = 2 * ch;
		if (param->split) {
			next = 2;
			w_off = 0;
		}
	} else if (bit == 24) {
		next = 3 * nCh;
		w_off = 3 * ch;
		if (param->split) {
			next = 3;
			w_off = 0;
		}
	} else if (bit == 32) {
		next = 4 * nCh;
		w_off = 4 * ch;
		if (param->split) {
			next = 4;
			w_off = 0;
		}
	} else if (bit == 64) {
		next = 8 * nCh;
		w_off = 8 * ch;
		if (param->split) {
			next = 8;
			w_off = 0;
		}
	}

	if (ch == 0) {
		nx = param->nx1;
	} else if (ch == 1) {
		nx = param->nx2;
	} else if (ch == 2) {
		nx = param->nx3;
	} else if (ch == 3) {
		nx = param->nx4;
	} else if (ch == 4) {
		nx = param->nx5;
	} else if (ch == 5) {
		nx = param->nx6;
	}

	ss[0] = ss[1] = ss[2] = 0;
	sd[0] = sd[1] = sd[2] = 0;

	buf = malloc(nSample * sizeof (SSIZE));
	if (buf == NULL) {
		param->errLine = __LINE__;
		return STATUS_MEM_ALLOC_ERR;
	}
	
	pastSample = 0;
	do {
		remain = fio_read(buf,sizeof (SSIZE),nSample,fp_r);
		for (i = 1;i < remain + 1 && pastSample < nSample;i++,pastSample++) {
			noise = normalNoise();

			ss[0] = buf[i-1];
			if (i < remain) {
				ss[1] = buf[i];
			} else {
				ss[1] = ss[0];
			}
			if (i + 1 < remain) {
				ss[2] = buf[i+1];
			} else {
				ss[2] = ss[0];
			}

			ss[0] = CLIP_NX(ss[0],nx);
			ss[1] = CLIP_NX(ss[1],nx);
			ss[2] = CLIP_NX(ss[2],nx);

			if (bit == 16) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x10000000000;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFFFF) {
						ss[1] += 0x10000000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFFFF) {
						ss[1] -= 0x10000000000;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x10000000000;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFFFF) {
						ss[0] += 0x10000000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFFFF) {
						ss[0] -= 0x10000000000;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 16);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,16);
				s_os = (short *)&buffer[w_off];
				*s_os = (short)s;
				w_off += next;
			} else if (bit == 24) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x100000000;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFF) {
						ss[1] += 0x100000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFF) {
						ss[1] -= 0x100000000;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x100000000;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7FFFFFFF) {
						ss[0] += 0x100000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFFFF) {
						ss[0] -= 0x100000000;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 24);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,24);
				buffer[w_off + 0] = (unsigned char)s;
				buffer[w_off + 1] = (unsigned char)(s >> 8);
				buffer[w_off + 2] = (unsigned char)(s >> 16);
				w_off += next;
			} else if (bit == 32) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x1000000;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7FFFFF) {
						ss[1] += 0x1000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFF) {
						ss[1] -= 0x1000000;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x1000000;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7FFFFF) {
						ss[0] += 0x1000000;
					}
				} else {
					if ((val * -1) > 0x7FFFFF) {
						ss[0] -= 0x1000000;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 32);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,32);
				s <<= (64 - 32);
				f_os = (float *)&buffer[w_off];
				*f_os = (float)s / 0x7FFFFFFFFFFFFFFF;
				w_off += next;
			} else if (bit == 64) {
				ss[1] += sd[1];
				s = ss[1];
				val = ss[1] % 0x100;
				ss[1] = ss[1] - val;
				if (val > 0) {
					if (val > 0x7F) {
						ss[1] += 0x100;
					}
				} else {
					if ((val * -1) > 0x7F) {
						ss[1] -= 0x100;
					}
				}
				sd[1]  = ss[1] - s;
				sd[0] += sd[1] / 2;
				sd[2] += sd[1] / 2;
				sd[1] = 0;
				ss[0] += sd[0];
				s = ss[0];
				val = ss[0] % 0x100;
				ss[0] = ss[0] - val;
				if (val > 0) {
					if (val > 0x7F) {
						ss[0] += 0x100;
					}
				} else {
					if ((val * -1) > 0x7F) {
						ss[0] -= 0x100;
					}
				}
				sd[0] = ss[0] - s;
				sd[1] = sd[0] / 2;

				sd[0] = sd[1];
				sd[1] = sd[2];
				sd[2] = 0;
				s = ss[0];
				s >>= ((64-8) - 48);

				if (param->ditherLv > 0) {
					noise *= (1 << (param->ditherLv - 1));
				} else {
					noise = 0;
				}
				s = CLIP_ADD(s,(SSIZE)noise);
				s = CLIP_MAX_N(s,48);
				s <<= (64 - 48);
				d_os = (double *)&buffer[w_off];
				*d_os = (double)s / 0x7FFFFFFFFFFFFFFF;
				w_off += next;
			}
		}
	} while (pastSample < nSample);
	free(buf);

	return STATUS_SUCCESS;
}

//---------------------------------------------------------------------------
// Function   : normalNoise
// Description: 正規乱数生成
//
double normalNoise(void)
/*
 * 正規乱数
 */
{
	double x1,x2;

	x1 = (double)rand() / RAND_MAX;
	x1 = 0.99999 * x1 + 0.00001;
	x2 = (double)rand() / RAND_MAX;
	return sqrt(-log(x1)) * cos(2.0 * 3.1415926 * x2) / 3.5;
}
//---------------------------------------------------------------------------
// Function   : UpdateBext
// Description: bext チャンク更新
// ---
//	bext	:bext構造体へのアドレス
//	inFmt	:入力ファイル音声形式情報
//	inFmt	:入力ファイル音声形式情報
//	outFmt	:出力ファイル音声形式情報
//	param	:パラメーター構造体
//	bwf_size:bextのバイト数
//
int UpdateBext(BROADCAST_EXT *bext,SOUNDFMT *inFmt,SOUNDFMT *outFmt,PARAM_INFO *param,long bwf_size)
{
	__int64 timeRef;
	double d_timeRef;
	char strbuf[512];
	strbuf[0] = '\0';

	// version
	bext->version = 1;

	// TimeRef
	timeRef = bext->timeReferenceHigh;timeRef <<= 32;
	timeRef |= bext->timeReferenceLow;
	if (timeRef != 0) {
		d_timeRef = (double)timeRef;
		d_timeRef /= inFmt->sample;
		d_timeRef *= outFmt->sample;
		timeRef = (__int64)d_timeRef;
		bext->timeReferenceHigh = (DWORD)(timeRef >> 32);
		bext->timeReferenceLow	= (DWORD)timeRef;
	}
	// Coding History
	if (strlen(bext->codingHistory) == 0) {
		sprintf(strbuf,"A=PCM,F=%d,W=%d,M=%s,T=original,\r\n",inFmt->sample,inFmt->bitsPerSample,inFmt->channel == 1 ? "mono" : "stereo");
		strcat(bext->codingHistory,strbuf);
	}
	sprintf(strbuf,"A=PCM,F=%d,W=%d,M=%s,T=Upconv0.6.x;",outFmt->sample,outFmt->bitsPerSample,outFmt->channel == 1 ? "mono" : "stereo");
	strcat(bext->codingHistory,strbuf);
	if (param->abe) {
		strcpy(strbuf,"ABE;");
		strcat(bext->codingHistory,strbuf);
	}
	sprintf(strbuf,"HFA%d",param->hfa);
	strcat(bext->codingHistory,strbuf);
	strcat(bext->codingHistory,",\r\n");
	return 0;
}

//---------------------------------------------------------------------------
// Function   : WriteData
// Description: データ書き込み
//
int WriteData(FIO *fp_w,int bit,SSIZE data)
{
	unsigned char buffer[3];
	float f;
	double d;
	int ret;
	
	ret = FALSE;
	switch (bit) {
		case 16:
			buffer[0] = (unsigned char)data;
			buffer[1] = (unsigned char)(data >> 8);
			fio_write(buffer,2,1,fp_w);
			if (!fp_w->error) {
				ret = TRUE;
			}
			break;
		case 24:
			buffer[0] = (unsigned char)data;
			buffer[1] = (unsigned char)(data >> 8);
			buffer[2] = (unsigned char)(data >> 16);
			fio_write(buffer,3,1,fp_w);
			if (!fp_w->error) {
				ret = TRUE;
			}
			break;
		case 32:
			f = (float)data / 0x7FFFFFFFFFFFFFFF;
			fio_write(&f,sizeof (float),1,fp_w);
			if (!fp_w->error) {
				ret = TRUE;
			}
			break;
		case 64:
			d = (double)( (__float128)data / (__float128)(0x7FFFFFFFFFFFFFFF) );
			fio_write(&d,sizeof (double),1,fp_w);
			if (!fp_w->error) {
				ret = TRUE;
			}
			break;
		default:
			break;
	}
	
	return ret;
}
