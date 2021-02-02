//---------------------------------------------------------------------------
/****************************************************************************/
/* upconv (C) 2007-2012 By 59414d41											*/
/*																			*/
/*																			*/
/****************************************************************************/

/*--- Log ------------------------------------------------------------------
 * Ver 0.10 <07/03/15> - 作成
 * Ver 0.20 <07/04/19> - DFT バージョン追加
 * Ver 0.30 <08/08/17> - 新方式に変更
 * Ver 0.40 <09/02/20> - fftw バージョンに変更しいったんFix
 * Ver 0.50 <09/05/18> - 高域の補間方法をいったんFix
 * Ver 0.60 <09/06/02> - BCB6とVC6に対応
 * Ver 0.61 <09/06/06> - メモリマップドファイルを使用しないようにした
 * Ver 0.70 <09/06/28> - 処理方法を変更
 * Ver 0.80 <09/07/10> - hfa2 の処理を改良しマルチコアに対応するために処理を分離
 * Ver 0.90 <09/09/24> - メモリ使用量削減
 * Ver 0.91 <09/09/27> - ダウンサンプリング時のバグ修正、ノーマライズ時のバグ修正
 * Ver 0.92 <09/09/29> - hfa2の補間方法を変更
 * Ver 0.93 <09/10/06> - hfa2の補間処理を変更
 * Ver 0.94 <09/11/01> - hfa3追加、バグ修正、パラメータファイルの採用
 * Ver 0.95 <09/11/11> - hfa3にノイズとのブレンド具合を指定可能にした
 * Ver 0.96 <09/11/15> - hfa補間時の周波数特性を修正
 * Ver 0.97 <09/11/22> - ビット拡張処理追加、hfa補間時の周波数特性を指定可能にした
 * Ver 0.98 <10/01/11> - デエンファシスに対応
 * Ver 0.99 <10/03/14> - hfa3の補間方法を変更
 * Ver 1.00 <10/03/14> - GCC(MinGW)に対応
 * Ver 1.01 <10/04/11> - OpenMPに対応
 * Ver 1.02 <10/07/26> - スピーカーごとの調整機能追加
 * Ver 1.03 <10/09/14> - hfc autoに対応
 * Ver 1.04 <10/11/02> - スピーカーごとの調整機能バグ修正
 * Ver 1.05 <10/12/27> - イコライザ機能修正
 * Ver 1.06 <11/01/07> - lfa 対応、ソースコードの整理
 * Ver 1.07 <11/10/01> - テンポラリファイル対応
 * Ver 1.08 <12/02/28> - fio 対応
 */

#define STR_COPYRIGHT	"upconv.exe (c) 2012 Ver 1.08 By 59414d41\n\n"
#define STR_USAGE		"upconv.exe paramfile outfile\n"

#define STATUS_SUCCESS			(0)		/* 正常終了 */
#define STATUS_MEM_ALLOC_ERR	(-4)	/* メモリー確保エラー */
#define STATUS_FILE_READ_ERR	(-5)	/* ファイルリードエラー */
#define STATUS_FILE_WRITE_ERR	(-6)	/* ファイルライトエラー */

#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("/tmp/upconv.log","a");									\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n",__FUNCTION__,__LINE__,s);			\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* #include <windows.h> */
#include <time.h>
#include "fileio.h"
#include "fftw3.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#define M_PI		(3.14159265358979323846)
#endif

#ifndef DWORD
#define DWORD		unsigned int
#endif

#ifndef TRUE
#define TRUE		(1)
#endif

#ifndef FALSE
#define FALSE		(0)
#endif

// 最大サンプリング数
#define MAX_SAMPLE_N	(1536000)

// サンプルを処理するデータ型
#define SSIZE	signed long long

// ノーマライズ情報
typedef struct {
	SSIZE	min;
	SSIZE	max;
	SSIZE	avg;
} NORM_INFO;

typedef struct {
	long cutOff[257];
} BITEXT_INFO;

typedef struct {
	int err;
	int errLine;
	int disRc;
	int disEq;
	long inSampleR;
	long outSampleR;
	long validOutSampleR;
	int hfa;
	int lfa;
	int overSamp;
	int adjBE;
	int abeNX;
	int abeFnLevel;
	int cutLowData;
	int smLowData;
	int adaptiveFilder;
	int w;
	int iw;
	int hiDither;
	int hfc_auto;
	long hfc;
	long lfc;
	long lpf;
	long nr;
	long fio;
	int nrLV;

	// hfa2,3 option
	int sig1Enb,sig2Enb,sig3Enb;
	int sig1AvgLineNx;
	int sig1Phase,sig2Phase,sig3Phase;
	int hfaNB;
	int hfaFast;
	int hfaWide;
	int hfaDiffMin;
	// 
	int deEmphasis;
	BITEXT_INFO beInfo;
	char tempPath[4];
	int thread;
	int downSample;
	int ana;
	int sp_ana;
	double *ana_pw;
	double *ana_avgpw;
	long   ana_avgcount;
	double *lfa_eq;
	double *sp_eq;
	double *eq;
	char *sp_path;
	char *nd_path;
	int r1_flag;
	int eq_flag;
	int post_abe;
	double hfa3_max;
	int pwAdj;
} PARAM_INFO;

typedef struct {
	double	power[65536*2];
	double	phase[65536*2];
	double	pw[65536*2];
	double	diff[65536*2];
	double	avg[65536*2];
	double	base[65536*2];
	double	baseToAdj[65536*2];
	int		pw_cnt[65536*2];
	int		nSample;
	int		validSamplingRate;
	long	samplingRate;
	int		do2Idx[360];
	double	*phaseX;
	double	*phaseY;
	int log;
} OVERTONE_INFO;

SSIZE *diskBuffer;
SSIZE *diskBuffer2;
NORM_INFO NormInfo;

/*--- Function Prototype ---------------------------------------------------*/
void SamplingRateConvert(char *rawFile,PARAM_INFO *param);
void fftFilter(DWORD inSample,DWORD outSample,FIO *fp,FIO *tmp,PARAM_INFO *param);
void fftFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,PARAM_INFO *param,int id);
void anaLfaParam(PARAM_INFO *param);
void anaHFC_AutoParam(PARAM_INFO *param);
void spAnalyze(DWORD inSample,FIO *fp,PARAM_INFO *param);
void adjBitExtension(DWORD inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void genNoise(long hfc,DWORD inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void genOverTone(long hfc,DWORD inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void genOverToneSub(long hfc,SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,OVERTONE_INFO *ovInfo,PARAM_INFO *param);
int bpFilter(long lfc,long hfc,DWORD inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void bpFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,long lfc,long hfc,PARAM_INFO *param);
void anaLfaParam(PARAM_INFO *param);

void noiseCut(long nfs,DWORD inSample,FIO *fio_r,FIO *fio_w,PARAM_INFO *param);
void anaOverToneHFA2(OVERTONE_INFO *ovInfo,PARAM_INFO *param);
void anaOverToneHFA3(OVERTONE_INFO *ovInfo,PARAM_INFO *param);
void outTempFile(FIO *fio,void *in,DWORD size,PARAM_INFO *param);
double normalNoise(void);
void adjPinkFilter(int mode,long fftSizeOut,fftw_complex *fftw_out2,PARAM_INFO *param);
void merageTempFile(char type,int normFlag,FIO *in1,FIO *in2,FIO *out,DWORD inSample,PARAM_INFO *param);
void copyToFFTW(fftw_complex *fftw,SSIZE *buf,long size);
void windowFFTW(fftw_complex *fftw,long size);
void cutFFTW(fftw_complex *fftw,long index,long size);
void *al_malloc(long size);
void *al_free(void *ptr);

//---------------------------------------------------------------------------
// Function   : main
// Description: 引数を処理し変換関数を呼び出す
//
//
int main(int argc, char *argv[])
{
	FILE *fp;
	FILE *fp_files;
	char tmppath[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	char pparam[512];
	char ddBuf[100];
	long paramFlag;
	char *p1,*p2;
	double dd;
	long i;
	static PARAM_INFO param;
	long temp,temp2,temp3;
	int retCode = 0;
	SSIZE max,avg;
	double persent;
	int c;

	do {
		memset(&param,0,sizeof (PARAM_INFO));

		diskBuffer	= (SSIZE *)calloc(4 * 1024 * 1024L,sizeof (SSIZE));
		diskBuffer2 = (SSIZE *)calloc(4 * 1024 * 1024L,sizeof (SSIZE));
		if (diskBuffer == NULL || diskBuffer2 == NULL) {
			param.err = STATUS_MEM_ALLOC_ERR;param.errLine = __LINE__;
			break;
		}
		memset(&NormInfo,0,sizeof (NORM_INFO));

		param.ana_pw 	= malloc(MAX_SAMPLE_N * sizeof (double));
		param.ana_avgpw	= malloc(MAX_SAMPLE_N * sizeof (double));
		param.lfa_eq	= malloc(MAX_SAMPLE_N * sizeof (double));
		param.sp_eq		= malloc(MAX_SAMPLE_N * sizeof (double));
		param.eq		= malloc(MAX_SAMPLE_N * sizeof (double));
		if (param.ana_pw == NULL || param.ana_avgpw == NULL || param.lfa_eq == NULL || param.sp_eq == NULL || param.eq == NULL) {
			param.err = STATUS_MEM_ALLOC_ERR;param.errLine = __LINE__;
			break;
		}

		for (i = 0;i < MAX_SAMPLE_N;i++) {
			param.ana_pw[i] = 0;
			param.ana_avgpw[i] = 0;
		}
		param.ana_avgcount = 0;
		param.err = STATUS_SUCCESS;
		param.disRc = 0;
		param.disEq = 0;
		param.hfc = -1;
		param.lfc = -1;
		param.lpf = -1;
		param.nr = -1;
		param.nrLV = 1;
		param.hfaNB = 0;
		param.thread = 1;
		param.sp_ana = 0;
		param.hiDither = 0;
		param.hfc_auto = 0;
		param.r1_flag = 0;
		param.eq_flag = 0;
		param.post_abe = 0;
		param.fio = -1;
		param.abeNX = 1;
		param.adaptiveFilder = 1;
		param.pwAdj = 0;
		paramFlag = 0;
		if (argc == 3) {
			// r1 ファイルかをチェックする(sp_ana)
			_splitpath(argv[2],drive,dir,fname,ext);
			if (strcmp(ext,".r1") == 0) {
				param.r1_flag = 1;
			}
			// eq
			for (i = 0;i < 192000;i++) {
				param.eq[i] = 1.00;
			}
			_splitpath(argv[0],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,"eq","dat");
			fp = fopen(tmppath,"r");
			if (fp) {
				for (i = 0;i < 192000;i++) {
					fgets(ddBuf,100,fp);
					if (sscanf(ddBuf,"%lf",&dd) == 1) {
						param.eq[i] = dd;
					}
				}
				fclose(fp);
			}

			// sp_ana パス
			_splitpath(argv[0],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,"adj","dat");
			param.sp_path = malloc(strlen(tmppath) + 1);
			if (param.sp_path != NULL) {
				strcpy(param.sp_path,tmppath);
			}
			fp = fopen(param.sp_path,"r");
			if (fp) {
				for (i = 0;i < 192000;i++) {
					fgets(ddBuf,100,fp);
					if (sscanf(ddBuf,"%lf",&dd) == 1) {
						param.sp_eq[i] = dd;
					}
				}
				fclose(fp);
			}

			// GenNoise 用のデータパス
			_splitpath(argv[0],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,"nd","dat");
			param.nd_path = malloc(strlen(tmppath) + 1);
			if (param.nd_path != NULL) {
				strcpy(param.nd_path,tmppath);
			}

			// ビット拡張テーブル
			_splitpath(argv[0],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,"bit_extend_table","");
			
			fp = fopen(tmppath,"r");
			if (fp) {
				for (i = 0;i < 256;i++) {
					if (fscanf(fp,"%d,",&c) == 1) {
						param.beInfo.cutOff[i] = c;
					}
				}
				fclose(fp);
			}
			
			_splitpath(argv[1],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"files");
			for (i = 0;i < 10;i++) {
				// ファイルオープン
				fp_files = fopen(tmppath,"a");
				if (fp_files != NULL) {
					strcpy(tmppath,argv[2]);
					strcat(tmppath,".tmp");
					fprintf(fp_files,"%s\n",tmppath);
					strcpy(tmppath,argv[2]);
					strcat(tmppath,".tmp2");
					fprintf(fp_files,"%s\n",tmppath);
					strcpy(tmppath,argv[2]);
					strcat(tmppath,".param");
					fprintf(fp_files,"%s\n",tmppath);
					fclose(fp_files);
					break;
				}
				usleep(3);
			}
			if (i == 10) {
				param.err = STATUS_FILE_READ_ERR;param.errLine = __LINE__;
				break;
			}

			// param ファイル
			_splitpath(argv[1],drive,dir,fname,ext);
			_makepath(tmppath,drive,dir,fname,"param");
			printf("paramfile:%s\n", tmppath);
			// ファイルオープン
			fp = fopen(tmppath,"r");
			if (fp == NULL) {
				param.err = STATUS_FILE_READ_ERR;param.errLine = __LINE__;
				break;
			}

			if (fgets(pparam,512,fp) == NULL) {
				param.err = STATUS_FILE_READ_ERR;param.errLine = __LINE__;
				break;
			}

			fclose(fp);
			p1 = pparam;
			p2 = strchr(p1,(int)' ');

			for (;p1 != NULL;) {
				if (p2 != NULL) {
					*p2 = '\0';
				}

				if (sscanf(p1,"-is:%ld",&temp) == 1) {
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
							param.inSampleR = temp;
							paramFlag |= 0x0001;
							break;
					}
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
							paramFlag |= 0x0002;
							param.outSampleR = temp;
							param.validOutSampleR = temp;
							break;
					}
				}
				if (sscanf(p1,"-ms:%ld",&temp) == 1) {
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
							paramFlag |= 0x0002;
							param.outSampleR = temp;
							param.validOutSampleR = temp;
							break;
					}
				}
				if (sscanf(p1,"-hfa:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
						case 2:
						case 3:
							paramFlag |= 0x0004;
							param.hfa = (int)temp;
							break;
					}
				}
				if (sscanf(p1,"-lfa:%ld",&temp) == 1) {
					switch (temp) {
						case 0:
						case 1:
							param.lfa = (int)temp;
							break;
					}
				}
				if (sscanf(p1,"-hfc:%ld",&temp) == 1) {
					if (temp >= 1000 && temp <= (384000 / 2)) {
						param.hfc = temp;
					}
				}
				if (sscanf(p1,"-lfc:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= (384000 / 2)) {
						param.lfc = temp;
					}
				}
				if (sscanf(p1,"-lpf:%ld",&temp) == 1) {
					if (temp >= 1000 && temp <= (384000 / 2)) {
						param.lpf = temp;
					}
				}
				if (sscanf(p1,"-nr:%ld",&temp) == 1) {
					if (temp >= 100 && temp <= (384000 / 2)) {
						param.nr = temp;
					}
				}
				if (sscanf(p1,"-sig1:%ld,%ld,%ld",&temp,&temp2,&temp3) == 3) {
					if ((temp == 0 || temp == 1) &&
						(temp2 >= 1 && temp2 <= 25) &&
						(temp3 >= -44 && temp3 <= 44)) {

						param.sig1Enb = (int)temp;
						param.sig1AvgLineNx = (int)temp2;
						param.sig1Phase = (int)temp3;

					}
				}
				if (sscanf(p1,"-sig2:%ld,%ld",&temp,&temp2) == 2) {
					if ((temp == 0 || temp == 1) &&
						(temp2 >= -44 && temp2 <= 44)) {

						param.sig2Enb = (int)temp;
						param.sig2Phase = (int)temp2;
					}
				}
				if (sscanf(p1,"-sig3:%ld,%ld",&temp,&temp2) == 2) {
					if ((temp == 0 || temp == 1) &&
						(temp2 >= -44 && temp2 <= 44)) {

						param.sig3Enb = (int)temp;
						param.sig3Phase = (int)temp2;
					}
				}
				if (sscanf(p1,"-hfaNB:%ld",&temp) == 1) {
					if (temp >= 0 && temp <= 100) {
						param.hfaNB = temp;
					}
				}
				if (sscanf(p1,"-hfaNB:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 5) {
						param.hfaDiffMin = temp;
					}
				}
				if (sscanf(p1,"-overSamp:%ld",&temp) == 1) {
					if (temp == 0 || temp == 2 || temp == 6) {
						param.overSamp = (int)temp;
					}
				}
				if (strcmp(p1,"-hdc") == 0) {
					param.hiDither = 1;
				}
				if (strcmp(p1,"-hfcAuto") == 0) {
					param.hfc_auto = 1;
				}
				if (sscanf(p1,"-thread:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 24) {
						param.thread = (int)temp;
					}
				}

				if (strcmp(p1,"-adjBE") == 0) {
					param.adjBE = 1;
				}
				if (strcmp(p1,"-cutLowData") == 0) {
					param.cutLowData = 1;
				}
				if (strcmp(p1,"-smLowData") == 0) {
					param.smLowData = 1;
				}
				if (sscanf(p1,"-cutDither:%ld",&temp) == 1) {
					if (temp == 1) {
						param.abeFnLevel = 2;
					} else if (temp == 2) {
						param.abeFnLevel = 5;
					} else if (temp == 3) {
						param.abeFnLevel = 10;
					} else if (temp == 4) {
						param.abeFnLevel = 15;
					} else if (temp == 5) {
						param.abeFnLevel = 20;
					}
				}
				if (strcmp(p1,"-postABE") == 0) {
					param.post_abe = 1;
				}
				if (strcmp(p1,"-hfaFast") == 0) {
					param.hfaFast = 1;
				}
				if (strcmp(p1,"-hfaWide") == 0) {
					param.hfaWide = 1;
				}
				if (sscanf(p1,"-iw:%ld",&temp) == 1) {
					if (temp == 16 || temp == 24 || temp == 32 || temp == 64) {
						if(temp == 64){
							temp = 56;
						}
						param.iw = (int)temp;
					}
				}
				if (sscanf(p1,"-w:%ld",&temp) == 1) {
					if (temp == 16 || temp == 24 || temp == 32 || temp == 64) {
						param.w = (int)temp;
					}
				}
				if (strcmp(p1,"-adjFreq") == 0) {
					param.ana = 1;
				}
				if (strcmp(p1,"-spAna") == 0) {
					// 192khz 固定で解析
					param.sp_ana = 1;
				}
				if (sscanf(p1,"-fio:%ld",&temp) == 1) {
					if (temp >= 30 && temp <= 16000) {
						param.fio = temp / 10;
					}
				}
				if (strcmp(p1,"-eq") == 0) {
					param.eq_flag = 1;
				}
				if (sscanf(p1,"-nrLV:%ld",&temp) == 1) {
					if (temp >= 1 && temp <= 5) {
						param.nrLV = (int)temp - 1;
					}
				}
				if (sscanf(p1,"-deEmphasis:%ld",&temp) == 1) {
					if (temp == 0 || temp == 1 || temp == 2) {
						param.deEmphasis = (int)temp;
					}
				}
				if (sscanf(p1,"-temp:%c%c%c",tmppath,tmppath+1,tmppath+2) == 3) {
					param.tempPath[0] = tmppath[0];
					param.tempPath[1] = tmppath[1];
					param.tempPath[2] = '\0';
					param.tempPath[3] = '\0';
				}

				if (p2 == NULL) {
					break;
				}
				p1 = p2 + 1;
				p2 = strchr(p1,(int)' ');
			}
			
			if (param.ana == 1) {
				param.hiDither	 = 0;
			}
			if (param.sp_ana == 1) {
				param.outSampleR = 192000;
				param.hfa		 = 0;
				param.overSamp	 = 0;
				param.hfc		 = -1;
				param.lfc		 = -1;
				param.lpf		 = -1;
				param.nr		 = -1;
				param.hiDither	 = 0;
				param.hfc_auto	 = 0;
				param.eq_flag	 = 0;
				param.post_abe	 = 0;
			}
#ifdef _OPENMP
	int nCpu;
	nCpu = param.thread;
	omp_set_num_threads(nCpu);
#endif
			if (param.hfc != -1) {
				param.hfc_auto = 0;
			}

			if (paramFlag == 0x0007) {
				SamplingRateConvert(argv[2],&param);
				if (param.err == STATUS_SUCCESS) {
					if (NormInfo.max < 0) {
						NormInfo.max *= -1;
					}
					if (NormInfo.min < 0) {
						NormInfo.min *= -1;
					}
					max = NormInfo.max;
					if (max < NormInfo.min) {
						max = NormInfo.min;
					}
					persent = (double)((__float128)max / (__float128)0x7FFFFFFFFFFFFF );
					avg = NormInfo.avg << 40;
					strcpy(tmppath,argv[2]);
					strcat(tmppath,".param");
					fp = fopen(tmppath,"wb");
					if (fp == NULL) {
						param.err = STATUS_FILE_WRITE_ERR;
						break;
					}
					if (fprintf(fp,"%.10lf,%llx\n",persent,avg) == EOF) {
						param.err = STATUS_FILE_WRITE_ERR;
						break;
					}
					fclose(fp);
				}
			}
		}
		if (argc != 3 || paramFlag != 0x0007) {
			printf(STR_COPYRIGHT);
			printf(STR_USAGE);
			exit(0);
		}
	} while (0);

	fprintf(stdout,"\n");
	fflush(stdout);
	fprintf(stdout,"End\n");
	fflush(stdout);

	if (param.err != STATUS_SUCCESS) {
		_splitpath(argv[1],drive,dir,fname,ext);
		_makepath(tmppath,drive,dir,fname,"err");
		fp = fopen(tmppath,"a");
		if (fp) {
			switch (param.err) {
				case STATUS_FILE_READ_ERR:
					fprintf(fp,"upconv[%d]:File read error.\n",param.errLine);
					break;
				case STATUS_FILE_WRITE_ERR:
					fprintf(fp,"upconv[%d]:File write error.\n",param.errLine);
					break;
				case STATUS_MEM_ALLOC_ERR:
					fprintf(fp,"upconv[%d]:Memory Allocation error.\n",param.errLine);
					break;
				default:
					fprintf(fp,"upconv[%d]:Other error.\n",param.errLine);
			}
			fclose(fp);
		}
	}
	exit(retCode);

	return 0;
}
//---------------------------------------------------------------------------
// Function   : SamplingRateConvert
// Description: サンプリングレート変換処理をする
// ---
//	rawFile	: RAWデータファイル名
//	param	: 変換パラメータ構造体
//
void SamplingRateConvert(char *rawFile,PARAM_INFO *param)
/*
 *	サンプリングレート変換
 */
{
	char outFile[_MAX_PATH];
	char tmpFile[_MAX_PATH];
	char drive[_MAX_DRIVE];
	char dir[_MAX_DIR];
	char fname[_MAX_FNAME];
	char ext[_MAX_EXT];
	FIO fp_r;
	FIO fp_w;
	FIO fp_r2;
	fio_size size;
	DWORD inSample,outSample;
	int hfc;
	double dd;
	DWORD svInSampleR;
	int svIw;
	int downSampleFlag = FALSE;
	DWORD svOutSampleR;

	memset(&fp_r,0,sizeof (FIO));
	memset(&fp_r2,0,sizeof (FIO));
	memset(&fp_w,0,sizeof (FIO));

	do {
		// テンポラリファイル名の作成
		if (strlen(param->tempPath) > 0) {
			_splitpath(rawFile,drive,dir,fname,ext);
			//strcpy(dir,"\\");
			_makepath(outFile,param->tempPath,dir,fname,ext);
		} else {
			strcpy(outFile,rawFile);
		}
		strcat(outFile,".tmp");
		strcpy(tmpFile,outFile);
		strcat(tmpFile,"2");
		
		// リード専用ファイルオープン関数(バッファリング機能付)
		fio_open(&fp_r,rawFile,FIO_MODE_R);
		if (fp_r.error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		if (param->overSamp == 2) {
			// 2倍のサンプリングレートで計算し後でダウンサンプルする。
			svOutSampleR = param->outSampleR;
			param->outSampleR *= 2;
		} else if (param->overSamp == 6) {
			// 768000のサンプリングレートで計算し後でダウンサンプルする。
			svOutSampleR = param->outSampleR;
			param->outSampleR = 768000;
		}

		outSample = 0;
		
		// ファイルサイズ取得
		fio_get_filesize(&fp_r,&size);
		if (fp_r.error) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}
		
		// 出力サンプル数を計算する。
		inSample = (DWORD)(size / sizeof (SSIZE));
		dd = inSample;
		dd /= param->inSampleR;
		dd *= param->outSampleR;
		outSample = (DWORD)dd;

		if (outSample == 0) {
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		// 量子化ビット拡張処理
		if (param->adjBE != 0) {
			fprintf(stdout,"ABE\n");
			fflush(stdout);

			// 出力用にファイルオープン
			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			fio_set_memory_limit(&fp_w,param->fio);
			adjBitExtension(inSample,&fp_r,&fp_w,param);
			if (param->err) {
				PRINT_LOG("ERROR");
				break;
			}
			fio_close(&fp_r);
			if (fp_r.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
		}

		param->cutLowData = 0;

		fprintf(stdout,"SRC\n");
		fflush(stdout);
		if (inSample > outSample) {
			downSampleFlag = TRUE;
		}

		fio_open(&fp_w,outFile,FIO_MODE_W);
		if (fp_w.error) {
				PRINT_LOG("ERROR");
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}

		// ファイルに出力するサイズを制限する(outSample数)
		fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
		fio_set_memory_limit(&fp_w,param->fio);

		// 音量調査、パワー調査用の fft
		// eq,lfc,hfc 等はしない
		param->disEq = 1;
		fftFilter(inSample,outSample,&fp_r,&fp_w,param);
		if (param->err) {
				PRINT_LOG("ERROR");
			break;
		}

		// 音量調査
		fio_flush(&fp_w);
		merageTempFile(' ',1,&fp_w,NULL,NULL,outSample,param);
		if (param->err) {
				PRINT_LOG("ERROR");
			break;
		}

		fio_close(&fp_w);

		// hfc auto 用パラメーター設定(hfc)
		if (param->hfc_auto == 1) {
			anaHFC_AutoParam(param);
		}

		// lfa 用パラメーター作成(lfa_eq)
		if (param->lfa != 0) {
			anaLfaParam(param);
		}

		fio_open(&fp_w,outFile,FIO_MODE_W);
		if (fp_w.error) {
			PRINT_LOG("ERROR");
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		// ファイルに出力するサイズを制限する(outSample数)
		fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
		fio_set_memory_limit(&fp_w,param->fio);

		// fft
		fftFilter(inSample,outSample,&fp_r,&fp_w,param);
		if (param->err) {
			PRINT_LOG("ERROR");
			break;
		}
		fio_close(&fp_r);
		if (fp_r.error) {
			PRINT_LOG("ERROR");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		// 生成したデータを読み込み用に設定する。
		fio_setmode_r(&fp_w,&fp_r,rawFile);
		if (fp_w.error) {
			PRINT_LOG("ERROR");
			param->err = STATUS_FILE_WRITE_ERR;
			break;
		}
		if (fp_r.error) {
			PRINT_LOG("ERROR");
			param->err = STATUS_FILE_READ_ERR;
			break;
		}

		param->hiDither = 0;
		param->lfa = 0;
		param->lpf = 0;
		param->ana = 0;
		param->eq_flag = 0;

		// スピーカー用周波数解析
		if (param->sp_ana == 1) {
			spAnalyze(outSample,&fp_r,param);
		}

		//
		// ノイズカット
		if (param->nr != -1) {
			fprintf(stdout,"NR\n");
			fflush(stdout);
			if (param->hfc != -1) {
				hfc = param->hfc;
			} else {
				if (downSampleFlag == FALSE) {
					hfc = param->inSampleR / 2;
				} else {
					hfc = param->outSampleR / 2;
				}
			}
			if (downSampleFlag == FALSE) {
				if (hfc > param->inSampleR / 2) {
					hfc = param->inSampleR / 2;
				}
			} else {
				if (hfc > param->outSampleR / 2) {
					hfc = param->outSampleR / 2;
				}
			}
			if (param->nr < hfc) {
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

				noiseCut(param->nr,outSample,&fp_r,&fp_w,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

				bpFilter(-1,hfc,outSample,&fp_r2,&fp_w,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_close(&fp_r2);
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

				merageTempFile('-',0,&fp_r,&fp_r2,&fp_w,outSample,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_close(&fp_r);
				fio_close(&fp_r2);
				if (fp_r.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				if (fp_r2.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_setmode_r(&fp_w,&fp_r,rawFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
			}
		}

		//
		// 高域補間処理
		if (param->hfa != 0) {
			fflush(stdout);
			if (param->hfc != -1) {
				hfc = param->hfc;
			} else {
				if (downSampleFlag == FALSE) {
					hfc = param->inSampleR / 2;
				} else {
					hfc = param->outSampleR / 2;
				}
			}
			if (downSampleFlag == FALSE) {
				if (hfc > param->inSampleR / 2) {
					hfc = param->inSampleR / 2;
				}
			} else {
				if (hfc > param->outSampleR / 2) {
					hfc = param->outSampleR / 2;
				}
			}
			if (param->hfa != 1) {
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

				fprintf(stdout,"HFA%d\n",param->hfa);
				fflush(stdout);
				genOverTone(hfc,outSample,&fp_r,&fp_w,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

//				bpFilter(hfc,-1,outSample,&fp_w,param);
				bpFilter(hfc,-1,outSample,&fp_r2,&fp_w,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_close(&fp_r2);
				fio_setmode_r(&fp_w,&fp_r2,tmpFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r2.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

				merageTempFile('+',0,&fp_r,&fp_r2,&fp_w,outSample,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_close(&fp_r);
				fio_close(&fp_r2);
				if (fp_r.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				if (fp_r2.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_setmode_r(&fp_w,&fp_r,rawFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
			}
			if (param->hfa == 1 || param->hfaNB > 0) {
				fprintf(stdout,"HFA1\n");
				fflush(stdout);
				fio_open(&fp_w,outFile,FIO_MODE_W);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				// ファイルに出力するサイズを制限する(outSample数)
				fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
				fio_set_memory_limit(&fp_w,param->fio);

				genNoise(hfc,outSample,&fp_r,&fp_w,param);
				if (param->err) {
					PRINT_LOG("ERROR");
					break;
				}
				fio_close(&fp_r);
				if (fp_r.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
				fio_setmode_r(&fp_w,&fp_r,rawFile);
				if (fp_w.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
				if (fp_r.error) {
					PRINT_LOG("ERROR");
					param->err = STATUS_FILE_READ_ERR;
					break;
				}
			}
		}
		if (param->post_abe == 1) {
			fprintf(stdout,"Post ABE\n");
			fflush(stdout);
			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			// ファイルに出力するサイズを制限する(outSample数)
			fio_set_maxsize(&fp_w,(fio_size)outSample * sizeof (SSIZE));
			fio_set_memory_limit(&fp_w,param->fio);

			svInSampleR = param->inSampleR;
			svIw		= param->iw;
			param->inSampleR = param->outSampleR;
			param->iw = param->w;
			if(param->iw == 64){
				param->iw = 56;
			}
			param->cutLowData = 0;
			param->abeFnLevel = 0;
			param->smLowData  = 0;
			param->adaptiveFilder = 0;
			param->abeNX = 0;
			adjBitExtension(outSample,&fp_r,&fp_w,param);
			if (param->err) {
				PRINT_LOG("ERROR");
				break;
			}
			param->inSampleR = svInSampleR;
			param->iw = svIw;

			fio_close(&fp_r);
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_READ_ERR;
				break;
			}

		}

		// オーバーサンプリング用のダウンサンプル処理
		if (param->overSamp != 0) {
			DWORD wkSample;
			memset(&NormInfo,0,sizeof (NORM_INFO));

			fio_open(&fp_w,outFile,FIO_MODE_W);
			if (fp_w.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			fio_set_memory_limit(&fp_w,param->fio);

			fprintf(stdout,"SRC(DS)\n");
			fflush(stdout);

			param->inSampleR  = param->outSampleR;
			param->outSampleR = svOutSampleR;
			
			dd = outSample;
			dd /= param->inSampleR;
			dd *= param->outSampleR;
			wkSample = (DWORD)dd;

			// ファイルに出力するサイズを制限する(outSample数)
			fio_set_maxsize(&fp_w,(fio_size)wkSample * sizeof (SSIZE));

			param->disEq = 1;
			//param->pwAdj = 1;
			fftFilter(outSample,wkSample,&fp_r,&fp_w,param);
			if (param->err) {
				PRINT_LOG("ERROR");
				break;
			}
			fio_close(&fp_r);
			if (fp_r.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			fio_setmode_r(&fp_w,&fp_r,rawFile);
			if (fp_w.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_WRITE_ERR;
				break;
			}
			if (fp_r.error) {
				PRINT_LOG("ERROR");
				param->err = STATUS_FILE_READ_ERR;
				break;
			}
			merageTempFile(' ',1,&fp_r,NULL,NULL,wkSample,param);
			if (param->err) {
				PRINT_LOG("ERROR");
				break;
			}
		}
		fio_close(&fp_r);
		//
		// 終わり
	} while (0);
}

//---------------------------------------------------------------------------
// Function   : fftFilter
// Description: FFT によるフィルタ処理
// ---
//	inBuffer	:入力データバッファ
//	inSample	:入力データのサンプル数(ch毎)
//	fp_r		:入力ファイル用構造体
//	fp_w		:テンポラリファイル用構造体
//	param		:変換パラメータ
//
void fftFilter(DWORD inSample,DWORD outSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	SSIZE level,level2;
	double nx;
	long wkMemSize;
	long inSampleR,outSampleR;
	long fftSizeIn,fftSizeOut,i,j,n;
	long wr;
	double persent,per;
	SSIZE *pIn[3],*pOut[3];
	SSIZE startInSample,nSample;
	DWORD outRemain;
	fftw_complex *fftw_in[3],*fftw_out[3];
	fftw_plan fftw_p[3],fftw_ip[3];
	SSIZE membyte;
	char s[50];

	fftw_in[0] = NULL;
	fftw_in[1] = NULL;
	fftw_in[2] = NULL;
	fftw_out[0] = NULL;
	fftw_out[1] = NULL;
	fftw_out[2] = NULL;

	fio_rewind(fp_r);

	inSampleR = param->inSampleR;
	outSampleR = param->outSampleR;

	fftSizeIn = inSampleR * 2;
	fftSizeOut = outSampleR * 2;
	if (param->disRc) {
		inSampleR = outSampleR;
		fftSizeIn = fftSizeOut;
	}

	wkMemSize = fftSizeIn;
	if (wkMemSize < fftSizeOut) {
		wkMemSize = fftSizeOut;
	}
	wkMemSize *= 2;

membyte = 0;
	// 入力用
	mem1 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem1 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
membyte += wkMemSize * sizeof (SSIZE);
	// 出力用
	mem2 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem2 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
membyte += wkMemSize * sizeof (SSIZE);
	// ワーク用(別スレッド用)
	mem3 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem3 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
membyte += wkMemSize * sizeof (SSIZE);
	mem4 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem4 == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
membyte += wkMemSize * sizeof (SSIZE);

	// 1
	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_in[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
membyte += wkMemSize * sizeof (fftw_complex);

	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_out[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
membyte += wkMemSize * sizeof (fftw_complex);
	fftw_p[0] = fftw_plan_dft_1d(fftSizeIn,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSizeOut,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		PRINT_LOG("ERROR");
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	if (param->overSamp != 32) {
		// 2
		fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_in[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
membyte += wkMemSize * sizeof (fftw_complex);
		fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_out[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
membyte += wkMemSize * sizeof (fftw_complex);
		fftw_p[1] = fftw_plan_dft_1d(fftSizeIn,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
		if (fftw_p[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_ip[1] = fftw_plan_dft_1d(fftSizeOut,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
		if (fftw_ip[1] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		// 3
		fftw_in[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_in[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
membyte += wkMemSize * sizeof (fftw_complex);
		fftw_out[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
		if (fftw_out[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
membyte += wkMemSize * sizeof (fftw_complex);
		fftw_p[2] = fftw_plan_dft_1d(fftSizeIn,fftw_in[2],fftw_out[2],FFTW_FORWARD,FFTW_ESTIMATE);
		if (fftw_p[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
		fftw_ip[2] = fftw_plan_dft_1d(fftSizeOut,fftw_out[2],fftw_in[2],FFTW_BACKWARD,FFTW_ESTIMATE);
		if (fftw_ip[2] == NULL) {
			PRINT_LOG("ERROR");
			param->err = STATUS_MEM_ALLOC_ERR;
			return;
		}
	}

	outRemain = outSample;

	pIn[0]	= &mem1[((fftSizeIn / 2) * 0)];
	pOut[0] = &mem2[((fftSizeOut / 2) * 0)];
	pIn[1]	= &mem1[((fftSizeIn / 2) * 1)];
	pOut[1] = &mem3[((fftSizeOut / 2) * 1)];
	pIn[2]	= &mem1[((fftSizeIn / 2) * 2)];
	pOut[2] = &mem4[((fftSizeOut / 2) * 2)];
	
	per = -1;
	for (startInSample = ((fftSizeIn + (fftSizeIn / 2)) * -1);startInSample < inSample + (fftSizeIn + fftSizeIn / 2);startInSample += fftSizeIn) {
		if (startInSample >= 0 && startInSample < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
		}

		memset(mem1,0,wkMemSize * sizeof (SSIZE));

		if (startInSample >= 0 && startInSample + (fftSizeIn * 2) < inSample) {
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			nSample = fftSizeIn * 2;
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSizeIn * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSizeIn * 2;
		}

		// 音のレベルを調査しておく
		if (param->pwAdj) {
			level = 0;
			for (i = fftSizeIn / 2,j = 0,n = 0;n < fftSizeIn;i++,j++,n++) {
				if (mem1[i] > 0) {
					level += mem1[i] >> (56 - 16);
					j++;
				}
			}
			if (j > 0) {
				level /= j;
			}
		}

		memset(mem2,0,wkMemSize * sizeof (SSIZE));
		memset(mem3,0,wkMemSize * sizeof (SSIZE));
		memset(mem4,0,wkMemSize * sizeof (SSIZE));

		if (param->overSamp != 32) {
			#pragma omp parallel
			{
				#pragma omp sections
				{
					#pragma omp section
					{
						// 1
						fftFilterSub(pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,0);
					}
					#pragma omp section
					{
						// 2
						fftFilterSub(pIn[1],pOut[1],fftw_in[1],fftw_out[1],fftw_p[1],fftw_ip[1],param,1);
					}
					#pragma omp section
					{
						// 3
						fftFilterSub(pIn[2],pOut[2],fftw_in[2],fftw_out[2],fftw_p[2],fftw_ip[2],param,2);
					}
				}
				#pragma omp for
				for (i = 0;i < wkMemSize;i++) {
					mem2[i] += mem3[i] + mem4[i];
				}
			}
		} else {
			// オーバーサンプリング(1536000)のときはメモリ消費が大きいのでfft用メモリは1つだけ確保し、並列処理はしない

			// 1
			fftFilterSub(pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,0);
			// 2
			fftFilterSub(pIn[1],pOut[1],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,1);
			// 3
			fftFilterSub(pIn[2],pOut[2],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],param,2);
			#pragma omp parallel for
			for (i = 0;i < wkMemSize;i++) {
				mem2[i] += mem3[i] + mem4[i];
			}
		}
		if (startInSample + fftSizeIn / 2 >= 0) {
			// 音のレベルを調整する
			if (param->pwAdj) {
				level2 = 0;
				for (i = fftSizeOut / 2,j = 0,n = 0;n < fftSizeOut;i++,j++,n++) {
					if (mem2[i] > 0) {
						level2 += mem2[i] >> (56 - 16);
						j++;
					}
				}
				if (j > 0) {
					level2 /= j;
				}
				if (level > 0 && level2 > 0) {
					nx = ((double)level / (double)level2);
					for (i = fftSizeOut / 2,n = 0;n < fftSizeOut;i++,n++) {
						mem2[i] = (SSIZE)((double)mem2[i] * nx);
					}
				}
			}

			if (outRemain >= fftSizeOut) {
				wr = fio_write(mem2 + (fftSizeOut / 2),sizeof (SSIZE),fftSizeOut,fp_w);
				if (wr != fftSizeOut) {
					char s[100];
					sprintf(s,"%ld:fio_write(%ld,%ld)",wr,sizeof (SSIZE),fftSizeOut);
					param->err = STATUS_FILE_WRITE_ERR;
					PRINT_LOG(s);
					return;
				}
			} else {
				wr = fio_write(mem2 + (fftSizeOut / 2),sizeof (SSIZE),outRemain,fp_w);
				if (wr != outRemain) {
					param->err = STATUS_FILE_WRITE_ERR;
					PRINT_LOG("ERROR");
					return;
				}
			}
			if (outRemain >= fftSizeOut) {
				outRemain -= fftSizeOut;
			} else {
				break;
			}
		}
	}
	param->disEq = 0;
	param->disRc = 0;

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);
	
	// 1
	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_ip[0]);
	fftw_free(fftw_in[0]);
	fftw_free(fftw_out[0]);


	if (param->overSamp != 32) {
		// 2
		fftw_destroy_plan(fftw_p[1]);
		fftw_destroy_plan(fftw_ip[1]);
		fftw_free(fftw_in[1]);
		fftw_free(fftw_out[1]);

		// 3
		fftw_destroy_plan(fftw_p[2]);
		fftw_destroy_plan(fftw_ip[2]);
		fftw_free(fftw_in[2]);
		fftw_free(fftw_out[2]);
	}
}
//---------------------------------------------------------------------------
// Function   : fftFilterSub
// Description: FFT によるフィルタ処理(サブ関数)
// ---
//	param		:変換パラメータ
//
void fftFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,PARAM_INFO *param,int id)
{
	long inSampleR,outSampleR;
	long wkSampleR;
	long fftSizeIn,fftSizeOut,i,n;
	long cutOff;
	long hfc,lfc;
	
	double p;
	long validIndex;
	long index15k,index19k,width1k;
	long h;
	double pw15k,pw19k;

	inSampleR = param->inSampleR;
	outSampleR = param->outSampleR;

	fftSizeIn = inSampleR * 2;
	fftSizeOut = outSampleR * 2;

	validIndex = ((double)fftSizeOut / param->outSampleR) * (inSampleR / 2);

	if (param->disRc) {
		inSampleR = outSampleR;
		fftSizeIn = fftSizeOut;
	}

	// FFT 初期設定
	copyToFFTW(fftw_in,pIn,fftSizeIn);

	windowFFTW(fftw_in,fftSizeIn);

	for (i = 0;i < fftSizeOut;i++) {
		fftw_out[i][0] = 0;
		fftw_out[i][1] = 0;
	}
	
	// FFT
	fftw_execute(fftw_p);

	if (param->disEq == 0) {
		if (param->deEmphasis != 0) {
			adjPinkFilter(3,fftSizeOut,fftw_out,param);
		}
	}

	// 高域削除
	if (inSampleR <= outSampleR) {
		wkSampleR = inSampleR;
	} else {
		wkSampleR = outSampleR;
	}
	hfc = wkSampleR / 2;
	if (param->disEq == 0 && param->hfc != -1) {
		hfc = param->hfc;
	}
	if (hfc > wkSampleR / 2) {
		hfc = wkSampleR / 2;
	}
	lfc = -1;
	if (param->disEq == 0 && param->lfc != -1) {
		if (param->lfc + 1000 < hfc) {
			lfc = param->lfc;
		}
	}

	// 周波数補正用データ作成
	if (param->disEq == 1 && id == 1) {
		for (i = 1;i < validIndex;i++) {
			p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
			if (p > 0) {
				p = sqrt(p);
			}
			if (param->ana_pw[i] < p) {
				param->ana_pw[i] = p;
			}
			param->ana_avgpw[i] += p;
		}
		param->ana_avgcount++;
	}

	if (param->disEq == 0) {
		// ハイパスディザのカット
		if (param->hiDither != 0) {
			index15k = ((double)(outSampleR * 2) / outSampleR) * 14000;
			width1k = ((double)(outSampleR * 2) / outSampleR) * 1000;
			pw15k = 0;
			for (i = index15k,n = 0;n < width1k;i++,n++) {
				p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
				if (p > 0) {
					p = sqrt(p);
				}
				pw15k += p;
			}
			if (n > 0) {
				pw15k /= n;
			}
			pw15k = pw15k - (pw15k / 10);
			index19k = ((double)(outSampleR * 2) / outSampleR) * 18000;
			pw19k = 0;
			for (i = index19k,n = 0;n < width1k;i++,n++) {
				p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
				if (p > 0) {
					p = sqrt(p);
				}
				pw19k += p;
			}
			if (n > 0) {
				pw19k /= n;
			}
			if (pw15k < pw19k) {
				for (i = index15k;i < outSampleR;i++) {
					pw15k = pw15k - (pw15k / 40000);
					p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
					if (p > 0) {
						p = sqrt(p);
						fftw_out[i][0] *= (pw15k / p);
						fftw_out[i][1] *= (pw15k / p);
					}
				}
			}
		}

		// lfa (低域調整)
		if (param->lfa == 1) {
			for (i = 1;i < validIndex;i++) {
				if (param->lfa_eq[i] != 0) {
					fftw_out[i][0] *= param->lfa_eq[i];
					fftw_out[i][1] *= param->lfa_eq[i];
				}
			}
		}

		if (param->lpf != -1) {
			adjPinkFilter(4,fftSizeOut,fftw_out,param);
		}

		// スピーカーごとの調整
		if (param->ana == 1 && param->sp_ana == 0) {
			for (i = 1;i < validIndex;i++) {
				h = ((double)param->outSampleR / fftSizeOut) * i;
				if (param->sp_eq[h] != 0) {
					fftw_out[i][0] *= (param->sp_eq[h] * 0.50);
					fftw_out[i][1] *= (param->sp_eq[h] * 0.50);
				}
			}
		}

		// イコライザ
		if (param->eq_flag == 1) {
			for (i = 1;i < validIndex;i++) {
				h = ((double)param->outSampleR / fftSizeOut) * i;
				if (param->eq[h] != 1) {
					fftw_out[i][0] *= param->eq[h];
					fftw_out[i][1] *= param->eq[h];
				}
			}
		}
	}

	// カットオフ(hfc)
	cutOff = ((double)fftSizeOut / outSampleR) * hfc;
	cutFFTW(fftw_out,cutOff,fftSizeOut);

	if (param->disEq == 0 && lfc != -1) {
		// カットオフ(lfc)
		cutOff = ((double)fftSizeOut / outSampleR) * lfc;
		for (i = 1;i < cutOff;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}

	// 半分のデータを復元
	for (i = 1;i < fftSizeOut / 2;i++) {
		fftw_out[fftSizeOut - i][0] = fftw_out[i][0];
		fftw_out[fftSizeOut - i][1] = fftw_out[i][1] * -1;
	}

	// invert FFT
	fftw_execute(fftw_ip);

	// 出力
	for (i = 0;i < fftSizeOut;i++) {
		pOut[i] = (SSIZE)(fftw_in[i][0] / fftSizeOut);
	}
}
//---------------------------------------------------------------------------
// Function   : anaLfaParam
// Description: LFA パラメーター作成
// ---
//	param		:変換パラメータ
//
void anaLfaParam(PARAM_INFO *param)
{
	long width_l;
	long index_b,index_l,index_h;
	long i,n;
	long fftSizeOut;
	long validIndex;
	double pw;
	double div;
	double rate;
	fftSizeOut = param->outSampleR * 2;
	validIndex = ((double)fftSizeOut / param->outSampleR) * (param->inSampleR / 2);

	for (i = 0;i < 192000;i++) {
		param->lfa_eq[i] = 0;
	}

#if 0
	// 20Hz ごとの平均をとる
	width_l = ((double)fftSizeOut / param->outSampleR) * 20;
	for (index_l = 0;index_l + width_l < validIndex;index_l += width_l) {
		pw = 0;
		for (i = index_l,n = 0;n < width_l;n++,i++) {
			pw += param->ana_pw[i];
		}
		if (n > 0) {
			pw /= n;
		}
		for (i = index_l,n = 0;n < width_l;n++,i++) {
			param->lfa_eq[i] = pw;
		}
	}

	// 20Hz から 30Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 15;
	index_h = ((double)fftSizeOut / param->outSampleR) * 30;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.05;
		}
	}
	// 30Hz から 40Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 30;
	index_h = ((double)fftSizeOut / param->outSampleR) * 40;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.07;
		}
	}
	// 40Hz から 50Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 40;
	index_h = ((double)fftSizeOut / param->outSampleR) * 50;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.09;
		}
	}
	// 50Hz から 60Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 50;
	index_h = ((double)fftSizeOut / param->outSampleR) * 60;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.11;
		}
	}
	// 60Hz から 80Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 60;
	index_h = ((double)fftSizeOut / param->outSampleR) * 80;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.09;
		}
	}
	// 80Hz から 100Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 80;
	index_h = ((double)fftSizeOut / param->outSampleR) * 100;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.07;
		}
	}
	// 100Hz から 120Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 100;
	index_h = ((double)fftSizeOut / param->outSampleR) * 120;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.05;
		}
	}
	// 120Hz から 150Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 120;
	index_h = ((double)fftSizeOut / param->outSampleR) * 150;
	index_b = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 1.03;
		}
	}
	for (;i < validIndex;i++) {
		param->lfa_eq[i] = 0;
	}

	// 10Hz から 20Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 0;
	index_h = ((double)fftSizeOut / param->outSampleR) * 20;
	div = index_h - index_l;
	div = (double)1.00 / div;
	rate = 1.00;
	for (i = index_h;i > index_l;i--) {
		param->lfa_eq[i] *= rate;
		if ((rate - div) > 0) {
			rate = rate - div;
			if (rate < 0.5) break;
		}
	}

	index_l = ((double)fftSizeOut / param->outSampleR) * 15;
	for (i = 1;i < index_l;i++) {
		param->lfa_eq[i] = 0;
	}
#else
	// 10Hz ごとの平均をとる
	width_l = ((double)fftSizeOut / param->outSampleR) * 10;
	for (index_l = 0;index_l + width_l < validIndex;index_l += width_l) {
		pw = 0;
		for (i = index_l,n = 0;n < width_l;n++,i++) {
			pw += param->ana_pw[i];
		}
		if (n > 0) {
			pw /= n;
		}
		for (i = index_l,n = 0;n < width_l;n++,i++) {
			param->lfa_eq[i] = pw;
		}
	}
	// 0Hz から 10Hz
	index_l = ((double)fftSizeOut / param->outSampleR) * 0;
	index_h = ((double)fftSizeOut / param->outSampleR) * 10;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 0;
	}

	// 10Hz から 20Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 10;
	index_h = ((double)fftSizeOut / param->outSampleR) * 20;
	index_b = ((double)fftSizeOut / param->outSampleR) * 50;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 0.80;
		}
	}
	// 20Hz から 30Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 20;
	index_h = ((double)fftSizeOut / param->outSampleR) * 30;
	index_b = ((double)fftSizeOut / param->outSampleR) * 50;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 0.90;
		}
	}
	// 30Hz から 40Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 30;
	index_h = ((double)fftSizeOut / param->outSampleR) * 40;
	index_b = ((double)fftSizeOut / param->outSampleR) * 55;
	for (i = index_l;i < index_h;i++) {
		if (param->lfa_eq[i] != 0) {
			param->lfa_eq[i] = param->lfa_eq[index_b] / param->lfa_eq[i];
			param->lfa_eq[i] *= 0.97;
		}
	}
	// 40Hz から 50Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 40;
	index_h = ((double)fftSizeOut / param->outSampleR) * 50;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.50;
	}
	// 50Hz から 60Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 50;
	index_h = ((double)fftSizeOut / param->outSampleR) * 60;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.45;
	}
	// 60Hz から 70Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 60;
	index_h = ((double)fftSizeOut / param->outSampleR) * 70;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.40;
	}
	// 70Hz から 80Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 70;
	index_h = ((double)fftSizeOut / param->outSampleR) * 80;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.40;
	}
	// 80Hz から 90Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 80;
	index_h = ((double)fftSizeOut / param->outSampleR) * 90;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.40;
	}
	// 90Hz から 100Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 90;
	index_h = ((double)fftSizeOut / param->outSampleR) * 100;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.40;
	}
	// 100Hz から 110Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 100;
	index_h = ((double)fftSizeOut / param->outSampleR) * 110;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.40;
	}
	// 110Hz から 120Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 110;
	index_h = ((double)fftSizeOut / param->outSampleR) * 120;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.30;
	}
	// 120Hz から 150Hz の調整
	index_l = ((double)fftSizeOut / param->outSampleR) * 120;
	index_h = ((double)fftSizeOut / param->outSampleR) * 150;
	for (i = index_l;i < index_h;i++) {
		param->lfa_eq[i] = 1.20;
	}
	for (;i < validIndex;i++) {
		param->lfa_eq[i] = 0;
	}
#endif
}
//---------------------------------------------------------------------------
// Function   : anaHFC_autoParam
// Description: HFC Auto パラメーター作成
// ---
//	param		:変換パラメータ
//
void anaHFC_AutoParam(PARAM_INFO *param)
{
	
	long i,n;
	static double avg[MAX_SAMPLE_N];
	static double pw[MAX_SAMPLE_N];
	double p;
	long fftSizeOut;
	long index;
	long width;
	long validIndex;
	

	for (i = 0;i < MAX_SAMPLE_N;i++) {
		avg[i] = 0;
		pw[i]  = 0;
	}

	fftSizeOut = param->outSampleR * 2;
	validIndex = ((double)fftSizeOut / param->outSampleR) * (param->inSampleR / 2);

	if (param->ana_avgcount > 0) {
		for (i = 0;i < validIndex;i++) {
			avg[i] = param->ana_avgpw[i] / param->ana_avgcount;
		}
	}

	// 200Hz ごとの平均をとる
	width = ((double)fftSizeOut / param->outSampleR) * 200;
	for (index = 0;index + width < validIndex;index += width) {
		p = 0;
		for (i = index,n = 0;n < width;n++,i++) {
			p += param->ana_pw[i];
		}
		if (n > 0) {
			p /= n;
		}
		for (i = index,n = 0;n < width;n++,i++) {
			pw[i] = p;
		}
	}

	for (i = validIndex - 1;i > 0;i--) {
		if (pw[i] < (double)30000000000000000) {
			pw[i] = 0;
		} else {
			break;
		}
	}

	for (i = validIndex - 1;i > 0;i--) {
		if (avg[i] < (double)12000000000000000) {
			pw[i] = 0;
		} else {
			break;
		}
	}

	for (i = validIndex - 1;i > 0;i--) {
		if (pw[i] != 0) {
			break;
		}
	}

	param->hfc = ((double)param->outSampleR / fftSizeOut) * (i);

}

//---------------------------------------------------------------------------
// Function   : adjBitExtension
// Description: ビット分解能を高める処理
// ---
//	inSample 	:処理するサンプル数
//	fp			:入力ファイル
//	tmpFp		:出力ファイル
//	param		:変換パラメータ
//
void adjBitExtension(DWORD inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	SSIZE level,level2;
	double nx;
	long wkMemSize;
	long inSampleR;
	long fftSize,i,j,k,n;
	long cnt;
	long cutOff;
	long hfc;
	long wr;
	double persent,per;
	SSIZE *pIn,*pIn2,*pOut;
	SSIZE d1,d2,d3,d4,d5,d6,d7,d8;
	SSIZE dd1,dd2,dd3,dd4;
	SSIZE startInSample,nSample;
	DWORD outRemain;
	double samp[9];
	double power;
	double shresh;
	int next;
	fftw_complex *fftw_in,*fftw_out;
	fftw_plan fftw_p,fftw_ip;
	int sm_ignore;
	SSIZE sm_avg;

	fio_rewind(fp_r);

	inSampleR = param->inSampleR;

	fftSize = inSampleR * 2;
	wkMemSize = fftSize * 2;

	mem1 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem4 = (SSIZE *)al_malloc(wkMemSize * sizeof (SSIZE));
	if (mem4 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_in == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * wkMemSize);
	if (fftw_out == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p = fftw_plan_dft_1d(fftSize,fftw_in,fftw_out,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip = fftw_plan_dft_1d(fftSize,fftw_out,fftw_in,FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	if (param->iw == 16) {
		shresh = (double)12592500000000000.0;	// 1k
		shresh = (double)4200000000000000.0;	// 1k
	} else if (param->iw == 24) {
		shresh = (double)85120000000000.0;	// 1k
	} else {
		shresh = -1;
	}
	outRemain = inSample;
	per = -1;
	for (startInSample = (((fftSize * 2) + (fftSize / 2)) * -1);startInSample < inSample + (fftSize * 1);startInSample += fftSize) {
		if (startInSample >= 0 && startInSample  < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
		}

		memset(mem1,0,wkMemSize * sizeof (SSIZE));
		nSample = fftSize * 2;

		if (startInSample >= 0 && startInSample + (fftSize * 1) < inSample) {
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			nSample = fftSize * 2;
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else if (startInSample + (fftSize * 1) >= 0 || startInSample < inSample) {
			mem0 = mem1;
			nSample = 0;
			if (startInSample < 0) {
				fio_seek(fp_r,0,SEEK_SET);
				if ((startInSample * -1) < fftSize * 2) {
					mem0 += (startInSample * -1);
					nSample = (fftSize * 2) + startInSample;
				}
			} else if (startInSample < inSample) {
				fio_seek(fp_r,startInSample,SEEK_SET);
			}
			if (nSample > 0 && startInSample < inSample && startInSample + nSample > inSample) {
				nSample = inSample - startInSample;
			}

			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}
		
		// 音のレベルを調査しておく
		if (param->abeNX) {
			level = 0;
			for (i = fftSize / 2,j = 0,n = 0;n < fftSize;i++,j++,n++) {
				if (mem1[i] > 0) {
					level += mem1[i] >> (56 - 16);
					j++;
				}
			}
			if (j > 0) {
				level /= j;
			}
		}
		if (param->abeFnLevel > 0) {
			// ディザキャンセラ
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			// ２値平均
			for (i = 0;i + 2 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
//				if (dd1 < 0) {
//					dd1 *= -1;
//				}
//				if (dd2 < 0) {
//					dd2 *= -1;
//				}
//				if (dd3 < 0) {
//					dd3 *= -1;
//				}
				dd1 /= param->abeFnLevel;
				dd2 /= param->abeFnLevel;
				dd3 /= param->abeFnLevel;
				sm_ignore = 1;
				if (dd1 + 1 == dd2 && dd2 == dd3 + 1) {
					sm_ignore = 0;
				}
				if (dd1 == dd2 + 1 && dd2 + 1 == dd3) {
					sm_ignore = 0;
				}

				if (sm_ignore == 0) {
					sm_avg = 0;
					sm_avg += pIn[i + 0];
					sm_avg += pIn[i + 1];
					sm_avg /= 2;
					pOut[i] = sm_avg;
				}
			}
			if (param->abeFnLevel > 1) {
				// 3値平均
				for (i = nSample - 1;i + 1 > 0;i--) {
					d1 = pIn[i - 0];
					d2 = pIn[i - 1];
					d3 = pIn[i - 2];
					dd1 = d1;
					dd2 = d2;
					dd3 = d3;
					dd1 >>= (56 - param->iw);
					dd2 >>= (56 - param->iw);
					dd3 >>= (56 - param->iw);
//					if (dd1 < 0) {
//						dd1 *= -1;
//					}
//					if (dd2 < 0) {
//						dd2 *= -1;
//					}
//					if (dd3 < 0) {
//						dd3 *= -1;
//					}
					dd1 /= param->abeFnLevel;
					dd2 /= param->abeFnLevel;
					dd3 /= param->abeFnLevel;
					sm_ignore = 1;
					if (dd1 + 1 == dd2 && dd2 == dd3 + 1) {
						sm_ignore = 0;
					}
					if (dd1 == dd2 + 1 && dd2 + 1 == dd3) {
						sm_ignore = 0;
					}

					if (sm_ignore == 0) {
						sm_avg = 0;
						sm_avg += pIn[i - 0];
						sm_avg += pIn[i - 1];
						sm_avg += pIn[i - 2];
						sm_avg /= 3;
						pOut[i] = sm_avg;
					}
				}
			}
			if (param->abeFnLevel > 3) {
				// 3値平均
				for (i = 0;i + 1 < nSample;i++) {
					d1 = pIn[i + 0];
					d2 = pIn[i + 1];
					d3 = pIn[i + 2];
					dd1 = d1;
					dd2 = d2;
					dd3 = d3;
					dd1 >>= (56 - param->iw);
					dd2 >>= (56 - param->iw);
					dd3 >>= (56 - param->iw);
//					if (dd1 < 0) {
//						dd1 *= -1;
//					}
//					if (dd2 < 0) {
//						dd2 *= -1;
//					}
//					if (dd3 < 0) {
//						dd3 *= -1;
//					}
					dd1 /= param->abeFnLevel;
					dd2 /= param->abeFnLevel;
					dd3 /= param->abeFnLevel;
					sm_ignore = 1;
					if (dd1 + 1 == dd2 && dd2 == dd3 + 1) {
						sm_ignore = 0;
					}
					if (dd1 == dd2 + 1 && dd2 + 1 == dd3) {
						sm_ignore = 0;
					}

					if (sm_ignore == 0) {
						sm_avg = 0;
						sm_avg += pIn[i + 0];
						sm_avg += pIn[i + 1];
						sm_avg += pIn[i + 2];
						sm_avg /= 3;
						pOut[i] = sm_avg;
					}
				}
			}
		}

		if (param->smLowData > 0) {
			// 2値同値でその左右隣が異なる値の調整
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = 0;i + 3 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				d4 = pIn[i + 3];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd4 = d4;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				dd4 >>= (56 - param->iw);
				if (dd1 != dd2 && dd2 == dd3 && dd3 != dd4) {
					sm_ignore = 0;
					if ((dd1 > dd2 && (dd1 - dd2) > 2) || (dd1 < dd2 && (dd2 - dd1) > 2)) {
						sm_ignore = 1;
					}
					if ((dd3 > dd4 && (dd3 - dd4) > 2) || (dd3 < dd4 && (dd4 - dd3) > 2)) {
						sm_ignore = 1;
					}
					if (sm_ignore == 0) {
						sm_avg = (d1 + d2 + d3) / 3;
						pOut[i + 1] = sm_avg;
						sm_avg = (d2 + d3 + d4) / 3;
						pOut[i + 2] = sm_avg;
						i++;
					}
				}
			}
			// 2値同値でその左右隣が異なる値の調整(逆順)
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = nSample - 1;i > 2;i--) {
				d1 = pIn[i - 0];
				d2 = pIn[i - 1];
				d3 = pIn[i - 2];
				d4 = pIn[i - 3];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd4 = d4;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				dd4 >>= (56 - param->iw);
				if (dd1 != dd2 && dd2 == dd3 && dd3 != dd4) {
					sm_ignore = 0;
					if ((dd1 > dd2 && (dd1 - dd2) > 2) || (dd1 < dd2 && (dd2 - dd1) > 2)) {
						sm_ignore = 1;
					}
					if ((dd3 > dd4 && (dd3 - dd4) > 2) || (dd3 < dd4 && (dd4 - dd3) > 2)) {
						sm_ignore = 1;
					}
					if (sm_ignore == 0) {
						sm_avg = (d1 + d2 + d3) / 3;
						pOut[i + 1] = sm_avg;
						sm_avg = (d2 + d3 + d4) / 3;
						pOut[i + 2] = sm_avg;
						i--;
					}
				}
			}
			// 山や谷の形の波形調整
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = 0;i + 2 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (d1 < d2 && d2 > d3 && (dd2 - dd1) <= 2 && (dd2 - dd3) <= 2) {
					// 山
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				} else if (d1 > d2 && d2 < d3 && (dd1 - dd2) <= 2 && (dd3 - dd2) <= 2) {
					// 谷
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				}
			}
			// 山や谷の形の波形調整(逆順)
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = nSample - 1;i > 1;i--) {
				d1 = pIn[i - 0];
				d2 = pIn[i - 1];
				d3 = pIn[i - 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (d1 < d2 && d2 > d3 && (dd2 - dd1) <= 2 && (dd2 - dd3) <= 2) {
					// 山
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				} else if (d1 > d2 && d2 < d3 && (dd1 - dd2) <= 2 && (dd3 - dd2) <= 2) {
					// 谷
					sm_avg = (((d1 + d3) / 2) + d2) / 2;
					pOut[i + 1] = sm_avg;
				}
			}
			// 同値以外の移動平均
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = 0;i + 2 < nSample;i++) {
				d1 = pIn[i + 0];
				d2 = pIn[i + 1];
				d3 = pIn[i + 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (dd1 < 0) {
					dd1 *= -1;
				}
				if (dd2 < 0) {
					dd2 *= -1;
				}
				if (dd3 < 0) {
					dd3 *= -1;
				}
				dd1 >>= 2;
				dd2 >>= 2;
				dd3 >>= 2;
				sm_ignore = 1;
				if (dd1 == dd2 && dd2 == dd3) {
					sm_ignore = 0;
				}
				if (sm_ignore == 0) {
					sm_avg = (d1 + d2 + d3) / 3;
					pOut[i + 1] = sm_avg;
				}
			}
			// 同値以外の移動平均(逆順)
			pIn  = (SSIZE *)mem1;
			pOut = (SSIZE *)mem1;
			for (i = nSample - 1;i > 1;i--) {
				d1 = pIn[i - 0];
				d2 = pIn[i - 1];
				d3 = pIn[i - 2];
				dd1 = d1;
				dd2 = d2;
				dd3 = d3;
				dd1 >>= (56 - param->iw);
				dd2 >>= (56 - param->iw);
				dd3 >>= (56 - param->iw);
				if (dd1 < 0) {
					dd1 *= -1;
				}
				if (dd2 < 0) {
					dd2 *= -1;
				}
				if (dd3 < 0) {
					dd3 *= -1;
				}
				dd1 >>= 3;
				dd2 >>= 3;
				dd3 >>= 3;
				sm_ignore = 1;
				if (dd1 == dd2 && dd2 == dd3) {
					sm_ignore = 0;
				}
				if (sm_ignore == 0) {
					sm_avg = (d1 + d2 + d3) / 3;
					pOut[i + 1] = sm_avg;
				}
			}
		}
#if 1
		// 適応型フィルター処理
		// 同値が続く場合に同値の個数に応じたfftで高域カットフィルターをする。
		// フィルター後は波形がなめらかになるように制御する。
		memset(mem2,0,wkMemSize * sizeof (SSIZE));
		memset(mem3,0,wkMemSize * sizeof (SSIZE));
		memset(mem4,0,wkMemSize * sizeof (SSIZE));
		if (param->adaptiveFilder == 1) {
			pIn = mem1;							// 入力波形
			pOut = mem2;						// 同値個数
			// mem2に同値が続く場合に同値の個数を記録していく。
			for (i = 0,j = 1,cnt = 1;j < nSample;j++) {
				d1 = pIn[i] >> (56 - param->iw);
				d2 = pIn[j] >> (56 - param->iw);
				if (d1 == d2 && cnt < 255) {
					cnt++;
				} else {
					for (k = i;k < j;k++) {
						pOut[k] = cnt;
					}
					i = j;
					cnt = 1;
				}
			}
			for (i = 0;i < nSample;i++) {
				if (pOut[i] >= 0 && pOut[i] < 2) {
					pOut[i] = 0;
				}
			}
			// 同値が3つ以上続くものに対して、fft をかける
			do {
				pIn = mem1;
				pIn2 = mem2;
				memset(mem3,0,wkMemSize * sizeof (SSIZE));
				cnt = 0;
				for (i = 0;i < nSample;i++) {
					if (pIn2[i] > 0) {
						cnt = pIn2[i];
						break;
					}
				}
				if (cnt == 0) {
					break;
				}
				for (n = 0;n < 3;n++) {
					// FFT 初期設定
					copyToFFTW(fftw_in,&pIn[((fftSize / 2) * n)],fftSize);
					
					// 窓関数
					windowFFTW(fftw_in,fftSize);

					// FFT
					fftw_execute(fftw_p);

					if (param->beInfo.cutOff[cnt] > 0) {
						// 高域削除
						hfc = inSampleR / param->beInfo.cutOff[cnt];
						cutOff = ((double)fftSize / inSampleR) * hfc;
						for (i = cutOff;i < fftSize;i++) {
							fftw_out[i][0] = 0;
							fftw_out[i][1] = 0;
						}
					}
					// 半分のデータを復元
					//#pragma omp parallel for
					for (i = 1;i < fftSize / 2;i++) {
						fftw_out[fftSize - i][0] = fftw_out[i][0];
						fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
					}

					// invert FFT
					fftw_execute(fftw_ip);

					// 出力
					pOut = (SSIZE *)&mem3[((fftSize / 2) * n)];
					//#pragma omp parallel for
					for (i = 0;i < fftSize;i++) {
						pOut[i] += fftw_in[i][0] / fftSize;
					}
				}
				pIn   = (SSIZE *)mem3;
				pIn2  = (SSIZE *)mem2;
				pOut  = (SSIZE *)mem4;
				for (i = 0;i < nSample;i++) {
					if (pIn2[i] > 0 && pIn2[i] == cnt) {
						pOut[i] = pIn[i];
						pIn2[i] *= -1;
					}
				}
			} while (1);
			pIn = (SSIZE *)mem4;
			pIn2 = (SSIZE *)mem2;
			pOut = (SSIZE *)mem1;
			for (i = 0;i < nSample;i++) {
				d1 = pIn[i];
				d1 >>= (56 - param->iw);
				d2 = pOut[i];
				d2 >>= (56 - param->iw);
				if (pIn2[i] < 0) {
					if (d1 <= d2 && (d2 - d1) <= 3) {
						pOut[i] = pIn[i];
					} else if (d1 > d2 && (d1 - d2) <= 3) {
						pOut[i] = pIn[i];
					}
				}
			}
		}
#endif
		if (param->cutLowData && shresh > 0) {
			// 閾値より低いパワーの音は削除する
			pIn = mem1;
			pOut = mem3;
			memset(mem3,0,wkMemSize * sizeof (SSIZE));
			for (n = 0;n < 3;n++) {
				// FFT 初期設定
				copyToFFTW(fftw_in,&pIn[((fftSize / 2) * n)],fftSize);

				// 窓関数
				windowFFTW(fftw_in,fftSize);

				// FFT
				fftw_execute(fftw_p);

				// 削除
				for (i = 1;i < fftSize / 2;i++) {
					power = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
					if (power > 0) {
						power = sqrt(power);
					}
					if (power < shresh) {
						fftw_out[i][0] = 0;
						fftw_out[i][1] = 0;
					}
				}

				// 半分のデータを復元
				#pragma omp parallel for
				for (i = 1;i < fftSize / 2;i++) {
					fftw_out[fftSize - i][0] = fftw_out[i][0];
					fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
				}

				// invert FFT
				fftw_execute(fftw_ip);

				// 出力
				pOut = (SSIZE *)mem3;
				pOut = &pOut[(fftSize / 2) * n];
				#pragma omp parallel for
				for (i = 0;i < fftSize;i++) {
					pOut[i] += fftw_in[i][0] / fftSize;
				}
			}
			memcpy(mem1,mem3,wkMemSize * sizeof (SSIZE));
		}

#if 1
		// 波形の調整処理
		// 上がりや、下がりの波形の量子化誤差を少なくする
		pIn  = (SSIZE *)mem1;
		pOut = (SSIZE *)mem1;
		for (i = 4;i + 3 < nSample;) {
			next = 1;
			d1 = pIn[i - 4];
			d2 = pIn[i - 3];
			d3 = pIn[i - 2];
			d4 = pIn[i - 1];
			d5 = pIn[i];
			d6 = pIn[i + 1];
			d7 = pIn[i + 2];
			d8 = pIn[i + 3];
			if ((d2 < d3 && d3 <= d4 && d4 < d5 && d5 <= d6 && d6 < d7) ||
						(d2 <= d3 && d3 < d4 && d4 <= d5 && d5 < d6 && d6 <= d7)) {
				// 上がり波形
				samp[1] = pIn[i - 2] - pIn[i - 3];
				samp[2] = pIn[i - 1] - pIn[i - 2];
				samp[3] = pIn[i] - pIn[i - 1];
				samp[4] = pIn[i + 1] - pIn[i];
				samp[5] = pIn[i + 2] - pIn[i + 1];
				for (j = 1;j < 5;j++) {
					for (k = j + 1;k < 6;k++) {
						if (samp[j] > samp[k]) {
							samp[8] = samp[j];
							samp[j] = samp[k];
							samp[k] = samp[8];
						}
					}
				}
				samp[2] = samp[2] + samp[3] + samp[4];
				samp[2] /= 3;
				d1 = pIn[i];
				d2 = pIn[i - 1] + (SSIZE)samp[2];
				d1 >>= (56 - param->iw);
				d2 >>= (56 - param->iw);
				if (d1 == d2) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 < d2 && (d2 - d1) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 > d2 && (d1 - d2) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				}

			} else if ((d2 > d3 && d3 >= d4 && d4 > d5 && d5 >= d6 && d6 > d7) ||
						(d2 >= d3 && d3 > d4 && d4 >= d5 && d5 > d6 && d6 >= d7)) {
				// 下がり波形
				samp[1] = pIn[i - 2] - pIn[i - 3];
				samp[2] = pIn[i - 1] - pIn[i - 2];
				samp[3] = pIn[i] - pIn[i - 1];
				samp[4] = pIn[i + 1] - pIn[i];
				samp[5] = pIn[i + 2] - pIn[i + 1];
				for (j = 1;j < 5;j++) {
					for (k = j + 1;k < 6;k++) {
						if (samp[j] > samp[k]) {
							samp[8] = samp[j];
							samp[j] = samp[k];
							samp[k] = samp[8];
						}
					}
				}
				samp[2] = samp[2] + samp[3] + samp[4];
				samp[2] /= 3;
				d1 = pIn[i];
				d2 = pIn[i - 1] + (SSIZE)samp[2];
				d1 >>= (56 - param->iw);
				d2 >>= (56 - param->iw);
				if (d1 == d2) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 < d2 && (d2 - d1) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				} else if (d1 > d2 && (d1 - d2) <= 3) {
					pOut[i] = pIn[i - 1] + (SSIZE)samp[2];
				}
			}
			i += next;
		}
#endif

		if (startInSample + (fftSize / 2) >= 0) {
#if 0
			// 音のレベルに変化がないか調査
			if (param->abeNX == 1) {
				level2 = 0;
				for (i = fftSize / 2,j = 0,n = 0;n < fftSize;i++,j++,n++) {
					if (mem1[i] > 0) {
						level2 += mem1[i] >> (56 - 16);
						j++;
					}
				}
				if (j > 0) {
					level2 /= j;
				}
				if (level > 0 && level2 > 0) {
					nx = ((double)level / (double)level2);
					for (i = fftSize / 2,n = 0;n < fftSize;i++,n++) {
						mem1[i] = (SSIZE)((double)mem1[i] * nx);
					}
				}
			}
#endif
			if (outRemain >= fftSize) {
				wr = fio_write(mem1 + (fftSize / 2),sizeof (SSIZE),fftSize,fp_w);
				if (wr != fftSize) {
					param->err = STATUS_FILE_WRITE_ERR;
					return;
				}
			} else {
				wr = fio_write(mem1 + (fftSize / 2),sizeof (SSIZE),outRemain,fp_w);
				if (wr != outRemain) {
					param->err = STATUS_FILE_WRITE_ERR;
					return;
				}
			}
			if (outRemain >= fftSize) {
				outRemain -= fftSize;
			} else {
				break;
			}
		}
	}
	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);

	fftw_destroy_plan(fftw_p);
	fftw_destroy_plan(fftw_ip);
	fftw_free(fftw_in);
	fftw_free(fftw_out);
}
//---------------------------------------------------------------------------
// Function   : genNoise
// Description: 失われた高域の再現処理(正規分布のノイズ付加)
// ---
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域にデータを追加する)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
void genNoise(long hfc,DWORD inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3;
	long outSampleR;
	long wkMemSize;
	long fftSize,i,j,n,nn;
	
	long lowIndex,highIndex;
	long nrIndex;
	double persent,per;
	double nx;
	double p,refPw,noisePw;
	SSIZE *pIn,*pIn2,*pOut;
	SSIZE startInSample,nSample;
	
	fftw_complex *fftw_in[2],*fftw_out[2],*fftw_in2;
	fftw_plan fftw_p[2],fftw_ip[2],fftw_p2;
	double hfaNB;
	
//	// ノイズ生成用データの読み込み
//	fpNoise = fopen(param->nd_path,"rb");
//	if (fpNoise == NULL) {
//		// データがなければ作成する
//		fpNoise = fopen(param->nd_path,"wb");
//		if (fpNoise) {
//			for (i = 0;i < (long)(MAX_SAMPLE_N) * 60;i++) {
//				nrData = normalNoise() * 0x100;
//				fwrite(&nrData,1,sizeof (long),fpNoise);
//			}
//			fflush(fpNoise);
//			fclose(fpNoise);
//		}
//		fpNoise = fopen(param->nd_path,"rb");
//	}
	outSampleR = param->outSampleR;
	fio_rewind(fp_r);

	fftSize = param->outSampleR * 2;

	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	mem1 = (SSIZE *)al_malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize);
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize);
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize / 100);
	if (fftw_in2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_p[0] = fftw_plan_dft_1d(fftSize,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p[1] = fftw_plan_dft_1d(fftSize,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSize,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSize,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p2 = fftw_plan_dft_1d(fftSize / 100,fftw_in2,fftw_in2,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	per = -1;
	nrIndex = 0;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + (fftSize + (fftSize / 2));startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
		}

		memset(mem1,0,wkMemSize);
		memset(mem3,0,wkMemSize);

		// mem2 には正規分布のノイズを格納する
		pIn2 = (SSIZE *)mem2;
		#pragma omp parallel for
		for (i = 0;i < fftSize * 2;i++) {
			pIn2[i] = normalNoise() * 0x100;
		}

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		pIn = (SSIZE *)mem1;
		for (n = 0;n < 3;n++) {
			lowIndex = ((double)fftSize / outSampleR) * (hfc - 2000);
			highIndex = ((double)fftSize / outSampleR) * hfc;
			// FFT 初期設定
			copyToFFTW(fftw_in[0],&pIn[((fftSize / 2) * n)],fftSize);
			windowFFTW(fftw_in[0],fftSize);

			// FFT
			fftw_execute(fftw_p[0]);

			// 元信号の高域のパワーを調べる
			refPw = 0;
			for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
				p = fftw_out[0][i][0] * fftw_out[0][i][0] + fftw_out[0][i][1] * fftw_out[0][i][1];
				if (p != 0) {
					p = sqrt(p);
				}
				refPw += p;
			}
			if (nn > 0) {
				refPw /= nn;
			}

			// 付加する信号
			copyToFFTW(fftw_in[1],&pIn2[((fftSize / 2) * n)],fftSize);
			windowFFTW(fftw_in[1],fftSize);

			// FFT
			fftw_execute(fftw_p[1]);
			// 1/f 信号にする
			adjPinkFilter(1,fftSize,fftw_out[1],param);

			// 付加する信号のパワーを調べる
			noisePw = 0;
			lowIndex = ((double)fftSize / outSampleR) * hfc;
			highIndex = ((double)fftSize / outSampleR) * (hfc + 2000);
			for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
				p = fftw_out[1][i][0] * fftw_out[1][i][0] + fftw_out[1][i][1] * fftw_out[1][i][1];
				if (p != 0) {
					p = sqrt(p);
				}
				noisePw += p;
			}
			if (nn > 0) {
				noisePw /= nn;
			}
			if (noisePw == 0) {
				noisePw = 1;
			}
			nx = refPw / noisePw;
			highIndex = ((double)fftSize / outSampleR) * hfc;

			#pragma omp parallel for
			for (i = highIndex;i < fftSize / 2;i++) {
				fftw_out[1][i][0] *= nx;
				fftw_out[1][i][1] *= nx;
			}

			if (param->hfa != 1 && param->hfaNB > 0) {
				hfaNB = param->hfaNB / 100.0;
				// hfa1 信号
				for (i = highIndex;i < fftSize / 2;i++) {
					fftw_out[1][i][0] *= hfaNB;
					fftw_out[1][i][1] *= hfaNB;
				}
				hfaNB = 1.0 - hfaNB;
				// hfa2 信号
				for (i = highIndex;i < fftSize / 2;i++) {
					fftw_out[0][i][0] *= hfaNB;
					fftw_out[0][i][1] *= hfaNB;
				}
				// 合成
				for (i = highIndex;i < fftSize / 2;i++) {
					fftw_out[1][i][0] += fftw_out[0][i][0];
					fftw_out[1][i][1] += fftw_out[0][i][1];
				}
			}

			// 低域カット
			highIndex = ((double)fftSize / outSampleR) * hfc;

			#pragma omp parallel for
			for (i = 1;i < highIndex;i++) {
				fftw_out[1][i][0] = 0;
				fftw_out[1][i][1] = 0;
			}

			// 半分のデータを復元
			#pragma omp parallel for
			for (i = 1;i < fftSize / 2;i++) {
				fftw_out[1][fftSize - i][0] = fftw_out[1][i][0];
				fftw_out[1][fftSize - i][1] = fftw_out[1][i][1] * -1;
			}
			// 直流部分除去
			fftw_out[1][0][0] = 0;
			fftw_out[1][0][1] = 0;

			// invert FFT
			fftw_execute(fftw_ip[1]);

			pOut = (SSIZE *)mem3;
			pOut = &pOut[(fftSize / 2) * n];
			#pragma omp parallel for
			for (i = 0;i < fftSize;i++) {
				pOut[i] += fftw_in[1][i][0] / fftSize;
			}
		}
		if (startInSample + fftSize / 2 >= 0) {
			if (param->hfa == 1) {
				pIn = (SSIZE *)mem1;
			} else {
				pIn = (SSIZE *)mem3;
			}
			pIn = pIn + (fftSize / 2);
			pOut = (SSIZE *)mem3;
			pOut = pOut + (fftSize / 2);
			for (i = 0;i < 100;i++) {
				// 元信号に音があるかを調べる
				for (j = 0;j < fftSize / 100;j++) {
					fftw_in2[j][0] = pIn[j];
					fftw_in2[j][1] = 0;
				}
				// 窓関数
				for (j = 0;j < ((fftSize / 100) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 100));
				}
				for (j = ((fftSize / 100) - 1) / 2;j < (fftSize / 100);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 100));
				}

				// FFT
				fftw_execute(fftw_p2);

				// 元信号の高域のパワーを調べる
				refPw = 0;
				lowIndex = (((double)fftSize / 100) / outSampleR) * (hfc - 2000);
				highIndex = (((double)fftSize / 100) / outSampleR) * hfc;
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					refPw += p;
				}
				if (nn > 0) {
					refPw /= nn;
				}
				// 付加信号の高域のパワーを調べる
				for (j = 0;j < fftSize / 100;j++) {
					fftw_in2[j][0] = pOut[j];
					fftw_in2[j][1] = 0;
				}
				// 窓関数
				for (j = 0;j < ((fftSize / 100) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 100));
				}
				for (j = ((fftSize / 100) - 1) / 2;j < (fftSize / 100);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 100));
				}

				// FFT
				fftw_execute(fftw_p2);

				noisePw = 0;
				lowIndex = (((double)fftSize / 100) / outSampleR) * hfc;
				highIndex = (((double)fftSize / 100) / outSampleR) * (hfc + 2000);
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					noisePw += p;
				}
				if (nn > 0) {
					noisePw /= nn;
				}
				if (refPw > 0) {
					if (param->hfa == 1) {
						for (j = 0;j < fftSize / 100;j++) {
							pOut[j] *= ((refPw / noisePw) * 0.70);
						}
					} else {
						for (j = 0;j < fftSize / 100;j++) {
							pOut[j] *= (refPw / noisePw);
						}
					}
				} else {
					for (j = 0;j < fftSize / 100;j++) {
						pOut[j] = 0;
					}
				}
				pIn  += fftSize / 100;
				pOut += fftSize / 100;
			}
			// 再び低域カット処理をする
			pIn = (SSIZE *)mem1;
			pOut = (SSIZE *)mem3;
			for (n = 0;n < 3;n++) {
				// FFT 初期設定
				for (i = 0;i < fftSize;i++) {
					fftw_in[1][i][0] = pOut[((fftSize / 2) * n) + i];
				   	fftw_in[1][i][1] = 0;
				}
				// 窓関数
				for (i = 0;i < (fftSize - 1) / 2;i++) {
					fftw_in[1][i][0] = fftw_in[1][i][0] * (2.0 * i / (double)fftSize);
				}
				for (i = (fftSize - 1) / 2;i < fftSize;i++) {
					fftw_in[1][i][0] = fftw_in[1][i][0] * (2.0 - 2.0 * i / (double)fftSize);
				}

				// FFT
				fftw_execute(fftw_p[1]);

				highIndex = ((double)fftSize / outSampleR) * (hfc - 1000);

				// 低域カット
				for (i = 1;i < highIndex;i++) {
					fftw_out[1][i][0] = 0;
					fftw_out[1][i][1] = 0;
				}

				// 半分のデータを復元
				for (i = 1;i < fftSize / 2;i++) {
					fftw_out[1][fftSize - i][0] = fftw_out[1][i][0];
					fftw_out[1][fftSize - i][1] = fftw_out[1][i][1] * -1;
				}

				// 直流部分除去
				fftw_out[1][0][0] = 0;
				fftw_out[1][0][1] = 0;

				// invert FFT
				fftw_execute(fftw_ip[1]);

				// 出力
				for (i = 0;i < fftSize;i++) {
					pOut[((fftSize / 2) * n) + i] += fftw_in[1][i][0] / fftSize;
				}
			}
			for (i = 0;i < fftSize;i++) {
				pOut[(fftSize / 2) + i] += pIn[(fftSize / 2) + i];
			}
			fio_seek(fp_w,(startInSample + (fftSize / 2)) * sizeof (SSIZE),SEEK_SET);
			outTempFile(fp_w,mem3 + fftSize / 2,fftSize,param);
			if (param->err) {
				break;
			}
		}
	}

	fio_flush(fp_w);

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_ip[0]);
	fftw_destroy_plan(fftw_ip[1]);
	fftw_destroy_plan(fftw_p2);
	fftw_free(fftw_in[0]);
	fftw_free(fftw_in[1]);
	fftw_free(fftw_out[0]);
	fftw_free(fftw_out[1]);
	fftw_free(fftw_in2);

}
//---------------------------------------------------------------------------
// Function   : genOverTone
// Description: 失われた高域の再現処理(倍音解析)
// ---
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域にデータを追加する)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
void genOverTone(long hfc,DWORD inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	long outSampleR;
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	long wkMemSize;
	long fftSize,i,j,k,nn;
	long lowIndex,highIndex;
	double persent,per;
	double nx;
	double p,refPw,ovTonePw;
	double avg;
	SSIZE *pIn[4],*pOut[4];
	SSIZE startInSample,nSample;
	SSIZE s[2];
	SSIZE p1,p2,p3;
	fftw_complex *fftw_in[3],*fftw_out[3],*fftw_in2;
	fftw_plan fftw_p[3],fftw_ip[3],fftw_p2;
	OVERTONE_INFO *ovInfo[3];
	double *phaseX,*phaseY;
	outSampleR = param->outSampleR;
	fio_rewind(fp_r);
	fio_rewind(fp_w);
param->hfa3_max = 0;
	fftSize = 4096;
	if (outSampleR == 44100 * 2 || outSampleR == 48000 * 2) {
		fftSize = outSampleR / 10;
		fftSize = 4096 * 2;
	}
	if (outSampleR == 32000 * 6 || outSampleR == 44100 * 4 || outSampleR == 48000 * 4) {
		fftSize = outSampleR / 10;
		fftSize = 8192 * 2;
	}
	if (outSampleR == 32000 * 12 || outSampleR == 44100 * 8 || outSampleR == 48000 * 8) {
		fftSize = outSampleR / 10;
		fftSize = 16384 * 2;
	}
	if (outSampleR == 32000 * 24 || outSampleR == 44100 * 16 || outSampleR == 48000 * 16) {
		fftSize = 32768 * 2;
	}
	if (outSampleR == 32000 * 48 || outSampleR == 44100 * 32 || outSampleR == 48000 * 32) {
		fftSize = 65536 * 2;
	}
//fftSize = 4096;
	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	mem1 = (SSIZE *)al_malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize);
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize);
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem4 = (SSIZE *)al_malloc(wkMemSize);
	if (mem4 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	phaseX = (double *)al_malloc(65536 * 2 * sizeof (double));
	if (phaseX == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	phaseY = (double *)al_malloc(65536 * 2 * sizeof (double));
	if (phaseY == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	ovInfo[0] = (OVERTONE_INFO *)malloc(sizeof (OVERTONE_INFO));
	if (ovInfo[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	ovInfo[1] = (OVERTONE_INFO *)malloc(sizeof (OVERTONE_INFO));
	if (ovInfo[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	ovInfo[2] = (OVERTONE_INFO *)malloc(sizeof (OVERTONE_INFO));
	if (ovInfo[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 1
	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_out[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 2
	fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_out[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	// 3
	fftw_in[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_out[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 1
	fftw_p[0] = fftw_plan_dft_1d(fftSize,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[0] = fftw_plan_dft_1d(fftSize,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 2
	fftw_p[1] = fftw_plan_dft_1d(fftSize,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSize,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	// 3
	fftw_p[2] = fftw_plan_dft_1d(fftSize,fftw_in[2],fftw_out[2],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip[2] = fftw_plan_dft_1d(fftSize,fftw_out[2],fftw_in[2],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[2] == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536 * 2);
	if (fftw_in2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p2 = fftw_plan_dft_1d(fftSize/12,fftw_in2,fftw_in2,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	avg = 0;
	for (i = 1;i < fftSize / 2;i++) {
		double t = ((2 * M_PI) / ((fftSize / 2) - 1)) * i;
		phaseX[i] = sin(t) * 1000000;
		phaseY[i] = cos(t) * 1000000;
		p = (phaseX[i] * phaseX[i]) + (phaseY[i] * phaseY[i]);
		if (p != 0) {
			p = sqrt(p);
		}
		avg += p;
	}
	if (avg != 0) {
		avg /= (fftSize / 2);
	}
	for (i = 1;i < fftSize / 2;i++) {
		p = (phaseX[i] * phaseX[i]) + (phaseY[i] * phaseY[i]);
		if (p != 0) {
			p = sqrt(p);
			nx = avg / p;
			phaseX[i] *= nx;
			phaseY[i] *= nx;
		}
	}

	memset(ovInfo[0],0,sizeof (OVERTONE_INFO));
	memset(ovInfo[1],0,sizeof (OVERTONE_INFO));
	memset(ovInfo[2],0,sizeof (OVERTONE_INFO));

	for (i = 0;i < 360;i++) {
		ovInfo[0]->do2Idx[i] = -1;
		ovInfo[1]->do2Idx[i] = -1;
		ovInfo[2]->do2Idx[i] = -1;
	}
	ovInfo[0]->nSample = fftSize / 2;
	ovInfo[0]->validSamplingRate = hfc;
	ovInfo[0]->samplingRate = outSampleR;
	ovInfo[0]->phaseX = phaseX;
	ovInfo[0]->phaseY = phaseY;

	ovInfo[1]->nSample = fftSize / 2;
	ovInfo[1]->validSamplingRate = hfc;
	ovInfo[1]->samplingRate = outSampleR;
	ovInfo[1]->phaseX = phaseX;
	ovInfo[1]->phaseY = phaseY;
	ovInfo[1]->log = 0;

	ovInfo[2]->nSample = fftSize / 2;
	ovInfo[2]->validSamplingRate = hfc;
	ovInfo[2]->samplingRate = outSampleR;
	ovInfo[2]->phaseX = phaseX;
	ovInfo[2]->phaseY = phaseY;

	pIn[0]	= &mem1[((fftSize / 2) * 0)];
	pOut[0] = &mem2[((fftSize / 2) * 0)];
	pIn[1]	= &mem1[((fftSize / 2) * 1)];
	pOut[1] = &mem3[((fftSize / 2) * 1)];
	pIn[2]	= &mem1[((fftSize / 2) * 2)];
	pOut[2] = &mem4[((fftSize / 2) * 2)];

	per = -1;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + ((fftSize / 2));startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
		}
		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		memset(mem2,0,wkMemSize);
		memset(mem3,0,wkMemSize);
		memset(mem4,0,wkMemSize);
		
		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					// 1
					genOverToneSub(hfc,pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],ovInfo[0],param);
				}
				#pragma omp section
				{
					// 2
					genOverToneSub(hfc,pIn[1],pOut[1],fftw_in[1],fftw_out[1],fftw_p[1],fftw_ip[1],ovInfo[1],param);
				}
				#pragma omp section
				{
					// 3
					genOverToneSub(hfc,pIn[2],pOut[2],fftw_in[2],fftw_out[2],fftw_p[2],fftw_ip[2],ovInfo[2],param);
				}
			}
		}

		for (i = 0;i < fftSize * 2;i++) {
			mem2[i] += mem3[i] + mem4[i];
		}

		if (startInSample + fftSize / 2 >= 0) {
			SSIZE *pIn2,*pOut2;
			// レベル調整
			pIn2  = (SSIZE *)mem1;
			pOut2 = (SSIZE *)mem2;
			pIn2  += (fftSize / 2);
			pOut2 += (fftSize / 2);
			for (i = 0;i < 12;i++) {
				for (j = 0;j < fftSize / 12;j++) {
					fftw_in2[j][0] = pIn2[j];
					fftw_in2[j][1] = 0;
				}
				pIn2 += fftSize / 12;
				// 窓関数
				for (j = 0;j < ((fftSize / 12) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 12));
				}
				for (j = ((fftSize / 12) - 1) / 2;j < (fftSize / 12);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 12));
				}

				// FFT
				fftw_execute(fftw_p2);

				// 元信号の高域のパワーを調べる
				refPw = 0;
				lowIndex = (((double)fftSize / 12) / outSampleR) * (hfc - 2000);
				highIndex = (((double)fftSize / 12) / outSampleR) * hfc;
				memset(ovInfo[0]->avg,0,sizeof (double));
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					ovInfo[0]->avg[nn] = p;
				}
				for (j = 0;j + 1 < nn;j++) {
					for (k = j;k < nn;k++) {
						if (ovInfo[0]->avg[j] > ovInfo[0]->avg[k]) {
							p = ovInfo[0]->avg[j];
							ovInfo[0]->avg[j] = ovInfo[0]->avg[k];
							ovInfo[0]->avg[k] = p;
						}
					}
				}
				for (j = 0;j + 2 < nn;j++) {
					refPw += ovInfo[0]->avg[j];
				}
				if (j > 0) {
					refPw /= j;
				}

				// 付加信号の高域のパワーを調べる
				for (j = 0;j < fftSize / 12;j++) {
					fftw_in2[j][0] = pOut2[j];
					fftw_in2[j][1] = 0;
				}
				// 窓関数
				for (j = 0;j < ((fftSize / 12) - 1) / 2;j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 * j / ((double)fftSize / 12));
				}
				for (j = ((fftSize / 12) - 1) / 2;j < (fftSize / 12);j++) {
					fftw_in2[j][0] = fftw_in2[j][0] * (2.0 - 2.0 * j / ((double)fftSize / 12));
				}
				// FFT
				fftw_execute(fftw_p2);

				ovTonePw = 0;
				lowIndex = (((double)fftSize / 12) / outSampleR) * hfc;
				highIndex = (((double)fftSize / 12) / outSampleR) * (hfc + 2000);
				memset(ovInfo[0]->avg,0,sizeof (double));
				for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
					p = fftw_in2[j][0] * fftw_in2[j][0] + fftw_in2[j][1] * fftw_in2[j][1];
					if (p != 0) {
						p = sqrt(p);
					}
					ovInfo[0]->avg[nn] = p;
				}
				for (j = 0;j + 1 < nn;j++) {
					for (k = j;k < nn;k++) {
						if (ovInfo[0]->avg[j] > ovInfo[0]->avg[k]) {
							p = ovInfo[0]->avg[j];
							ovInfo[0]->avg[j] = ovInfo[0]->avg[k];
							ovInfo[0]->avg[k] = p;
						}
					}
				}
				for (j = 0;j + 3 < nn;j++) {
					ovTonePw += ovInfo[0]->avg[j];
				}
				if (j > 0) {
					ovTonePw /= j;
				}
				if (ovTonePw > 0) {
					for (j = 0;j < fftSize / 12;j++) {
						pOut2[j] *= ((refPw / ovTonePw));
						pOut2[j] *= 0.95;
					}
				} else {
					for (j = 0;j < fftSize / 12;j++) {
						pOut2[j] = 0;
					}
				}
				pOut2 += fftSize / 12;
			}

#if 1
			// パルス除去
			pIn2  = (SSIZE *)mem2;
			for (i = 1;i + 1 < fftSize;i++) {
				p1 = pIn2[i - 1];
				p2 = pIn2[i];
				p3 = pIn2[i + 1];
				p1 >>= (56 - param->w);
				p2 >>= (56 - param->w);
				p3 >>= (56 - param->w);
				if (p1 < p2 && p2 > p3) {
					if (p2 - ((p1 + p3) / 2) >= 4) {
						pIn2[i] = ((pIn2[i - 1] + pIn2[i + 1]) / 2);
					}
				} else if (p1 > p2 && p2 < p3) {
					if (((p1 + p3) / 2) - p2 >= 4) {
						pIn2[i] = ((pIn2[i - 1] + pIn2[i + 1]) / 2);
					}
				}
			}
#endif
			outTempFile(fp_w,mem2 + fftSize / 2,fftSize,param);
			if (param->err) {
				break;
			}
		}
	}

	fio_flush(fp_w);
	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);

	al_free(phaseX);
	al_free(phaseY);

	free(ovInfo[0]);
	free(ovInfo[1]);
	free(ovInfo[2]);

	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_ip[0]);

	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_ip[1]);

	fftw_destroy_plan(fftw_p[2]);
	fftw_destroy_plan(fftw_ip[2]);

	fftw_destroy_plan(fftw_p2);

	fftw_free(fftw_in[0]);
	fftw_free(fftw_out[0]);

	fftw_free(fftw_in[1]);
	fftw_free(fftw_out[1]);

	fftw_free(fftw_in[2]);
	fftw_free(fftw_out[2]);

	fftw_free(fftw_in2);
}

//---------------------------------------------------------------------------
// Function   : genOverToneSub
// Description: 失われた高域の再現処理のサブ関数(倍音解析)
// ---
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域にデータを追加する)
//	pIn			:入力バッファ
//	pOut		:出力バッファ
//	fftw_in		:FFTW 入力
//	fftw_out	:FFTW 出力
//	fftw_p		:FFTW プラン
//	fftw_ip		:FFTW プラン
//	ovInfo		:高域生成用構造体
//	param		:変換パラメータ
//
void genOverToneSub(long hfc,SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	long outSampleR;
	long fftSize,i,j,n,nn;
	long lowIndex,highIndex;
	double nx;
	double p,p2,refPw,ovTonePw;
	double phaseTemp;
	int d,d2;
	int overToneNotFound;

	outSampleR = param->outSampleR;
	fftSize = 4096;
	if (outSampleR == 44100 * 2 || outSampleR == 48000 * 2) {
		fftSize = outSampleR / 10;
		fftSize = 4096 * 2;
	}
	if (outSampleR == 32000 * 6 || outSampleR == 44100 * 4 || outSampleR == 48000 * 4) {
		fftSize = outSampleR / 10;
		fftSize = 8192 * 2;
	}
	if (outSampleR == 32000 * 12 || outSampleR == 44100 * 8 || outSampleR == 48000 * 8) {
		fftSize = outSampleR / 10;
		fftSize = 16384 * 2;
	}
	if (outSampleR == 32000 * 24 || outSampleR == 44100 * 16 || outSampleR == 48000 * 16) {
		fftSize = 32768 * 2;
	}
	if (outSampleR == 32000 * 48 || outSampleR == 44100 * 32 || outSampleR == 48000 * 32) {
		fftSize = 65536 * 2;
	}

//fftSize = 4096;
	// FFT 初期設定
	copyToFFTW(fftw_in,pIn,fftSize);

	// 窓関数
	windowFFTW(fftw_in,fftSize);

	// FFT
	fftw_execute(fftw_p);

	// 元信号の高域のパワーを調べる
	memset(ovInfo->pw_cnt,0,65536*2 * sizeof (int));
	refPw = 0;
	lowIndex = ((double)fftSize / outSampleR) * (hfc - 2000);
	highIndex = ((double)fftSize / outSampleR) * hfc;
#if 0
	for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
		p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
		if (p != 0) {
			p = sqrt(p);
		}
		refPw += p;
	}
#endif
	for (i = lowIndex;i < highIndex;i++) {
		p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
		if (p != 0) {
			p = sqrt(p);
			p = log10(p) * 20;
			p /= 3;
		}
		for (j = lowIndex;j < highIndex;j++) {
			p2 = fftw_out[j][0] * fftw_out[j][0] + fftw_out[j][1] * fftw_out[j][1];
			if (p2 != 0) {
				p2 = sqrt(p2);
				p2 = log10(p2) * 20;
				p2 /= 3;
			}
			if (p > 0 || p2 > 0) {
				if ((long)p == (long)p2) {
					ovInfo->pw_cnt[i]++;
				}
			}
		}
	}
	for (i = lowIndex + 1,nn = lowIndex;i < highIndex;i++) {
		if (ovInfo->pw_cnt[i] > ovInfo->pw_cnt[nn]) {
			nn = i;
		}
	}
	refPw = fftw_out[nn][0] * fftw_out[nn][0] + fftw_out[nn][1] * fftw_out[nn][1];
	if (refPw > 0) {
		refPw = sqrt(refPw);
	}

	//
	// 倍音解析
	memset(ovInfo->power,0,65536*2 * sizeof (double));
	memset(ovInfo->phase,0,65536*2 * sizeof (double));
	memset(ovInfo->pw,0,65536*2 * sizeof (double));
	memset(ovInfo->avg,0,65536*2 * sizeof (double));
	memset(ovInfo->diff,0,65536*2 * sizeof (double));
	memset(ovInfo->base,0,65536*2 * sizeof (double));
	memset(ovInfo->baseToAdj,0,65536*2 * sizeof (double));

	for (i = 1;i < fftSize / 2;i++) {
		p = (fftw_out[i][0] * fftw_out[i][0]) + (fftw_out[i][1] * fftw_out[i][1]);
		if (p != 0) {
			p = sqrt(p);
		}
		if (fftw_out[i][0] != 0 && fftw_out[i][1] != 0) {
			phaseTemp = atan2(fftw_out[i][1],fftw_out[i][0]) * 180 / M_PI;
			phaseTemp += 180;
		} else {
			phaseTemp = 0;
		}
		ovInfo->phase[i] = phaseTemp;
		ovInfo->power[i] = p;
	}
#if 0
if (param->r1_flag && ovInfo->log) {
	FILE *fp;
	fp = fopen("d:\\fft_p.csv","a");
	if (fp) {
		fprintf(fp,"\n");
		for (i = 1;i < fftSize / 2;i++) {
			fprintf(fp,",%lf",ovInfo->power[i]);
		}
		fclose(fp);
	}
	fp = fopen("d:\\fft_h.csv","a");
	if (fp) {
		fprintf(fp,"\n");
		for (i = 1;i < fftSize / 2;i++) {
			fprintf(fp,",%lf",ovInfo->phase[i]);
		}
		fclose(fp);
	}
}
#endif

#if 0
if (ovInfo->log) {
	FILE *ofp;
	ofp = fopen("d:\\fft.csv","a");
	if (ofp) {
		for (i = 300;i < 750;i++) {
			fprintf(ofp,"%lf,",ovInfo->power[i]);
		}
		fprintf(ofp,"\n");
		fclose(ofp);
	}
}
#endif
	if (param->hfa == 2) {
		anaOverToneHFA2(ovInfo,param);
	} else {
		anaOverToneHFA3(ovInfo,param);
	}

	//
	// 信号生成
	for (i = 1;i < fftSize / 2;i++) {
		if (ovInfo->pw[i] > 0) {
			d = ovInfo->phase[i];
			if (d < 0 || d >= 360) {
				d = 0;
			}
			if (ovInfo->do2Idx[d] == -1) {
				for (j = 1;j < fftSize / 2;j++) {
					if (ovInfo->phaseX[j] != 0 && ovInfo->phaseY[j] != 0) {
						phaseTemp = atan2(ovInfo->phaseY[j],ovInfo->phaseX[j]) * 180 / M_PI;
						phaseTemp += 180;
					} else {
						phaseTemp = 0;
					}
					d2 = phaseTemp;
					if (d == d2) {
						break;
					}
				}
			} else {
				j = ovInfo->do2Idx[d];
			}
			if (j < fftSize / 2) {
				// 位相が一致
				fftw_out[i][0] = ovInfo->phaseX[j];
				fftw_out[i][1] = ovInfo->phaseY[j];
				ovInfo->do2Idx[d] = j;
			} else {
				ovInfo->pw[i] = 0;
			}
		}
	}
	overToneNotFound = 1;
	for (i = 1;i < fftSize / 2;i++) {
		if (ovInfo->pw[i] > 0) {
			fftw_out[i][0] *= ovInfo->pw[i];
			fftw_out[i][1] *= ovInfo->pw[i];
			overToneNotFound = 0;
		} else {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}

	memset(ovInfo->pw_cnt,0,65536*2 * sizeof (int));
	lowIndex = ((double)fftSize / outSampleR) * (hfc);
	highIndex = ((double)fftSize / outSampleR) * (hfc + 2000);
	for (i = lowIndex;i < highIndex;i++) {
		p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
		if (p != 0) {
			p = sqrt(p);
			p = log10(p) * 20;
			p /= 3;
		}

		for (j = lowIndex;j < highIndex;j++) {
			p2 = fftw_out[j][0] * fftw_out[j][0] + fftw_out[j][1] * fftw_out[j][1];
			if (p2 != 0) {
				p2 = sqrt(p2);
				p2 = log10(p2) * 20;
				p2 /= 3;
			}
			if (p > 0 || p2 > 0) {
				if ((long)p == (long)p2) {
					ovInfo->pw_cnt[i]++;
				}
			}
		}
	}
	for (i = lowIndex + 1,nn = lowIndex;i < highIndex;i++) {
		if (ovInfo->pw_cnt[i] > ovInfo->pw_cnt[nn]) {
			nn = i;
		}
	}
	ovTonePw = fftw_out[nn][0] * fftw_out[nn][0] + fftw_out[nn][1] * fftw_out[nn][1];
	if (ovTonePw > 0) {
		ovTonePw = sqrt(ovTonePw);
	} else {
		ovTonePw = 1;
	}
	if (overToneNotFound == 0) {
		// 付加する信号のパワーを調べる
#if 0
		ovTonePw = 0;
		for (i = lowIndex,nn = 0;i < highIndex;i++,nn++) {
			p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
			if (p != 0) {
				p = sqrt(p);
			}
			ovTonePw += p;
		}
		if (nn > 0) {
			ovTonePw /= nn;
		}
		if (ovTonePw == 0) {
			ovTonePw = 1;
		}
#endif
		nx = refPw / ovTonePw;
		nx *= 0.9;

		for (i = 1;i < fftSize / 2;i++) {
			fftw_out[i][0] *= nx;
			fftw_out[i][1] *= nx;
		}
		adjPinkFilter(2,fftSize,fftw_out,param);

		lowIndex = ((double)fftSize / outSampleR) * (hfc - 2000);
		highIndex = ((double)fftSize / outSampleR) * hfc;

		// 低域カット
		for (i = 1;i <= highIndex;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
		// 半分のデータを復元
		for (i = 1;i < fftSize / 2;i++) {
			fftw_out[fftSize - i][0] = fftw_out[i][0];
			fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
		}

		// 直流成分を削除
		fftw_out[0][0] = 0;
		fftw_out[0][1] = 0;

		// invert FFT
		fftw_execute(fftw_ip);

		// 窓関数
		for (i = 0;i < (fftSize - 1) / 2;i++) {
			fftw_in[i][0] = fftw_in[i][0] * (2.0 * i / (double)fftSize);
		}
		for (i = (fftSize - 1) / 2;i < fftSize;i++) {
			fftw_in[i][0] = fftw_in[i][0] * (2.0 - 2.0 * i / (double)fftSize);
		}
	} else {
		// 無音生成
		for (i = 0;i < fftSize;i++) {
			fftw_in[i][0] = 0;
		}
	}
	// 出力
	for (i = 0;i < fftSize;i++) {
		pOut[i] += (fftw_in[i][0] / fftSize);
	}
}

//---------------------------------------------------------------------------
// Function   : anaOverToneHFA2
// Description: 倍音解析
//				window幅ごとにpowerの平均値をとりpowerが大きいものを抽出する
// ---
//	ovInfo :倍音生成用構造体
//	param  :パラメータ
//
void anaOverToneHFA2(OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	DWORD ofs,window,width,swidth;
	DWORD n,i;
	double avg,maxAvg;
	double avgLine;
	int maxOfs,maxWin;
	
	
	int lowIndex,highIndex;
	int pha;
	
	
	long nSample;
	long lowHz,wid;
	long skipCnt;
	nSample = ovInfo->nSample;
	//
	// 計算対象のインデックスを求める
	for (i = 1;i < ovInfo->nSample;i++) {
		ovInfo->pw[i] = 0;
	}

	lowHz = 9000;
	wid   = 2500;
	if (lowHz + wid > ovInfo->validSamplingRate) {
		lowHz = ovInfo->validSamplingRate - wid;
	}
	if (lowHz < 4000) {
		// 高域の情報がないので解析しない
		return;
	}

	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 200;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
	highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	if (swidth == 0) {
		swidth = 1;
	}

	// パワーの平均を計算する
	//
	avg = 0;
	for (i = lowIndex,n = 0;i < highIndex;i++) {
		avg += ovInfo->power[i];
		n++;
	}
	if (n > 0) {
		avg /= n;
	}
	avgLine = avg;
	ofs = lowIndex;
	//
	// powerが強いものを採用する(sig1)
	if (param->sig1Enb == 1) {
		window = swidth;
		maxOfs = ofs;
		maxWin = swidth;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			maxAvg = -1;
			skipCnt = 0;
			for (window = swidth;window < width;window++) {
				skipCnt++;
				if (param->hfaFast && (skipCnt % 4) != 0) {
					continue;
				}
				avg = 0;
				for (i = ofs,n = 0;i < highIndex;i += window) {
					avg += ovInfo->power[i];
					n++;
				}
				if (n > 0) {
					avg /= n;
				}

				if (maxAvg == -1 || maxAvg < avg) {
					maxAvg = avg;
					maxOfs = ofs;
					maxWin = window;
				}

			}
			if (avgLine * param->sig1AvgLineNx < maxAvg) {
				pha = ovInfo->phase[maxOfs];
				for (i = maxOfs,n = 1;i < nSample;i += maxWin,n++) {
					if (ovInfo->pw[i] < maxAvg / n) {
						ovInfo->pw[i] = maxAvg / n;
						ovInfo->phase[i] = pha;

						pha += param->sig1Phase;
						if (pha >= 360) {
							pha -= 360;
						}
						if (pha < 0) {
							pha += 360;
						}
					}
				}
			}
		}
	}

	//
	// window間隔の前後のpowerの差の累積が少ないものを採用する(sig2)
	if (param->sig2Enb == 1) {
		for (ofs = lowIndex;ofs < lowIndex + width;ofs++) {
			maxAvg = -1;
			skipCnt = 0;
			for (window = swidth;window < width;window++) {
				skipCnt++;
				if (param->hfaFast && (skipCnt % 4) != 0) {
					continue;
				}
				avg = 0;
				for (i = ofs,n = 0;i < highIndex;i += window) {
					if (1 + window < i) {
						if (ovInfo->power[i - window] <= ovInfo->power[i]) {
							avg += ovInfo->power[i] - ovInfo->power[i - window];
						} else {
							avg += ovInfo->power[i] - ovInfo->power[i - window];
						}
						n++;
					}
				}
				if (n > 0) {
					avg /= n;
				}

				if (maxAvg == -1 || maxAvg > avg) {
					maxAvg = avg;
					maxOfs = ofs;
					maxWin = window;
				}

			}
			avg = 0;
			for (i = maxOfs,n = 0;i < highIndex;i += maxWin) {
				avg += ovInfo->power[i];
				n++;
			}
			if (n > 0) {
				avg /= n;
			}
			maxAvg = avg;
			pha = ovInfo->phase[maxOfs];
			for (i = maxOfs,n = 1;i < nSample;i += maxWin,n++) {
				if (ovInfo->pw[i] > maxAvg) {
					ovInfo->pw[i] = maxAvg;
					ovInfo->phase[i] = pha;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha -= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
				} else if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = maxAvg / 10;
					ovInfo->phase[i] = pha;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha -= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
				}
			}
		}
	}
}
#if 1	// 0.7 の HFA3
void anaOverToneHFA3(OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	DWORD ofs,window,width,swidth;
	DWORD validIndex;
	DWORD n,i,j,k;
	double avg;
	double avgRef;
	double avgLine;
	double diff,diff0,diff1,diff2,diff3,diff4,diff5,diffP;
	double refPw[8];
	double avgPw,avgPw2,avgPw3;
	double tmpAvgPw,tmpAvgPw2,tmpAvgPw3;
	int    avgPwNX,avgPwNX2,avgPwNX3;
	long skipCnt;
	double tbl_hfaDiffMin[5] = {0.84,1.04,1.24,1.74,2.14};

	// 予測との最小パラメータ保存
	int minWin;
	int minType;
	int max_i;
	double minDiff;
	int nn;
	int odd;
	double hz;
	DWORD  baseOfs;
	double tmpPw,tmpPw2;
	int lowIndex,highIndex;
	int lowRange,highRange;
	int minWidth;
	int pha;
	double phaRand;
	int phaTmp = 0;
	long nSample;
	long lowHz,wid;
	double areaAvg;
	nSample = ovInfo->nSample;

	//
	// 初期化
	for (i = 1;i < ovInfo->nSample;i++) {
		ovInfo->pw[i] = 0;
		ovInfo->diff[i] = -1;
		ovInfo->baseToAdj[i] = 0;
	}
	swidth = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1600;
#if 1
	for (i = 1;i < ovInfo->nSample;i+= swidth) {
		avg = 0;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			avg += ovInfo->power[j];
		}
		avg /= n;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			ovInfo->baseToAdj[j] = ovInfo->power[j] - avg;
			ovInfo->power[j] -= ovInfo->baseToAdj[j];
		}
	}
#endif
	if (ovInfo->validSamplingRate < 8000) {
		// 高域の情報がないので解析しない
		return;
	}

	lowHz	= 8000;
	wid		= 3000;
	if (lowHz + wid  >= ovInfo->validSamplingRate) {
		wid = ovInfo->validSamplingRate - lowHz;
	}

	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
	minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
	if (ovInfo->validSamplingRate < 16000) {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	} else {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
	}
	if (swidth == 0) {
		swidth = 1;
	}

	avg = 0;
	for (i = lowIndex,n = 0;i < highIndex;i++) {
		avg += ovInfo->power[i];
		n++;
	}
	if (n > 0) {
		avg /= n;
	}
	avgLine = avg;

	if (param->sig2Enb == 1) {
		// 前後のwindowで振幅の差が少ない音声の補間
		window = width;
		minWin = window;
		minType = 0;

		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = diff1 = diff2 = diff3 = diff4 = diff5 = diffP = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				odd = 1;	// 奇数倍のみ倍音がある
				refPw[0] = -1;
				refPw[1] = -1;
				refPw[2] = -1;
				refPw[3] = -1;
				refPw[4] = -1;
				refPw[5] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						refPw[1] = ovInfo->power[i] * n;
						refPw[2] = ovInfo->power[i] * odd;
						refPw[3] = ovInfo->power[i] * (odd * odd);
						refPw[4] = ovInfo->power[i];
						refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					}
					// 前後のパワーの差の計算
					if (i - window >= ofs) {
						if (ovInfo->power[i - window] >= ovInfo->power[i]) {
							diff = ovInfo->power[i - window] - ovInfo->power[i];
						} else {
							diff = ovInfo->power[i] - ovInfo->power[i - window];
						}
					}
					diffP += (diff * tbl_hfaDiffMin[param->hfaDiffMin - 1]);

					avgPw += ovInfo->power[i];
					avgPwNX++;
					if ((avgPwNX & 0x01) == 0) {
						avgPw2 += ovInfo->power[i];
						avgPwNX2++;
					}
					if ((avgPwNX % 3) == 0) {
						avgPw3 += ovInfo->power[i];
						avgPwNX3++;
					}
					
					// 1/f(振幅が1/fになっているもの)
					diff = refPw[0] / hz;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff0 += diff;

					// 鋸波(nの逆数で小さくなる)
					diff = refPw[1] / n;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff1 += diff;

					// 短形波(奇数倍音,nの逆数で小さくなる)
					diff = refPw[2] / odd;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff2 += diff;

					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					diff = refPw[3] / (odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff3 += diff;

					// パルス(n番目の倍音でもパワーは同じ)
					diff = refPw[4];
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff4 += diff;

					// その他(もっとパワーが小さくなるパターン)
					diff = refPw[5] / (odd * odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff5 += diff;

					nn++;
				}

				diff0 += diffP * 2;
				diff1 += diffP * 2;
				diff2 += diffP * 2;
				diff3 += diffP * 2;
				diff4 += diffP * 2;
				diff5 += diffP * 2;

				if (nn > 0) {
					diff0 /= nn;
					diff1 /= nn;
					diff2 /= nn;
					diff3 /= nn;
					diff4 /= nn;
					diff5 /= nn;
					diffP /= nn;
//					if (refPw[4] > 0) {
//						diff0 /= refPw[4];
//						diff1 /= refPw[4];
//						diff2 /= refPw[4];
//						diff3 /= refPw[4];
//						diff4 /= refPw[4];
//						diff5 /= refPw[4];
//						diffP /= refPw[4];
//					}
				}

//				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
//					tmpAvgPw  = avgPw / avgPwNX;
//					tmpAvgPw2 = avgPw2 / avgPwNX2;
//					tmpAvgPw3 = avgPw3 / avgPwNX3;
//					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {//						continue;
//					}
//				}

				if (minDiff == -1 || minDiff > diffP) {
					minDiff = diffP;
					minWin = window;
					minType = 0;
				}
				if (minDiff > diff1) {
					minDiff = diff1;
					minWin = window;
					minType = 1;
				}
				if (minDiff > diff2) {
					minDiff = diff2;
					minWin = window;
					minType = 2;
				}
				if (minDiff > diff3) {
					minDiff = diff3;
					minWin = window;
					minType = 3;
				}
				if (minDiff > diff4) {
					minDiff = diff4;
					minWin = window;
					minType = 4;
				}
				if (minDiff > diff5) {
					minDiff = diff5;
					minWin = window;
					minType = 5;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			refPw[4] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}

			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					refPw[1] = ovInfo->power[i] * n;
					refPw[2] = ovInfo->power[i] * odd;
					refPw[3] = ovInfo->power[i] * (odd * odd);
					refPw[4] = ovInfo->power[i];
					refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (minType == 0) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[0] / hz;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (tmpPw2 > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 1) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[1] / n;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (tmpPw2 > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 2) {
					// 短形波(奇数倍音,nの逆数で小さくなる)
					tmpPw = refPw[2] / odd;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
						if (tmpPw2 > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 3) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
						if (tmpPw2 > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 4) {
					// パルス(n番目の倍音でもパワーは同じ)
					tmpPw = refPw[4];
					phaRand = rand() * 6;
					phaRand -= 3;
					pha += phaRand;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.25;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.25;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (tmpPw2 > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.25;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 5) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					tmpPw2 = ovInfo->pw[i];
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (tmpPw2 > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.39;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				}
			}
		}
	}

#if 1
	if (param->sig1Enb == 1) {
		// powerが強いもの優先で補間する
		lowHz	= 9500;
		wid		= 2500;
		if (lowHz + wid  >= ovInfo->validSamplingRate) {
			wid = ovInfo->validSamplingRate - lowHz;
		}

		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
		width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
		minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
		if (ovInfo->validSamplingRate < 16000) {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		} else {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
		}
		if (swidth == 0) {
			swidth = 1;
		}

		areaAvg = 0;
		for (i = lowIndex,n = 0;n < width;i++,n++) {
			areaAvg += ovInfo->power[i];
		}
		if (n > 0) {
			areaAvg /= n;
		}

		window = width;
		minWin = window;
		minType = 0;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			if (ovInfo->power[ofs] < areaAvg) {
				continue;
			}
			
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}


					if (i - window >= ofs) {
						if (ovInfo->power[i - window] < ovInfo->power[i]) {
							avgPw /= 75;
						}
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					if (ovInfo->power[i] > areaAvg) {
						avgPw += ovInfo->power[i];
					} else {
						avgPw /= 75;
					}
					avgPwNX++;
					nn++;
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}
				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff < avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが強いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig1Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = tmpPw * 0.7;
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}
#endif

	if (param->sig3Enb == 1) {
//	if (0) {
		// powerが弱いものを補間する
		lowHz	= 8000;
		wid		= 4000;
		if (lowHz + wid  >= ovInfo->validSamplingRate) {
			wid = ovInfo->validSamplingRate - lowHz;
		}

		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
		width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
		minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
		if (ovInfo->validSamplingRate < 16000) {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		} else {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
		}
		if (swidth == 0) {
			swidth = 1;
		}

		window = width;
		minWin = swidth;
		
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					avgPw += ovInfo->power[i];
					avgPwNX++;
					nn++;
				}

				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff > avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが弱いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig3Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = tmpPw * 0.14;
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}

	if (param->sig2Enb == 1 && param->hfaWide) {
		lowHz	= 7000;
		wid		= 4000;
		if (lowHz + wid  >= ovInfo->validSamplingRate) {
			wid = ovInfo->validSamplingRate - lowHz;
		}

		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
		width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
		minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
		if (ovInfo->validSamplingRate < 16000) {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		} else {
			highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
		}
		if (swidth == 0) {
			swidth = 1;
		}

		// 前後のwindowで振幅の差が少ない音声の補間
		window = width;
		minWin = window;
		minType = 0;

		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = diff1 = diff2 = diff3 = diff4 = diff5 = diffP = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				odd = 1;	// 奇数倍のみ倍音がある
				refPw[0] = -1;
				refPw[1] = -1;
				refPw[2] = -1;
				refPw[3] = -1;
				refPw[4] = -1;
				refPw[5] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						refPw[1] = ovInfo->power[i] * n;
						refPw[2] = ovInfo->power[i] * odd;
						refPw[3] = ovInfo->power[i] * (odd * odd);
						refPw[4] = ovInfo->power[i];
						refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					}
					// 前後のパワーの差の計算
					if (i - window >= ofs) {
						if (ovInfo->power[i - window] >= ovInfo->power[i]) {
							diff = ovInfo->power[i - window] - ovInfo->power[i];
						} else {
							diff = ovInfo->power[i] - ovInfo->power[i - window];
						}
					}
					diffP += (diff * 1.48);

					avgPw += ovInfo->power[i];
					avgPwNX++;
					if ((avgPwNX & 0x01) == 0) {
						avgPw2 += ovInfo->power[i];
						avgPwNX2++;
					}
					if ((avgPwNX % 3) == 0) {
						avgPw3 += ovInfo->power[i];
						avgPwNX3++;
					}
					
					// 1/f(振幅が1/fになっているもの)
					diff = refPw[0] / hz;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff0 += diff;

					// 鋸波(nの逆数で小さくなる)
					diff = refPw[1] / n;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff1 += diff;

					// 短形波(奇数倍音,nの逆数で小さくなる)
					diff = refPw[2] / odd;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff2 += diff;

					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					diff = refPw[3] / (odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff3 += diff;

					// パルス(n番目の倍音でもパワーは同じ)
					diff = refPw[4];
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff4 += diff;

					// その他(もっとパワーが小さくなるパターン)
					diff = refPw[5] / (odd * odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff5 += diff;

					nn++;
				}

				diff0 += diffP * 4;
				diff1 += diffP * 4;
				diff2 += diffP * 4;
				diff3 += diffP * 4;
				diff4 += diffP * 4;
				diff5 += diffP * 4;

				if (nn > 0) {
					diff0 /= nn;
					diff1 /= nn;
					diff2 /= nn;
					diff3 /= nn;
					diff4 /= nn;
					diff5 /= nn;
					diffP /= nn;
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}

				if (minDiff == -1 || minDiff > diffP) {
					minDiff = diffP;
					minWin = window;
					minType = 0;
				}
				if (minDiff > diff1) {
					minDiff = diff1;
					minWin = window;
					minType = 1;
				}
				if (minDiff > diff2) {
					minDiff = diff2;
					minWin = window;
					minType = 2;
				}
				if (minDiff > diff3) {
					minDiff = diff3;
					minWin = window;
					minType = 3;
				}
				if (minDiff > diff4) {
					minDiff = diff4;
					minWin = window;
					minType = 4;
				}
				if (minDiff > diff5) {
					minDiff = diff5;
					minWin = window;
					minType = 5;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			refPw[4] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}

			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					refPw[1] = ovInfo->power[i] * n;
					refPw[2] = ovInfo->power[i] * odd;
					refPw[3] = ovInfo->power[i] * (odd * odd);
					refPw[4] = ovInfo->power[i];
					refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (minType == 0) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[0] / hz;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.87 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 1) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[1] / n;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.87 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 2) {
					// 短形波(奇数倍音,nの逆数で小さくなる)
					tmpPw = refPw[2] / odd;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.87 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 3) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					phaTmp = pha;
					if (n & 0x01) {
						phaTmp = pha + 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.87 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 4) {
					// パルス(n番目の倍音でもパワーは同じ)
					tmpPw = refPw[4];
					phaRand = rand() * 6;
					phaRand -= 3;
					pha += phaRand;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.87 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				} else if (minType == 5) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] != 0 && ovInfo->pw[i] * 0.87 > tmpPw) {
						ovInfo->pw[i] = tmpPw * 0.32;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					}
				}
			}
		}
	}

//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

#if 1
	//
	// 位相の調整
	for (i = 1;i < validIndex;i++) {
		ovInfo->diff[i] = -1;
	}
	window = width;
	for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
		minDiff = -1;
		skipCnt = 0;
		for (window = swidth;window < width;window++) {
			if (window < minWidth) {
				if ((ofs - lowIndex) > window * 1) {
					continue;
				}
			} else {
				if ((ofs - lowIndex) > window) {
					continue;
				}
			}
			skipCnt++;
			if (param->hfaFast && (skipCnt % 8) != 0) {
				continue;
			}
			// 位相を調べる
			diffP = 0;
			baseOfs = ofs - ((ofs / window) * window);
			if (baseOfs == 0) {
				baseOfs = window;
			}
			n = 1;
			refPw[0] = -1;
			for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}

				if (refPw[0] == -1) {
					refPw[0] = ovInfo->phase[i];
				}
				// 前後の位相の差の計算
				if (i - window >= ofs) {
					if (refPw[0] <= ovInfo->phase[i]) {
						diff = ovInfo->phase[i] - refPw[0];
					} else {
						diff = refPw[0] - ovInfo->phase[i];
					}
				}
				diffP += diff;
				nn++;
			}

			if (nn > 0) {
				diffP /= nn;
			}

			if (minDiff == -1 || minDiff > diffP) {
				minDiff = diffP;
				minWin = window;
			}
		}

		// 一番予測誤差が少なかったものを採用する。

		baseOfs = ofs - ((ofs / minWin) * minWin);
		if (baseOfs == 0) {
			baseOfs = minWin;
		}

		pha = ovInfo->phase[baseOfs];
		n = 1;		// 奇数偶数倍すべてに倍音がある

		refPw[0] = -1;

		for (i = baseOfs;i < validIndex;i += minWin,n++) {
			hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
			if (hz < lowHz) {
				continue;
			}
			if (refPw[0] == -1) {
				refPw[0] = ovInfo->phase[i];
			}

			phaRand = rand() * 2;
			phaRand -= 1;
			pha = refPw[0];
			//pha += phaRand;
			if (pha >= 360) {
				pha %= 360;
			}
			if (pha < 0) {
				pha += 360;
			}
			if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
				ovInfo->phase[i] = pha;
				ovInfo->diff[i] = minDiff;
			}
		}
	}
#endif
	// 補間されていない箇所の対応
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	if (param->sig3Phase >= -6 || param->sig3Phase <= 6) {
		// 位相の修正
		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->phase[i];
					}
					
					diff = refPw[0];
					if (diff >= ovInfo->phase[i]) {
						diff = diff - ovInfo->phase[i];
					} else {
						diff = ovInfo->phase[i] - diff;
					}
					diff0 += diff;
					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
				}

				if (minDiff == -1 || minDiff > diff0) {
					minDiff = diff0;
					minWin = window;
					minType = 0;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
					pha += 360;
					pha += 360;
					pha %= 360;
				}
				if (ovInfo->pw[i] != 0) {
					if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
						ovInfo->phase[i] = pha;
					}
				}
			}
		}
	}
#if 1
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 5000);
	highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 2000);
	i = 1;
	do {
		for (j = i,n = 0,k = lowIndex;n < highIndex - lowIndex && j < ovInfo->nSample;j++,n++) {
			ovInfo->pw[j] += ovInfo->baseToAdj[k];
			k++;
			if (k > highIndex) {
				k = lowIndex;
			}
		}
		i += n;
	} while (i < ovInfo->nSample);
#endif
#if 0
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	highIndex = nSample;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 5000;
	i = lowIndex;
	do {
		avg = 0;
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			avg += ovInfo->pw[j];
		}
		if (n > 0) {
			avg /= n;
			avg = log10(avg) * 20;
		}
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			tmpPw = ovInfo->pw[j];
			if (tmpPw > 0) {
				tmpPw = log10(tmpPw) * 20;
			}
			if (tmpPw + 15 < avg) {
				ovInfo->pw[j] *= 13;
			} else if (tmpPw > avg + 15) {
				ovInfo->pw[j] /= 13;
			}
		}
		i = j;
	} while (i < ovInfo->nSample);
#endif
//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			for (i = 0;i < 4096;i++) {
//				fprintf(ofp,"%lf,%lf\n",ovInfo->pw[i],ovInfo->phase[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//		ovInfo->log = 0;
//	}

#if 1
	if (param->hfaWide) {
		// 解析情報を平均化し、急激な変化を抑える
		for (i = baseOfs;i + 1< nSample;i++) {
			hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
			if (hz < lowHz) {
				continue;
			}
			if (ovInfo->pw[i] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
#endif
	if (1) {
		double nx;
		// 付加した広域で強すぎるものがあれば弱める。
		lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
		highIndex = nSample;
		swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 5000;
		i = lowIndex;
		nx = 7.7;
		if (ovInfo->validSamplingRate < 12000) {
			nx = 4.2;
		}
		do {
			avg = 0;
			for (j = i,n = 0;n < swidth && j < ovInfo->nSample;j++,n++) {
				avg += ovInfo->pw[j];
			}
			if (n > 0) {
				avg /= n;
			}
			for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
				if (avg * nx < ovInfo->pw[j]) {
					ovInfo->pw[j] *= avg * nx / ovInfo->pw[j];
				}
			}
			i = j;
		} while (i < ovInfo->nSample);
	}
if (0) {
	if (ovInfo->pw[i] > 0 && ovInfo->pw[i + 1] > 0) {
		ovInfo->pw[i] = (ovInfo->pw[i] + ovInfo->pw[i+1]) / 2;
		ovInfo->phase[i] = (ovInfo->phase[i] + ovInfo->phase[i+1]) / 2;
		if (ovInfo->phase[i] >= 360) {
			pha = ovInfo->phase[i];
			ovInfo->phase[i] = (pha % 360);
		}
		if (ovInfo->phase[i] < 0) {
			ovInfo->phase[i] += 360;
		}
	}
}
}
#else
void anaOverToneHFA3(OVERTONE_INFO *ovInfo,PARAM_INFO *param)
{
	DWORD ofs,window,width,swidth;
	DWORD n,i,j,k,lll;
	int ii,jj;
	DWORD validN;
	double avg,minAvg,maxAvg,peekPw,diffPw,diffAvg;
	double keyPw;
	double avgLine;
	double diff,diff0,diff1,diff2,diff3,diff4,diff5,diff6,diff7,diffP;
	double refPw[8];
	double avgPw,avgPw2,avgPw3;
	double tmpAvgPw,tmpAvgPw2,tmpAvgPw3;
	int    avgPwNX,avgPwNX2,avgPwNX3;
	double maxPwAvg;
	long skipCnt;

	// 予測との最小パラメータ保存
	int minWin;
	int minType;
	int max_i;
	double minDiff;

	int nd,nn;
	int odd;
	double ndLv;
	double hz,hz2;
	DWORD  baseOfs;
	int minOfs,maxOfs,maxWin,maxType;
	double tmpPw;
	double nx,nx2;
	int tmpOfs,tmpWin;
	int lowIndex,highIndex;
	int minWidth;
	int miOffset;
	int pha;
	double phaRand;
	int phaTmp = 0;
	long nSample;
	long lowHz,wid;
	double areaAvg;
	nSample = ovInfo->nSample;

	//
	// 初期化
	for (i = 1;i < ovInfo->nSample;i++) {
		ovInfo->pw[i] = 0;
		ovInfo->diff[i] = -1;
		ovInfo->baseToAdj[i] = 0;
	}
	swidth = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 2000;
	//swidth = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 400;
#if 0
	for (i = 1;i < ovInfo->nSample;i+= swidth) {
		avg = 0;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			avg += ovInfo->power[j];
		}
		avg /= n;
		for (j = i,n = 0;n < swidth && j < ovInfo->nSample;n++) {
			ovInfo->baseToAdj[j] = ovInfo->power[j] - avg;
			ovInfo->power[j] -= ovInfo->baseToAdj[j];
		}
	}
#endif
	lowHz	= 6000;
	wid		= 3500;
	if (param->hfaWide) {
		lowHz	= 4500;
		wid		= 5500;
	}
	if (lowHz + wid + 2000 >= ovInfo->validSamplingRate) {
		wid = 2000;
		lowHz = 4500;
	}

	if (lowHz < 4000) {
		// 高域の情報がないので解析しない
		return;
	}

	swidth	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 180;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * wid;
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * lowHz;
	minWidth  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 1000;
	if (ovInfo->validSamplingRate < 16000) {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	} else {
		highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 16000;
	}
	if (swidth == 0) {
		swidth = 1;
	}

	avg = 0;
	for (i = lowIndex,n = 0;i < highIndex;i++) {
		avg += ovInfo->power[i];
		n++;
	}
	if (n > 0) {
		avg /= n;
	}
	avgLine = avg;

	if (param->sig2Enb == 1) {
		// 前後のwindowで振幅の差が少ない音声の補間
		window = width;
		minWin = window;
		minType = 0;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = diff1 = diff2 = diff3 = diff4 = diff5 = diffP = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				odd = 1;	// 奇数倍のみ倍音がある
				refPw[0] = -1;
				refPw[1] = -1;
				refPw[2] = -1;
				refPw[3] = -1;
				refPw[4] = -1;
				refPw[5] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						refPw[1] = ovInfo->power[i] * n;
						refPw[2] = ovInfo->power[i] * odd;
						refPw[3] = ovInfo->power[i] * (odd * odd);
						refPw[4] = ovInfo->power[i];
						refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					}
					// 前後のパワーの差の計算
					if (i - window >= ofs) {
						if (ovInfo->power[i - window] >= ovInfo->power[i]) {
							diff = ovInfo->power[i - window] - ovInfo->power[i];
						} else {
							diff = ovInfo->power[i] - ovInfo->power[i - window];
						}
					}
					diffP += diff;

					avgPw += ovInfo->power[i];
					avgPwNX++;
					if ((avgPwNX & 0x01) == 0) {
						avgPw2 += ovInfo->power[i];
						avgPwNX2++;
					}
					if ((avgPwNX % 3) == 0) {
						avgPw3 += ovInfo->power[i];
						avgPwNX3++;
					}
					
					// 1/f(振幅が1/fになっているもの)
					diff = refPw[0] / hz;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff0 += diff;

					// 鋸波(nの逆数で小さくなる)
					diff = refPw[1] / n;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff1 += diff;

					// 短形波(奇数倍音,nの逆数で小さくなる)
					diff = refPw[2] / odd;
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff2 += diff;

					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					diff = refPw[3] / (odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff3 += diff;

					// パルス(n番目の倍音でもパワーは同じ)
					diff = refPw[4];
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff4 += diff;

					// その他(もっとパワーが小さくなるパターン)
					diff = refPw[5] / (odd * odd * odd);
					if (diff >= ovInfo->power[i]) {
						diff = diff - ovInfo->power[i];
					} else {
						diff = ovInfo->power[i] - diff;
					}
					diff5 += diff;

					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
					diff1 /= nn;
					diff2 /= nn;
					diff3 /= nn;
					diff4 /= nn;
					diff5 /= nn;
					diffP /= nn;
				}

				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
					tmpAvgPw  = avgPw / avgPwNX;
					tmpAvgPw2 = avgPw2 / avgPwNX2;
					tmpAvgPw3 = avgPw3 / avgPwNX3;
					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {						continue;
					}
				}

				if (minDiff == -1 || minDiff > diffP) {
					minDiff = diffP;
					minWin = window;
					minType = 0;
				}
				if (minDiff > diff1) {
					minDiff = diff1;
					minWin = window;
					minType = 1;
				}
				if (minDiff > diff2) {
					minDiff = diff2;
					minWin = window;
					minType = 2;
				}
				if (minDiff > diff3) {
					minDiff = diff3;
					minWin = window;
					minType = 3;
				}
				if (minDiff > diff4) {
					minDiff = diff4;
					minWin = window;
					minType = 4;
				}
				if (minDiff > diff5) {
					minDiff = diff5;
					minWin = window;
					minType = 5;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			refPw[4] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					refPw[1] = ovInfo->power[i] * n;
					refPw[2] = ovInfo->power[i] * odd;
					refPw[3] = ovInfo->power[i] * (odd * odd);
					refPw[4] = ovInfo->power[i];
					refPw[5] = ovInfo->power[i] * (odd * odd * odd);
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (minType == 0) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[0] / hz;
					phaRand = 1;
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 1) {
					// 1/f(振幅が1/fになっているもの)
					tmpPw = refPw[1] / n;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 2) {
					// 短形波(奇数倍音,nの逆数で小さくなる)
					tmpPw = refPw[2] / odd;
					phaRand = 1;
					pha = ovInfo->phase[max_i];
					phaTmp = pha + param->sig2Phase;
					if (n & 0x01) {
						phaTmp += 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 3) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					phaTmp = pha + param->sig2Phase;
					if (n & 0x01) {
						phaTmp += 180;
					}
					if (phaTmp >= 360) {
						phaTmp %= 360;
					}
					if (phaTmp < 0) {
						phaTmp += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = phaTmp;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = phaTmp;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 4) {
					// パルス(n番目の倍音でもパワーは同じ)
					tmpPw = refPw[4];
					phaRand = rand() * 6;
					phaRand -= 3;
					pha += phaRand;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				} else if (minType == 5) {
					// 三角波(奇数倍音,n^2の逆数で小さくなる)
					tmpPw = refPw[3] / (odd * odd * odd);
					phaRand = 1;
//					pha = ovInfo->phase[max_i];
					pha += param->sig2Phase;
					if (pha >= 360) {
						pha %= 360;
					}
					if (pha < 0) {
						pha += 360;
					}
					if (ovInfo->pw[i] == 0) {
						ovInfo->pw[i] = tmpPw * 0.42;
						ovInfo->phase[i] = pha;
						ovInfo->diff[i] = minDiff;
					} else {
						if (ovInfo->diff[i] != -1 && ovInfo->diff[i] > minDiff) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
						if (ovInfo->pw[i] > tmpPw) {
							ovInfo->pw[i] = tmpPw * 0.42;
							ovInfo->phase[i] = pha;
							ovInfo->diff[i] = minDiff;
						}
					}
				}
			}
		}
	}
//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

#if 0
	if (param->sig3Enb == 1 && (param->sig3Phase < -6 || param->sig3Phase > 6)) {
		// 位相の修正
		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->phase[i];
					}
					
					diff = refPw[0];
					if (diff >= ovInfo->phase[i]) {
						diff = diff - ovInfo->phase[i];
					} else {
						diff = ovInfo->phase[i] - diff;
					}
					diff0 += diff;
					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
				}

				if (minDiff == -1 || minDiff > diff0) {
					minDiff = diff0;
					minWin = window;
					minType = 0;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
					pha += 360;
					pha += 360;
					pha %= 360;
				}
				if (ovInfo->pw[i] != 0) {
					if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
						ovInfo->phase[i] = pha;
					}
				}
			}
		}
	}
#endif
	if (param->sig1Enb == 1) {
		// powerが強いもの優先で補間する

		areaAvg = 0;
		for (i = lowIndex,n = 0;n < width;i++,n++) {
			areaAvg += ovInfo->power[i];
		}
		if (n > 0) {
			areaAvg /= n;
		}

		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			if (ovInfo->power[ofs] < areaAvg) {
				continue;
			}
			
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}


					if (i - window >= ofs) {
						if (ovInfo->power[i - window] < ovInfo->power[i]) {
							avgPw /= 75;
						}
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					if (ovInfo->power[i] > areaAvg) {
						avgPw += ovInfo->power[i];
					} else {
						avgPw /= 75;
					}
					avgPwNX++;
					nn++;
				}

//				if (avgPwNX > 0 && avgPwNX2 > 0 && avgPwNX3 > 0) {
//					tmpAvgPw  = avgPw / avgPwNX;
//					tmpAvgPw2 = avgPw2 / avgPwNX2;
//					tmpAvgPw3 = avgPw3 / avgPwNX3;
//					if ((tmpAvgPw  - (tmpAvgPw / 10)) > tmpAvgPw2 || tmpAvgPw + (tmpAvgPw / 10) < tmpAvgPw2 || (tmpAvgPw2 - (tmpAvgPw2 / 10)) > tmpAvgPw3 || tmpAvgPw2 + (tmpAvgPw2 / 10) < tmpAvgPw3) {//						continue;
//					}
//				}
				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff < avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが強いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig1Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = (tmpPw * 0.70);
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}
//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			fprintf(ofp,"\n");
//			for (i = 300;i < 750;i++) {
//				fprintf(ofp,"%lf,",ovInfo->pw[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//	}

	if (param->sig3Enb == 1) {
		// powerが弱いものと補間値がないものを補間する

		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}
				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				avgPw  = 0;
				avgPw2 = 0;
				avgPw3 = 0;
				avgPwNX  = 0;
				avgPwNX2 = 0;
				avgPwNX3 = 0;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->power[i] * hz;
						max_i = i;
					}
					avgPw += ovInfo->power[i];
					avgPwNX++;
					nn++;
				}

				if (avgPwNX > 0) {
					avgPw /= avgPwNX;
				}
				if (minDiff == -1 || minDiff > avgPw) {
					minDiff = avgPw;
					minWin = window;
					minType = 0;
				}
			}

			// 一番累積のパワーが弱いものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある
			odd = 1;
			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin,n++,odd+=2) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					refPw[0] = ovInfo->power[i] * hz;
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				// 1/f(振幅が1/fになっているもの)
				tmpPw = refPw[0] / hz;
				phaRand = 1;
				pha += param->sig3Phase;
				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
				}
				if (ovInfo->pw[i] == 0) {
					ovInfo->pw[i] = (tmpPw * 0.35);
					ovInfo->phase[i] = pha;
					ovInfo->diff[i] = -1;
				}
			}
		}
	}

	// 補間されていない箇所の対応
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	for (i = baseOfs;i + 1< nSample;i++) {
		hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
		if (hz < lowHz) {
			continue;
		}
		if (ovInfo->pw[i] == 0) {
			if (i - 2 >= 0 && i + 2 < nSample && ovInfo->pw[i - 2] > 0 && ovInfo->pw[i + 2] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-2] + ovInfo->pw[i+2]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-2] + ovInfo->phase[i+2]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
			if (ovInfo->pw[i - 1] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i-1] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i-1] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					ovInfo->phase[i] -= 360;
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
	if (param->sig3Phase >= -6 || param->sig3Phase <= 6) {
		// 位相の修正
		window = width;
		for (ofs = lowIndex;ofs < lowIndex + window;ofs++) {
			minDiff = -1;
			skipCnt = 0;
			for (window = swidth;window < (width / 2);window++) {
				if (window < minWidth) {
					if ((ofs - lowIndex) > window * 1) {
						continue;
					}
				} else {
					if ((ofs - lowIndex) > window) {
						continue;
					}
				}
				skipCnt++;
				if (param->hfaFast && (skipCnt % 8) != 0) {
					continue;
				}

				// スペクトル成分の遷移を調べる
				diff0 = 0;
				baseOfs = ofs - ((ofs / window) * window);
				if (baseOfs == 0) {
					baseOfs = window;
				}
				n = 1;		// 奇数偶数倍すべてに倍音がある
				refPw[0] = -1;
				for (i = baseOfs,nn = 0;i < highIndex; i+= window,n++) {
					hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
					if (hz < lowHz) {
						continue;
					}

					if (refPw[0] == -1) {
						refPw[0] = ovInfo->phase[i];
					}
					
					diff = refPw[0];
					if (diff >= ovInfo->phase[i]) {
						diff = diff - ovInfo->phase[i];
					} else {
						diff = ovInfo->phase[i] - diff;
					}
					diff0 += diff;
					nn++;
				}

				if (nn > 0) {
					diff0 /= nn;
				}

				if (minDiff == -1 || minDiff > diff0) {
					minDiff = diff0;
					minWin = window;
					minType = 0;
				}
			}

			// 一番予測誤差が少なかったものを採用する。

			baseOfs = ofs - ((ofs / minWin) * minWin);
			if (baseOfs == 0) {
				baseOfs = minWin;
			}

			pha = ovInfo->phase[baseOfs];
			n = 1;		// 奇数偶数倍すべてに倍音がある

			refPw[0] = -1;
			if (minWin == swidth || minWin == width - 1) {
				continue;
			}
			for (i = baseOfs;i < nSample;i += minWin) {
				hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
				if (hz < lowHz) {
					continue;
				}
				if (refPw[0] == -1) {
					max_i = i;
					pha = ovInfo->phase[max_i];
				}

				if (pha >= 360) {
					pha %= 360;
				}
				if (pha < 0) {
					pha += 360;
					pha += 360;
					pha += 360;
					pha %= 360;
				}
				if (ovInfo->pw[i] != 0) {
					if (ovInfo->diff[i] == -1 || ovInfo->diff[i] > minDiff) {
						ovInfo->phase[i] = pha;
					}
				}
			}
		}
	}
#if 0
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 5000);
	highIndex = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * (ovInfo->validSamplingRate - 2000);
	i = 1;
	do {
		for (j = i,n = 0,k = lowIndex;n < highIndex - lowIndex && j < ovInfo->nSample;j++,n++) {
			ovInfo->pw[j] += ovInfo->baseToAdj[k];
			k++;
			if (k > highIndex) {
				k = lowIndex;
			}
		}
		i += n;
	} while (i < ovInfo->nSample);
#endif
#if 0
	lowIndex  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * ovInfo->validSamplingRate;
	highIndex = nSample;
	width	  = ((double)ovInfo->nSample / (ovInfo->samplingRate / 2)) * 5000;
	i = lowIndex;
	do {
		avg = 0;
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			avg += ovInfo->pw[j];
		}
		if (n > 0) {
			avg /= n;
			avg = log10(avg) * 20;
		}
		for (j = i,n = 0;n < width && j < ovInfo->nSample;j++,n++) {
			tmpPw = ovInfo->pw[j];
			if (tmpPw > 0) {
				tmpPw = log10(tmpPw) * 20;
			}
			if (tmpPw + 15 < avg) {
				ovInfo->pw[j] *= 13;
			} else if (tmpPw > avg + 15) {
				ovInfo->pw[j] /= 13;
			}
		}
		i = j;
	} while (i < ovInfo->nSample);
#endif

#if 1
	if (param->hfaWide) {
		// 解析情報を平均化し、急激な変化を抑える
		for (i = baseOfs;i + 1< nSample;i++) {
			hz = ((ovInfo->samplingRate / 2) / (double)ovInfo->nSample) * i;
			if (hz < lowHz) {
				continue;
			}
			if (ovInfo->pw[i] > 0 && ovInfo->pw[i + 1] > 0) {
				ovInfo->pw[i] = (ovInfo->pw[i] + ovInfo->pw[i+1]) / 2;
				ovInfo->phase[i] = (ovInfo->phase[i] + ovInfo->phase[i+1]) / 2;
				if (ovInfo->phase[i] >= 360) {
					pha = ovInfo->phase[i];
					ovInfo->phase[i] = (pha % 360);
				}
				if (ovInfo->phase[i] < 0) {
					ovInfo->phase[i] += 360;
				}
			}
		}
	}
#endif

//	if (ovInfo->log) {
//		FILE *ofp;
//		ofp = fopen("d:\\fft.csv","a");
//		if (ofp) {
//			for (i = 0;i < 4096;i++) {
//				fprintf(ofp,"%lf,%lf\n",ovInfo->pw[i],ovInfo->phase[i]);
//			}
//			fprintf(ofp,"\n");
//			fflush(ofp);
//			fclose(ofp);
//		}
//		ovInfo->log = 0;
//	}
}
#endif
//---------------------------------------------------------------------------
// Function   : noiseCut
// Description: ノイズカット処理
// ---
//	nfs		 	:ノイズカットオフ周波数(この周波数以上の領域のノイズをカットする)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
void noiseCut(long nfs,DWORD inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3;
	long hfc;
	long outSampleR;
	long wkMemSize;
	long fftSize,i,j,n,nn,h;
	
	long lowIndex,highIndex;
	long hz,idx;
	double persent,per;
	
	double p;
	static double refPw[1000];
	double *pw;
	double cutLV[5]={0.99,1.06,1.12,1.18,1.30};
	SSIZE *pIn,*pOut;
	SSIZE startInSample,nSample;
	fftw_complex *fftw_in,*fftw_out;
	fftw_plan fftw_p,fftw_ip;

   	outSampleR = param->outSampleR;
	if (param->hfc != -1) {
		hfc = param->hfc;
	} else {
		hfc = param->inSampleR / 2;
	}
	if (hfc > param->inSampleR / 2) {
		hfc = param->inSampleR / 2;
	}

	fio_rewind(fp_r);
	fio_rewind(fp_w);

	if ((outSampleR == 32000) || (outSampleR == 44100) || (outSampleR == 48000)) {
		fftSize = 4096;
	}
	if ((outSampleR == 44100 * 2) || (outSampleR == 48000 * 2)) {
		fftSize = 8192;
	}
	if ((outSampleR == 32000 * 6) || (outSampleR == 44100 * 4) || (outSampleR == 48000 * 4)) {
		fftSize = 16384;
	}
	if ((outSampleR == 32000 * 12) || (outSampleR == 44100 * 8) || (outSampleR == 48000 * 8)) {
		fftSize = 32768;
	}
	if ((outSampleR == 32000 * 24) || (outSampleR == 44100 * 16) || (outSampleR == 48000 * 16)) {
		fftSize = 65536;
	}
	if ((outSampleR == 32000 * 48) || (outSampleR == 44100 * 32) || (outSampleR == 48000 * 32)) {
		fftSize = 65536*2;
	}

	wkMemSize = (65536 * 32) * sizeof (SSIZE);

	mem1 = (SSIZE *)malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem2 = (SSIZE *)malloc(wkMemSize);
	if (mem2 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	mem3 = (SSIZE *)malloc(wkMemSize);
	if (mem3 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	pw = (double *)malloc(65536*2 * sizeof (double));
	if (pw == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536*2);
	if (fftw_in == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 65536*2);
	if (fftw_out == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_p = fftw_plan_dft_1d(fftSize,fftw_in,fftw_out,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	fftw_ip = fftw_plan_dft_1d(fftSize,fftw_out,fftw_in,FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	per = -1;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + (fftSize + (fftSize / 2));startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
//				Sleep(1);
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
		}

		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 8) < inSample) {
			nSample = fftSize * 8;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 8;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 8;
		}

		memset(mem2,0,wkMemSize);
		memset(mem3,0,wkMemSize);
		memset(pw,0,65536 * sizeof (double));

		pIn = (SSIZE *)mem1;
		for (n = 0;n < 12;n++) {
			// FFT 初期設定
			for (i = 0;i < fftSize;i++) {
				fftw_in[i][0] = pIn[((fftSize / 2) * n) + i];
				fftw_in[i][1] = 0;
			}
			// 窓関数
			for (i = 0;i < (fftSize - 1) / 2;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 * i / (double)fftSize);
			}
			for (i = (fftSize - 1) / 2;i < fftSize;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 - 2.0 * i / (double)fftSize);
			}

			// FFT
			fftw_execute(fftw_p);

			// 元信号のパワーを累積する
			for (i = 1;i < fftSize / 2;i++) {
				p = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
				if (p != 0) {
					p = sqrt(p);
				}
				pw[i] += p;
			}
		}
		for (i = 0;i < fftSize / 2;i++) {
			pw[i] /= 12;
		}

		for (i = 0,h = 0;h < hfc - 1000;i++,h += 1000) {
			// 1khz 範囲のパワーを調べる
			lowIndex  = ((double)fftSize / outSampleR) * h;
			highIndex = ((double)fftSize / outSampleR) * (h + 1000);
			refPw[i] = 0;
			for (j = lowIndex,nn = 0;j < highIndex;j++,nn++) {
				refPw[i] += pw[j];
			}
			if (nn > 0) {
				refPw[i] /= nn;
			}
		}

		for (n = 0;n < 3;n++) {
			// FFT 初期設定
			for (i = 0;i < fftSize;i++) {
				fftw_in[i][0] = pIn[((fftSize / 2) * n) + i];
				fftw_in[i][1] = 0;
			}
			// 窓関数
			for (i = 0;i < (fftSize - 1) / 2;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 * i / (double)fftSize);
			}
			for (i = (fftSize - 1) / 2;i < fftSize;i++) {
				fftw_in[i][0] = fftw_in[i][0] * (2.0 - 2.0 * i / (double)fftSize);
			}

			// FFT
			fftw_execute(fftw_p);

			// 閾値より大きい音はカットする
			for (i = 1;i < fftSize / 2;i++) {
				hz = ((double)outSampleR / fftSize) * i;
				if (hz >= 100 && hz >= param->nr) {
					idx /= 1000;
					if (pw[i] >= refPw[idx] * cutLV[param->nrLV]) {
						fftw_out[i][0] = 0;
						fftw_out[i][1] = 0;
					} else {
						fftw_out[i][0] = (double)fftw_out[i][0] * 0.90;
						fftw_out[i][1] = (double)fftw_out[i][1] * 0.90;
					}
				} else {
					fftw_out[i][0] = 0;
					fftw_out[i][1] = 0;
				}
			}

			// 半分のデータを復元
			for (i = 1;i < fftSize / 2;i++) {
				fftw_out[fftSize - i][0] = fftw_out[i][0];
				fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
			}

			// invert FFT
			fftw_execute(fftw_ip);

			// 出力
			pOut = (SSIZE *)mem2;
			for (i = 0;i < fftSize;i++) {
				pOut[((fftSize / 2) * n) + i] += fftw_in[i][0] / fftSize;
			}
		}
		if (startInSample + fftSize / 2 >= 0) {
			outTempFile(fp_w,mem2 + fftSize / 2,fftSize,param);
		}
	}
	free(mem1);
	free(mem2);
	free(mem3);
	free(pw);
	fftw_destroy_plan(fftw_p);
	fftw_destroy_plan(fftw_ip);
	fftw_free(fftw_in);
	fftw_free(fftw_out);

}
//---------------------------------------------------------------------------
// Function   : bpFilter
// Description: 指定周波数をカットする
// ---
//	lfc		 	:低域のカットオフ周波数(この周波数以下の領域をカットする)
//	hfc		 	:高域のカットオフ周波数(この周波数以上の領域をカットする)
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	fp_w		:出力ファイル
//	param		:変換パラメータ
//
int bpFilter(long lfc,long hfc,DWORD inSample,FIO *fp_r,FIO *fp_w,PARAM_INFO *param)
{
	SSIZE *mem0,*mem1,*mem2,*mem3,*mem4;
	long outSampleR;
	long wkMemSize;
	long fftSize,i;
	
	
	
	
	SSIZE *pIn[3],*pOut[3];
	SSIZE startInSample,nSample;
	fftw_complex *fftw_in[3],*fftw_out[3];
	fftw_plan fftw_p[3],fftw_ip[3];

	outSampleR = param->outSampleR;
	fio_rewind(fp_r);
	fio_rewind(fp_w);

	fftSize = param->outSampleR * 2;

	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	mem1 = (SSIZE *)al_malloc(wkMemSize);
	if (mem1 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	mem2 = (SSIZE *)al_malloc(wkMemSize);
	if (mem2 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	mem3 = (SSIZE *)al_malloc(wkMemSize);
	if (mem3 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	mem4 = (SSIZE *)al_malloc(wkMemSize);
	if (mem4 == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_in[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_in[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_in[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_out[0] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_out[1] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_out[2] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_out[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_p[0] = fftw_plan_dft_1d(fftSize,fftw_in[0],fftw_out[0],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_p[1] = fftw_plan_dft_1d(fftSize,fftw_in[1],fftw_out[1],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_p[2] = fftw_plan_dft_1d(fftSize,fftw_in[2],fftw_out[2],FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}

	fftw_ip[0] = fftw_plan_dft_1d(fftSize,fftw_out[0],fftw_in[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[0] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_ip[1] = fftw_plan_dft_1d(fftSize,fftw_out[1],fftw_in[1],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[1] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}
	fftw_ip[2] = fftw_plan_dft_1d(fftSize,fftw_out[2],fftw_in[2],FFTW_BACKWARD,FFTW_ESTIMATE);
	if (fftw_ip[2] == NULL) {
		return STATUS_MEM_ALLOC_ERR;
	}


	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + fftSize + (fftSize / 2);startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
//			Sleep(0);
		}

		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		pIn[0]	= &mem1[((fftSize / 2) * 0)];
		pOut[0] = &mem2[((fftSize / 2) * 0)];
		pIn[1]	= &mem1[((fftSize / 2) * 1)];
		pOut[1] = &mem3[((fftSize / 2) * 1)];
		pIn[2]	= &mem1[((fftSize / 2) * 2)];
		pOut[2] = &mem4[((fftSize / 2) * 2)];

		memset(mem2,0,wkMemSize);
		memset(mem3,0,wkMemSize);
		memset(mem4,0,wkMemSize);

		#pragma omp parallel
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					// 1
					bpFilterSub(pIn[0],pOut[0],fftw_in[0],fftw_out[0],fftw_p[0],fftw_ip[0],lfc,hfc,param);
				}
				#pragma omp section
				{
					// 2
					bpFilterSub(pIn[1],pOut[1],fftw_in[1],fftw_out[1],fftw_p[1],fftw_ip[1],lfc,hfc,param);
				}
				#pragma omp section
				{
					// 3
					bpFilterSub(pIn[2],pOut[2],fftw_in[2],fftw_out[2],fftw_p[2],fftw_ip[2],lfc,hfc,param);
				}
			}
			#pragma omp for
			for (i = 0;i < nSample;i++) {
				mem2[i] += mem3[i] + mem4[i];
			}
		}

		if (startInSample + fftSize / 2 >= 0) {
			fio_seek(fp_w,(startInSample + (fftSize / 2)) * sizeof (SSIZE),SEEK_SET);
			outTempFile(fp_w,mem2 + fftSize / 2,fftSize,param);
			if (param->err) {
				break;
			}
		}
	}

	fio_flush(fp_w);

	al_free(mem1);
	al_free(mem2);
	al_free(mem3);
	al_free(mem4);

	fftw_destroy_plan(fftw_p[0]);
	fftw_destroy_plan(fftw_p[1]);
	fftw_destroy_plan(fftw_p[2]);

	fftw_destroy_plan(fftw_ip[0]);
	fftw_destroy_plan(fftw_ip[1]);
	fftw_destroy_plan(fftw_ip[2]);

	fftw_free(fftw_in[0]);
	fftw_free(fftw_in[1]);
	fftw_free(fftw_in[2]);

	fftw_free(fftw_out[0]);
	fftw_free(fftw_out[1]);
	fftw_free(fftw_out[2]);

	return STATUS_SUCCESS;
}
//---------------------------------------------------------------------------
// Function   : bpFilterSub
// Description: 指定周波数をカットする
// ---
//	param		:変換パラメータ
//
void bpFilterSub(SSIZE *pIn,SSIZE *pOut,fftw_complex *fftw_in,fftw_complex *fftw_out,fftw_plan fftw_p,fftw_plan fftw_ip,long lfc,long hfc,PARAM_INFO *param)
{
	long fftSize;
	long lowIndex,highIndex;
	long outSampleR;
	long i;

	fftSize = param->outSampleR * 2;
	outSampleR = param->outSampleR;
	
	// FFT 初期設定
	copyToFFTW(fftw_in,pIn,fftSize);

	// 窓関数
	windowFFTW(fftw_in,fftSize);

	// FFT
	fftw_execute(fftw_p);

	// 元信号の高域のパワーを調べる
	if (lfc != -1) {
		lowIndex = ((double)fftSize / outSampleR) * (lfc);
		
		// 低域カット
		for (i = 1;i < lowIndex;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}
	if (hfc != -1) {
		highIndex = ((double)fftSize / outSampleR) * (hfc);
		for (i = highIndex;i < fftSize / 2;i++) {
			fftw_out[i][0] = 0;
			fftw_out[i][1] = 0;
		}
	}
	// 半分のデータを復元
	for (i = 1;i < fftSize / 2;i++) {
		fftw_out[fftSize - i][0] = fftw_out[i][0];
		fftw_out[fftSize - i][1] = fftw_out[i][1] * -1;
	}

	// invert FFT
	fftw_execute(fftw_ip);

	// 出力
	for (i = 0;i < fftSize;i++) {
		pOut[i] += fftw_in[i][0] / fftSize;
	}
}
//---------------------------------------------------------------------------
// Function   : spAnalyze
// Description: スピーカーの周波数特性に応じて調整値パラメーターを出力する
// ---
//	inSample 	:処理するサンプル数(ch毎)
//	fp_r		:入力ファイル
//	param		:変換パラメータ
//
void spAnalyze(DWORD inSample,FIO *fp_r,PARAM_INFO *param)
{
	SSIZE *mem1,*mem0;
	long outSampleR;
	long inSampleR;
	long wkMemSize;
	long fftSize,i,n,hz;
	
	
	long validIndex,adjIndex,adjWidth;
	double persent,per;
	double div,step;
	long cnt;
	double p;
	static double adjData[192000];
	SSIZE *pIn;
	SSIZE startInSample,nSample;
	fftw_complex *fftw_in;
	fftw_plan fftw_p;
	double *adjFrom,*adjTo,*adjNx;
	FILE *ofp;
	outSampleR = param->outSampleR;
	inSampleR = param->inSampleR;

	fio_rewind(fp_r);

	fftSize = param->outSampleR * 2;

	wkMemSize = (fftSize * 2) * sizeof (SSIZE);

	validIndex = ((double)fftSize / param->outSampleR) * (inSampleR / 2);
	adjIndex   = ((double)fftSize / param->outSampleR) * (2000);
	adjWidth   = ((double)fftSize / param->outSampleR) * (1000);
	if (validIndex < adjIndex) {
		adjIndex = validIndex;
	}

	adjFrom = (double *)malloc(sizeof (double) * fftSize);
	adjTo	= (double *)malloc(sizeof (double) * fftSize);
	adjNx	= (double *)malloc(sizeof (double) * fftSize);
	if (adjFrom == NULL || adjTo == NULL || adjNx == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}
	for (i = 0;i < fftSize;i++) {
		adjFrom[i] = 0;
		adjTo[i] = 0;
		adjNx[i] = 0;
	}
	
	mem1 = (SSIZE *)malloc(wkMemSize);
	if (mem1 == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * fftSize);
	if (fftw_in == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	fftw_p = fftw_plan_dft_1d(fftSize,fftw_in,fftw_in,FFTW_FORWARD,FFTW_ESTIMATE);
	if (fftw_p == NULL) {
		param->err = STATUS_MEM_ALLOC_ERR;
		return;
	}

	per = -1;
	cnt = 0;
	for (startInSample = ((fftSize + (fftSize / 2)) * -1);startInSample < inSample + fftSize + (fftSize / 2);startInSample += fftSize) {
		if (startInSample >= 0 && startInSample < inSample) {
			persent = ((double)startInSample / inSample);
			persent *= 100;
			if (persent != per) {
				fprintf(stdout,"%d%%\n",(int)persent);
				fflush(stdout);
			}
			per = persent;
//			Sleep(1);
		}

		memset(mem1,0,wkMemSize);

		if (startInSample >= 0 && startInSample + (fftSize * 2) < inSample) {
			nSample = fftSize * 2;
			fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			fio_read(mem1,sizeof (SSIZE),nSample,fp_r);
		} else {
			mem0 = mem1;
			nSample = fftSize * 2;
			if (startInSample >= 0) {
				fio_seek(fp_r,startInSample * sizeof (SSIZE),SEEK_SET);
			} else {
				fio_seek(fp_r,0,SEEK_SET);
				mem0 += (startInSample * -1);
				if (nSample > startInSample * -1) {
					nSample -= startInSample * -1;
				} else {
					nSample = 0;
				}
			}

			if (startInSample >= inSample) {
				nSample = 0;
			} else {
				if (nSample != 0) {
					if (nSample > inSample - startInSample) {
						nSample = inSample - startInSample;
					}
				}
			}
			if (nSample > 0) {
				fio_read(mem0,sizeof (SSIZE),nSample,fp_r);
			}
			nSample = fftSize * 2;
		}

		pIn = (SSIZE *)mem1;
		for (n = 0;n < 3;n++) {
			// FFT 初期設定
			copyToFFTW(fftw_in,&pIn[((fftSize / 2) * n)],fftSize);

			// 窓関数
			windowFFTW(fftw_in,fftSize);

			// FFT
			fftw_execute(fftw_p);
			
			for (i = 1;i < validIndex;i++) {
				adjTo[i] = fftw_in[i][0] * fftw_in[i][0] + fftw_in[i][1] * fftw_in[i][1];
				if (adjTo[i] != 0) {
					adjTo[i] = sqrt(adjTo[i]);
				}
			}
			
			// 128 サイズの移動平均
			for (i = 1;i + 128 < validIndex;i++) {
				p = 0;
				for (n = 0;n < 128;n++) {
					p += adjTo[i + n];
				}
				if (p > 0) {
					p /= 128;
				}
				adjFrom[i] = p;
			}
			for (;i < fftSize / 2;i++) {
				adjFrom[i] = p;
			}

			p = 0;
			for (i = 1;i < 101;i++) {
				p += adjFrom[i];
			}
			p /= 100;
			
			adjTo[0] = 0;
			for (i = 1;i < fftSize / 2;i++) {
				adjTo[i] = p;
			}
			for (i = 1;i < fftSize / 2;i++) {
				if (adjFrom[i] != 0) {
					adjNx[i] += (adjTo[i] / adjFrom[i]);
				}
			}
			cnt++;
		}
	}
	if (cnt > 0) {
		for (i = 0;i < fftSize / 2;i++) {
			adjNx[i] /= cnt;
		}
	}
	
	p = 0;
	for (i = adjIndex,n = 0;n < adjWidth;adjIndex++,n++) {
		p += adjNx[i];
	}
	if (n > 0) {
		p /= n;
	}
	for (i = 1;i < adjIndex;i++) {
		adjNx[i] = p;
	}
	div = 1.0;
	step = 1.0 / 192000;
	for (i = adjIndex;i < fftSize / 2;i++) {
		adjNx[i] *= div;
		if (div - step > 0) {
			div -= step;
		} else {
			div = 0.08;
		}
	}

	if (param->r1_flag == 1) {
		unlink(param->sp_path);
		ofp = fopen(param->sp_path,"w");
		if (ofp) {
			for (i = 0;i < fftSize / 2;i++) {
				hz = ((double)param->outSampleR / fftSize) * i;
				adjData[hz] = adjNx[i];
			}
			for (i = 1;i < 192000;i++) {
				fprintf(ofp,"%lf\n",adjData[i]);
			}
			fclose(ofp);
		}
	}
	free(adjFrom);
	free(adjTo);
	free(adjNx);

	free(mem1);

	fftw_destroy_plan(fftw_p);
	fftw_free(fftw_in);

}

//---------------------------------------------------------------------------
// Function   : outTempFile
// Description: データをテンポラリファイルへ出力する
// ---
//	fp_w	:出力ファイル
//	in		:データのアドレス
//	size	:データー数
//	param	:パラメーター
//
void outTempFile(FIO *fp_w,void *in,DWORD size,PARAM_INFO *param)
{
	fio_size r;

	r = fio_write(in,1,size * sizeof (SSIZE),fp_w);
	if (r != size * sizeof (SSIZE)) {
		param->err = STATUS_FILE_WRITE_ERR;
	}
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
// Function   : adjPinkFilter
// Description: 1/f 特性にするフィルター
// ---
//	mode	  :モード(0,1,2,3)
//	fftSizeOut:FFT数
//	fftw_out2 :FFTW OUT 変数
//	param	  :変換パラメータ
//
void adjPinkFilter(int mode,long fftSizeOut,fftw_complex *fftw_out2,PARAM_INFO *param)
{
	long i;
	long startIdx,endIdx;
	double hz,div,step;
	
	long cutOff;
	long outSampleR;

	outSampleR = param->outSampleR;
	if (mode == 4 && param->lpf != -1) {
		startIdx = ((double)fftSizeOut / outSampleR) * param->lpf;
		endIdx	 = ((double)fftSizeOut / outSampleR) * (param->lpf + 20000);
		step = 1.00 / ((endIdx - startIdx) * 1.55);

		div = 1;
		for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
			if (div - step > 0) {
		   		div -= step;
			} else {
		   		div = 0.01;
			}
		}
		for (;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
		}
		div = 1;
		for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
			if (div - step > 0) {
		   		div -= step;
			} else {
		   		div = 0.01;
			}
		}
		for (;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
		}
		return;
	}

	if (mode == 1) {
		// 1/f 特性にするフィルター(hfa1)
		for (i = 1;i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			if (hz > 0) {
				fftw_out2[i][0] /= hz;
				fftw_out2[i][1] /= hz;
			}
		}
	}
	if (mode != 3) {
		// hfa1、hfa2、hfa3用の高域補間時の周波数調整
		if (param->hfa != 0 && param->hfc >= 8000 && param->hfc <= 23000) {
//			if (param->hfc > 13000) {
				startIdx = ((double)fftSizeOut / outSampleR) * 13000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 22000;
				step = 1.00 / ((endIdx - startIdx) * 1.20);
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				startIdx = ((double)fftSizeOut / outSampleR) * 17000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 45000;
				step = 1.00 / ((endIdx - startIdx) * 1.20);
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				startIdx = ((double)fftSizeOut / outSampleR) * 40000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 96000;
				step = 1.00 / ((endIdx - startIdx) * 1.40);
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				startIdx = ((double)fftSizeOut / outSampleR) * 41000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 120000;
				step = 1.00 / ((endIdx - startIdx) * 1.20);
				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
#if 0
			} else {
				startIdx = ((double)fftSizeOut / outSampleR) * 10000;
				endIdx	 = ((double)fftSizeOut / outSampleR) * 20000;

				step = 1.00 / ((endIdx - startIdx) * 1.68);

				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}

				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
				div = 1;
				for (i = startIdx;i < endIdx;i++) {
					hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
					if (hz >= (param->hfc + 1500)) {
						fftw_out2[i][0] *= div;
						fftw_out2[i][1] *= div;
						if (div - step > 0) {
					   		div -= step;
						} else {
					   		div = 0.01;
						}
					}
				}
				for (;i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
				}
			}
#endif
		}

		if (mode == 2) {
			// 独自のローパスフィルター
			cutOff = 40000;
			if (param->lpf > 1 && cutOff > param->lpf) {
				cutOff = param->lpf;
			}
			if ((outSampleR / 2) >= cutOff) {
				startIdx  = ((double)fftSizeOut / outSampleR) * cutOff;
				if (param->overSamp == 0) {
					endIdx = fftSizeOut / 2;
				} else {
					endIdx = fftSizeOut / 4;
				}
				step = 1.00 / ((endIdx - startIdx) * 1.30);

				div = 1;
				for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
					fftw_out2[i][0] *= div;
					fftw_out2[i][1] *= div;
					if (div - step > 0) {
						div -= step;
					} else {
						div = 0.01;
					}
				}
			}
		}
	} else {
		// デエンファシス用の処理
		if (param->deEmphasis == 1) {
			startIdx = ((double)fftSizeOut / outSampleR) * 3180;
			endIdx	 = ((double)fftSizeOut / outSampleR) * 10600;
		} else {
			startIdx = ((double)fftSizeOut / outSampleR) * 2100;
			endIdx	 = ((double)fftSizeOut / outSampleR) * 9520;
		}
		step = 1.00 / ((endIdx - startIdx) * 1.75);

		div = 1;
		for (i = startIdx;i < endIdx && i < fftSizeOut / 2;i++) {
			hz = (((double)(outSampleR / 2)) / (fftSizeOut / 2)) * i;
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
			if (div - step > 0) {
		   		div -= step;
			} else {
		   		div = 0.01;
			}
		}
		for (;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] *= div;
			fftw_out2[i][1] *= div;
		}
	}

	if (param->overSamp == 0) {
		for (i = (fftSizeOut / 2) - 5;i < fftSizeOut / 2;i++) {
			fftw_out2[i][0] = 0;
			fftw_out2[i][1] = 0;
		}
	} else {
		for (i = (fftSizeOut / 4) - 5;i < fftSizeOut / 4;i++) {
			fftw_out2[i][0] = 0;
			fftw_out2[i][1] = 0;
		}
	}
}
//---------------------------------------------------------------------------
// Function   : merageTempFile
// Description: 出力結果のファイルをマージする
// ---
//	type	 :マージの種類
//	normFlag :ノーマライズ用変数更新フラグ
//	fp_r 	 :入力ファイル1
//	fp_r2	 :入力ファイル2
//	fp_w	 :出力ファイル
//	inSample :サンプル数
//	param	 :パラメーター
//
void merageTempFile(char type,int normFlag,FIO *fp_r,FIO *fp_r2,FIO *fp_w,DWORD inSample,PARAM_INFO *param)
{
	
	SSIZE min,max;
	SSIZE maxLv,maxLv2;
	DWORD remainSample;
	DWORD ns;
	long i;
	fio_size remain1,remain2;
	fio_size wr_n;

	fio_rewind(fp_r);

	if (fp_r2 != NULL) {
		fio_rewind(fp_r2);
	}

	if (fp_w != NULL) {
		fio_rewind(fp_w);
	}

	min = max = 0;
	ns	= 0;
	maxLv = 0;
	maxLv2 = 0;
	remainSample = inSample;

	do {
		if (type == '+') {
//			Sleep(1);
			remain1 = fio_read(diskBuffer,sizeof (SSIZE),1 * 1024 * 1024,fp_r);
			if (fp_r2 != NULL) {
				remain2 = fio_read(diskBuffer2,sizeof (SSIZE),1 * 1024 * 1024,fp_r2);
			}
			if (remain1 == 0 || remain2 == 0) {
				break;
			}
			for (i = 0;i < remain1;i++) {
				if (diskBuffer[i] != 0) {
					diskBuffer[i] += diskBuffer2[i];
					if (diskBuffer[i] < min) {
						min = diskBuffer[i];
					}
					if (diskBuffer[i] > max) {
						max = diskBuffer[i];
					}
					if (diskBuffer[i] >> 40 > 0) {
						maxLv2 += diskBuffer[i] >> 40;
						ns++;
						if (maxLv2 >= 0x1000000000000) {
							maxLv2 /= ns;
							if (maxLv > 0) {
								maxLv = (maxLv + maxLv2) / 2;
							} else {
								maxLv = maxLv2;
							}
							maxLv2 = 0;
							ns = 0;
						}
					}
				}
			}
			remain1 =  remain1 < remainSample ? remain1 : remainSample;
			if (fp_w != NULL) {
//				fio_seek(fp_w,pos,SEEK_SET);
				wr_n = fio_write(diskBuffer,sizeof (SSIZE),remain1,fp_w);
				if (wr_n != remain1) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
			}
			remainSample -= remain1;
		} else if (type == '-') {
//			Sleep(1);
			remain1 = fio_read(diskBuffer,sizeof (SSIZE),1 * 1024 * 1024,fp_r);
			if (fp_r2 != NULL) {
				remain2 = fio_read(diskBuffer2,sizeof (SSIZE),1 * 1024 * 1024,fp_r2);
			}
			if (remain1 == 0) {
				break;
			}
			for (i = 0;i < remain1;i++) {
				if (diskBuffer[i] != 0) {
					diskBuffer[i] -= diskBuffer2[i];
					if (diskBuffer[i] < min) {
						min = diskBuffer[i];
					}
					if (diskBuffer[i] > max) {
						max = diskBuffer[i];
					}
					if (diskBuffer[i] >> 40 > 0) {
						maxLv2 += diskBuffer[i] >> 40;
						ns++;
						if (maxLv2 >= 0x1000000000000) {
							maxLv2 /= ns;
							if (maxLv > 0) {
								maxLv = (maxLv + maxLv2) / 2;
							} else {
								maxLv = maxLv2;
							}
							maxLv2 = 0;
							ns = 0;
						}
					}
				}
			}
			remain1 =  remain1 < remainSample ? remain1 : remainSample;
			if (fp_w != NULL) {
//				fio_seek(fp_w,pos,SEEK_SET);
				wr_n = fio_write(diskBuffer,sizeof (SSIZE),remain1,fp_w);
				if (wr_n != remain1) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
			}
			remainSample -= remain1;
		} else if (type == ' ') {
			//Sleep(1);
			remain1 = fio_read(diskBuffer,sizeof (SSIZE),1 * 1024 * 1024,fp_r);
			if (remain1 == 0) {
				break;
			}
			for (i = 0;i < remain1;i++) {
				if (diskBuffer[i] < min) {
					min = diskBuffer[i];
				}
				if (diskBuffer[i] > max) {
					max = diskBuffer[i];
				}
				if (diskBuffer[i] >> 40 > 0) {
					maxLv2 += diskBuffer[i] >> 40;
					ns++;
					if (maxLv2 >= 0x1000000000000) {
						maxLv2 /= ns;
						if (maxLv > 0) {
							maxLv = (maxLv + maxLv2) / 2;
						} else {
							maxLv = maxLv2;
						}
						ns = 0;
					}
				}
			}
			remain1 =  remain1 < remainSample ? remain1 : remainSample;
			if (fp_w != NULL) {
//				fio_seek(fp_w,pos,SEEK_SET);
				wr_n = fio_write(diskBuffer,sizeof (SSIZE),remain1,fp_w);
				if (wr_n != remain1) {
					param->err = STATUS_FILE_WRITE_ERR;
					break;
				}
			}
			remainSample -= remain1;
		}
	} while (remain1 == 1 * 1024 * 1024L && remainSample > 0);

	if (remainSample > 0 && param->err == STATUS_SUCCESS) {
		param->err = STATUS_FILE_READ_ERR;param->errLine = __LINE__;
	}

	if (param->err) {
		param->err = STATUS_FILE_READ_ERR;
		// エラー終了
		return;
	}
	
	if (fp_w != NULL) {
		fio_flush(fp_w);
	}
	if (normFlag == 1) {
		if (max > NormInfo.max) {
			NormInfo.max = max;
		}
		if (min < NormInfo.min) {
			NormInfo.min = min;
		}
		if (ns > 0) {
			maxLv2 /= ns;
		}
		if (maxLv > 0) {
			maxLv = (maxLv + maxLv2) / 2;
		} else {
			maxLv = maxLv2;
		}
		NormInfo.avg = maxLv;
	}
}
//---------------------------------------------------------------------------
// Function   : copyToFFTW
// Description: fftw用配列に値をコピーする
// ---
//	
//
void copyToFFTW(fftw_complex *fftw,SSIZE *buf,long size)
{
	long i;
	
	for (i = 0;i + 64 < size;i+=64) {
		fftw[i + 0][0] = buf[0];
		fftw[i + 0][1] = 0;
		fftw[i + 1][0] = buf[1];
		fftw[i + 1][1] = 0;
		fftw[i + 2][0] = buf[2];
		fftw[i + 2][1] = 0;
		fftw[i + 3][0] = buf[3];
		fftw[i + 3][1] = 0;
		fftw[i + 4][0] = buf[4];
		fftw[i + 4][1] = 0;
		fftw[i + 5][0] = buf[5];
		fftw[i + 5][1] = 0;
		fftw[i + 6][0] = buf[6];
		fftw[i + 6][1] = 0;
		fftw[i + 7][0] = buf[7];
		fftw[i + 7][1] = 0;
		fftw[i + 8][0] = buf[8];
		fftw[i + 8][1] = 0;
		fftw[i + 9][0] = buf[9];
		fftw[i + 9][1] = 0;
		fftw[i + 10][0] = buf[10];
		fftw[i + 10][1] = 0;
		fftw[i + 11][0] = buf[11];
		fftw[i + 11][1] = 0;
		fftw[i + 12][0] = buf[12];
		fftw[i + 12][1] = 0;
		fftw[i + 13][0] = buf[13];
		fftw[i + 13][1] = 0;
		fftw[i + 14][0] = buf[14];
		fftw[i + 14][1] = 0;
		fftw[i + 15][0] = buf[15];
		fftw[i + 15][1] = 0;
		fftw[i + 16][0] = buf[16];
		fftw[i + 16][1] = 0;
		fftw[i + 17][0] = buf[17];
		fftw[i + 17][1] = 0;
		fftw[i + 18][0] = buf[18];
		fftw[i + 18][1] = 0;
		fftw[i + 19][0] = buf[19];
		fftw[i + 19][1] = 0;
		fftw[i + 20][0] = buf[20];
		fftw[i + 20][1] = 0;
		fftw[i + 21][0] = buf[21];
		fftw[i + 21][1] = 0;
		fftw[i + 22][0] = buf[22];
		fftw[i + 22][1] = 0;
		fftw[i + 23][0] = buf[23];
		fftw[i + 23][1] = 0;
		fftw[i + 24][0] = buf[24];
		fftw[i + 24][1] = 0;
		fftw[i + 25][0] = buf[25];
		fftw[i + 25][1] = 0;
		fftw[i + 26][0] = buf[26];
		fftw[i + 26][1] = 0;
		fftw[i + 27][0] = buf[27];
		fftw[i + 27][1] = 0;
		fftw[i + 28][0] = buf[28];
		fftw[i + 28][1] = 0;
		fftw[i + 29][0] = buf[29];
		fftw[i + 29][1] = 0;
		fftw[i + 30][0] = buf[30];
		fftw[i + 30][1] = 0;
		fftw[i + 31][0] = buf[31];
		fftw[i + 31][1] = 0;
		fftw[i + 32][0] = buf[32];
		fftw[i + 32][1] = 0;
		fftw[i + 33][0] = buf[33];
		fftw[i + 33][1] = 0;
		fftw[i + 34][0] = buf[34];
		fftw[i + 34][1] = 0;
		fftw[i + 35][0] = buf[35];
		fftw[i + 35][1] = 0;
		fftw[i + 36][0] = buf[36];
		fftw[i + 36][1] = 0;
		fftw[i + 37][0] = buf[37];
		fftw[i + 37][1] = 0;
		fftw[i + 38][0] = buf[38];
		fftw[i + 38][1] = 0;
		fftw[i + 39][0] = buf[39];
		fftw[i + 39][1] = 0;
		fftw[i + 40][0] = buf[40];
		fftw[i + 40][1] = 0;
		fftw[i + 41][0] = buf[41];
		fftw[i + 41][1] = 0;
		fftw[i + 42][0] = buf[42];
		fftw[i + 42][1] = 0;
		fftw[i + 43][0] = buf[43];
		fftw[i + 43][1] = 0;
		fftw[i + 44][0] = buf[44];
		fftw[i + 44][1] = 0;
		fftw[i + 45][0] = buf[45];
		fftw[i + 45][1] = 0;
		fftw[i + 46][0] = buf[46];
		fftw[i + 46][1] = 0;
		fftw[i + 47][0] = buf[47];
		fftw[i + 47][1] = 0;
		fftw[i + 48][0] = buf[48];
		fftw[i + 48][1] = 0;
		fftw[i + 49][0] = buf[49];
		fftw[i + 49][1] = 0;
		fftw[i + 50][0] = buf[50];
		fftw[i + 50][1] = 0;
		fftw[i + 51][0] = buf[51];
		fftw[i + 51][1] = 0;
		fftw[i + 52][0] = buf[52];
		fftw[i + 52][1] = 0;
		fftw[i + 53][0] = buf[53];
		fftw[i + 53][1] = 0;
		fftw[i + 54][0] = buf[54];
		fftw[i + 54][1] = 0;
		fftw[i + 55][0] = buf[55];
		fftw[i + 55][1] = 0;
		fftw[i + 56][0] = buf[56];
		fftw[i + 56][1] = 0;
		fftw[i + 57][0] = buf[57];
		fftw[i + 57][1] = 0;
		fftw[i + 58][0] = buf[58];
		fftw[i + 58][1] = 0;
		fftw[i + 59][0] = buf[59];
		fftw[i + 59][1] = 0;
		fftw[i + 60][0] = buf[60];
		fftw[i + 60][1] = 0;
		fftw[i + 61][0] = buf[61];
		fftw[i + 61][1] = 0;
		fftw[i + 62][0] = buf[62];
		fftw[i + 62][1] = 0;
		fftw[i + 63][0] = buf[63];
		fftw[i + 63][1] = 0;
		buf += 64;
	}
	for (;i < size;i++) {
		fftw[i][0] = *buf++;
		fftw[i][1] = 0;
	}
}
//---------------------------------------------------------------------------
// Function   : windowFFTW
// Description: FFTW用Window関数
// ---
//	
//
void windowFFTW(fftw_complex *fftw,long size)
{
	long i;

	// ウインドウサイズ毎に定数化する
	switch (size) {
		case (4096 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 1));
			}
			#pragma omp parallel for
			for (i = ((4096 * 1) - 1) / 2;i < (4096 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 1));
			}
			break;
		case (4096 * 2):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 2) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 2));
			}
			#pragma omp parallel for
			for (i = ((4096 * 2) - 1) / 2;i < (4096 * 2);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 2));
			}
			break;
		case (4096 * 4):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 4) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 4));
			}
			#pragma omp parallel for
			for (i = ((4096 * 4) - 1) / 2;i < (4096 * 4);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 4));
			}
			break;
		case (4096 * 8):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 8) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 8));
			}
			#pragma omp parallel for
			for (i = ((4096 * 8) - 1) / 2;i < (4096 * 8);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 8));
			}
			break;
		case (4096 * 16):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 16) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 16));
			}
			#pragma omp parallel for
			for (i = ((4096 * 16) - 1) / 2;i < (4096 * 16);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 16));
			}
			break;
		case (4096 * 32):
			#pragma omp parallel for
			for (i = 0;i < ((4096 * 32) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(4096 * 32));
			}
			#pragma omp parallel for
			for (i = ((4096 * 32) - 1) / 2;i < (4096 * 32);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(4096 * 32));
			}
			break;
		case (32000 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((32000 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(32000 * 1));
			}
			#pragma omp parallel for
			for (i = ((32000 * 1) - 1) / 2;i < (32000 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(32000 * 1));
			}
			break;
		case (44100 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 1));
			}
			#pragma omp parallel for
			for (i = ((44100 * 1) - 1) / 2;i < (44100 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 1));
			}
			break;
		case (48000 * 1):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 1) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 1));
			}
			#pragma omp parallel for
			for (i = ((48000 * 1) - 1) / 2;i < (48000 * 1);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 1));
			}
			break;
		case (44100 * 2):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 2) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 2));
			}
			#pragma omp parallel for
			for (i = ((44100 * 2) - 1) / 2;i < (44100 * 2);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 2));
			}
			break;
		case (48000 * 2):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 2) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 2));
			}
			#pragma omp parallel for
			for (i = ((48000 * 2) - 1) / 2;i < (48000 * 2);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 2));
			}
			break;
		case (44100 * 4):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 4) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 4));
			}
			#pragma omp parallel for
			for (i = ((44100 * 4) - 1) / 2;i < (44100 * 4);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 4));
			}
			break;
		case (48000 * 4):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 4) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 4));
			}
			#pragma omp parallel for
			for (i = ((48000 * 4) - 1) / 2;i < (48000 * 4);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 4));
			}
			break;
		case (44100 * 8):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 8) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 8));
			}
			#pragma omp parallel for
			for (i = ((44100 * 8) - 1) / 2;i < (44100 * 8);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 8));
			}
			break;
		case (48000 * 8):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 8) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 8));
			}
			#pragma omp parallel for
			for (i = ((48000 * 8) - 1) / 2;i < (48000 * 8);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 8));
			}
			break;
		case (44100 * 16):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 16) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 16));
			}
			#pragma omp parallel for
			for (i = ((44100 * 16) - 1) / 2;i < (44100 * 16);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 16));
			}
			break;
		case (48000 * 16):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 16) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 16));
			}
			#pragma omp parallel for
			for (i = ((48000 * 16) - 1) / 2;i < (48000 * 16);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 16));
			}
			break;
		case (44100 * 32):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 32) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 32));
			}
			#pragma omp parallel for
			for (i = ((44100 * 32) - 1) / 2;i < (44100 * 32);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 32));
			}
			break;
		case (48000 * 32):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 32) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 32));
			}
			#pragma omp parallel for
			for (i = ((48000 * 32) - 1) / 2;i < (48000 * 32);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 32));
			}
			break;
		case (44100 * 64):
			#pragma omp parallel for
			for (i = 0;i < ((44100 * 64) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(44100 * 64));
			}
			#pragma omp parallel for
			for (i = ((44100 * 64) - 1) / 2;i < (44100 * 64);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(44100 * 64));
			}
			break;
		case (48000 * 64):
			#pragma omp parallel for
			for (i = 0;i < ((48000 * 64) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(48000 * 64));
			}
			#pragma omp parallel for
			for (i = ((48000 * 64) - 1) / 2;i < (48000 * 64);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(48000 * 64));
			}
			break;
		default:
			for (i = 0;i < ((size) - 1) / 2;i++) {
				fftw[i][0] = fftw[i][0] * (2.0 * i / (double)(size));
			}
			for (i = ((size) - 1) / 2;i < (size);i++) {
				fftw[i][0] = fftw[i][0] * (2.0 - 2.0 * i / (double)(size));
			}

			break;
	}
}
//---------------------------------------------------------------------------
// Function   : cutFFTW
// Description: FFTW用カットオフ関数
// ---
//	
//
void cutFFTW(fftw_complex *fftw,long index,long size)
{
	long i;

	// 64 個ずつ
	for (i = index;i + 64 < size;i+= 64) {
		fftw[i + 0][0] = 0;
		fftw[i + 0][1] = 0;
		fftw[i + 1][0] = 0;
		fftw[i + 1][1] = 0;
		fftw[i + 2][0] = 0;
		fftw[i + 2][1] = 0;
		fftw[i + 3][0] = 0;
		fftw[i + 3][1] = 0;
		fftw[i + 4][0] = 0;
		fftw[i + 4][1] = 0;
		fftw[i + 5][0] = 0;
		fftw[i + 5][1] = 0;
		fftw[i + 6][0] = 0;
		fftw[i + 6][1] = 0;
		fftw[i + 7][0] = 0;
		fftw[i + 7][1] = 0;
		fftw[i + 8][0] = 0;
		fftw[i + 8][1] = 0;
		fftw[i + 9][0] = 0;
		fftw[i + 9][1] = 0;
		fftw[i + 10][0] = 0;
		fftw[i + 10][1] = 0;
		fftw[i + 11][0] = 0;
		fftw[i + 11][1] = 0;
		fftw[i + 12][0] = 0;
		fftw[i + 12][1] = 0;
		fftw[i + 13][0] = 0;
		fftw[i + 13][1] = 0;
		fftw[i + 14][0] = 0;
		fftw[i + 14][1] = 0;
		fftw[i + 15][0] = 0;
		fftw[i + 15][1] = 0;
		fftw[i + 16][0] = 0;
		fftw[i + 16][1] = 0;
		fftw[i + 17][0] = 0;
		fftw[i + 17][1] = 0;
		fftw[i + 18][0] = 0;
		fftw[i + 18][1] = 0;
		fftw[i + 19][0] = 0;
		fftw[i + 19][1] = 0;
		fftw[i + 20][0] = 0;
		fftw[i + 20][1] = 0;
		fftw[i + 21][0] = 0;
		fftw[i + 21][1] = 0;
		fftw[i + 22][0] = 0;
		fftw[i + 22][1] = 0;
		fftw[i + 23][0] = 0;
		fftw[i + 23][1] = 0;
		fftw[i + 24][0] = 0;
		fftw[i + 24][1] = 0;
		fftw[i + 25][0] = 0;
		fftw[i + 25][1] = 0;
		fftw[i + 26][0] = 0;
		fftw[i + 26][1] = 0;
		fftw[i + 27][0] = 0;
		fftw[i + 27][1] = 0;
		fftw[i + 28][0] = 0;
		fftw[i + 28][1] = 0;
		fftw[i + 29][0] = 0;
		fftw[i + 29][1] = 0;
		fftw[i + 30][0] = 0;
		fftw[i + 30][1] = 0;
		fftw[i + 31][0] = 0;
		fftw[i + 31][1] = 0;
		fftw[i + 32][0] = 0;
		fftw[i + 32][1] = 0;
		fftw[i + 33][0] = 0;
		fftw[i + 33][1] = 0;
		fftw[i + 34][0] = 0;
		fftw[i + 34][1] = 0;
		fftw[i + 35][0] = 0;
		fftw[i + 35][1] = 0;
		fftw[i + 36][0] = 0;
		fftw[i + 36][1] = 0;
		fftw[i + 37][0] = 0;
		fftw[i + 37][1] = 0;
		fftw[i + 38][0] = 0;
		fftw[i + 38][1] = 0;
		fftw[i + 39][0] = 0;
		fftw[i + 39][1] = 0;
		fftw[i + 40][0] = 0;
		fftw[i + 40][1] = 0;
		fftw[i + 41][0] = 0;
		fftw[i + 41][1] = 0;
		fftw[i + 42][0] = 0;
		fftw[i + 42][1] = 0;
		fftw[i + 43][0] = 0;
		fftw[i + 43][1] = 0;
		fftw[i + 44][0] = 0;
		fftw[i + 44][1] = 0;
		fftw[i + 45][0] = 0;
		fftw[i + 45][1] = 0;
		fftw[i + 46][0] = 0;
		fftw[i + 46][1] = 0;
		fftw[i + 47][0] = 0;
		fftw[i + 47][1] = 0;
		fftw[i + 48][0] = 0;
		fftw[i + 48][1] = 0;
		fftw[i + 49][0] = 0;
		fftw[i + 49][1] = 0;
		fftw[i + 50][0] = 0;
		fftw[i + 50][1] = 0;
		fftw[i + 51][0] = 0;
		fftw[i + 51][1] = 0;
		fftw[i + 52][0] = 0;
		fftw[i + 52][1] = 0;
		fftw[i + 53][0] = 0;
		fftw[i + 53][1] = 0;
		fftw[i + 54][0] = 0;
		fftw[i + 54][1] = 0;
		fftw[i + 55][0] = 0;
		fftw[i + 55][1] = 0;
		fftw[i + 56][0] = 0;
		fftw[i + 56][1] = 0;
		fftw[i + 57][0] = 0;
		fftw[i + 57][1] = 0;
		fftw[i + 58][0] = 0;
		fftw[i + 58][1] = 0;
		fftw[i + 59][0] = 0;
		fftw[i + 59][1] = 0;
		fftw[i + 60][0] = 0;
		fftw[i + 60][1] = 0;
		fftw[i + 61][0] = 0;
		fftw[i + 61][1] = 0;
		fftw[i + 62][0] = 0;
		fftw[i + 62][1] = 0;
		fftw[i + 63][0] = 0;
		fftw[i + 63][1] = 0;
	}
	// 残り
	for (;i < size;i++) {
		fftw[i + 0][0] = 0;
		fftw[i + 0][1] = 0;
	}
}
//---------------------------------------------------------------------------
// Function   : al_malloc
// Description: 16バイト境界対応malloc関数
// ---
// 返すポインタの16バイト前にmallocで確保したメモリ領域のアドレスを入れる
//
void *al_malloc(long size)
{
	void *ptr;
	void *ret;
	int align;

	ptr = malloc(size + 32);
	if (ptr) {
		ret = ptr;
		align = (int)ptr % 16;
		if (align != 0) {
			align = 16 - align;
			ret = (char *)ret + align;
		}
		*((SSIZE *)ret) = (SSIZE)ptr;

		ret = (char *)ret + 16;
	} else {
		ret = NULL;
	}
	return ret;
}
//---------------------------------------------------------------------------
// Function   : al_free
// Description: 16バイト境界対応free関数
// ---
// 
//
void *al_free(void *ptr)
{
	void *p;
	
	if (ptr) {
		p = (char *)ptr - 16;
		p = (void *)(*((SSIZE *)p));
		free(p);
	}
}
