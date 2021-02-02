// ---------------------------------------------------------------------------
/** ************************************************************************* */
/* PLG_AUDIO_IO (C) 2008-2012 By 59414d41 */
/* */
/* */
/** ************************************************************************* */

/* --- Log -------------------------------------------------------------------
 * Ver 0.10 <09/06/02> - VC6用のDLLとして作成
 * Ver 0.50 <10/10/23> - MP3,FLAC,WAVPACKファイル対応(情報取得のみ対応)
 *					   - BWF 対応(保存のみ)
 *					   - RF64 対応(保存のみ)
 * Ver 0.60 <10/11/25> - Multi Channel 対応
 * Ver 0.61 <10/12/30> - Extensible 対応
 * Ver 0.70 <11/07/24> - APIを整理
 * Ver 0.80 <12/02/11> - DLL をやめて、各プログラムへリンクするように修正
 *						 mmsystemを使用しないように修正
 *						 RF64 ファイルの読み込みに対応
 */

// ---------------------------------------------------------------------------

#define _LARGEFILE_SOURCE				// for GCC
#define _FILE_OFFSET_BITS 64			// for GCC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "PLG_AudioIO.h"
#define strnicmp strncasecmp
#include <stdint.h>
#ifdef __GUNC__
#include <sys/fcntl.h>
#include <sys/stat.h>
#else
#include <fcntl.h>
#include <sys/stat.h>
//#include <sys/io.h>
#include<sys/types.h>
#endif
#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("/tmp/lg_audio.log","a");								\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n","PLG",__LINE__,s);					\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif

unsigned char tempBuffer[1 * 1024 * 1024L];

#ifndef TRUE
#define TRUE	(1)
#endif
#ifndef FALSE
#define FALSE	(0)
#endif

// Wav ファイル
#pragma pack(push, 1)
typedef struct {
	char id[4];
	unsigned int size;
	unsigned char data[1];
}WAVE_HEADER;

// Wave file metadata(LIST)
typedef struct {
	BYTE h_list[4];
	unsigned int s_list;
	BYTE h_info[4];
	BYTE h_isrc[4];
	unsigned int s_isrc;
	BYTE *p_isrc;
	BYTE h_icrd[4];
	unsigned int s_icrd;
	BYTE *p_icrd;
	BYTE h_iprd[4];
	unsigned int s_iprd;
	BYTE *p_iprd;
	BYTE h_ieng[4];
	unsigned int s_ieng;
	BYTE *p_ieng;
	BYTE h_inam[4];
	unsigned int s_inam;
	BYTE *p_inam;
	BYTE h_icop[4];
	unsigned int s_icop;
	BYTE *p_icop;
	BYTE h_ignr[4];
	unsigned int s_ignr;
	BYTE *p_ignr;
	BYTE h_isft[4];
	unsigned int s_isft;
	BYTE *p_isft;
	BYTE h_icmt[4];
	unsigned int s_icmt;
	BYTE *p_icmt;
	BYTE h_iart[4];
	unsigned int s_iart;
	BYTE *p_iart;
}META_LIST;
#pragma pack(pop)

// サンプルを処理するデータ型
#define SSIZE	long
#define UI64	unsigned long
#define int64_t long

//
// DSF ファイルフォーマット仕様書を参照
#pragma pack(push, 1)
typedef struct {
	char id[4];
	UI64 chunk_size;
	UI64 file_size;
	UI64 ptr_meta;
}DSF;

typedef struct {
	char id[4];
	UI64 chunk_size;
	unsigned int fmt_version;
	unsigned int fmt_id;
	unsigned int channel_type;
	unsigned int channel_count;
	unsigned int sampling;
	unsigned int sample_bit_count;
	UI64 sample_count;
	unsigned int block_size;
	unsigned int reserved;
}DSF_FMT;

typedef struct {
	char id[4];
	UI64 chunk_size;
}DSF_DATA;
#pragma pack(pop)

int infoWavFile(char*filename, SOUNDFMT * pFmt, DWORD * pInSample,FILEINFO * pFileinfo);
int infoRF64File(char*filename, SOUNDFMT * pFmt, unsigned int * pInSample,FILEINFO * pFileinfo);
int infoMP3File(char*filename, SOUNDFMT * pFmt, unsigned int * pInSample,FILEINFO * pFileinfo);
int infoFlacFile(char*filename, SOUNDFMT * pFmt, unsigned int * pInSample,FILEINFO * pFileinfo);
int infoWavPackFile(char*filename, SOUNDFMT * pFmt, unsigned int * pInSample,FILEINFO * pFileinfo);
int infoDsfFile(char*filename, SOUNDFMT * pFmt, unsigned int * pInSample,FILEINFO * pFileinfo);
int findWavChunk(FILE * fp, char*id, long*size);
static int64_t ftell_local(FILE * fp);
static int fseek_local(FILE * fp, int64_t offset, int origin);

int infoWavFile(char *filename, SOUNDFMT *pFmt, DWORD *pInSample,FILEINFO *pFileinfo)
/*
 * WavFile の情報取得(2G以下のwavファイルのみ対応)
 * 2019.04.09 サイズチェックは行っていない。内部的には4Gまで許容。
 */ 
{
	WAVE_HEADER *wh;
	WFEX wfex;
	int64_t seekPtr;
	unsigned int pcmSize;
	unsigned int chunkSize;
	char subTypePCM[16] = {
		0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,
		0x00, 0x38, 0x9b, 0x71
	};
	char subTypeFloat[16] = {
		0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,
		0x00, 0x38, 0x9b, 0x71
	};
	FILE *fp;
	long rs;
	int fmt_chunk, data_chunk;

	int fileIsWAV = FALSE;
	fmt_chunk = 0;
	data_chunk = 0;
	wh = (WAVE_HEADER*)tempBuffer;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// WAV
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RIFF", 4) == 0 || memcmp(tempBuffer, "riff",4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE", 4) == 0 || memcmp(tempBuffer + 8,"wave", 4) == 0)) {
			break;
		}
		// 2019.04.09 サイズチェックは無視
		/*if (wh->size > 2147483647 || wh->size + 8 > 2147483647) {
			// size over
			break;
		}*/

		// fmt チャンクのサーチ
		if (findWavChunk(fp, "fmt ", &chunkSize)) {
			fmt_chunk = 1;
			rs = fread(&wfex, 1, sizeof(WFEX), fp);
			if (rs < 16) {
				break;
			}
			if (wfex.Format.wFormatTag == WF_PCM) {
				if ((wfex.Format.wBitsPerSample == 16 || wfex.Format.wBitsPerSample == 24|| wfex.Format.wBitsPerSample == 32 || wfex.Format.wBitsPerSample == 64) && wfex.Format.nChannels <= 2) {
					pFmt->sample = wfex.Format.nSamplesPerSec;
					pFmt->channel = (unsigned char)wfex.Format.nChannels;
					pFmt->bitsPerSample = (unsigned char)wfex.Format.wBitsPerSample;
					pFmt->type = 1;
					fileIsWAV = TRUE;
				}
			} else if (wfex.Format.wFormatTag == WF_IEEE_FLOAT) {
				if ((wfex.Format.wBitsPerSample == 32 || wfex.Format.wBitsPerSample == 64) && wfex.Format.nChannels <= 2) {
					pFmt->sample = wfex.Format.nSamplesPerSec;
					pFmt->channel = (unsigned char)wfex.Format.nChannels;
					pFmt->bitsPerSample = (unsigned char)wfex.Format.wBitsPerSample;
					pFmt->type = 3;
					fileIsWAV = TRUE;
				}
			} else if (wfex.Format.wFormatTag == 0xFFFE && chunkSize >= sizeof(WFEX)) {
				if ((memcmp(wfex.subFormat, subTypePCM,16) == 0 && wfex.Format.wBitsPerSample == 16 || wfex.Format.wBitsPerSample == 24) ||
					(memcmp(wfex.subFormat, subTypeFloat,16) == 0 && wfex.Format.wBitsPerSample == 32 || wfex.Format.wBitsPerSample == 64)) {
					if (wfex.Format.nChannels <= 6) {
						pFmt->sample = wfex.Format.nSamplesPerSec;
						pFmt->channel = (unsigned char)wfex.Format.nChannels;
						pFmt->bitsPerSample = (unsigned char)wfex.Format.wBitsPerSample;
						pFmt->type = 3;
						fileIsWAV = TRUE;
					}
				}
			}
		}
		if (findWavChunk(fp, "data", &chunkSize)) {
			data_chunk = 1;
			seekPtr = ftell_local(fp);
			if (pFileinfo != NULL) {
				pFileinfo->dataOffset = seekPtr;
				pFileinfo->dataSize = wh->size;
			}
			pcmSize = wh->size;
			*pInSample = (pcmSize / pFmt->channel) / (wfex.Format.wBitsPerSample / 8);
		}
		// 17.5.6 fix && -> ||
		if (fmt_chunk || data_chunk) {
			pFmt->fmt[0] = 'w';
			pFmt->fmt[1] = 'a';
			pFmt->fmt[2] = 'v';
			pFmt->fmt[3] = '\0';
			fileIsWAV = TRUE;
			fclose(fp);
			fp = NULL;
		} else {
			fileIsWAV = FALSE;
			fclose(fp);
			fp = NULL;
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsWAV;
}

int infoRF64File(char *filename, SOUNDFMT *pFmt, unsigned int *pInSample,FILEINFO *pFileinfo)
/*
 * RF64ファイルの情報取得
 */
{
	WAVE_HEADER *wh;
	RF64_DS64 ds64;
	RF64_FMT ds64_fmt;
	int64_t seekPtr;
	int64_t pcmSize;
	unsigned int chunkSize;
	FILE *fp;
	long rs;
	int fmt_chunk, data_chunk, ds64_chunk;

	int fileIsRF64 = FALSE;
	fmt_chunk = 0;
	ds64_chunk = 0;
	data_chunk = 0;
	pcmSize = 0;
	wh = (WAVE_HEADER*)tempBuffer;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// RF64
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RF64", 4) == 0 || memcmp(tempBuffer, "rf64",4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE", 4) == 0 || memcmp(tempBuffer + 8,"wave", 4) == 0)) {
			break;
		}
		if (wh->size != 0xFFFFFFFF) {
			break;
		}
		// ds64 チャンクのサーチ
		if (findWavChunk(fp, "ds64", &chunkSize)) {
			if (fseek_local(fp, -8, SEEK_CUR)) {
				break;
			}
			ds64_chunk = 1;
			rs = fread(&ds64, sizeof(RF64_DS64), 1, fp);
			if (rs != 1) {
				break;
			}
			if (ds64.riffSizeLow == 0 && ds64.riffSizeHigh == 0) {
				break;
			}
			if (ds64.dataSizeLow == 0 && ds64.dataSizeHigh == 0) {
				break;
			}
			if (ds64.sampleCountLow == 0 && ds64.sampleCountHigh == 0) {
				break;
			}
			pcmSize = ((int64_t)(ds64.dataSizeHigh) << 32) + (int64_t)(ds64.dataSizeLow);
			if (pFileinfo != NULL) {
				pFileinfo->dataSize = pcmSize;
			}
		}
		// fmt チャンクのサーチ
		if (findWavChunk(fp, "fmt ", &chunkSize)) {
			if (fseek_local(fp, -8, SEEK_CUR)) {
				break;
			}
			fmt_chunk = 1;
			rs = fread(&ds64_fmt, sizeof(RF64_FMT), 1, fp);
			if (rs != 1) {
				break;
			}
			if (ds64_fmt.formatType == WF_PCM) {
				if((ds64_fmt.bitsPerSample == 16 || ds64_fmt.bitsPerSample == 24) && ds64_fmt.channelCount <= 2) {
					pFmt->sample = ds64_fmt.sampleRate;
					pFmt->channel = (unsigned char)ds64_fmt.channelCount;
					pFmt->bitsPerSample = (unsigned char)ds64_fmt.bitsPerSample;
					fileIsRF64 = TRUE;
				}
			}
			else if (ds64_fmt.formatType == WF_IEEE_FLOAT) {
				if((ds64_fmt.bitsPerSample == 32 || ds64_fmt.bitsPerSample == 64) && ds64_fmt.channelCount <= 2) {
					pFmt->sample = ds64_fmt.sampleRate;
					pFmt->channel = (unsigned char)ds64_fmt.channelCount;
					pFmt->bitsPerSample = (unsigned char)ds64_fmt.bitsPerSample;
					fileIsRF64 = TRUE;
				}
			}
		}
		if (findWavChunk(fp, "data", &chunkSize)) {
			data_chunk = 1;
			seekPtr = ftell_local(fp);
			if (pFileinfo != NULL) {
				pFileinfo->dataOffset = seekPtr;
			}
		}
		if (fmt_chunk && data_chunk && ds64_chunk) {
			pFmt->fmt[0] = 'r';
			pFmt->fmt[1] = 'f';
			pFmt->fmt[2] = '6';
			pFmt->fmt[3] = '4';
			pFmt->fmt[4] = '\0';
			PRINT_LOG("");
			*pInSample = ((pcmSize / pFmt->channel) / (pFmt->bitsPerSample / 8)) <= 0xFFFFFFFF ? ((pcmSize / pFmt->channel) / (pFmt->bitsPerSample / 8)) : 0;
			{
				char ppp[50];
				sprintf(ppp, "inSample:%ld", *pInSample);
				PRINT_LOG(ppp);
			}
			if (*pInSample > 0) {
				// サンプル数を32bitに制限する(本ソフトの仕様)
				fileIsRF64 = TRUE;
			}
			fclose(fp);
			fp = NULL;
		} else {
			fileIsRF64 = FALSE;
			fclose(fp);
			fp = NULL;
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsRF64;
}

int infoMP3File(char *filename, SOUNDFMT *pFmt, unsigned int *pInSample,
	FILEINFO *pFileinfo)
/*
 * MP3 Fileの情報取得
 */ {
	FILE *fp;
	long rs;
	int64_t seekPtr;
	int fileIsMP3 = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// MP3
		rs = fread(tempBuffer, 1, 10, fp);
		if (rs != 10) {
			break;
		}
		if (tempBuffer[0] == 'I' && tempBuffer[1] == 'D' && tempBuffer[2] == '3') {
			// ID3 v2 Tag
			seekPtr = (tempBuffer[6] & 0x7F);
			seekPtr <<= 7;
			seekPtr |= (tempBuffer[7] & 0x7F);
			seekPtr <<= 7;
			seekPtr |= (tempBuffer[8] & 0x7F);
			seekPtr <<= 7;
			seekPtr |= (tempBuffer[9] & 0x7F);
			seekPtr += 10;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(tempBuffer, 1, 4, fp);
			if (rs != 4) {
				break;
			}
		}
		if (tempBuffer[0] == 0xFF && (tempBuffer[1] & 0xF0) == 0xF0 && (tempBuffer[1] & 0x08) == 0x08) {
			int layer = (tempBuffer[1] >> 1) & 0x03;
			int sf = (tempBuffer[2] >> 2) & 0x03;
			long st[3] = {
				44100, 48000, 32000
			};
			int ch = (tempBuffer[3] >> 6) & 0x03;
			if (sf >= 0 && sf < 3) {
				*pInSample = 0;
				pFmt->sample = st[sf];
				pFmt->channel = 2;
				if (ch == 3) {
					pFmt->channel = 1;
				}
				pFmt->bitsPerSample = 16;
				pFmt->fmt[0] = 'm';
				pFmt->fmt[1] = 'p';
				pFmt->fmt[2] = '3';
				pFmt->fmt[3] = '\0';
				fileIsMP3 = TRUE;
				fclose(fp);
				fp = NULL;
			}
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsMP3;
}

int infoFlacFile(char *filename, SOUNDFMT *pFmt, unsigned int *pInSample,FILEINFO *pFileinfo)
/*
 * Flac Fileの情報取得
 */ {
	FILE *fp;
	long rs;
	int fileIsFLAC = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// Flac
		rs = fread(tempBuffer, 1, 10, fp);
		if (rs != 10) {
			break;
		}
		if (tempBuffer[0] == 'f' && tempBuffer[1] == 'L' && tempBuffer[2] == 'a' && tempBuffer[3] == 'C') {
			if (fseek_local(fp, 4, SEEK_SET)) {
				break;
			}
			rs = fread(tempBuffer, 1, 38, fp);
			if (rs != 38) {
				break;
			}
			if ((tempBuffer[0] & 0x7F) == 0x00) {
				int sr;
				int ch;
				int bit;
				sr = tempBuffer[14];
				sr <<= 8;
				sr |= tempBuffer[15];
				sr <<= 4;
				sr |= (tempBuffer[16] >> 4) & 0x0F;
				ch = (tempBuffer[16] >> 1) & 0x07;
				ch++;
				bit = (tempBuffer[16] & 0x01) << 4;
				bit |= (tempBuffer[17] >> 4) & 0x0F;
				bit++;
				*pInSample = 0;
				pFmt->sample = sr;
				pFmt->channel = ch;
				pFmt->bitsPerSample = bit;
				pFmt->fmt[0] = 'f';
				pFmt->fmt[1] = 'l';
				pFmt->fmt[2] = 'a';
				pFmt->fmt[3] = 'c';
				pFmt->fmt[4] = '\0';
				fileIsFLAC = TRUE;
				fclose(fp);
				fp = NULL;
			}
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsFLAC;
}

int infoWavPackFile(char *filename, SOUNDFMT *pFmt, unsigned int *pInSample,
	FILEINFO *pFileinfo)
/*
 * WavPack File の情報取得
 */ {
	FILE *fp;
	long rs;
	int fileIsWV = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// WavPack
		rs = fread(tempBuffer, 1, 10, fp);
		if (rs != 10) {
			break;
		}
		if (tempBuffer[0] == 'w' && tempBuffer[1] == 'v' && tempBuffer[2] == 'p' && tempBuffer[3] == 'k') {
			if (fseek_local(fp, 24, SEEK_SET)) {
				break;
			}
			rs = fread(tempBuffer, 1, 4, fp);
			if (rs == 4) {
				int bit[4] = {
					8, 16, 24, 32
				};
				long st[16] = {6000, 8000, 9600, 11025, 12000, 16000, 22050, 24000, 32000, 44100, 48000, 64000, 88200, 96000, 192000, 0};
				int stIndex;
				int bitIndex;
				*pInSample = 0;
				bitIndex = (tempBuffer[0] & 0x03);
				stIndex = ((tempBuffer[2] >> 7) & 0x01) | ((unsigned short)tempBuffer[3] << 1);
				if (bitIndex != 0 && stIndex >= 8 && stIndex <= 14 && stIndex != 11) {
					pFmt->sample = st[stIndex];
					pFmt->channel = ((tempBuffer[0] >> 2) & 1) == 0 ? 2 : 1;
					pFmt->bitsPerSample = bit[bitIndex];
					pFmt->fmt[0] = 'w';
					pFmt->fmt[1] = 'a';
					pFmt->fmt[2] = 'v';
					pFmt->fmt[3] = 'p';
					pFmt->fmt[4] = '\0';
					fileIsWV = TRUE;
					fclose(fp);
					fp = NULL;
				}
			}
		}
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsWV;
}

int infoDsfFile(char *filename, SOUNDFMT *pFmt, unsigned int *pInSample,
	FILEINFO *pFileinfo)
/*
 * DSF File の情報取得
 */ {
	FILE *fp;
	long rs;
	DSF dsf;
	DSF_FMT dsf_fmt;
	DSF_DATA dsf_data;
	int fileIsDSF = FALSE;

	do {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		// DSF
		rs = fread(&dsf, 1, sizeof(DSF), fp);
		if (rs != sizeof(DSF)) {
			break;
		}
		if (memcmp(dsf.id, "DSD ", 4)) {
			break;
		}
		if (dsf.chunk_size < 28) {
			break;
		}
		if (dsf.file_size < (28 + 52 + 12)) {
			break;
		}

		if (fseek_local(fp, dsf.chunk_size, SEEK_SET)) {
			break;
		}
		rs = fread(&dsf_fmt, 1, sizeof(DSF_FMT), fp);
		if (rs != sizeof(DSF_FMT)) {
			break;
		}
		if (memcmp(dsf_fmt.id, "fmt ", 4)) {
			break;
		}
		if (dsf_fmt.chunk_size < 52) {
			break;
		}
		if (dsf_fmt.fmt_version != 1) {
			break;
		}
		if (dsf_fmt.fmt_id != 0) {
			break;
		}
		if (!(dsf_fmt.channel_type == 1 || dsf_fmt.channel_type == 2)) {
			break;
		}
		if (!(dsf_fmt.channel_count == 1 || dsf_fmt.channel_count == 2)) {
			break;
		}
		if (!(dsf_fmt.sampling == 2822400 || dsf_fmt.sampling == 2822400 * 2)) {
			break;
		}
		if (dsf_fmt.sample_bit_count != 1) {
			break;
		}
		if (dsf_fmt.sample_count < 1) {
			break;
		}
		if (dsf_fmt.block_size != 4096) {
			break;
		}
		if (fseek_local(fp, dsf.chunk_size + dsf_fmt.chunk_size, SEEK_SET)) {
			break;
		}
		rs = fread(&dsf_data, 1, sizeof(DSF_DATA), fp);
		if (rs != sizeof(DSF_DATA)) {
			break;
		}
		if (dsf_data.chunk_size <= 12) {
			break;
		}
		*pInSample = 0;
		pFmt->sample = dsf_fmt.sampling;
		pFmt->channel = dsf_fmt.channel_count;
		pFmt->bitsPerSample = dsf_fmt.sample_bit_count;
		pFmt->fmt[0] = 'd';
		pFmt->fmt[1] = 's';
		pFmt->fmt[2] = 'f';
		pFmt->fmt[3] = '\0';
		fileIsDSF = TRUE;
		fclose(fp);
		fp = NULL;
	}
	while (0);

	if (fp != NULL) {
		fclose(fp);
	}

	return fileIsDSF;

}

/* */
// ---------------------------------------------------------------------------
/*
 * Function    : PLG_InfoAudioData
 * Description : Wave ファイルの情報を得る
 *				 サポートするフォーマットは16,20,24,32(float),64(float) のファイル(それ以外はサポートしない)
 *				 RF64ファイル
 *				 MP3ファイル
 *				 Flacファイル
 *				 Wavpackファイル
 *				 DSDファイル
 * ---
 *	filename   : ファイル名
 *	pFmt	   : 音声ファイル情報を格納する構造体のポインタ
 *	pInSample  : 1Ch あたりのサンプル数を格納する変数のポインタ
 *	pFileinfo  : 曲情報を格納する構造体のポインタ
 *
 */
int PLG_InfoAudioData(char *filename, SOUNDFMT *pFmt, unsigned int *pInSample,
	FILEINFO *pFileinfo)
/*
 * Audio Data の情報取得
 */ {
	int fd;
	struct stat st;
	struct tm *s_time;
	int retValue;

	memset(pFmt, 0, sizeof(SOUNDFMT));
	if (pFileinfo != NULL) {
		memset(pFileinfo, 0, sizeof(FILEINFO));
	}
	fd = open(filename, O_RDONLY);
	if (fd >= 0) {
		if (!fstat(fd, &st)) {
			s_time = localtime(&st.st_mtime);
			sprintf(pFmt->date, "%4d-%02d-%02d", 1900 + s_time->tm_year,s_time->tm_mon, s_time->tm_mday);
			sprintf(pFmt->time, "%02d:%02d:%02d", s_time->tm_hour,s_time->tm_min, s_time->tm_sec);
		}
		close(fd);
	}
	retValue = STATUS_FILE_READ_ERR;

	do {
		retValue = STATUS_UNKNOWN_FORMAT;
		if (infoWavFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoRF64File(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoMP3File(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoFlacFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoWavPackFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}
		if (infoDsfFile(filename, pFmt, pInSample, pFileinfo)) {
			retValue = STATUS_SUCCESS;
			break;
		}

	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_MakeHeaderWAV
 * Description : Wav ファイルのヘッダを作成する
 * ---
 *	pInFmt		: Wave ファイル情報を格納する構造体のポインタ(元)
 *	pOutFmt 	: Wave ファイル情報を格納する構造体のポインタ(先)
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	size		: bufferのサイズ
 *
 */
int PLG_MakeHeaderWAV(SOUNDFMT *pInFmt, SOUNDFMT *pOutFmt, char *buffer,
	long size, long *header_size)
/*
 * Wav ヘッダ作成
 */ {
	int retValue;
	char subTypePCM[16] = {0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,0x00, 0x38, 0x9b, 0x71};
	char subTypeFloat[16] = {0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x80, 0x00, 0x00, 0xaa,0x00, 0x38, 0x9b, 0x71};
	WAVE_HEADER *wh;
	WFE *wfe;
	WFEX *wfex;
	long ws;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (pOutFmt->channel <= 2) {
			if (size <= 20+sizeof(WFE)) {
				break;
			}
		} else {
			if (size <= 20+sizeof(WFEX)) {
				break;
			}
		}
		ws = 0;
		wh = (WAVE_HEADER*)buffer;
		wh->id[0] = 'R';
		wh->id[1] = 'I';
		wh->id[2] = 'F';
		wh->id[3] = 'F';
		wh->size = 0;
		wh->data[0] = 'W';
		wh->data[1] = 'A';
		wh->data[2] = 'V';
		wh->data[3] = 'E';
		//
		wh->data[4] = 'f';
		wh->data[5] = 'm';
		wh->data[6] = 't';
		wh->data[7] = ' ';
		if (pOutFmt->channel <= 2) {
			wh->data[8] = (unsigned char)((sizeof(WFE)	- 2));
			wh->data[9] = (unsigned char)((sizeof(WFE)	- 2) >> 8);
			wh->data[10] = (unsigned char)((sizeof(WFE) - 2) >> 16);
			wh->data[11] = (unsigned char)((sizeof(WFE) - 2) >> 24);
			ws += 20;
			wfe = (WFE*) & buffer[ws];
			memset(wfe, 0, sizeof(WFE));
			/* Set WaveformatEx */
			if (pOutFmt->bitsPerSample == 16 || pOutFmt->bitsPerSample == 24) {
				wfe->wFormatTag = WF_PCM;
			} else if (pOutFmt->bitsPerSample == 32 || pOutFmt->bitsPerSample == 64) {
				wfe->wFormatTag = WF_IEEE_FLOAT;
			}
			wfe->nChannels = pOutFmt->channel;
			wfe->nSamplesPerSec = pOutFmt->sample;
			wfe->wBitsPerSample = pOutFmt->bitsPerSample;
			wfe->nBlockAlign = pOutFmt->channel * pOutFmt->bitsPerSample / 8;
			wfe->nAvgBytesPerSec = pOutFmt->sample * wfe->nBlockAlign;
//			wfe->cbSize = 0;
			ws += (sizeof(WFE) - 2);
		} else {
			wh->data[8] = (unsigned char)(sizeof(WFEX));
			wh->data[9] = (unsigned char)(sizeof(WFEX) >> 8);
			wh->data[10] = (unsigned char)(sizeof(WFEX) >> 16);
			wh->data[11] = (unsigned char)(sizeof(WFEX) >> 24);
			ws += 20;
			wfex = (WFEX*) & buffer[ws];
			memset(wfex, 0, sizeof(WFEX));
			if (pOutFmt->bitsPerSample == 16 || pOutFmt->bitsPerSample == 24 ) {
				memcpy(wfex->subFormat, subTypePCM, 16);
			} else if (pOutFmt->bitsPerSample == 32 || pOutFmt->bitsPerSample == 64) {
				memcpy(wfex->subFormat, subTypeFloat, 16);
			}
			wfex->Format.wFormatTag = 0xFFFE;
			wfex->Format.nChannels = pOutFmt->channel;
			wfex->Format.nSamplesPerSec = pOutFmt->sample;
			wfex->Format.wBitsPerSample = pOutFmt->bitsPerSample;
			wfex->Format.nBlockAlign = wfex->Format.nChannels * pOutFmt->bitsPerSample / 8;
			wfex->Format.nAvgBytesPerSec = pOutFmt->sample * wfex->Format.nBlockAlign;
			wfex->Format.cbSize = 22;
			wfex->Samples.wValidBitsPerSample = pOutFmt->bitsPerSample;

			ws += sizeof(WFEX);
		}
		wh = (WAVE_HEADER*) & buffer[ws];
		wh->id[0] = 'd';
		wh->id[1] = 'a';
		wh->id[2] = 't';
		wh->id[3] = 'a';
		wh->size = 0;
		ws += 8;
		*header_size = ws;
		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_UpdateHeaderWAV
 * Description : Wav ファイルのヘッダを更新する
 * ---
 *	filesize	: ファイルサイズ
 *	datasize	: 音声データサイズ
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	header_size	: bufferのサイズ
 *
 */
int PLG_UpdateHeaderWAV(SOUNDFMT *pOutFmt, int64_t filesize, int64_t datasize,char *buffer, long header_size)
/*
 * Wav ヘッダ更新
 */
{
	int retValue;
	WAVE_HEADER *wh;
	WFE *wfe;
	WFEX *wfex;
	long ws;

	retValue = STATUS_MEM_ALLOC_ERR;
	do {
		ws = 0;
		if (pOutFmt->channel <= 2) {
			if (header_size < 20+(sizeof(WFE) - 2)) {
				break;
			}
		} else {
			if (header_size < 20+sizeof(WFEX)) {
				break;
			}
		}
		wh = (WAVE_HEADER*)buffer;
		wh->size = filesize - 8;
		//
		ws += 12 + 8;
		if (pOutFmt->channel <= 2) {
			ws += (sizeof(WFE) - 2);
		} else {
			ws += sizeof(WFEX);
		}
		wh = (WAVE_HEADER*) & buffer[ws];
		wh->size = datasize;
		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_MakeHeaderRF64
 * Description : RF64 ファイルのヘッダを作成する
 * ---
 *	pInFmt		: Wave ファイル情報を格納する構造体のポインタ(元)
 *	pOutFmt 	: Wave ファイル情報を格納する構造体のポインタ(先)
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	size		: bufferのサイズ
 *
 */
int PLG_MakeHeaderRF64(SOUNDFMT *pInFmt, SOUNDFMT *pOutFmt, char *buffer,long size, long *header_size)
/*
 * Wav ヘッダ作成
 */
{
	int retValue;
	RF64 *rf64;
	RF64_DS64 *ds64;
	RF64_FMT *fmt;
	RF64_DATA *data;
	int64_t ws;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (size < 12+sizeof(RF64_DS64)+sizeof(RF64_FMT) + (sizeof(RF64_DATA) - 1)) {
			break;
		}
		ws = 0;

		rf64 = (RF64*)buffer;
		/* RF64 */
		rf64->chunkId[0] = 'R';
		rf64->chunkId[1] = 'F';
		rf64->chunkId[2] = '6';
		rf64->chunkId[3] = '4';
		rf64->chunkSize = 0xFFFFFFFF;
		rf64->type[0] = 'W';
		rf64->type[1] = 'A';
		rf64->type[2] = 'V';
		rf64->type[3] = 'E';

		/* ds64 */
		ws += sizeof(RF64);
		ds64 = (RF64_DS64*) & buffer[ws];
		memset(ds64, 0, sizeof(RF64_DS64));
		ds64->chunkId[0] = 'd';
		ds64->chunkId[1] = 's';
		ds64->chunkId[2] = '6';
		ds64->chunkId[3] = '4';
		ds64->chunkSize = sizeof(RF64_DS64) - 8;

		/* fmt */
		ws += sizeof(RF64_DS64);
		fmt = (RF64_FMT*) & buffer[ws];
		fmt->chunkId[0] = 'f';
		fmt->chunkId[1] = 'm';
		fmt->chunkId[2] = 't';
		fmt->chunkId[3] = ' ';
		fmt->chunkSize = sizeof(RF64_FMT) - 8;
		if (pOutFmt->bitsPerSample == 16 || pOutFmt->bitsPerSample == 24) {
			fmt->formatType = WF_PCM;
		}
		else if (pOutFmt->bitsPerSample == 32 || pOutFmt->bitsPerSample == 64) {
			fmt->formatType = WF_IEEE_FLOAT;
		}
		fmt->channelCount = pOutFmt->channel;
		fmt->sampleRate = pOutFmt->sample;
		fmt->bitsPerSample = pOutFmt->bitsPerSample;
		fmt->blockAlignment = pOutFmt->channel * pOutFmt->bitsPerSample / 8;
		fmt->bytesPerSecond = pOutFmt->sample * fmt->blockAlignment;

		/* data */
		ws += sizeof(RF64_FMT);
		data = (RF64_DATA*) & buffer[ws];
		data->chunkId[0] = 'd';
		data->chunkId[1] = 'a';
		data->chunkId[2] = 't';
		data->chunkId[3] = 'a';
		data->chunkSize = 0xFFFFFFFF;
		ws += 8;
		*header_size = ws;
		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_UpdateHeaderRF64
 * Description : Wav ファイルのヘッダを更新する
 * ---
 *	filesize	: ファイルサイズ
 *	datasize	: 音声データサイズ
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	header_size	: bufferのサイズ
 *
 */
int PLG_UpdateHeaderRF64(SOUNDFMT *pOutFmt, int64_t filesize, int64_t datasize,char *buffer, long header_size)
/*
 * Wav ヘッダ更新
 */
{
	int retValue;
	RF64_DS64 *ds64;
	unsigned int dwLow, dwHigh;
	int64_t ws;
	int64_t samplecount;

	retValue = STATUS_MEM_ALLOC_ERR;

	do {
		if (header_size < 12+sizeof(RF64_DS64)+sizeof(RF64_FMT) + (sizeof(RF64_DATA) - 1)) {
			break;
		}
		ws = 0;

		/* ds64 */
		ws += sizeof(RF64);
		ds64 = (RF64_DS64*) & buffer[ws];
		dwLow = (unsigned int)(filesize - 8);
		dwHigh = (unsigned int)((filesize - 8) >> 32);
		ds64->riffSizeLow = dwLow;
		ds64->riffSizeHigh = dwHigh;
		dwLow = (unsigned int)datasize;
		dwHigh = (unsigned int)(datasize >> 32);
		ds64->dataSizeLow = dwLow;
		ds64->dataSizeHigh = dwHigh;
		samplecount = datasize;
		samplecount /= (pOutFmt->channel * (pOutFmt->bitsPerSample / 8));
		ds64->sampleCountLow = (unsigned int)(samplecount);
		ds64->sampleCountHigh = (unsigned int)((samplecount >> 32));

		retValue = STATUS_SUCCESS;
	}
	while (0);

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_GetExtraChunk
 * Description : Wav ファイルのfmt,data以外のチャンクを取得する
 * ---
 *	filesize	: ファイルサイズ
 *	datasize	: 音声データサイズ
 *	buffer 		: 作成したヘッダを格納するバッファ
 *	header_size	: bufferのサイズ
 *
 */
int PLG_GetExtraChunk(char *filename, int nChunk, char **infoChunk,long *infoSize)
{
	int retValue;
	char *p1, *p2;
	BYTE header[12];
	BYTE *ovcom;
	long ovcom_count;
	FILE *fp;
	SOUNDFMT inFmt;
	unsigned int inSample;
	FILEINFO fileinfo;
	WAVE_HEADER *wh;
	META_LIST meta_list;
	unsigned int rs;
	int64_t seekPtr;
	unsigned int wk_long;
	int i, next, end;
	int err;
	retValue = STATUS_FILE_READ_ERR;

	if (*infoChunk) {
		free(*infoChunk);
		*infoChunk = NULL;
	}
	*infoSize = 0;
	wh = (WAVE_HEADER*)header;
	ovcom = NULL;
	do {
		memset(&meta_list, 0, sizeof(META_LIST));

		fp = NULL;
		if (infoWavFile(filename, &inFmt, &inSample, &fileinfo) == TRUE) {
			fp = fopen(filename, "rb");
			if (fp == NULL) {
				break;
			}
			// ヘッダ読み込み
			rs = fread(wh, 1, 12, fp);
			if (rs != 12) {
				break;
			}
			seekPtr = 12;
			// 指定チャンクまでスキップ
			next = 0;
			for (i = 0; i < nChunk; ) {
				err = 0;
				if (next) {
					seekPtr += wh->size + 8;
					next = 0;
				}
				if (fseek_local(fp, seekPtr, SEEK_SET)) {
					err = 1;
					break;
				}
				rs = fread(wh, 1, 8, fp);
				if (rs != 8) {
					err = 1;
					break;
				}
				if (!((wh->id[0] == 'f' && wh->id[1] == 'm' && wh->id[2] == 't' && wh->id[3] == ' ') ||
					(wh->id[0] == 'd' && wh->id[1] == 'a' && wh->id[2] == 't' && wh->id[3] == 'a'))) {
					i++;
				} else {
					err = 1;
				}
				next = 1;
			}
			if (err == 0 && wh->size > 0) {
				*infoChunk = (char*)malloc(wh->size + 8);
				if (*infoChunk) {
					if (fseek_local(fp, seekPtr, SEEK_SET) == 0) {
						rs = fread(*infoChunk, 1, wh->size + 8, fp);
						if (rs == wh->size + 8) {
							*infoSize = wh->size + 8;
						}
					}
					if (*infoSize != wh->size + 8) {
						free(*infoChunk);
						*infoSize = 0;
					}
				}
			}
		} else if (infoFlacFile(filename, &inFmt, &inSample, &fileinfo) == TRUE) {
			if (nChunk > 1) {
				break;
			}
			memcpy(meta_list.h_list, "LIST", 4);
			meta_list.s_list = 0;
			memcpy(meta_list.h_info, "INFO", 4);
			memcpy(meta_list.h_isrc, "ISRC", 4);
			meta_list.s_isrc = 0;
			meta_list.p_isrc = NULL;
			memcpy(meta_list.h_icrd, "ICRD", 4);
			meta_list.s_icrd = 0;
			meta_list.p_icrd = NULL;
			memcpy(meta_list.h_iprd, "IPRD", 4);
			meta_list.s_iprd = 0;
			meta_list.p_iprd = NULL;
			memcpy(meta_list.h_ieng, "IENG", 4);
			meta_list.s_ieng = 0;
			meta_list.p_ieng = NULL;
			memcpy(meta_list.h_inam, "INAM", 4);
			meta_list.s_inam = 0;
			meta_list.p_inam = NULL;
			memcpy(meta_list.h_icop, "ICOP", 4);
			meta_list.s_icop = 0;
			meta_list.p_icop = NULL;
			memcpy(meta_list.h_ignr, "IGNR", 4);
			meta_list.s_ignr = 0;
			meta_list.p_ignr = NULL;
			memcpy(meta_list.h_isft, "ISFT", 4);
			meta_list.s_isft = 0;
			meta_list.p_isft = NULL;
			memcpy(meta_list.h_icmt, "ICMT", 4);
			meta_list.s_icmt = 0;
			meta_list.p_icmt = NULL;
			memcpy(meta_list.h_iart, "IART", 4);
			meta_list.s_iart = 0;
			meta_list.p_iart = NULL;

			fp = fopen(filename, "rb");
			if (fp == NULL) {
				break;
			}
			seekPtr = 4;
			for (end = 0; end == 0; ) {
				if (fseek_local(fp, seekPtr, SEEK_SET)) {
					break;
				}
				// ヘッダ読み込み
				rs = fread(header, 1, 4, fp);
				if (rs != 4) {
					break;
				}
				if ((header[0] & 0x80)) {
					end = 1;
				}
				wk_long = header[1] << 16 | header[2] << 8 | header[3];
				seekPtr += wk_long + 4;
				if ((header[0] & 0x7F) == 0x04) {
					end = 1;
					// vendor_length
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					wk_long = header[3] << 24 | header[2] << 16 | header[1]
						<< 8 | header[0];
					if (fseek_local(fp, wk_long, SEEK_CUR)) {
						break;
					}
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					ovcom_count = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					for (i = 0; i < ovcom_count; i++) {
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
						ovcom = (BYTE*)malloc(wk_long + 2);
						if (ovcom) {
							memset(ovcom, 0, wk_long + 2);
							rs = fread(ovcom, 1, wk_long, fp);
							if (wk_long != rs) {
								break;
							}
							p1 = (char*)ovcom;
							p2 = strchr((const char*)ovcom, '=');
							if (p2 != NULL) {
								*p2 = '\0';
								p2++;
							}
							if (strnicmp(p1, "title", 5) == 0 && strlen(p2) > 0) {
								if (meta_list.p_inam != NULL) {
									free(meta_list.p_inam);
									meta_list.p_inam = NULL;
									meta_list.s_inam = 0;
								}
								meta_list.p_inam = (BYTE*)malloc(strlen(p2));
								if (meta_list.p_inam != NULL) {
									meta_list.s_inam = strlen(p2);
									memcpy(meta_list.p_inam, p2,
										meta_list.s_inam);
								}
							} else if (strnicmp(p1, "album", 5) == 0 && strlen(p2) > 0) {
								if (meta_list.p_iprd != NULL) {
									free(meta_list.p_iprd);
									meta_list.p_iprd = NULL;
									meta_list.s_iprd = 0;
								}
								meta_list.p_iprd = (BYTE*)malloc(strlen(p2));
								if (meta_list.p_iprd != NULL) {
									meta_list.s_iprd = strlen(p2);
									memcpy(meta_list.p_iprd, p2,
										meta_list.s_iprd);
								}
							} else if (strnicmp(p1, "artist", 6) == 0 && strlen(p2) > 0) {
								if (meta_list.p_iart != NULL) {
									free(meta_list.p_iart);
									meta_list.p_iart = NULL;
									meta_list.s_iart = 0;
								}
								meta_list.p_iart = (BYTE*)malloc(strlen(p2));
								if (meta_list.p_iart != NULL) {
									meta_list.s_iart = strlen(p2);
									memcpy(meta_list.p_iart, p2,
										meta_list.s_iart);
								}
							} else if (strnicmp(p1,"copyright", 9) == 0 && strlen(p2) > 0) {
								if (meta_list.p_icop != NULL) {
									free(meta_list.p_icop);
									meta_list.p_icop = NULL;
									meta_list.s_icop = 0;
								}
								meta_list.p_icop = (BYTE*)malloc(strlen(p2));
								if (meta_list.p_icop != NULL) {
									meta_list.s_icop = strlen(p2);
									memcpy(meta_list.p_icop, p2,
										meta_list.s_icop);
								}

							} else if (strnicmp(p1,"description", 11) == 0 && strlen(p2) > 0) {
								if (meta_list.p_icmt != NULL) {
									free(meta_list.p_icmt);
									meta_list.p_icmt = NULL;
									meta_list.s_icmt = 0;
								}
								meta_list.p_icmt = (BYTE*)malloc(strlen(p2));
								if (meta_list.p_icmt != NULL) {
									meta_list.s_icmt = strlen(p2);
									memcpy(meta_list.p_icmt, p2,
										meta_list.s_icmt);
								}
							} else if (strnicmp(p1, "date", 4) == 0 && strlen(p2) > 0) {
								if (meta_list.p_icrd != NULL) {
									free(meta_list.p_icrd);
									meta_list.p_icrd = NULL;
									meta_list.s_icrd = 0;
								}
								meta_list.p_icrd = (BYTE*)malloc(strlen(p2));
								if (meta_list.p_icrd != NULL) {
									meta_list.s_icrd = strlen(p2);
									memcpy(meta_list.p_icrd, p2,
										meta_list.s_icrd);
								}
							}
							free(ovcom);
							ovcom = NULL;
						}
					}
					if (i == ovcom_count) {
						// サイズ計算
						wk_long = 0;
						if (meta_list.s_isrc > 0) {
							wk_long += meta_list.s_isrc + 8;
						}
						if (meta_list.s_icrd > 0) {
							wk_long += meta_list.s_icrd + 8;
						}
						if (meta_list.s_iprd > 0) {
							wk_long += meta_list.s_iprd + 8;
						}
						if (meta_list.s_ieng > 0) {
							wk_long += meta_list.s_ieng + 8;
						}
						if (meta_list.s_inam > 0) {
							wk_long += meta_list.s_inam + 8;
						}
						if (meta_list.s_icop > 0) {
							wk_long += meta_list.s_icop + 8;
						}
						if (meta_list.s_ignr > 0) {
							wk_long += meta_list.s_ignr + 8;
						}
						if (meta_list.s_isft > 0) {
							wk_long += meta_list.s_isft + 8;
						}
						if (meta_list.s_icmt > 0) {
							wk_long += meta_list.s_icmt + 8;
						}
						if (meta_list.s_iart > 0) {
							wk_long += meta_list.s_iart + 8;
						}
						meta_list.s_list = wk_long + 4;
						if (meta_list.s_list > 4) {
							*infoChunk = (char*)malloc
								(meta_list.s_list + 8 + 16);
							if (*infoChunk != NULL) {
								p1 = *infoChunk;
								memcpy(p1, meta_list.h_list, 12);
								p1 += 12;
								if (meta_list.s_isrc > 0) {
									memcpy(p1, meta_list.h_isrc, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_isrc,
										meta_list.s_isrc);
									p1 += meta_list.s_isrc;
									if (meta_list.s_isrc & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_icrd > 0) {
									memcpy(p1, meta_list.h_icrd, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_icrd,
										meta_list.s_icrd);
									p1 += meta_list.s_icrd;
									if (meta_list.s_icrd & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_iprd > 0) {
									memcpy(p1, meta_list.h_iprd, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_iprd,
										meta_list.s_iprd);
									p1 += meta_list.s_iprd;
									if (meta_list.s_iprd & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_ieng > 0) {
									memcpy(p1, meta_list.h_ieng, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_ieng,
										meta_list.s_ieng);
									p1 += meta_list.s_ieng;
									if (meta_list.s_ieng & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_inam > 0) {
									memcpy(p1, meta_list.h_inam, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_inam,
										meta_list.s_inam);
									p1 += meta_list.s_inam;
									if (meta_list.s_inam & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_icop > 0) {
									memcpy(p1, meta_list.h_icop, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_icop,
										meta_list.s_icop);
									p1 += meta_list.s_icop;
									if (meta_list.s_icop & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_ignr > 0) {
									memcpy(p1, meta_list.h_ignr, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_ignr,
										meta_list.s_ignr);
									p1 += meta_list.s_ignr;
									if (meta_list.s_ignr & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_isft > 0) {
									memcpy(p1, meta_list.h_isft, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_isft,
										meta_list.s_isft);
									p1 += meta_list.s_isft;
									if (meta_list.s_isft & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_icmt > 0) {
									memcpy(p1, meta_list.h_icmt, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_icmt,
										meta_list.s_icmt);
									p1 += meta_list.s_icmt;
									if (meta_list.s_icmt & 0x01) {
										*p1++ = '\0';
									}
								}
								if (meta_list.s_iart > 0) {
									memcpy(p1, meta_list.h_iart, 8);
									p1 += 8;
									memcpy(p1, meta_list.p_iart,
										meta_list.s_iart);
									p1 += meta_list.s_iart;
									if (meta_list.s_iart & 0x01) {
										*p1++ = '\0';
									}
								}
								*infoSize = (p1 - *infoChunk) - 8;
								meta_list.s_list = *infoSize;
							}
						}
						break;
					}
					break;
				}
			}
		}
	}
	while (0);

	if (ovcom != NULL) {
		free(ovcom);
	}

	if (meta_list.p_isrc != NULL)
		free(meta_list.p_isrc);
	if (meta_list.p_icrd != NULL)
		free(meta_list.p_icrd);
	if (meta_list.p_iprd != NULL)
		free(meta_list.p_iprd);
	if (meta_list.p_ieng != NULL)
		free(meta_list.p_ieng);
	if (meta_list.p_inam != NULL)
		free(meta_list.p_inam);
	if (meta_list.p_icop != NULL)
		free(meta_list.p_icop);
	if (meta_list.p_ignr != NULL)
		free(meta_list.p_ignr);
	if (meta_list.p_isft != NULL)
		free(meta_list.p_isft);
	if (meta_list.p_icmt != NULL)
		free(meta_list.p_icmt);
	if (meta_list.p_iart != NULL)
		free(meta_list.p_iart);

	if (fp != NULL) {
		fclose(fp);
	}

	return retValue;
}

// ---------------------------------------------------------------------------
/*
 * Function    : PLG_GetChunkInfo
 * Description : Wave ファイルのChunk情報を返す
 * ---
 *	filename   : ファイル名
 *	info	   : 情報を格納するバッファ
 *	infoSize   : info のバイト数
 *
 */
int PLG_GetChunkInfo(char *filename, char *info, int infoSize)
/*
 * チャンク情報の取得
 */
{
	FILE *fp;
	WFEX wfex;
	BROADCAST_EXT *bext;
	SOUNDFMT inFmt;
	unsigned int inSample;
	int retValue;
	long chunkSize;
	unsigned int rs;
	unsigned int dw_temp;
	int remInfoSize;
	char work[300];
	char *listData;
	long listSize;
	char *data;
	char *p, *cr;
	long len;
	int i;
	int64_t offset;
	int64_t timeRef;
	int key_writed;
	int64_t seekPtr;
	int end;
	BYTE header[40];
	long wk_long;
	long ovcom_count;
	BYTE *ovcom;
	char *p1, *p2;
	DSF dsf;
	DSF_FMT dsf_fmt;
	DSF_DATA dsf_data;
	RF64_DS64 ds64;

	retValue = STATUS_FILE_READ_ERR;
	data = NULL;
	listData = NULL;
	bext = NULL;
	ovcom = NULL;
	end = 0;
	do {
		remInfoSize = infoSize;
		memset(info, 0, infoSize);

		fp = fopen(filename, "rb");
		if (fp == NULL) {
			break;
		}
		retValue = STATUS_UNKNOWN_FORMAT;

		if (infoWavFile(filename, &inFmt, &inSample, NULL)) {
			// Wave File
			strcpy(info, "AWAVE\n"); // Add "Wave"
			sprintf(work, "IOffset\t-\n");
			strcat(info, work);
			sprintf(work, "ISize\t-\n");
			strcat(info, work);

			if (findWavChunk(fp, "fmt ", &chunkSize)) {
				strcat(info, "D\n"); // Down Tree
				strcat(info, "Afmt\n"); // Add "fmt"
				offset = ftell_local(fp);
				rs = fread(&wfex, 1, chunkSize, fp);
				if (rs != chunkSize) {
					break;
				}
				sprintf(work, "IID\tfmt \n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", rs + 8);
				strcat(info, work);

				sprintf(work, "IwFormatTag\t%04x", wfex.Format.wFormatTag);
				strcat(info, work);
				if (wfex.Format.wFormatTag == WF_PCM) {
					strcat(info, "(WAVE_FORMAT_PCM)\n");
				} else if (wfex.Format.wFormatTag == WF_IEEE_FLOAT) {
					strcat(info, "(WAVE_FORMAT_IEEE_FLOAT)\n");
				} else if (wfex.Format.wFormatTag == 0xFFFE) {
					strcat(info, "(WAVE_FORMAT_EXTENSIBLE)\n");
				} else {
					strcat(info, "\n");
				}
				sprintf(work, "InChannels\t%d\n", wfex.Format.nChannels);
				strcat(info, work);
				sprintf(work, "InSamplesPerSec\t%ld\n",wfex.Format.nSamplesPerSec);
				strcat(info, work);
				sprintf(work, "IwBitsPerSample\t%d\n",wfex.Format.wBitsPerSample);
				strcat(info, work);
				sprintf(work, "InBlockAlign\t%d\n", wfex.Format.nBlockAlign);
				strcat(info, work);
				sprintf(work, "InAvgBytesPerSec\t%ld\n",wfex.Format.nAvgBytesPerSec);
				strcat(info, work);
				if (wfex.Format.wFormatTag == 0xFFFE) {
					sprintf(work, "IwValidBitsPerSample\t%d\n",wfex.Samples.wValidBitsPerSample);
					strcat(info, work);
					sprintf(work, "IwSamplesPerBlock\t%d\n",wfex.Samples.wSamplesPerBlock);
					strcat(info, work);
					sprintf(work, "IdwChannelMask\t%X\n", wfex.dwChannelMask);
					strcat(info, work);
				}
			}
			if (findWavChunk(fp, "data", &chunkSize)) {
				strcat(info, "Adata\n"); // Add "data"
				offset = ftell_local(fp);
				sprintf(work, "IID\tdata\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);

				if (chunkSize > 0L) {
					dw_temp = chunkSize;
					dw_temp /= (wfex.Format.wBitsPerSample / 8);
					dw_temp /= (wfex.Format.nChannels);
					sprintf(work, "ISampleCount\t%ld\n", dw_temp);
					strcat(info, work);
				}
			}

			retValue = STATUS_SUCCESS;

			// List(info)
			if (findWavChunk(fp, "LIST", &chunkSize)) {
				if (chunkSize > 20) {
					offset = ftell_local(fp); listSize = chunkSize;
					listData = malloc(chunkSize);
					if (listData != NULL) {
						rs = fread(listData, 1, listSize, fp);
						if (rs != listSize) {
							free(listData);
							listData = NULL;
						} else {
							if (memcmp(listData, "INFO", 4) != 0) {
								free(listData);
								listData = NULL;
							}
						}
					}
				}
				if (listData != NULL) {
					strcat(info, "ALIST(INFO)\n");
					// Add "INFO"
					sprintf(work, "IID\tLIST(INFO)\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%ld\n", listSize + 8);
					strcat(info, work); p = listData + 4; listSize -= 4;
					while (listSize >= 8) {
						strcat(info, "I");
						memset(work, 0, 5);
						strncpy(work, p, 4);
						strcat(info, work);
						strcat(info, "\t");
						p += 4;
						len = *(long*)p;
						p += 4;
						memset(work, 0, 256);
						strncpy(work, p, len < 256 ? len : 255);
						if (len & 1) {
							len++;
						}
						strcat(info, work);
						strcat(info, "\n");
						p += len;
						if (listSize <= len + 8) {
							break;
						}
						listSize -= (len + 8);
					}
					free(listData);
					listData = NULL;
				}
			}
			if (findWavChunk(fp, "bext", &chunkSize)) {
				// BWF
				strcat(info, "Abext\n"); // Add "bext"
				offset = ftell_local(fp);
				bext = (BROADCAST_EXT*)malloc(chunkSize);
				if (bext != NULL) {
					rs = fread(bext, 1, chunkSize, fp);
					if (rs != chunkSize) {
						free(bext); 
						bext = NULL;
					}
				}
				if (bext != NULL) {
					sprintf(work, "IID\tbext\n");
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset - 8);
					strcat(info, work);
					sprintf(work, "ISize\t%ld\n", chunkSize + 8);
					strcat(info, work); strcat(info, "IDescription\t");
					strcpy(work, "-");
					if (strlen(bext->description) > 0) {
						memset(work, 0, 257);
						strncpy(work, bext->description, 256);
					}
					strcat(info, work);
					strcat(info, "\n");
					strcat(info, "IOriginator\t");
					strcpy(work, "-");
					if (strlen(bext->originator) > 0) {
						memset(work, 0, 257);
						strncpy(work, bext->originator, 32);
					}
					strcat(info, work); strcat(info, "\n");

					strcat(info, "IOriginatorReference\t");
					strcpy(work, "-");
					if (strlen(bext->originatorReference) > 0) {
						memset(work, 0, 257);
						strncpy(work, bext->originatorReference, 32);
					}
					strcat(info, work); 
					strcat(info, "\n");

					strcat(info, "IOriginationDate\t"); 
					strcpy(work, "-");
					if (strlen(bext->originationDate) > 0) {
						memset(work, 0, 257);
						strncpy(work, bext->originationDate, 10);
					}
					strcat(info, work);
					strcat(info, "\n");

					strcat(info, "IOriginationTime\t");
					strcpy(work, "-");
					if (strlen(bext->originationTime) > 0) {
						memset(work, 0, 257);
						strncpy(work, bext->originationTime, 8);
					}
					strcat(info, work); 
					strcat(info, "\n");

					strcat(info, "ITimeReference\t");
					timeRef = bext->timeReferenceHigh; timeRef <<= 32;
					timeRef |= bext->timeReferenceLow;
					timeRef /= wfex.Format.nSamplesPerSec;
					sprintf(work, "%d:%02d:%02d", (int)(timeRef / 3600),(int)((timeRef / 60) % 60), (int)(timeRef % 60));
					strcat(info, work); strcat(info, "\n");

					strcat(info, "IVersion\t");
					sprintf(work, "%d", bext->version);
					strcat(info, work); strcat(info, "\n");

					strcat(info, "IUMID\t"); strcpy(work, "-");
					if (strlen((char*)bext->UMID) > 0) {
						memset(work, 0, 257);
						strncpy(work, (char*)bext->UMID, 64);
					}
					strcat(info, work); strcat(info, "\n");

					strcat(info, "ICodingHistory\t"); strcpy(work, "-");
					key_writed = 1;
					if (strlen(bext->codingHistory) > 0) {
						p = bext->codingHistory;
						while (strlen(p) > 0) {
							if (key_writed == 0) {
								strcat(info, "\nI \t");
							}
							cr = strchr(p, '\n');
							if (cr) {
								*cr = '\0';
							}
							memset(work, 0, 257);
							strncpy(work, p, 256);
							strcat(info, work);
							if (cr) {
								p = (cr + 1);
							}
							key_writed = 0;
						}
					}
					if (key_writed == 1) {
						strcat(info, work);
					}
					strcat(info, "\n"); 
					free(bext); 
					bext = NULL;
				}
			}
			if (findWavChunk(fp, "qlty", &chunkSize)) {
				// BWF(qlty)
				strcat(info, "Aqlty\n"); // Add "qlty"
				offset = ftell_local(fp);
				sprintf(work, "IID\tqlty\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "levl", &chunkSize)) {
				// BWF(qlty)
				strcat(info, "Alevl\n"); // Add "levl"
				offset = ftell_local(fp);
				sprintf(work, "IID\tlevl\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "link", &chunkSize)) {
				// BWF(link)
				strcat(info, "Alink\n"); // Add "link"
				offset = ftell_local(fp);
				sprintf(work, "IID\tlink\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "axml", &chunkSize)) {
				// BWF(axml)
				strcat(info, "Aaxml\n"); // Add "axml"
				offset = ftell_local(fp);
				sprintf(work, "IID\taxml\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "cue ", &chunkSize)) {
				// BWFJ(cue)
				strcat(info, "Acue\n"); // Add "cue"
				offset = ftell_local(fp);
				sprintf(work, "IID\tcue \n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "plst", &chunkSize)) {
				// BWFJ(plst)
				strcat(info, "Aplst\n"); // Add "plst"
				offset = ftell_local(fp);
				sprintf(work, "IID\tplst\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
			if (findWavChunk(fp, "adtl", &chunkSize)) {
				// BWFJ Associated data chunk(adtl)
				strcat(info, "ALIST(adtl)\n"); // Add "adtl"
				offset = ftell_local(fp);
				sprintf(work, "IID\tadtl\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset - 8);
				strcat(info, work);
				sprintf(work, "ISize\t%ld\n", chunkSize + 8);
				strcat(info, work);
			}
		}
		if (infoFlacFile(filename, &inFmt, &inSample, NULL)) {
			strcpy(info, "AFLAC\n"); // Add "Flac"
			sprintf(work, "IOffset\t-\n");
			strcat(info, work);
			sprintf(work, "ISize\t-\n");
			strcat(info, work);
			seekPtr = 4;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(header, 1, 38, fp);
			if (rs != 38) {
				break;
			}
			if ((header[0] & 0x80) != 0x00) {
				end = 1;
			}
			if ((header[0] & 0x7F) == 0x00) {
				//
				long sr;
				long ch;
				long bit;
				unsigned int dw;
				long ld;
				int64_t total;
				strcat(info, "D\n"); // Down Tree
				strcat(info, "ASTREAMINFO\n"); // Add "STREAMINFO"
				offset = 4;
				sprintf(work, "IBlock type\t%d\n", header[0] & 0x7F);
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work);
				ld = header[1];
				ld <<= 8;ld |= header[2];
				ld <<= 8;ld |= header[3];
				sprintf(work, "ISize\t%ld\n", ld + 4); strcat(info, work);

				ld = header[4];
				ld <<= 8;ld |= header[5];

				sprintf(work, "Iminimum block size\t%ld\n", ld);
				strcat(info, work);

				ld = header[6]; ld <<= 8; ld |= header[7];

				sprintf(work, "Imaximum block size\t%ld\n", ld);
				strcat(info, work);

				ld = header[8]; ld <<= 8; ld |= header[9]; ld <<= 8;
				ld |= header[10];

				sprintf(work, "Iminimum frame size\t%ld\n", ld);
				strcat(info, work);

				ld = header[11]; ld <<= 8; ld |= header[12]; ld <<= 8;
				ld |= header[13];

				sprintf(work, "Imaximum frame size\t%ld\n", ld);
				strcat(info, work);

				sr = header[14]; sr <<= 8; sr |= header[15]; sr <<= 4;
				sr |= (header[16] >> 4) & 0x0F;

				sprintf(work, "Isample rate\t%ld\n", sr);
				strcat(info, work);

				ch = (header[16] >> 1) & 0x07; ch++;

				sprintf(work, "Ichannels\t%ld\n", ch); strcat(info, work);

				bit = (header[16] & 0x01) << 4;
				bit |= (header[17] >> 4) & 0x0F; bit++;

				sprintf(work, "Ibits per sample\t%ld\n", bit);
				strcat(info, work);

				total = (header[17] & 0x0F); total <<= 8;
				total |= header[18]; total <<= 8; total |= header[19];
				total <<= 8; total |= header[20];

#ifdef __GUNC__
				sprintf(work, "Itotal samples\t%lld\n", total);
#else
				sprintf(work, "Itotal samples\t%I64d\n", total);
#endif
				strcat(info, work);

				strcat(info, "IMD5\t"); 
				memset(work, 0, 300); 
				p1 = work;
				for (i = 0; i < 16; i++) {
					sprintf(p1, "%02x", header[21 + i]);
					p1 = work + strlen(work);
				}
				strcat(p1, "\n"); 
				strcat(info, work);
				retValue = STATUS_SUCCESS;
			}

			seekPtr = 4;
			for (; end == 0; ) {
				if (fseek_local(fp, seekPtr, SEEK_SET)) {
					break;
				}
				// ヘッダ読み込み
				rs = fread(header, 1, 4, fp);
				if (rs != 4) {
					break;
				}
				if ((header[0] & 0x80)) {
					end = 1;
				}
				wk_long = header[1] << 16 | header[2] << 8 | header[3];
				offset = seekPtr; seekPtr += wk_long + 4;
				if (!((header[0] & 0x7F) == 0x00 || (header[0] & 0x7F) == 0x04)) {
					long ld;
					wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					if ((header[0] & 0x7F) == 1) {
						strcat(info, "APADDING\n"); // Add "PADDING"
					} else if ((header[0] & 0x7F) == 2) {
						strcat(info, "AAPPLICATION\n"); // Add "APPLICATION"
					} else if ((header[0] & 0x7F) == 3) {
						strcat(info, "ASEEKTABLE\n"); // Add "APPLICATION"
					} else if ((header[0] & 0x7F) == 5) {
						strcat(info, "ACUESHEET\n"); // Add "CUESHEET"
					} else if ((header[0] & 0x7F) == 6) {
						strcat(info, "APICTURE\n"); // Add "PICTURE"
					}
					strcat(info, work);
					sprintf(work, "IBlock type\t%d\n", header[0] & 0x7F);
					strcat(info, work);
					sprintf(work, "IOffset\t%08lXh\n", offset);
					strcat(info, work);
					sprintf(work, "ISize\t%ld\n", wk_long);
					strcat(info, work);
				}
				if ((header[0] & 0x7F) == 0x04) {
					// VORBIS_COMMENT
					// end = 1;
					// vendor_length
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					if (fseek_local(fp, wk_long, SEEK_CUR)) {
						break;
					}
					rs = fread(header, 1, 4, fp);
					if (rs != 4) {
						break;
					}
					ovcom_count = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
					if (ovcom_count > 0){
						long ld;
						strcat(info, "AVORBIS_COMMENT\n");
						// Add "VORBIS_COMMENT"
						strcat(info, work);
						sprintf(work, "IBlock type\t%d\n",header[0] & 0x7F); strcat(info, work);
						sprintf(work, "IOffset\t%08lXh\n", offset);
						strcat(info, work);
						sprintf(work, "ISize\t%ld\n", wk_long);
						strcat(info, work);
					}

					for (i = 0; i < ovcom_count; i++) {
						rs = fread(header, 1, 4, fp);
						if (rs != 4) {
							break;
						}
						wk_long = header[3] << 24 | header[2] << 16 | header[1] << 8 | header[0];
						ovcom = (BYTE*)malloc(wk_long + 2);
						if (ovcom) {
							memset(ovcom, 0, wk_long + 2);
							rs = fread(ovcom, 1, wk_long, fp);
							if (wk_long != rs) {
								break;
							}
							p1 = (char*)ovcom;
							p2 = strchr((const char*)ovcom, '=');
							if (p2 != NULL) {
								*p2 = '\0';
								p2++;
							}
							if (p1 != NULL && p2 != NULL && strlen(p1) > 0 && strlen(p2) > 0) {
								sprintf(work, "I%s\t%s\n", p1, p2);
								strcat(info, work);
							}
							free(ovcom); 
							ovcom = NULL;
						}
					}
				}
			}
		}
		if (infoDsfFile(filename, &inFmt, &inSample, NULL)) {
			strcpy(info, "ADSF\n"); // Add "DSF"
			sprintf(work, "IOffset\t-\n"); strcat(info, work);
			sprintf(work, "ISize\t-\n"); strcat(info, work);

			seekPtr = 0;
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			// DSF
			rs = fread(&dsf, 1, sizeof(DSF), fp);
			if (rs != sizeof(DSF)) {
				break;
			}
			strcat(info, "D\n"); // Down Tree
			strcat(info, "AHeader\n"); // Add "header"
			sprintf(work, "IOffset\t%08lXh\n", 0); strcat(info, work);
			sprintf(work, "ISize\t%ld\n", sizeof(DSF)); strcat(info, work);
#ifdef __GUNC__
			sprintf(work, "Ichunk_size\t%lld\n", dsf.chunk_size);
#else
			sprintf(work, "Ichunk_size\t%I64d\n", dsf.chunk_size);
#endif
			strcat(info, work);
#ifdef __GNUC__
			sprintf(work, "Ifile_size\t%lld\n", dsf.file_size);
#else
			sprintf(work, "Ifile_size\t%I64d\n", dsf.file_size);
#endif
			strcat(info, work);
#ifdef __GUNC__
			sprintf(work, "Iptr_meta\t%ld\n", dsf.ptr_meta);
#else
			sprintf(work, "Iptr_meta\t%I64d\n", dsf.ptr_meta);
#endif
			strcat(info, work);

			if (fseek_local(fp, dsf.chunk_size, SEEK_SET)) {
				break;
			}
			rs = fread(&dsf_fmt, 1, sizeof(DSF_FMT), fp);
			if (rs != sizeof(DSF_FMT)) {
				break;
			}
			strcat(info, "Afmt\n"); // Add "fmt"
			strcat(info, work);
			sprintf(work, "IOffset\t%08lXh\n", dsf.chunk_size);
			strcat(info, work);
			sprintf(work, "ISize\t%ld\n", sizeof(DSF_FMT));
			strcat(info, work);

			sprintf(work, "Ifmt_version\t%ld\n", dsf_fmt.fmt_version);
			strcat(info, work);
			sprintf(work, "Ifmt_id\t%ld\n", dsf_fmt.fmt_id);
			strcat(info, work);
			sprintf(work, "Ichannel_type\t%ld\n", dsf_fmt.channel_type);
			strcat(info, work);
			sprintf(work, "Ichannel_count\t%ld\n", dsf_fmt.channel_count);
			strcat(info, work);
			sprintf(work, "Isampling\t%ld\n", dsf_fmt.sampling);
			strcat(info, work);
			sprintf(work, "Isample_bit_count\t%ld\n",dsf_fmt.sample_bit_count); strcat(info, work);

#ifdef __GNUC__
			sprintf(work, "Isample_count\t%lld\n", dsf_fmt.sample_count);
#else
			sprintf(work, "Isample_count\t%I64d\n", dsf_fmt.sample_count);
#endif
			strcat(info, work);

			sprintf(work, "Iblock_size\t%ld\n", dsf_fmt.block_size);
			strcat(info, work);

			if (fseek_local(fp, dsf.chunk_size + dsf_fmt.chunk_size,SEEK_SET)) {
				break;
			}
			rs = fread(&dsf_data, 1, sizeof(DSF_DATA), fp);
			if (rs != sizeof(DSF_DATA)) {
				break;
			}

			strcat(info, "Adata\n"); // Add "data"
			strcat(info, work);
			sprintf(work, "IOffset\t%08lXh\n",dsf.chunk_size + dsf_fmt.chunk_size); strcat(info, work);
			sprintf(work, "ISize\t%ld\n", dsf_data.chunk_size + 12);
			strcat(info, work); retValue = STATUS_SUCCESS;
		}
		if (infoRF64File(filename, &inFmt, &inSample, NULL)) {
			// Wave File
			strcpy(info, "ARF64\n"); // Add "RF64"
			sprintf(work, "IOffset\t-\n"); strcat(info, work);
			sprintf(work, "ISize\t-\n"); strcat(info, work);

			if (findWavChunk(fp, "ds64", &chunkSize)) {
				int64_t qw; 
				strcat(info, "D\n"); // Down Tree
				strcat(info, "Ads64\n"); // Add "fmt"
				offset = ftell_local(fp);
				if (fseek_local(fp,-8,SEEK_CUR)) {
					break;
				}
				rs = fread(&ds64, 1, chunkSize + 8, fp);
				if (rs != chunkSize + 8) {
					break;
				}
				sprintf(work, "IID\tds64\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work); sprintf(work, "ISize\t%ld\n", rs);
				strcat(info, work); qw = ds64.riffSizeHigh; qw <<= 32;
				qw |= ds64.riffSizeLow;
#ifdef __GNUC__
				sprintf(work, "IriffSize\t%lld\n", qw);
#else
				sprintf(work, "IriffSize\t%I64d\n", qw);
#endif
				strcat(info, work);

				qw = ds64.dataSizeHigh; qw <<= 32; qw |= ds64.dataSizeLow;
#ifdef __GNUC__
				sprintf(work, "IdataSize\t%lld\n", qw);
#else
				sprintf(work, "IdataSize\t%I64d\n", qw);
#endif
				strcat(info, work);

				qw = ds64.sampleCountHigh; qw <<= 32;
				qw |= ds64.sampleCountLow;
#ifdef __GNUC__
				sprintf(work, "IsampleCount\t%lld\n", qw);
#else
				sprintf(work, "IsampleCount\t%I64d\n", qw);
#endif
				strcat(info, work);
				sprintf(work, "ItableLength\t%ld\n", ds64.tableLength);
				strcat(info, work);
			}
			if (findWavChunk(fp, "fmt ", &chunkSize)) {
				strcat(info, "Afmt\n"); // Add "fmt"
				offset = ftell_local(fp);
				rs = fread(&wfex, 1, chunkSize, fp);
				if (rs != chunkSize) {
					break;
				}
				sprintf(work, "IID\tfmt \n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work); sprintf(work, "ISize\t%ld\n", rs);
				strcat(info, work);

				sprintf(work, "IwFormatTag\t%04x", wfex.Format.wFormatTag);
				strcat(info, work);
				if (wfex.Format.wFormatTag == WF_PCM) {
					strcat(info, "(WAVE_FORMAT_PCM)\n");
				} else if (wfex.Format.wFormatTag == WF_IEEE_FLOAT) {
					strcat(info, "(WAVE_FORMAT_IEEE_FLOAT)\n");
				} else if (wfex.Format.wFormatTag == 0xFFFE) {
					strcat(info, "(WAVE_FORMAT_EXTENSIBLE)\n");
				} else {
					strcat(info, "\n");
				}
				sprintf(work, "InChannels\t%d\n", wfex.Format.nChannels);
				strcat(info, work);
				sprintf(work, "InSamplesPerSec\t%ld\n",wfex.Format.nSamplesPerSec); strcat(info, work);
				sprintf(work, "IwBitsPerSample\t%d\n",wfex.Format.wBitsPerSample); strcat(info, work);
				sprintf(work, "InBlockAlign\t%d\n",wfex.Format.nBlockAlign); strcat(info, work);
				sprintf(work, "InAvgBytesPerSec\t%ld\n",wfex.Format.nAvgBytesPerSec); strcat(info, work);
				if (wfex.Format.wFormatTag == 0xFFFE) {
					sprintf(work, "IwValidBitsPerSample\t%d\n",wfex.Samples.wValidBitsPerSample);
					strcat(info, work);
					sprintf(work, "IwSamplesPerBlock\t%d\n",wfex.Samples.wSamplesPerBlock);
					strcat(info, work);
					sprintf(work, "IdwChannelMask\t%X\n",wfex.dwChannelMask);
					strcat(info, work);
				}
			}
			if (findWavChunk(fp, "data", &chunkSize)) {
				int64_t qw;
				strcat(info, "Adata\n"); // Add "data"
				offset = ftell_local(fp);
				sprintf(work, "IID\tdata\n");
				strcat(info, work);
				sprintf(work, "IOffset\t%08lXh\n", offset);
				strcat(info, work); qw = ds64.dataSizeHigh; qw <<= 32;
				qw |= ds64.dataSizeLow;
#ifdef __GNUC__
				sprintf(work, "ISize\t%lld\n", qw);
#else
				sprintf(work, "ISize\t%I64d\n", qw);
#endif
				strcat(info, work);

				if (qw > 0L) {
					qw /= (wfex.Format.wBitsPerSample / 8);
					qw /= (wfex.Format.nChannels);
#ifdef __GNUC__
					sprintf(work, "ISampleCount\t%lld\n", qw);
#else
					sprintf(work, "ISampleCount\t%I64d\n", qw);
#endif
					strcat(info, work);
				}
			}

			retValue = STATUS_SUCCESS;
		}
	} while (0);

	if (listData != NULL) {
		free(listData);
	}
	if (bext != NULL) {
		free(bext);
	}
	if (fp != NULL) {
		fclose(fp);
	}

	return retValue; 
}
int findWavChunk(FILE * fp, char*id, long*size)
{
	WAVE_HEADER *wh;
	long rs;
	int64_t seekPtr;
	int found;
	wh = (WAVE_HEADER*)tempBuffer;

	found = 0;
	do {
		// WAV header
		if (fseek_local(fp, 0, SEEK_SET)) {
			break;
		}
		rs = fread(tempBuffer, 1, 12, fp);
		if (rs != 12) {
			break;
		}
		if (!(memcmp(tempBuffer, "RIFF", 4) == 0 || memcmp(tempBuffer,"riff", 4) == 0 || 
			memcmp(tempBuffer, "RF64",4) == 0 || memcmp(tempBuffer, "rf64", 4) == 0)) {
			break;
		}
		if (!(memcmp(tempBuffer + 8, "WAVE",4) == 0 || memcmp(tempBuffer + 8, "wave", 4) == 0)) {
			break;
		}
		seekPtr = 12;
		// 指定チャンクのサーチ
		for (; ; ) {
			if (fseek_local(fp, seekPtr, SEEK_SET)) {
				break;
			}
			rs = fread(wh, 8, 1, fp);
			if (rs != 1) {
				break;
			}
			if (memcmp(wh->id, id, 4) == 0 && wh->size > 0) {
				found = 1;
				*size = wh->size;
				break;
			}
			if (wh->size & 1) {
				// 奇数なら00のパディングがある
				seekPtr += wh->size + 9;
			} else {
				seekPtr += wh->size + 8;
			}
		}
	}
	while (0);
	return found;
}
static int64_t ftell_local(FILE * fp) {
	int64_t pos;
#ifdef __GNUC__
	pos = ftello64(fp);
#else
	pos = _ftelli64(fp);
#endif
	return pos;
}
static int fseek_local(FILE * fp, int64_t offset,int origin) {
	int ret;
#ifdef __GNUC__
	ret = fseeko64(fp, offset, origin);
#else
	ret = _fseeki64(fp, offset, origin);
#endif
	return ret;
}
