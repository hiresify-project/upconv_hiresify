//---------------------------------------------------------------------------
/****************************************************************************/
/* PLG_AUDIO_IO (C) 2008-2012 By 59414d41									*/
/*																			*/
/*																			*/
/****************************************************************************/

#ifndef PLG_AUDIOIO_H
#define PLG_AUDIOIO_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef unsigned char 	BYTE;
typedef unsigned short	WORD;
typedef unsigned int	DWORD;

typedef struct {
	DWORD	sample;
	BYTE	channel;
	BYTE	bitsPerSample;
	BYTE  type;
	char	fmt[16];
	char	date[10];
	char	time[8];
} SOUNDFMT;

typedef struct {
	int64_t	dataOffset;
	int64_t	dataSize;
} FILEINFO;

/* broadcast_audio_extension typedef struct */
#pragma pack(push, 1)
typedef struct broadcast_audio_extension {
	char description[256];			/* Description */
	char originator[32];			/* originator */
	char originatorReference[32];	/* Reference of the originator */
	char originationDate[10];		/* yyyy-mm-dd */
	char originationTime[8];		/* hh-mm-ss */
	DWORD timeReferenceLow;			/* First sample count since midnight low word */
	DWORD timeReferenceHigh;		/* First sample count since midnight, high word */
	WORD version;					/* Version of the BWF; unsigned binary number */
	BYTE UMID[64];					/* Binary byte 0 of SMPTE UMID */
	BYTE reserved[190];				/* 190 bytes, reserved for future use, set to “NULL”*/
	char codingHistory[1];			/* ASCII : History coding */
} BROADCAST_EXT;
#pragma pack(pop)

/* RF64 */
#pragma pack(push, 1)
typedef struct rf64 {
	char chunkId[4];				/* RF64 */
	DWORD chunkSize;				/* 0xFFFFFFFF */
	char type[4];					/* WAVE */
} RF64;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct rf64_ds64 {
	char	chunkId[4];				/* ds64 */
	DWORD	chunkSize;				/* size of the ds64 chunk */
	DWORD	riffSizeLow;			/* size of RF64 block */
	DWORD	riffSizeHigh;			/* size of RF64 block */
	DWORD	dataSizeLow;			/* size of data chunk */
	DWORD	dataSizeHigh;			/* size of data chunk */
	DWORD	sampleCountLow;			/* sample count of fact chunk */
	DWORD	sampleCountHigh;		/* sample count of fact chunk */
	DWORD	tableLength;			/* number of valid entries in table array */
} RF64_DS64;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct rf64_fmt {
	char	chunkId[4];				/* fmt */
	DWORD	chunkSize;				/* size of the fmt chunk */
	WORD	formatType;				/* WAVE_FORMAT_PCM */
	WORD	channelCount;			/* 1 = mono, 2 = stereo */
	DWORD	sampleRate;				/* sampling rate */
	DWORD	bytesPerSecond;			/* data rate */
	WORD	blockAlignment;			
	WORD	bitsPerSample;			/* bits per sample */
} RF64_FMT;
#pragma pack(pop)

#pragma pack(push, 1)
typedef struct rf64_data {
	char	chunkId[4];				/* fmt */
	DWORD	chunkSize;				/* 0xFFFFFFFF */
	char	data[1];				/* data */
} RF64_DATA;
#pragma pack(pop)

/* wave format */
#define WF_PCM			(0x0001)
#define WF_IEEE_FLOAT	(0x0003)
#define WF_EXTENSIBLE	(0xFFFE)

#pragma pack(push, 1)
typedef struct {
	WORD	wFormatTag;
	WORD	nChannels;
	DWORD	nSamplesPerSec;
	DWORD	nAvgBytesPerSec;
	WORD	nBlockAlign;
	WORD	wBitsPerSample;
	WORD	cbSize;
} WFE;

typedef struct {
	WFE		Format;
	union {
		WORD wValidBitsPerSample;
		WORD wSamplesPerBlock;
		WORD wReserved;
	} Samples;
	DWORD	 dwChannelMask;
	char	 subFormat[16];
} WFEX;
#pragma pack(pop)

#define SPEAKER_FRONT_LEFT				(0x1)
#define SPEAKER_FRONT_RIGHT				(0x2)
#define SPEAKER_FRONT_CENTER			(0x4)
#define SPEAKER_LOW_FREQUENCY			(0x8)
#define SPEAKER_BACK_LEFT				(0x10)
#define SPEAKER_BACK_RIGHT				(0x20)

/* Return Status */
#define STATUS_SUCCESS			(0)		/* 正常終了 */
#define STATUS_NOT_IMPLEMENT	(-1)	/* その機能はインプリメントされていない */
#define STATUS_UNKNOWN_FORMAT	(-2)	/* 未知のフォーマット */
#define STATUS_DATA_ERR			(-3)	/* データが壊れている */
#define STATUS_MEM_ALLOC_ERR	(-4)	/* メモリーが確保出来ない */
#define STATUS_FILE_READ_ERR	(-5)	/* ファイルリードエラー */
#define STATUS_FILE_WRITE_ERR	(-6)	/* ファイルライトエラー */
#define STATUS_PLUGIN_ERR		(-15)	/* 内部エラー */

/* Function Prototype */
int PLG_InfoAudioData(char * filename,SOUNDFMT *pFmt,DWORD *pAudioSize,FILEINFO *pFileinfo);
int PLG_MakeHeaderWAV(SOUNDFMT *pInFmt,SOUNDFMT *pOutFmt,char *header,long size,long *header_size);
int PLG_UpdateHeaderWAV(SOUNDFMT *pOutFmt,int64_t filesize,int64_t datasize,char *header,long header_size);
int PLG_MakeHeaderRF64(SOUNDFMT *pInFmt,SOUNDFMT *pOutFmt,char *header,long size,long *header_size);
int PLG_UpdateHeaderRF64(SOUNDFMT *pOutFmt,int64_t filesize,int64_t datasize,char *header,long header_size);
int PLG_GetExtraChunk(char * filename,int nChunk,char **infoChunk,long *infoSize);
int PLG_GetChunkInfo(char * filename,char *info,int infoSize);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
