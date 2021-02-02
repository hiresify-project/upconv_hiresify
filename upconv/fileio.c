/****************************************************************************/
/* fileio (C) 2011-2012 By 59414d41											*/
/* �t�@�C�����o��IO															*/
/*																			*/
/****************************************************************************/

/*--- Log ------------------------------------------------------------------
 * Ver 0.80 <11/07/24> - �V�K�쐬
 *
 * �e�v���O�����ƃ����N���Ďg�p����
 */

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fileio.h"

#define P_BUF_SIZE_N	(3200)

unsigned char fio_ibuffer[10 * 1024 * 1024L];

typedef struct fio_data {
	unsigned char *buf;
	int	r_flag;
	int w_flag;
	int block_no;
	struct fio_data *next;
} FIO_DATA;

void fio_store_buffer(FIO *fio,FIO_DATA *p_data);
void fio_copy_file(FIO *fio_r,FILE *ofp);
void fio_clear_buffer(FIO *fio);
void fio_fill_buffer(FIO *fio);
FIO_DATA *fio_alloc(FIO *fio);
FIO_DATA *fio_free(FIO *fio,FIO_DATA *p_data);
static __int64 ftell_local(FILE *fp);
static int fseek_local(FILE *fp,__int64 offset,int origin);

#if 0
#define	PRINT_LOG(s)	do {																	\
							FILE *log;															\
							log = fopen("d:\\fileio.log","a");									\
							if (log) {															\
								fprintf(log,"%s [%d] %s\n",__FUNCTION__,__LINE__,s);			\
								fclose(log);													\
							}																	\
						} while (0)
#else
#define	PRINT_LOG(s)	//
#endif

//---------------------------------------------------------------------------
// Function   : fio_open
// Description: �t�@�C���I�[�v������
// ---
//	fio			: fileio �\����
//	filename	: �t�@�C����
//	mode		: �t�@�C���I�[�v�����[�h(R,W)
//
void fio_open(FIO *fio,char *filename,int mode)
{
	FIO_DATA *p_data,*p_before;
	int i;
	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	memset(fio,0,sizeof (FIO));

	if (filename == NULL || !(mode == FIO_MODE_R || mode == FIO_MODE_W)) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}

	if (strlen(filename) + 1 > _MAX_PATH) {
		// �t�@�C������������
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	fio->magic = 5915;
	fio->mode  = mode;
	fio->block_size = 10 * 1000 * 1000L;	// 1block 10MB
	fio->data_offset = 0;
	fio->block_count = 0;
	fio->block_max = 5;
	fio->data_maxsize = 0;
	fio->error = 0;

	if (mode == FIO_MODE_R) {
		fio->fp = fopen(filename,"rb");
		if (fio->fp == NULL) {
			fio->error = FIO_ERR_READ;
			return;
		}
		fio_get_filesize(fio,&fio->data_size);
		if (fio->error) {
			PRINT_LOG("ERROR");
			return;
		}
	} else {
		fio->fp = fopen(filename,"wb+");
		if (fio->fp == NULL) {
			fio->error = FIO_ERR_WRITE;
			PRINT_LOG("ERROR");
			return;
		}
		fio->data_size = 0;
	}
	
	strcpy(fio->fname,filename);
	p_before = NULL;
	for (i = 0;i < fio->block_max;i++) {
		p_data = fio_alloc(fio);
		if (fio->error) {
			return;
		}
		p_data->next = p_before;
		p_before = p_data;
		fio->p_data = p_data;
		fio->block_count++;
	}
	return;
}
//---------------------------------------------------------------------------
// Function   : fio_close
// Description: �t�@�C���N���[�Y����
// ---
//	fio			: fileio �\����
//
void fio_close(FIO *fio)
{
	FIO_DATA *p_data;
	FIO_DATA *p_old;
	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}
	
	if (fio->mode == FIO_MODE_W) {
		if (fio->fp == NULL) {
			fio->error = FIO_ERR_WRITE;
			PRINT_LOG("ERROR");
			return;
		}
		// �o�b�t�@�����O����Ă���f�[�^���������݁B
PRINT_LOG("store:all");
		fio_store_buffer(fio,NULL);
		if (fio->error) {
			PRINT_LOG("ERROR");
			return;
		}
	}
	
	// fio data cleanup
	p_data = fio->p_data;
	while (p_data != NULL) {
		p_data = fio_free(fio,p_data);
	}
	fio->p_data = NULL;
	
	if (fio->fp != NULL) {
		fclose(fio->fp);
	}
	fio->fp = NULL;
	fio->error = 0;
	fio->magic = 0;
}
//---------------------------------------------------------------------------
// Function   : fio_get_filesize
// Description: �t�@�C���T�C�Y�擾
// ---
//	fio			: fileio �\����
//	p_size		: �T�C�Y��Ԃ����߂̃A�h���X
//
void fio_get_filesize(FIO *fio,fio_size *p_size)
{
	fio_size	size;
	int error;

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915 || p_size == NULL) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->error) {
		return;
	}

	error = 1;
	if (fseek_local(fio->fp,0,SEEK_END) == 0) {
		size = ftell_local(fio->fp);
		if (fseek_local(fio->fp,0,SEEK_SET) == 0) {
			*p_size = size;
			error = 0;
		}
	}

	if (error) {
		fio->error = FIO_ERR_READ;
		PRINT_LOG("ERROR");
	}
}
//---------------------------------------------------------------------------
// Function   : fio_get_datasize
// Description: �f�[�^�T�C�Y�擾
// ---
//	fio			: fileio �\����
//	p_size		: �T�C�Y��Ԃ����߂̃A�h���X
//
void fio_get_datasize(FIO *fio,fio_size *p_size)
{
	fio_size	size;
	int error;

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915 || p_size == NULL) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->error) {
		return;
	}
	*p_size = fio->data_size;
}
//---------------------------------------------------------------------------
// Function   : fio_setmode_r
// Description: �t�@�C�����������݂���ǂݍ��ݐ�p�ɐݒ�
// ---
//	fio_w		: fileio �\����(W)
//	fio_r		: fileio �\����(R)
//	filename	: �V�����t�@�C����
//
void fio_setmode_r(FIO *fio_w,FIO *fio_r,char *filename)
{
	FILE *fp;
	int i;
	char s[100];
	FIO_DATA *p_data;
	FIO_DATA *p_old;

	if (fio_r == NULL || fio_w == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio_w->magic != 5915) {
		fio_w->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}

	if (fio_r->error) {
		return;
	}

	if (fio_w->error) {
		return;
	}

	if (filename != NULL) {
		if (strlen(filename) + 1 > _MAX_PATH) {
			// �t�@�C������������
			fio_w->error = FIO_ERR_OTHER;
			PRINT_LOG("ERROR");
			return;
		}
	}

	if (fio_w->data_size > 0) {
PRINT_LOG("store:all");
		fio_store_buffer(fio_w,NULL);
		if (fio_w->error) {
			PRINT_LOG("ERROR");
			return;
		}
	}

	fp = NULL;
	if (filename != NULL) {
		if (fio_w->fp == NULL) {
			fio_w->error = FIO_ERR_WRITE;
			PRINT_LOG("ERROR");
			return;
		}
		fclose(fio_w->fp);
		remove(filename);
		rename(fio_w->fname,filename);
#if 0
		fp = fopen(filename,"wb");
		if (fp == NULL) {
			fio_w->error = FIO_ERR_WRITE;
			PRINT_LOG("ERROR");
			return;
		}
		if (fio_w->error) {
			PRINT_LOG("ERROR");
			return;
		}
		fio_copy_file(fio_w,fp);
		if (fio_w->error) {
			PRINT_LOG("ERROR");
			return;
		}
		remove(fio_w->fname);
		fclose(fp);
#endif
		fp = fopen(filename,"rb");
		if (fp == NULL) {
			fio_r->error = FIO_ERR_READ;
			PRINT_LOG("ERROR");
			return;
		}
	} else {
		fclose(fio_w->fp);
		fp = fopen(fio_w->fname,"rb");
		if (fp == NULL) {
			fio_r->error = FIO_ERR_READ;
			PRINT_LOG("ERROR");
			return;
		}
	}

	memset(fio_r,0,sizeof (FIO));
	fio_r->fp = fp;
	fio_r->magic = 5915;
	fio_r->mode = FIO_MODE_R;
	fio_r->error = 0;
	fio_r->p_data = fio_w->p_data;
	fio_r->block_size = fio_w->block_size;
	fio_r->data_size = fio_w->data_size;
	if (fio_w->data_maxsize > 0) {
		fio_r->data_size = fio_w->data_maxsize;
	}
	fio_r->data_maxsize = fio_w->data_maxsize;
	fio_r->data_offset = 0;
	fio_r->block_count = fio_w->block_count;
	fio_r->block_max = 5;
	fio_r->block_count = 5;

	// �u���b�N������
	p_data = fio_r->p_data;
	for (i = 1;p_data != NULL;i++) {
		if (i <= fio_r->block_count) {
			p_data->r_flag = 0;
			p_data->w_flag = 0;
			p_data->block_no = -1;
			if (i == fio_r->block_count) {
				p_old = p_data->next;
				p_data->next = NULL;
				p_data = p_old;
			} else {
				p_data = p_data->next;
			}
		} else {
			p_data = fio_free(fio_r,p_data);
		}
	}
	

	if (filename != NULL) {
		strcpy(fio_r->fname,filename);
	}
	memset(fio_w,0,sizeof (FIO));
}
//---------------------------------------------------------------------------
// Function   : fio_set_maxsize
// Description: �������݃f�[�^�T�C�Y�̍ő�l��ݒ�
// ---
//	fio			: fileio �\����
//	max			: �ő�l
//
void fio_set_maxsize(FIO *fio,fio_size max)
{

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915 || max < 0) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->error) {
		return;
	}
	
	fio->data_maxsize = max;
}
//---------------------------------------------------------------------------
// Function   : fio_set_memory_limit
// Description: �������m�ۂ̍ő�l��ݒ�
// ---
//	fio			: fileio �\����
//	max			: �ő�l
//
void fio_set_memory_limit(FIO *fio,int max)
{
	FIO_DATA *p_data,*p_before;
	int i;
	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915 || max < 0 || max > P_BUF_SIZE_N) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->error) {
		return;
	}

	if (max < 5) {
		p_data = fio->p_data;
		while (p_data != NULL) {
			p_data = fio_free(fio,p_data);
		}
		fio->p_data = NULL;
		fio->block_count = 0;
	}

	fio->block_max = max;

	p_before = fio->p_data;
	for (i = fio->block_count;i < fio->block_max;i++) {
		p_data = fio_alloc(fio);
		if (fio->error) {
			break;
		}
		p_data->next = p_before;
		p_before = p_data;
		fio->p_data = p_data;
		fio->block_count++;
	}
	
}
//---------------------------------------------------------------------------
// Function   : fio_seek
// Description: �t�@�C���|�C���^�ݒ�
// ---
//	fio			: fileio �\����
//	offset		: �I�t�Z�b�g
//	orign		: �I�t�Z�b�g�w��p�����[�^
//
void fio_seek(FIO *fio,fio_size offset,int orign)
{
	int error = 1;
	char s[100];

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}

	if (orign == SEEK_SET) {
		if (offset < fio->data_size && offset >= 0) {
			sprintf(s,"seek:%ld",offset);
			//PRINT_LOG(s);
			fio->data_offset = offset;
			error = 0;
		} else if (offset >= 0) {
			sprintf(s,"seek:%ld",offset);
			//PRINT_LOG(s);
			fio->data_offset = offset;
			fio->data_size	 = offset;
			error = 0;
		}
	}
	if (orign == SEEK_CUR) {
		if ((fio->data_offset + offset) >= 0) {
			fio->data_offset = fio->data_offset + offset;
			sprintf(s,"seek:%ld",fio->data_offset);
			//PRINT_LOG(s);
			if (fio->data_size < fio->data_offset) {
				fio->data_size = fio->data_offset;
			}
			error = 0;
		}
	}
	if (orign == SEEK_END) {
		if (fio->data_size + offset >= 0) {
			fio->data_offset = fio->data_offset + offset;
			sprintf(s,"seek:%ld",fio->data_offset);
			//PRINT_LOG(s);
			if (fio->data_size < fio->data_offset) {
				fio->data_size = fio->data_offset;
			}
			error = 0;
		}
	}
	
	if (error) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
	}
}
//---------------------------------------------------------------------------
// Function   : fio_read
// Description: �t�@�C���ǂݍ���
// ---
//	buf			: �ǂݍ��݃f�[�^�i�[�A�h���X
//	size		: �ǂݍ��݃o�C�g��
//	n			: �ǂݍ��݌�
//	fio			: fileio�\����
//
fio_size fio_read(void *buf,fio_size size,fio_size n,FIO *fio)
{
	fio_size	read_size;
	fio_size	rd;
	fio_size	block;
	fio_size	remain;
	fio_size	tmp_offset;
	fio_size	tmp_size;
	int 		error;
	int			i;
	char		s[100];
	unsigned char *p;
	FIO_DATA *p_data,*p_before;

	memset(buf,0,size * n);

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return 0;
	}

	if (fio->magic != 5915 || buf == NULL || size < 0 || n < 0) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return 0;
	}

	if (fio->error) {
		return 0;
	}
	
	remain = size * n;
	if (remain > fio->data_size - fio->data_offset) {
		remain = fio->data_size - fio->data_offset;
	}
	if (fio->data_maxsize > 0 && remain > fio->data_maxsize - fio->data_offset) {
		remain = fio->data_maxsize - fio->data_offset;
	}

	if (remain == 0) {
		PRINT_LOG("zero");
		return 0;
	}

	read_size = 0;
	for (;remain > 0;) {
		// �ǂݍ��ރf�[�^�̃u���b�N�����߂�
		block = fio->data_offset;
		block /= fio->block_size;

		p_before = NULL;
		p_data = fio->p_data;
		while (p_data != NULL) {
			if (p_data->block_no == block) {
//PRINT_LOG("find block:ok");
				break;
			}
			if (p_data->next == NULL) {
PRINT_LOG("find block:none");
				break;
			}
			p_before = p_data;
			p_data = p_data->next;
		}

		if (p_data == NULL) {
			fio->error = FIO_ERR_OTHER;
			PRINT_LOG("ERROR");
			return 0;
		}

		if (p_data->block_no != block) {
			// �o�b�t�@�ɂȂ�
			if (p_before != NULL) {
				// �V�����f�[�^��擪�Ɏ����Ă���
				p_before->next = p_data->next;
				p_data->next = fio->p_data;
				fio->p_data = p_data;
			}
			if (fio->mode == FIO_MODE_R) {
				// �o�b�t�@���N���A���ĐV���ɓǂݍ��݂��Ȃ����B
				fio_clear_buffer(fio);
				fio_fill_buffer(fio);
			}
			p_data->r_flag = 0;
			p_data->w_flag = 0;
			p_data->block_no = -1;
			error = 1;
			tmp_offset = (fio->data_offset / fio->block_size) * fio->block_size;
			tmp_size   = fio->data_size - tmp_offset;
			if (fio->data_maxsize > 0 && fio->data_size > fio->data_maxsize) {
				tmp_size   = fio->data_maxsize - tmp_offset;
			}
			if (tmp_size > fio->block_size) {
				tmp_size = fio->block_size;
			}
			if (fseek_local(fio->fp,tmp_offset,SEEK_SET) == 0) {
				error = 2;
				rd = fread(p_data->buf,1,tmp_size,fio->fp);
				if (rd == tmp_size) {
					p_data->r_flag = 1;
					p_data->w_flag = 0;
					p_data->block_no = block;
					error = 0;
				}
			}
			if (error) {
				sprintf(s,"data_size:%ld\n",fio->data_size);
				PRINT_LOG(s);
				sprintf(s,"data_maxsize:%ld\n",fio->data_maxsize);
				PRINT_LOG(s);
				sprintf(s,"tmp_offset:%ld\n",tmp_offset);
				PRINT_LOG(s);
				sprintf(s,"tmp_size:%ld\n",tmp_size);
				PRINT_LOG(s);
				sprintf(s,"rd:%ld\n",rd);
				PRINT_LOG(s);
				PRINT_LOG("read zero");
				sprintf(s,"error:%d",error);
				PRINT_LOG(s);
				return 0;
			}
		}

		if (p_data == NULL || p_data->buf == NULL || p_data->block_no != block || p_data->r_flag != 1) {
			fio->error = FIO_ERR_OTHER;
			PRINT_LOG("ERROR");
			return 0;
		}

		tmp_offset = fio->data_offset % fio->block_size;
		tmp_size   = fio->block_size - tmp_offset;
		if (remain < tmp_size) {
			tmp_size = remain;
		}

		memcpy(buf,p_data->buf + tmp_offset,tmp_size);
		remain -= tmp_size;
		buf	+= tmp_size;
		fio->data_offset += tmp_size;
		read_size += tmp_size;
	}

	if (read_size == 0) {
		PRINT_LOG("zero");
	}
	return (read_size / size);
}
//---------------------------------------------------------------------------
// Function   : fio_write
// Description: �t�@�C���ǂݍ���
// ---
//	buf			: �������݃f�[�^�i�[�A�h���X
//	size		: �������݃o�C�g��
//	n			: �������݌�
//	fio			: fileio�\����
//
fio_size fio_write(void *buf,fio_size size,fio_size n,FIO *fio)
{
	fio_size	write_size;
	fio_size	block;
	fio_size	remain;
	fio_size	tmp_offset;
	fio_size	tmp_size;
	fio_size	rd;
	int			alloc_error;
	int			i;
	char		s[50];
	unsigned char *p;
	int log;
	int error;
	FIO_DATA *p_data,*p_before;

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return 0;
	}

	if (fio->magic != 5915 || buf == NULL || size < 0 || n < 0) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return 0;
	}
	
	if (fio->error) {
		return 0;
	}
	
	if (fio->mode != FIO_MODE_W) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return 0;
	}

	log = 1;

	write_size = 0;
	remain = size * n;
	
	if (remain == 0) {
		return 0;
	}

	for (;remain > 0;) {
		block = fio->data_offset;
		block /= fio->block_size;

		p_before = NULL;
		p_data = fio->p_data;
		while (p_data != NULL) {
			if (p_data->block_no == block) {
//PRINT_LOG("find block:ok");
				break;
			}
			if (p_data->next == NULL) {
PRINT_LOG("find block:none");
				break;
			}
			p_before = p_data;
			p_data = p_data->next;
		}

		if (p_data == NULL) {
			fio->error = FIO_ERR_OTHER;
			PRINT_LOG("ERROR");
			return 0;
		}

		if (p_data->block_no != block) {
			if (p_before != NULL) {
				// �V�����f�[�^��擪�Ɏ����Ă���
				p_before->next = p_data->next;
				p_data->next = fio->p_data;
				fio->p_data = p_data;
			}
PRINT_LOG("store:1");
			fio_store_buffer(fio,p_data);
			if (fio->error) {
				PRINT_LOG("ERROR");
				return 0;
			}
			p_data->r_flag = 0;
			p_data->w_flag = 0;
			p_data->block_no = -1;
			error = 0;
			tmp_offset = (fio->data_offset / fio->block_size) * fio->block_size;
			tmp_size   = fio->data_size;
			if (fio->data_maxsize > 0 && fio->data_size > fio->data_maxsize) {
				tmp_size   = fio->data_maxsize;
			}
			if (tmp_offset < tmp_size) {
				tmp_size = tmp_size - tmp_offset;
				if (tmp_size > fio->block_size) {
					tmp_size = fio->block_size;
				}
				if (fseek_local(fio->fp,tmp_offset,SEEK_SET) == 0) {
					error = 2;
					rd = fread(p_data->buf,1,tmp_size,fio->fp);
					if (rd == tmp_size) {
PRINT_LOG("block:read");
						p_data->r_flag = 1;
						p_data->w_flag = 0;
						p_data->block_no = block;
						error = 0;
					}
				}
				if (error) {
					sprintf(s,"data_size:%ld\n",fio->data_size);
					PRINT_LOG(s);
					sprintf(s,"data_maxsize:%ld\n",fio->data_maxsize);
					PRINT_LOG(s);
					sprintf(s,"tmp_offset:%ld\n",tmp_offset);
					PRINT_LOG(s);
					sprintf(s,"tmp_size:%ld\n",tmp_size);
					PRINT_LOG(s);
					sprintf(s,"rd:%ld\n",rd);
					PRINT_LOG(s);
					PRINT_LOG("read zero");
					sprintf(s,"error:%d",error);
					PRINT_LOG(s);
					return 0;
				}
			} else {
				// �ő�̏������݃T�C�Y�𒴂��Ă���ꍇ�͂����Ő���I������B
				if (fio->data_maxsize > 0 && fio->data_offset > fio->data_maxsize) {
					return n;
				}
			}
		}
		tmp_offset = fio->data_offset % fio->block_size;
		tmp_size   = fio->block_size - tmp_offset;
		if (remain < tmp_size) {
			tmp_size = remain;
		}
		memcpy(p_data->buf + tmp_offset,buf,tmp_size);
		p_data->w_flag = 1;
		p_data->r_flag = 1;
		p_data->block_no = block;
		remain -= tmp_size;
		buf	+= tmp_size;
		fio->data_offset += tmp_size;
		if (fio->data_offset > fio->data_size) {
			fio->data_size	 = fio->data_offset;
		}
		write_size += tmp_size;
	}
	return (write_size / size);
}
//---------------------------------------------------------------------------
// Function   : fio_flush
// Description: �t�@�C���������݃t���b�V��
// ---
//	fio			: fileio�\����
//
void fio_flush(FIO *fio)
{

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}

	if (fio->mode == FIO_MODE_W) {
PRINT_LOG("store:all");
		fio_store_buffer(fio,NULL);
	}
}
//---------------------------------------------------------------------------
// Function   : fio_rewind
// Description: �t�@�C���|�C���^�[�����Z�b�g
// ---
//	fio			: fileio�\����
//
void fio_rewind(FIO *fio)
{

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}

	fio->data_offset = 0;
}
//---------------------------------------------------------------------------
// Function   : fio_tell
// Description: �t�@�C���|�C���^�擾
// ---
//	fio			: fileio�\����
//
fio_size fio_tell(FIO *fio)
{

	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}

	return fio->data_offset;
}
//---------------------------------------------------------------------------
// Function   : fio_store_buffer
// Description: �t�@�C���o�b�t�@��������
// ---
//	fio			: fileio�\����
//	p_data		: �������݃f�[�^
//
void fio_store_buffer(FIO *fio,FIO_DATA *p_data)
{
	fio_size remain;
	fio_size data_size;
	fio_size offset;
	fio_size r_off;
	fio_size wr_n;
	FIO_DATA **fio_all = NULL;
	FIO_DATA *wk_data,*wk_tmp;
	int i,j,c;
	int all;
	char s[50];
	
	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->magic != 5915 || fio->mode != FIO_MODE_W) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->fp == NULL) {
		fio->error = FIO_ERR_WRITE;
		PRINT_LOG("ERROR");
		return;
	}

	if (fio->error) {
		return;
	}

	all = 0;
	if (p_data == NULL) {
		all = 1;
		fio_all = (FIO_DATA **)malloc(fio->block_count * sizeof (FIO_DATA *));
		if (fio_all != NULL) {
PRINT_LOG("fio_all");
			for (i = 0;i < fio->block_count;i++) {
				fio_all[i] = NULL;
			}
			for (i = 0,wk_data = fio->p_data;i < fio->block_count && wk_data != NULL;i++,wk_data = wk_data->next) {
				if (wk_data->block_no != -1 && wk_data->w_flag == 1) {
					for (j = 0;j < fio->block_count;j++) {
						if (fio_all[j] == NULL) {
							break;
						}
					}
sprintf(s,"j:%d,block_no:%d",j,wk_data->block_no);
PRINT_LOG(s);
					fio_all[j] = wk_data;
					c = j;
					if (c > 0) {
						for (--j;j >= 0;j--,c--) {
							if (fio_all[j]->block_no > fio_all[c]->block_no) {
sprintf(s,"j:%d  block_no:%d,c:%d block_no:%d",j,fio_all[j]->block_no,c,fio_all[c]->block_no);
PRINT_LOG(s);
								wk_tmp = fio_all[j];
								fio_all[j] = fio_all[c];
								fio_all[c] = wk_tmp;
							} else {
								break;
							}
						}
					}
				}
			}
		}
		if (fio_all) {
PRINT_LOG("c = 0");
			c = 0;
			p_data = fio_all[c];
		}
	}

	if (fio->data_size > 0) {
		while (p_data != NULL) {
sprintf(s,"block_no:%d",p_data->block_no);
PRINT_LOG(s);
			if (p_data->block_no != -1 && p_data->w_flag == 1) {
				offset = p_data->block_no * (fio_size)fio->block_size;
				data_size = fio->data_size;
				if (fio->data_maxsize > 0 && data_size > fio->data_maxsize) {
					data_size = fio->data_maxsize;
				}
				if (offset < data_size) {
					remain = data_size - offset;
					wr_n = 0;
					fflush(fio->fp);
					if (fseek_local(fio->fp,offset,SEEK_SET) == 0) {
						r_off = ftell_local(fio->fp);
						if (offset != r_off) {
							char s[100];
							sprintf(s,"ERROR fseeko64 bad:%lld,%lld",offset,r_off);
							PRINT_LOG(s);
							fio->error = FIO_ERR_WRITE;
							return;
						}
						if (remain >= fio->block_size) {
							wr_n = fwrite(p_data->buf,fio->block_size,1,fio->fp);
						} else {
							wr_n = fwrite(p_data->buf,remain,1,fio->fp);
						}
					} else {
						r_off = ftell_local(fio->fp);
						sprintf(s,"SEEKERR offset:%lld,r_off:%lld,data_size:%lld",offset,r_off,data_size);
						PRINT_LOG(s);
						fio->error = FIO_ERR_WRITE;
						return;
					}
					if (wr_n == 0) {
						char s[100];
						sprintf(s,"i:%d,wr_n:%lld,remain:%lld,block_size:%ld",i,wr_n,remain,fio->block_size);
						fio->error = FIO_ERR_WRITE;
						PRINT_LOG(s);
						return;
					}
				}
				p_data->w_flag = 0;
			}
			if (all == 0) {
				break;
			} else {
				if (fio_all) {
					if (c + 1 < fio->block_count) {
						c++;
						p_data = fio_all[c];
					} else {
						break;
					}
				} else {
					p_data = p_data->next;
				}
			}
		}
	}
	if (fio_all) {
		free(fio_all);
	}
	PRINT_LOG("fio_store_buffer end");
}
//---------------------------------------------------------------------------
// Function   : fio_copy_file
// Description: �t�@�C���R�s�[
// ---
//	fio			: fileio�\����
//
void fio_copy_file(FIO *fio_r,FILE *ofp)
{
	fio_size remain;
	size_t	 n;
	char	s[100];

	remain = fio_r->data_size;
	if (fio_r->data_maxsize > 0 && fio_r->data_size > fio_r->data_maxsize) {
		remain = fio_r->data_maxsize;
	}
	
	fseek_local(fio_r->fp,0,SEEK_SET);

	while (remain) {
		if (remain >= 10 * 1024 * 1024L) {
			n = fread(fio_ibuffer,10 * 1024 * 1024L,1,fio_r->fp);
		} else {
			n = fread(fio_ibuffer,remain,1,fio_r->fp);
		}
		if (n != 1) {
			fio_r->error = FIO_ERR_READ;
			return;
		}

		if (remain >= 10 * 1024 * 1024L) {
			sprintf(s,"copy_file:%ld\n",1 * 1024 * 1024L);
		} else {
			sprintf(s,"copy_file:%ld\n",remain);
		}

		if (remain >= 10 * 1024 * 1024L) {
			n = fwrite(fio_ibuffer,10 * 1024 * 1024L,1,ofp);
			remain -= 10 * 1024 * 1024L;
		} else {
			n = fwrite(fio_ibuffer,remain,1,ofp);
			remain = 0;
		}
		if (n != 1) {
			fio_r->error = FIO_ERR_WRITE;
			return;
		}
	}
}
//---------------------------------------------------------------------------
// Function   : fio_clear_buffer
// Description: �o�b�t�@�N���A
// ---
//	fio			: fileio�\����
//
void fio_clear_buffer(FIO *fio)
{
	FIO_DATA *p_data;
	FIO_DATA *p_old;
	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}
	
	if (fio->mode == FIO_MODE_W) {
		return;
	}
	// fio data cleanup
	p_data = fio->p_data;
	while (p_data != NULL) {
		p_data->r_flag = 0;
		p_data->w_flag = 0;
		p_data->block_no = -1;
		p_data = p_data->next;
	}
}
//---------------------------------------------------------------------------
// Function   : fio_fill_buffer
// Description: �o�b�t�@��ǂ�
// ---
//	fio			: fileio�\����
//
void fio_fill_buffer(FIO *fio)
{
	FIO_DATA *p_data;
	fio_size offset;
	fio_size remain;
	int block;
	long rs;
	int error;

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return;
	}

	// check magic no
	if (fio->magic != 5915) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return;
	}
	
	if (fio->error) {
		return;
	}
	
	if (fio->mode == FIO_MODE_W) {
		return;
	}
	
	offset = fio->data_offset;
	offset = (offset / fio->block_size) * fio->block_size;
	remain = fio->data_size - offset;
	if (remain == 0) {
		return;
	}

	p_data = fio->p_data;
	while (p_data != NULL) {
		remain = fio->data_size - offset;
		if (fio->data_maxsize > 0 && fio->data_size > fio->data_maxsize) {
			if (offset >= fio->data_maxsize) {
				return;
			}
			remain = fio->data_maxsize - offset;
		}
		if (remain > fio->block_size) {
			remain = fio->block_size;
		}
		error = 1;
		if (fseek_local(fio->fp,offset,SEEK_SET) == 0) {
			rs = fread(p_data->buf,1,remain,fio->fp);
			if (rs == remain) {
				p_data->r_flag = 1;
				p_data->w_flag = 0;
				p_data->block_no = offset / fio->block_size;
				error = 0;
			}
		}
		if (error) {
			fio->error = FIO_ERR_READ;
			return;
		}
		p_data = p_data->next;
		offset += fio->block_size;
		if (offset >= fio->data_size) {
			return;
		}
	}
}
//---------------------------------------------------------------------------
// Function   : fio_alloc
// Description: fio_data �������m��
// ---
//	fio_data  : fio_data�\����
//
FIO_DATA *fio_alloc(FIO *fio)
{
	FIO_DATA *p_data;

	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return NULL;
	}

	p_data = malloc(sizeof (FIO_DATA));
	if (p_data == NULL) {
		fio->error = FIO_ERR_MEMORY;
		PRINT_LOG("ERROR");
		return NULL;
	}
	memset(p_data,0,sizeof (FIO_DATA));
	p_data->r_flag = 0;
	p_data->w_flag = 0;
	p_data->block_no = -1;
	p_data->next = NULL;
	p_data->buf = malloc(fio->block_size);
	if (p_data->buf == NULL) {
		free(p_data);
		fio->error = FIO_ERR_MEMORY;
		PRINT_LOG("ERROR");
		return NULL;
	}
	memset(p_data->buf,0,fio->block_size);
	return p_data;
}

//---------------------------------------------------------------------------
// Function   : fio_free
// Description: fio_data �������J��
// ---
//	fio_data  : fio_data�\����
//
FIO_DATA *fio_free(FIO *fio,FIO_DATA *p_data)
{
	FIO_DATA *next;

	next = NULL;
	if (fio == NULL) {
		PRINT_LOG("ERROR");
		return NULL;
	}
	if (p_data == NULL) {
		fio->error = FIO_ERR_OTHER;
		PRINT_LOG("ERROR");
		return NULL;
	}

	next = p_data->next;
	if (p_data->buf != NULL) {
		free(p_data->buf);
		p_data->buf = NULL;
	}
	free(p_data);
	return next;
}
static __int64 ftell_local(FILE *fp)
{
	__int64 pos;
#ifdef __GNUC__
	pos = ftello64(fp);
#else
	pos = _ftelli64(fp);
#endif
	return pos;
}
static int fseek_local(FILE *fp,__int64 offset,int origin)
{
	int ret;
#ifdef __GNUC__
	ret = fseeko64(fp,offset,origin);
#else
	ret = _fseeki64(fp,offset,origin);
#endif
	return ret;
}
