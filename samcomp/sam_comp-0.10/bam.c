/*
 * Handles reading (only for now) from a BAM format file.
 * We use zlib to read the data in and fetch data one read at a time.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <fcntl.h>
#include <zlib.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#include "bam.h"
#include "os.h"

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/* Macros to store integers of various sizes in little endian byte order.
 * The value is put in the location pointed to by ucp, which should be 
 * an unsigned char pointer.  ucp is incremented by the size of the
 * stored value. */

#define STORE_UINT16(ucp, val)			\
    *(ucp)++ = ((uint16_t) val)      & 0xff;	\
    *(ucp)++ = ((uint16_t) val >> 8) & 0xff;

#define STORE_UINT32(ucp, val)			\
    *(ucp)++ = ((uint32_t) (val))       & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >>  8) & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >> 16) & 0xff;	\
    *(ucp)++ = ((uint32_t) (val) >> 24) & 0xff;

#define STORE_UINT64(ucp, val)			\
    *(ucp)++ = ((uint64_t) (val))       & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >>  8) & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >> 16) & 0xff;	\
    *(ucp)++ = ((uint64_t) (val) >> 24) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 32) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 40) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 48) & 0xff; \
    *(ucp)++ = ((uint64_t) (val) >> 56) & 0xff;

static int bam_more_input(bam_file_t *b);
static int bam_more_output(bam_file_t *b);

/*
 * Allow for unaligned memory access. This is used in cigar string processing
 * as the packed BAM data struct has cigar after read name instead of before.
 */
#define ALLOW_UAC

/*
 * Reads len bytes from fp into data.
 *
 * Returns the number of bytes read.
 *         0 for eof.
 *        -1 for failure.
 */
static int bam_read(bam_file_t *b, void *data, size_t len) {
    int nb = 0, n;
    unsigned char *cdata = (unsigned char *)data;

    while (len) {
	/* Consume any available uncompressed output */
	if (b->out_sz) {
	    size_t l = MIN(b->out_sz, len);
	    memcpy(cdata, b->out_p, l);
	    b->out_p += l;
	    b->out_sz -= l;
	    cdata += l;
	    len -= l;
	    nb += l;

	    if (!len)
		return nb;
	}

	if (!b->gzip) {
	    /* Already uncompressed, so easy to deal with */
	    if (!b->in_sz)
		if (-1 == bam_more_input(b))
		    return nb ? nb : 0;
		    
	    b->out_p  = b->in_p;
	    b->out_sz = b->in_sz;
	    b->in_sz  = 0;
	    continue;
	}

	n = bam_more_output(b);
	if (n == -1)
	    return -1;
	if (n == 0)
	    return nb;
    }

    return nb;
}

/*
 * Reads a line of text of unknown length.
 * 'str' is both input and output. If *str == NULL then memory is allocated
 * for the line. If *str != NULL then it is expected to point to an existing
 * block of memory that we can write into and realloc as required.
 *
 * Similarly *len is both input and output. It is expected to hold the
 * current allocated size of *str. It is modified if we realloc it.
 *
 * Lines have the \n removed and will be null terminated.
 *
 * Returns actual line length used (note note the same as *len) on success
 *        -1 on failure
 */
int bam_get_line(bam_file_t *b, unsigned char **str, size_t *len) {
    unsigned char *buf = *str;
    size_t used_l = 0;
    size_t alloc_l = *len;
    int next_condition;

    while (b->out_sz || bam_more_output(b) > 0) {
	int tmp;
	unsigned char *from = b->out_p;
	unsigned char *to   = &buf[used_l];

	/*
	 * Next condition is the number of loop iterations before something
	 * has to be done - either getting more uncompressed output or
	 * resizing the buffer. We don't care which, but it allows us to
	 * have just one check per loop instead of two. Once out of the loop
	 * we can then afford to determine which case is and deal with it.
	 */
	tmp = next_condition = MIN(b->out_sz, alloc_l-used_l);

	while (next_condition-- > 0) { /* these 3 lines are 50% of SAM cpu */
	    if (*from != '\n') {
		*to++ = *from++;
	    } else {
		b->out_p = from;
		used_l = to-buf;
		b->out_p++;
		buf[used_l] = 0;
		*str = buf;
		*len = alloc_l;
		b->out_sz -= (tmp - next_condition);
		return used_l;
	    }
	}

	used_l = to-buf;
	b->out_p = from;
	b->out_sz -= tmp;

	if (used_l >= alloc_l) {
	    alloc_l = alloc_l ? alloc_l * 2 : 1024;
	    if (NULL == (buf = (unsigned char *)realloc(buf, alloc_l)))
		return -1;
	}
    }

    return b->out_sz ? -1 : 0;
}

int load_bam_header(bam_file_t *b) {
    char magic[4];
    uint32_t i;

    if (4 != bam_read(b, magic, 4))
	return -1;
    if (memcmp(magic, "BAM\x01",4) != 0)
	return -1;
    if (4 != bam_read(b, &b->header_len, 4))
	return -1;
    b->header_len = le_int4(b->header_len);
    if (!(b->header = (char *)malloc(b->header_len)))
	return -1;
    if (b->header_len != (uint32_t)bam_read(b, b->header, b->header_len))
	return -1;
    if (4 != bam_read(b, &b->nref, 4))
	return -1;
    b->nref = le_int4(b->nref);

    b->ref = (bam_ref_t *)malloc(b->nref * sizeof(*b->ref));
    for (i = 0; i < b->nref; i++) {
	uint32_t nlen;
	if (4 != bam_read(b, &nlen, 4))
	    return -1;
	nlen = le_int4(nlen);

	b->ref[i].name = (char *)calloc(1, nlen);
	if (nlen != (uint32_t)bam_read(b, b->ref[i].name, nlen))
	    return -1;
	if (4 != bam_read(b, &b->ref[i].len, 4))
	    return -1;
	b->ref[i].len = le_int4(b->ref[i].len);
    }

    for (i = 0; i < b->header_len; i++) {
	if (b->header[i] == '\n')
	    b->line++;
    }

    return 0;
}

int load_sam_header(bam_file_t *b) {
    unsigned char *str = NULL;
    size_t alloc = 0;
    int header_pos = 0, len;

    b->header = NULL;
    b->header_len = 0;
    
    b->ref  = NULL;
    b->nref = 0;

    while ((b->out_sz > 0 || bam_more_output(b) > 0) && *b->out_p == '@') {
	b->line++;
	if ((len = bam_get_line(b, &str, &alloc)) == -1)
	    return -1;

	/* Add header lines to b->header */
	while (header_pos + len + 1>= (int)b->header_len) {
	    b->header_len += 8192;
	    b->header = (char *)realloc(b->header, b->header_len);
	    if (!b->header)
		return -1;
	}
	memcpy(&b->header[header_pos], str, len);
	b->header[header_pos+len] = '\n';
	header_pos += len+1;

	/* Also parse any reference lines while reading the header */
	if (str[1] == 'S' && str[2] == 'Q') {
	    int rlen = -1;
	    char *rname = NULL;
	    unsigned char *cp = str+3;
	    HashData hd;

	    //printf("line=%ld/%ld/'%s'\n", (long)len, (long)strlen(str), str);

	    if (*cp != '\t') {
		fprintf(stderr, "Malformed header line '%s'\n", str);
		return -1;
	    }
	    cp++;

	    /* Tokenise line into key/value pairs */
	    while (*cp) {
		char *key = (char *)cp;
		char *val;

		if (!key[0] || !key[1] || key[2] != ':') {
		    fprintf(stderr, "Malformed header line '%s'\n", str);
		    return -1;
		}
		key[2] = '\0';

		cp = (unsigned char *)(val = key+3);
		while (*cp && *cp != '\t')
		    cp++;

		if (*cp)
		    *cp++ = 0;

		if (0 == strcmp(key, "LN")) {
		    rlen = atoi(val);
		} else if (0 == strcmp(key, "SN")) {
		    rname = val;
		}
	    }

	    if (!rname || rlen == -1) {
		fprintf(stderr, "No SN or LN value in @SQ line\n");
		return -1;
	    }

	    b->ref = (bam_ref_t *)realloc(b->ref, (b->nref+1)*sizeof(*b->ref));
	    if (!b->ref)
		return -1;
	    b->ref[b->nref].len  = rlen;
	    b->ref[b->nref].name = strdup(rname);

	    hd.i = b->nref;
	    HashTableAdd(b->ref_hash, b->ref[b->nref].name, 0, hd, NULL);

	    b->nref++;
	}
    }

    /* Blank header if none supplied */
    if (!b->header)
	b->header = strdup("");

    if (header_pos)
	b->header[header_pos] = '\0';
    b->header_len = header_pos;

    return 0;
}

/*
 * Extracts @RG records from the header and places them in a hash table.
 * Returns 0 on success
 *        -1 on failure
 */
static int parse_header(bam_file_t *b) {
    int lno = 0;
    int ntabs;
    tag_list_t *tags;
    size_t i, j;

    if (!b->header)
	return -1;

    b->rg_hash = HashTableCreate(4, HASH_FUNC_HSIEH |
				    HASH_DYNAMIC_SIZE |
				    HASH_NONVOLATILE_KEYS);
    if (!b->rg_hash)
	return -1;

    for (i = 0; i < b->header_len; i++) {
	char *id = NULL;
	int id_len;
	HashData hd;
	int i_start = i, i_len;

	for (j = i; b->header[j] && j < b->header_len; j++)
	    if (b->header[j] == '\n')
		break;
	i_len = j - i_start;
	lno++;

	if (b->header[i] != '@') {
	    fprintf(stderr, "Header line does not start with '@' at line %d:\n"
		    "\"%.*s\"\n",
		    lno, i_len, &b->header[i_start]);
            i = i_start + i_len;
            continue;
	}

	if (b->header[i+1] != 'R' && b->header[i+2] != 'G') {
	    i = i_start + i_len;
	    continue;
	}

	/* Tokenise, not fast but doesn't matter */
	for (ntabs = 0, j = i; j < b->header_len && b->header[j] != '\n'; j++)
	    if (b->header[j] == '\t')
		ntabs++;

	if (ntabs == 0) {
	    fprintf(stderr, "Missing tab in header line %d:\n\"%.*s\"\n",
		    lno, i_len, &b->header[i_start]);
            i = i_start + i_len;
            continue;
	}

	if (NULL == (tags = (tag_list_t *)malloc((ntabs+1) * sizeof(*tags))))
	    return -1;
	
	while (i < b->header_len && b->header[i] != '\n') {
	    if (b->header[i] == '\t')
		break;
	    i++;
	}

	ntabs = 0;
	while (i < b->header_len && b->header[i] != '\n') {
	    if (b->header[i] == '\t') {
		int key;
		char *value;

		if (!b->header[i+1] || !b->header[i+2] ||
		    b->header[i+3] != ':') {
		    fprintf(stderr,
			    "Malformed key:value pair in header line %d:\n"
			    "\"%.*s\"\n",
			    lno, i_len, &b->header[i_start]);
		    i = i_start + i_len;
		    continue;
		}

		key = (b->header[i+1]<<8) + b->header[i+2]; 
		i+=4;
		
		value = &b->header[i];
		while (i < b->header_len && b->header[i] != '\n') {
		    if (b->header[i] == '\t')
			break;
		    i++;
		}

		if (key == (('I'<<8) | 'D')) {
		    id = value;
		    id_len = &b->header[i] - value;
		} else {
		    tags[ntabs].value  = value;
		    tags[ntabs].key    = key;
		    tags[ntabs].length = &b->header[i] - value;
		    ntabs++;
		}
	    }
	}
	tags[ntabs].value = NULL;

	if (id) {
	    hd.p = tags;
	    HashTableAdd(b->rg_hash, id, id_len, hd, NULL);
	} else {
	    fprintf(stderr, "No ID record in @RG line\n");
	}

	if (0) {
	    printf("RG ID=%.*s\n", id_len, id);
	    tag_list_t *t = tags;
	    while (t->value) {
		printf("%c%c -> '%.*s'\n",
		       t->key >> 8, t->key & 0xff,
		       t->length,
		       t->value);
		t++;
	    }
	}
    }

    return 0;
}

tag_list_t *bam_find_rg(bam_file_t *b, char *id) {
    HashItem *hi = HashTableSearch(b->rg_hash, id, 0);
    return hi ? (tag_list_t *)hi->data.p : NULL;
}

#ifndef O_BINARY
#    define O_BINARY 0
#endif

bam_file_t *bam_open(const char *fn, const char *mode) {
    bam_file_t *b = (bam_file_t *)calloc(1, sizeof *b);
    
    if (*mode != 'r')
	return NULL;

    if (!b)
	return NULL;

    if (strcmp(mode, "rb") == 0) {
	if (-1 == (b->fd = open(fn, O_RDONLY | O_BINARY, 0)))
	    goto error;
    } else {
	if (-1 == (b->fd = open(fn, O_RDONLY, 0)))
	    goto error;
    }
    b->in_p     = b->in;
    b->in_sz    = 0;
    b->out_p    = b->out;
    b->out_sz   = 0;
    b->next_len = -1;
    b->bs       = NULL;
    b->bs_size  = 0;
    b->z_finish = 1;
    b->bgzf     = 0;
    b->no_aux   = 0;
    b->line     = 0;

    b->ref_hash = HashTableCreate(16, HASH_FUNC_HSIEH |
				      HASH_DYNAMIC_SIZE |
				      HASH_NONVOLATILE_KEYS);

    /* Load first block so we can check */
    bam_more_input(b);
    if (b->in_sz < 2)
	return NULL;
    if (b->in_p[0] == 31 && b->in_p[1] == 139)
	b->gzip = 1;
    else
	b->gzip = 0;

    if (b->gzip) {
	/* Set up zlib */
	b->s.zalloc    = NULL;
	b->s.zfree     = NULL;
	b->s.opaque    = NULL;
	inflateInit2(&b->s, -15);
    }

    /* Load header */
    if (strcmp(mode, "rb") == 0) {
	if (-1 == load_bam_header(b))
	    goto error;
	b->bam = 1;
    } else {
	if (-1 == load_sam_header(b))
	    goto error;
	b->bam = 0;
    }

    /* Parse RG records in header */
    if (-1 == parse_header(b))
	return NULL;

    return b;

 error:
    if (b) {
	if (b->header)
	    free(b->header);
	free(b);
    }

    return NULL;
}

void bam_close(bam_file_t *b) {
    if (!b)
	return;

    if (b->bs)
	free(b->bs);
    if (b->header)
	free(b->header);
    if (b->ref) {
	size_t i;
	for (i = 0; i < b->nref; i++)
	    free(b->ref[i].name);
	free(b->ref);
    }

    if (b->gzip)
	inflateEnd(&b->s);

    if (b->ref_hash)
	HashTableDestroy(b->ref_hash, 0);

    if (b->rg_hash)
	HashTableDestroy(b->rg_hash, 1);

    free(b);
}

/*
 * Loads more data into the input (compressed) buffer.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
static int bam_more_input(bam_file_t *b) {
    size_t l;

    if (b->in != b->in_p) {
	memmove(b->in, b->in_p, b->in_sz);
	b->in_p = b->in;
    }

    l = read(b->fd, &b->in[b->in_sz], Z_BUFF_SIZE - b->in_sz);
    if (l <= 0)
	return -1;
    
    b->in_sz += l;
    return 0;
}

/*
 * Converts compressed input to the uncompressed output buffer
 *
 * Returns number of additional output bytes on success
 *         0 on eof
 *        -1 on failure.
 */
static int bam_more_output(bam_file_t *b) {
    int err = Z_OK;
    unsigned char *bgzf;
    size_t xlen;
    size_t bsize;

    assert(b->out_sz == 0);

    if (!b->gzip) {
	/* Already uncompressed, so easy to deal with */
	if (!b->in_sz)
	    if (-1 == bam_more_input(b))
		return 0;
		    
	b->out_p  = b->in_p;
	b->out_sz = b->in_sz;
	b->in_sz  = 0;
	return b->out_sz;
    }

    /* Uncompress another BGZF block */
    /* BGZF header */
    if (b->in_sz < 18) {
	if (-1 == bam_more_input(b))
	    return 0;
	if (b->in_sz < 18)
	    return -1;
	if (b->in_sz == 0)
	    return 0; /* eof */
    }
	
    if (b->z_finish) {
	/*
	 * BGZF header is gzip + extra fields.
	 */
	bgzf = b->in_p;
	b->in_p += 10; b->in_sz -= 10;

	if (bgzf[0] != 31 || bgzf[1] != 139)
	    return -1; /* magic number failure */
	if ((bgzf[3] & 4) == 4) {
	    /* has extra fields, eg BGZF */
	    xlen = bgzf[10] + bgzf[11]*256;
	    b->in_p += 2; b->in_sz -= 2;
	} else {
	    xlen = 0;
	}
    } else {
	/* Continuing with an existing data stream */
	xlen = 0;
    }

    /* BGZF */
    if (xlen == 6) {
	b->bgzf = 1;
	b->in_p += 6; b->in_sz -= 6;

	if (bgzf[12] != 'B' || bgzf[13] != 'C' ||
	    bgzf[14] !=  2  || bgzf[15] !=  0)
	    return -1;
	bsize = bgzf[16] + bgzf[17]*256;
	bsize -= 6+19;

	/* Inflate */
	if (b->in_sz < bsize + 8) {
	    do {
		if (bam_more_input(b) == -1)
		    return -1; /* Truncated */
	    } while (b->in_sz < bsize + 8);
	}
	b->s.avail_in  = bsize;
	b->s.next_in   = b->in_p;
	b->s.avail_out = Z_BUFF_SIZE;
	b->s.next_out  = b->out;
	b->s.total_out = 0;
	    
	inflateReset(&b->s);
	err = inflate(&b->s, Z_FINISH);
	if (err != Z_STREAM_END) {
	    fprintf(stderr, "Inflate returned error code %d\n", err);
	    return -1;
	}
	b->z_finish = 1;

	b->in_p   += bsize + 8; /* crc & isize */
	b->in_sz  -= bsize + 8;
	b->out_sz  = b->s.total_out;
	b->out_p   = b->out;
    } else {
	/* Some other gzip variant, but possibly still having xlen */
	while (xlen) {
	    int d = MIN(b->in_sz, xlen);
	    xlen     -= d;
	    b->in_p  += d;
	    b->in_sz -= d;
	    if (b->in_sz == 0)
		bam_more_input(b);
	    if (b->in_sz == 0)
		return -1; /* truncated file */
	}
	    
	b->s.avail_in  = b->in_sz;
	b->s.next_in   = b->in_p;
	b->s.avail_out = Z_BUFF_SIZE;
	b->s.next_out  = b->out;
	b->s.total_out = 0;
	if (b->z_finish)
	    inflateReset(&b->s);
	    
	do {
	    err = inflate(&b->s, Z_BLOCK);
	    //printf("err=%d\n", err);

	    if (err == Z_OK || err == Z_STREAM_END) {
		b->in_p  += b->in_sz - b->s.avail_in;
		b->in_sz  = b->s.avail_in;
		b->out_sz = b->s.total_out;
		b->out_p  = b->out;
	    }

	    if (err == Z_STREAM_END) {
		b->z_finish = 1;
		/* Consume (ignore) CRC & ISIZE */
		if (b->in_sz < 8)
		    bam_more_input(b);

		if (b->in_sz < 8)
		    return -1; /* truncated file */

		b->in_sz -= 8;
		b->in_p  += 8;
	    } else {
		b->z_finish = 0;
	    }
	} while (err != Z_OK && err != Z_STREAM_END);
    }

    return b->out_sz;
}

/*
 * Decodes the next line of SAM into a bam_seq_t struct.
 */
int sam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    static unsigned char *str = NULL;
    static size_t alloc_l = 0;
    int used_l, n, sign;
    unsigned char *cpf, *cpt, *cp;
    int cigar_len;
    bam_seq_t *bs;
    HashItem *hi;
    static int lookup[256] = {
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 00 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 20 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 30 */
	0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0, /* 40 */
	0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0, /* 50 */
	0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0, /* 60 */
	0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0, /* 70 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 80 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 90 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* a0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* b0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* c0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* d0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* e0 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};/* f0 */

    /* Fetch a single line */
    if ((used_l = bam_get_line(b, &str, &alloc_l)) <= 0) {
	return used_l;
    }

    used_l *= 4; // FIXME

    /* Over sized memory, for worst case? FIXME: cigar can break this! */
    if (!*bsp || used_l + sizeof(*bs) > (*bsp)->alloc) {
	if (!(*bsp = (bam_seq_t *)realloc(*bsp, used_l + sizeof(*bs))))
	    return -1;
	(*bsp)->alloc = used_l + sizeof(*bs);
	(*bsp)->blk_size = 0; /* compute later */
    }

    bs = *bsp;
    bs->flag_nc = 0;
    bs->bin_mq_nl = 0;
    
    /* Decode line */
    cpf = str;
    cpt = (unsigned char *)&bs->data;
    
    /* Name */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	*cpt++ = *cpf++;
    *cpt++ = 0;
    if (!*cpf++) return -1;
    bam_set_name_len(bs, cpf-cp);

    /* flag */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bam_set_flag(bs, n);

    /* ref */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	cpf++;
    if (*cp == '*') {
	/* Unmapped */
	bs->ref = -1;
    } else {
	hi = HashTableSearch(b->ref_hash, (char *)cp, cpf-cp);
	if (!hi) {
	    HashData hd;

	    //fprintf(stderr, "Reference seq %.*s unknown\n", (int)(cpf-cp), cp);

	    /* Fabricate it instead */
	    b->ref = (bam_ref_t *)realloc(b->ref, (b->nref+1)*sizeof(*b->ref));
	    if (!b->ref)
		return -1;
	    b->ref[b->nref].len  = 0; /* Unknown value */
	    b->ref[b->nref].name = (char *)malloc(cpf-cp+1);
	    memcpy(b->ref[b->nref].name, cp, cpf-cp);
	    b->ref[b->nref].name[cpf-cp] = 0;

	    hd.i = b->nref;
	    hi = HashTableAdd(b->ref_hash, b->ref[b->nref].name, 0, hd, NULL);
	    b->nref++;
	}
	bs->ref = hi->data.i;
    }
    cpf++;

    /* Pos */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->pos = n-1;

    /* map qual */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bam_set_map_qual(bs, n);

    /* cigar */
    n = 0;
    cigar_len = 0;
    if (*cpf == '*') {
	cpf++;
    } else {
	//while (*cpf && *cpf != '\t') {
	while (*cpf > '\t') {
	    if (isdigit(*cpf)) {
		n = n*10 + *cpf++ - '0';
	    } else {
		unsigned char op;
		union {
		    unsigned char c[4];
		    uint32_t i;
		} c4i;

		switch (*cpf++) {
		case 'M': op=BAM_CMATCH;         break;
		case 'I': op=BAM_CINS;           break;
		case 'D': op=BAM_CDEL;           break;
		case 'N': op=BAM_CREF_SKIP;      break;
		case 'S': op=BAM_CSOFT_CLIP;     break;
		case 'H': op=BAM_CHARD_CLIP;     break;
		case 'P': op=BAM_CPAD;           break;
		case '=': op=BAM_CBASE_MATCH;    break;
		case 'X': op=BAM_CBASE_MISMATCH; break;
		default:
		    fprintf(stderr, "Unknown cigar opcode '%c'\n", cpf[-1]);
		    return -1;
		}

		c4i.i = (n << 4) | op;
		*cpt++ = c4i.c[0];
		*cpt++ = c4i.c[1];
		*cpt++ = c4i.c[2];
		*cpt++ = c4i.c[3];

		n = 0;
		cigar_len++;
	    }
	}
    }
    bam_set_cigar_len(bs, cigar_len);
    if (!*cpf++) return -1;
    
    /* mate ref name */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	cpf++;
    if (!*cpf++) return -1;
    bs->mate_ref = 0; // FIXME

    /* mate pos */
    n = 0;
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->mate_pos = n-1;

    /* insert size */
    n = 0;
    if (*cpf == '-') {
	sign = -1;
	cpf++;
    } else {
	sign = 1;
    }
    //while (*cpf && *cpf != '\t')
    while (*cpf > '\t')
	n = n*10 + *cpf++-'0';
    if (!*cpf++) return -1;
    bs->ins_size = n*sign;

    /* seq */
    cp = cpf;
    n = 0;
    //while (*cpf && *cpf != '\t') {
    if (cpf[0] == '*' && cpf[1] == '\t') {
	cpf++;
	bs->len = 0;
    } else {
	while (*cpf > '\t') {
	    if (n == 0) {
		*cpt = lookup[*cpf]<<4;
		n = 1;
	    } else {
		n = 0;
		*cpt++ |= lookup[*cpf];
	    }
	    cpf++;
	}
	if (n == 1)
	    cpt++;
	bs->len = cpf-cp;
    }
    if (!*cpf++) return -1;

    /* qual */
    cp = cpf;
    //while (*cpf && *cpf != '\t')
    if (cpf[0] == '*' && (cpf[1] == '\0' || cpf[1] == '\t')) {
	/* no qual */
	memset(cpt, '\0', bs->len);
	cpt += bs->len;
	cpf++;
    } else {
	while (*cpf > '\t')
	    *cpt++ = *cpf++ - '!';
    }

    assert((char *)cpt == (char *)(bam_aux(bs)));

    if (!*cpf++ || b->no_aux) goto skip_aux;

    /* aux */
    while (*cpf) {
	unsigned char *key = cpf, *value;
	if (!(key[0] && key[1] && key[2] == ':' && key[3] && key[4] == ':'))
	    return -1;
	cpf += 5;

	value = cpf;
	while (*cpf && *cpf != '\t')
	    cpf++;

	*cpt++ = key[0];
	*cpt++ = key[1];

	switch(key[3]) {
	case 'A':
	    *cpt++ = 'A';
	    *cpt++ = *value;
	    break;

	case 'i':
	    n = atoi((char *)value);
	    if (n >= 0) {
		if (n < 256) {
		    *cpt++ = 'C';
		    *cpt++ = n;
		} else if (n < 65536) {
		    *cpt++ = 'S';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		} else {
		    *cpt++ = 'I';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		    *cpt++ = (n >>16) & 0xff;
		    *cpt++ = (n >>24) & 0xff;
		}
	    } else {
		if (n >= -128 && n < 128) {
		    *cpt++ = 'c';
		    *cpt++ = n;
		} else if (n >= -32768 && n < 32768) {
		    *cpt++ = 's';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		} else {
		    *cpt++ = 'i';
		    *cpt++ = n & 0xff;
		    *cpt++ = (n >> 8) & 0xff;
		    *cpt++ = (n >>16) & 0xff;
		    *cpt++ = (n >>24) & 0xff;
		}
	    }
	    break;

	case 'f': {
	    union {
		float f;
		int i;
	    } u;
	    u.f = atof((char *)value);
	    *cpt++ = 'f';
	    *cpt++ = (u.i) & 0xff;
	    *cpt++ = (u.i >> 8) & 0xff;
	    *cpt++ = (u.i >>16) & 0xff;
	    *cpt++ = (u.i >>24) & 0xff;
	    break;
	}

	case 'Z':
	    *cpt++ = 'Z';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'H':
	    *cpt++ = 'Z';
	    while (value != cpf)
		*cpt++=*value++;
	    *cpt++ = 0;
	    break;

	case 'B': {
	    char subtype = *value++;
	    unsigned char *sz;
	    int count = 0;

	    *cpt++ = 'B';
	    *cpt++ = subtype;
	    sz = cpt; cpt += 4; /* Fill out later */

	    while (*value == ',') {
		value++;
		switch (subtype) {
		case 'c': case 'C':
		    *cpt++ = strtol((char *)value, (char **)&value, 10);
		    break;

		case 's': case 'S':
		    n = strtol((char *)value, (char **)&value, 10);
		    STORE_UINT16(cpt, n);
		    break;
		    
		case 'i': case 'I':
		    n = strtoll((char *)value, (char **)&value, 10);
		    STORE_UINT32(cpt, n);
		    break;
		    
		case 'f': {
		    union {
			float f;
			int i;
		    } u;

		    u.f = strtod((char *)value, (char **)&value);
		    STORE_UINT32(cpt, u.i);
		    break;
		}
		}
		count++;
	    }
	    if (value != cpf) {
		fprintf(stderr, "Malformed %c%c:B:... auxiliary field\n",
			key[0], key[1]);
		value = cpf;
	    }
	    STORE_UINT32(sz, count);
	    break;
	}

	default:
	    fprintf(stderr, "Unknown aux format code '%c'\n", key[3]);
	    break;
	}

	if (*cpf == '\t')
	    cpf++;
    }

 skip_aux:
    *cpt++ = 0;

    bs->blk_size = (unsigned char *)cpt - (unsigned char *)&bs->ref;
    if (bs->blk_size >= bs->alloc)
	abort();

    return 1;
}

/*
 * Fills out the next bam_seq_t struct.
 * bs must be non-null, but *bs may be NULL or an existing bam_seq_t pointer.
 * This function will alloc and/or grow the memory accordingly, allowing for
 * efficient reuse.
 *
 * Returns 1 on success
 *         0 on eof
 *        -1 on error
 */
#ifdef ALLOW_UAC
int bam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    int32_t blk_size, blk_ret;
    bam_seq_t *bs;

    b->line++;

    if (!b->bam)
	return sam_next_seq(b, bsp);

    if (b->next_len > 0) {
	blk_size = b->next_len;
    } else {
	if (4 != bam_read(b, &blk_size, 4))
	    return 0;
	blk_size = le_int4(blk_size);
    }

    if (!*bsp || blk_size > (*bsp)->blk_size) {
	/* 20 extra is for bs->alloc to bs->cigar_len plus next_len */
	if (!(*bsp = (bam_seq_t *)realloc(*bsp, blk_size+20)))
	    return -1;
	(*bsp)->alloc = blk_size+20;
	(*bsp)->blk_size = blk_size;
    }
    bs = *bsp;
    
    if ((blk_ret = bam_read(b, &bs->ref, blk_size+4)) == 0)
	return 0;

    if (blk_size+4 != blk_ret) {
	if (blk_size != blk_ret) {
	    return -1;
	} else {
	    b->next_len = 0;
	    ((char *)(&bs->ref))[blk_size] = 0;
	}
    } else {
	memcpy(&b->next_len, &((char *)(&bs->ref))[blk_size], 4);
	((char *)(&bs->ref))[blk_size] = 0;
    }
    b->next_len = le_int4(b->next_len);

    bs->blk_size  = blk_size;
    bs->ref       = le_int4(bs->ref);
    bs->pos       = le_int4(bs->pos);
    bs->bin_mq_nl = le_int4(bs->bin_mq_nl);
    bs->flag_nc   = le_int4(bs->flag_nc);
    bs->len       = le_int4(bs->len);
    bs->mate_ref  = le_int4(bs->mate_ref);
    bs->mate_pos  = le_int4(bs->mate_pos);
    bs->ins_size  = le_int4(bs->ins_size);

    /* Unpack flag_nc into separate flag & cigar_len */
    bs->cigar_len = bs->flag_nc & 0xffff;
    
    if (10 == be_int4(10)) {
	int i, cigar_len = bam_cigar_len(bs);
	uint32_t *cigar = bam_cigar(bs);
	for (i = 0; i < cigar_len; i++) {
	    cigar[i] = le_int4(cigar[i]);
	}
    }

    return 1;
}

#else

int bam_next_seq(bam_file_t *b, bam_seq_t **bsp) {
    int32_t blk_size, blk_ret;
    bam_seq_t *bs;

    if (!b->bam)
	return sam_next_seq(b, bsp);

    if (b->next_len > 0) {
	blk_size = b->next_len;
    } else {
	if (4 != bam_read(b, &blk_size, 4))
	    return 0;
	blk_size = le_int4(blk_size);
    }

    if (!*bsp || blk_size > (*bsp)->blk_size) {
	if (!(*bsp = realloc(*bsp, blk_size+24)))
	    return -1;
	(*bsp)->alloc = blk_size+24;
	(*bsp)->blk_size = blk_size;
    }
    bs = *bsp;

    /* The fixed-sized fields */
    if ((blk_ret = bam_read(b, &bs->ref, 32)) == 0)
	return 0;

    if (blk_ret != 32)
	return -1;

    //bs->blk_size  = blk_size;
    bs->ref       = le_int4(bs->ref);
    bs->pos       = le_int4(bs->pos);
    bs->bin_mq_nl = le_int4(bs->bin_mq_nl);
    bs->flag_nc   = le_int4(bs->flag_nc);
    bs->len       = le_int4(bs->len);
    bs->mate_ref  = le_int4(bs->mate_ref);
    bs->mate_pos  = le_int4(bs->mate_pos);
    bs->ins_size  = le_int4(bs->ins_size);

    /* Unpack flag_nc into separate flag & cigar_len */
    bs->cigar_len = bs->flag_nc & 0xffff;

    /* Name */
    if (bam_read(b, &bs->data, bam_name_len(bs)) != bam_name_len(bs))
	return -1;

    /* Pad name out to end on a word-aligned boundary */
    blk_ret = blk_size - 32 - bam_name_len(bs);
    bam_set_name_len(bs, round4(bam_name_len(bs)));

    /* The remainder, word aligned */
    if (bam_read(b, &bs->data + bam_name_len(bs), blk_ret) != blk_ret)
	return -1;

    if (10 == be_int4(10)) {
	int i, cigar_len = bam_cigar_len(bs);
	uint32_t *cigar = bam_cigar(bs);
	for (i = 0; i < cigar_len; i++) {
	    bs[i] = le_int4(bs[i]);
	}
    }

    return 1;
}
#endif

/*
 * Looks for aux field 'key' and returns the value.
 * Returns NULL if not found.
 */
char *bam_aux_find(bam_seq_t *b, char *key) {
    char *cp = bam_aux(b);
    static int type_size[256] = {
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  4,  0,  4,  0,  0,  0,  0,  0,  7,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  4, 11,  0,  7,  0,  0,  7,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

    while (*cp) {
	int sz;

	//printf("%c%c:%c:?\n", cp[0], cp[1], cp[2]);

	if (cp[0] == key[0] && cp[1] == key[1])
	    return cp+3;
	
	if ((sz = type_size[(unsigned char)cp[2]])) {
	    /* Fixed length fields */
	    cp += sz;
	} else {
	    /* Variable length, null terminated */
	    if (cp[2] == 'Z' || cp[2] == 'H') {
		cp += 3;
		while (*cp++)
		    ;
	    } else {
		fprintf(stderr, "Unknown aux type code %c(%d) in seq %s\n",
			cp[2], cp[2], bam_name(b));
		return NULL;
	    }
	}
    }

    return NULL;
}


/*
 * An iterator on bam aux fields. NB: This code is not reentrant or multi-
 * thread capable. The values returned are valid until the next call to
 * this function.
 * key:  points to an array of 2 characters (eg "RG", "NM")
 * type: points to an address of 1 character (eg 'Z', 'i')
 * val:  points to an address of a bam_aux_t union.
 *
 * Pass in *iter_handle as NULL to initialise the search and then
 * pass in the modified value on each subsequent call to continue the search.
 *
 * Returns 0 if the next value is valid, setting key, type and val.
 *        -1 when no more found.
 */
int bam_aux_iter(bam_seq_t *b, char **iter_handle,
		 char *key, char *type, bam_aux_t *val) {
    char *s;

    if (!iter_handle || !*iter_handle) {
	s = (char *)bam_aux(b);
    } else {
	s = *iter_handle;
    }

    /* We null terminate our aux list for ease */
    if (s[0] == 0)
	return -1;

    key[0] = s[0];
    key[1] = s[1];
    
    switch (s[2]) {
    case 'A':
	if (type) *type = 'A';
	if (val) val->i = *(s+3);
	s+=4;
	break;

    case 'C':
	if (type) *type = 'i';
	if (val) val->i = *(uint8_t *)(s+3);
	s+=4;
	break;

    case 'c':
	if (type) *type = 'i';
	if (val) val->i = *(int8_t *)(s+3);
	s+=4;
	break;

    case 'S':
	if (type) *type = 'i';
	if (val) {
	    uint16_t tmp;
	    ((char *)&tmp)[0] = s[3];
	    ((char *)&tmp)[1] = s[4];
	    val->i = tmp;
	}
	s+=5;
	break;

    case 's':
	if (type) *type = 'i';
	if (val) {
	    int16_t tmp;
	    ((char *)&tmp)[0] = s[3];
	    ((char *)&tmp)[1] = s[4];
	    val->i = tmp;
	}
	s+=5;
	break;

    case 'I':
	if (type) *type = 'i';
	if (val) {
	    uint32_t tmp;
	    ((char *)&tmp)[0] = s[3];
	    ((char *)&tmp)[1] = s[4];
	    ((char *)&tmp)[2] = s[5];
	    ((char *)&tmp)[3] = s[6];
	    val->i = tmp;
	}
	s+=7;
	break;

    case 'i':
	if (type) *type = 'i';
	if (val) {
	    int32_t tmp;
	    ((char *)&tmp)[0] = s[3];
	    ((char *)&tmp)[1] = s[4];
	    ((char *)&tmp)[2] = s[5];
	    ((char *)&tmp)[3] = s[6];
	    val->i = tmp;
	}
	s+=7;
	break;

    case 'f':
	if (type) *type = 'f';
	if (val) memcpy(&val->f, s+3, 4);
	s+=7;
	break;

    case 'd':
	if (type) *type = 'd';
	if (val) memcpy(&val->d, s+3, 8);
	s+=11;
	break;

    case 'Z': case 'H':
	if (type) *type = s[2];
	s+=3;
	if (val) val->s = s;
	while (*s++);
	break;

    default:
	fprintf(stderr, "Unknown aux type '%c'\n", s[2]);
	return -1;
    }

    if (iter_handle)
	*iter_handle = s;

    return 0;
}


#ifdef TEST_BAM_MAIN

//#define RG_COUNT

#ifdef RG_COUNT
#include "io_lib/hash_table.h"
#endif

int main(int argc, char **argv) {
    size_t i;
    bam_file_t *b;
    bam_seq_t *s;
    int nr = 0;
    int nm[2];
    char *mode = "rb";

#ifdef RG_COUNT
    HashTable *h = HashTableCreate(4, HASH_FUNC_HSIEH | HASH_DYNAMIC_SIZE);
#endif

    if (argc != 2) {
	fprintf(stderr, "Usage: %s bam_file\n", argv[0]);
	return 1;
    }
    
    {
	size_t l = strlen(argv[1]);
	if (l >= 4 && argv[1][l-3] == 's')
	    mode = "r";
    }

    b = bam_open(argv[1], mode);
    if (!b) {
	fprintf(stderr, "Failed to open bam file\n");
	return 1;
    }

    //b->no_aux = 1;

    nm[0] = nm[1] = 0;
    s = NULL;
    while (bam_next_seq(b, &s) > 0) {
	nr++;
	nm[(bam_flag(s) & BAM_FUNMAP) > 0]++;
#ifdef RG_COUNT
	//printf("%s RG='%s'\n", bam_name(s), bam_aux_find(s, "RG"));
	{
	    char *rg = bam_aux_find(s, "RG");
	    HashItem *hi;
	    if (!rg)
		rg = "unknown";

	    hi = HashTableSearch(h, rg, 0);
	    if (!hi) {
		HashData hd;
		hd.p = calloc(2,sizeof(int));
		hi = HashTableAdd(h, rg, 0, hd, NULL);
	    }
	    ((int *)hi->data.p)[(bam_flag(s) & BAM_FUNMAP) > 0]++;
	}
#endif
    }

    printf("Mapped      = %d\n",nm[0]);
    printf("Unmpped     = %d\n",nm[1]);
    printf("Total reads = %d\n",nr);

#ifdef RG_COUNT
    {
	HashIter *iter = HashTableIterCreate();
	HashItem *hi;
	
	printf("Mapped\tUnmap\tRead-group\n");
	while (hi = HashTableIterNext(h, iter)) {
	    printf("%d\t%d\t%.*s\n",
		   ((int *)hi->data.p)[0],
		   ((int *)hi->data.p)[1],
		   hi->key_len, hi->key);
	}
	HashTableIterDestroy(iter);
	HashTableDestroy(h, 1);
    }
#endif

    bam_close(b);

    return 0;
}

#endif /* TEST_MAIN */
