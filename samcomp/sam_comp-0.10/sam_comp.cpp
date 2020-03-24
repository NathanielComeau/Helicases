/*
 * TODO:
 * - Add header with reference name and checksum.
 * - Add adler32 checksum of fields (names, seqs, qual) to detect encoding
 *   and decoding bugs..
 */

/*
Ideas:

Use offset into read (note strandedness) as a context for likelihood to differ
to consensus. (model_diff_cons, model_diff_ref)

More accurately we could compute a model for error rates along the
read and combine with depth for an indication of how likely the entire
column is to match with no errors. (model_col_type)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>

#include "bam.h"
#include "bit_vec.h"

/*
 * Define this if your reference sequences are in fasta .fa format.
 * Ours have been 2-bit encoded as a .fb format for speed and memory
 * efficiency.
 */
//#define FASTA_REFS

/* Maximum quality value */
#define MAX_QUAL 96

/* Maximum name length */
#define MAX_NAME 256

/*
 * Roughly maximum sequence. Going over this harms compression, but it
 * doesn't break the compressor or decompressor.
 */
#define MAX_SEQ  102400

/* Maximum cigar string length */
#define MAX_CIG  100000

typedef unsigned char uc;
enum format_t { SAM = 0, BAM = 1, FASTQ = 2, FQ1 = 3 };

/* Range Coder */
//static block_t *blk;
//static unsigned char *out_buf;
//static unsigned char *in_buf;
//static int in_ind, out_ind;
//static size_t out_sz;

/* Faster I/O methods instead of getchar/putchar */
#define F_BLK 1024*1024
static char f_blk[F_BLK];
static int f_off = 0;
static int f_sz = 0;

int f_getchar() {
    if (f_off >= f_sz) {
	f_sz = read(0, f_blk, F_BLK);
	if (f_sz <= 0)
	    return EOF;
    }

    return f_blk[f_off++];
}

int f_putchar(int c) {
    if (f_off >= F_BLK) {
	if (f_off != write(1, f_blk, f_off))
	    return EOF;
	f_off = 0;
    }

    f_blk[f_off++] = c;
    return 0;
}

int f_flush() {
    return (f_off == write(1, f_blk, f_off)) ? 0 : EOF;
}

//#define InpSrcByte() (in_buf[in_ind++])
//#define OutTgtByte(x) (out_buf[out_ind++] = (x))

#define InpSrcByte() getchar()
//#define InpSrcByte() f_getchar()

#define OutTgtByte(x) putchar(x)
//#define OutTgtByte(x) f_putchar(x)

#define ABS(a)   ((a)>0?(a):-(a))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))


//#include "clrf_1.cdr"
//#include "clr.cdr"
#include "clr_io.cdr"
//#include "rc.h"
static RangeCoder rc;

int lookup[256];
void lookup_init(void) {
    /* ACGTN* */
    memset(lookup, 4, 256);
    lookup[(uc)'A'] = lookup[(uc)'a'] = 0;
    lookup[(uc)'C'] = lookup[(uc)'c'] = 1;
    lookup[(uc)'G'] = lookup[(uc)'g'] = 2;
    lookup[(uc)'T'] = lookup[(uc)'t'] = 3;
    lookup[(uc)'N'] = lookup[(uc)'n'] = 4;
    lookup[(uc)'*'] = 5;
}

/* Models */
//#define ORDER_0_CODER SIMPLE_MODEL
//#include "order0_coder.inc"
#include "simple_model.h"

#include "base_coder.inc"

/* Sequence length */
SIMPLE_MODEL<256> model_len1;
SIMPLE_MODEL<256> model_len2;
SIMPLE_MODEL<2> model_same_len;

typedef uint16_t base_t;
//typedef uint8_t base_t;

/* Consensus initialisation model */
/* 125 is 2-3% faster and .5% larger */
#define MODEL_CONS_SZ 78125
//#define MODEL_CONS_SZ 125

static SIMPLE_MODEL<2> model_cons_type[256*128];  // consensus & ref match?
static SIMPLE_MODEL<8> model_base[36*6];      // consensus base ^ ref
static BASE_CODER<base_t> model_cons_base[MODEL_CONS_SZ];    // consensus base

/* Name level 2 */
#define MAX_TOK 1000
static SIMPLE_MODEL<10> model_name_type[MAX_TOK];
static SIMPLE_MODEL<256> model_name_alpha_len[MAX_TOK];
static SIMPLE_MODEL<256> model_name_alpha[MAX_TOK];
static SIMPLE_MODEL<256> model_name_zero[MAX_TOK];
static SIMPLE_MODEL<256> model_name_digit0[MAX_TOK];
static SIMPLE_MODEL<256> model_name_digit1[MAX_TOK];
static SIMPLE_MODEL<256> model_name_digit2[MAX_TOK];
static SIMPLE_MODEL<256> model_name_digit3[MAX_TOK];
static SIMPLE_MODEL<256> model_name_ddelta[MAX_TOK];
static SIMPLE_MODEL<256> model_name_char[MAX_TOK];

/* Reference name */
static SIMPLE_MODEL<2> model_same_ref;
static SIMPLE_MODEL<256> model_ref_name[256];

/* Positional models */
static SIMPLE_MODEL<256> model_pos_val1;
static SIMPLE_MODEL<256> model_pos_val2;
static SIMPLE_MODEL<256> model_pos_val3;
static SIMPLE_MODEL<256> model_pos_val4;

/* Flag, split into bytes */
static SIMPLE_MODEL<256> model_flag1;
static SIMPLE_MODEL<256> model_flag2;

/* Mapping quality */
static SIMPLE_MODEL<256> model_mapQ;

/* CIGAR models */
static SIMPLE_MODEL<2> model_cigar_allM;
static SIMPLE_MODEL<16>  model_cigar_op[16];
static SIMPLE_MODEL<256> model_cigar_len[16];

/* Differences to a known column */
#define MODSZ MAX_SEQ
static BASE_CODER<base_t> model_ref_pos[MODSZ];
static BASE_CODER<base_t> model_ins_pos[MODSZ];
static BASE_CODER<base_t> model_sclip_pos[MODSZ];

#define QBITS 12
#define QSIZE (1<<QBITS)
SIMPLE_MODEL<MAX_QUAL> model_qual[QSIZE*16*16];

#define NS 8
static BASE_CODER<base_t> model_seq8[1<<(2*NS)];

static void encode_len(int len) {
    static int last_len = 0;
    if (len != last_len) {
	model_same_len.encodeSymbol(&rc, 0);
	model_len1.encodeSymbol(&rc, len & 0xff);
	model_len2.encodeSymbol(&rc, (len >> 8) & 0xff);
	last_len = len;
    } else {
	model_same_len.encodeSymbol(&rc, 1);
    }
}

static int decode_len() {
    static int last_len = 0;
    if (model_same_len.decodeSymbol(&rc)) {
	return last_len;
    } else {
	int l1 = model_len1.decodeSymbol(&rc);
	int l2 = model_len2.decodeSymbol(&rc);
	last_len = l1 + (l2 << 8);
	return last_len;
    }
}

static void encode_rname(bam_file_t *fp, int ref) {
    static int last_ref = -1;
    static unsigned char last_char = 0;

    if (ref != last_ref) {
	char *rname = fp->ref[ref].name;
	size_t l = strlen(rname);

	//fprintf(stderr, "Encoding for %s\n", rname);
	model_same_ref.encodeSymbol(&rc, 0);
	    
	for (size_t i = 0; i <= l; i++) {
	    model_ref_name[last_char].encodeSymbol(&rc, rname[i]);
	    last_char = rname[i];
	}

	last_ref = ref;
    } else {
	model_same_ref.encodeSymbol(&rc, 1);
    }
}

static char *decode_rname(int *diff) {
    static unsigned char last_char = 0;
    static char rname[1024];

    if (model_same_ref.decodeSymbol(&rc)) {
	*diff = 0;
	return rname;
    }

    *diff = 1;
    int i = 0;
    do {
	rname[i] = model_ref_name[last_char].decodeSymbol(&rc);
	last_char = rname[i];
    } while (rname[i++]);

    return rname;
}

static void encode_pos(uint32_t pos) {
    /*
     * Variable sized integers. Seems easiest encoding is just
     * to do each byte at a time and let it accumulate stats.
     *
     * Ideally we'd use the previous byte as context, as val2 being 0
     * usually implies val3/4 are 0 too, but in reality the difference
     * is negligible.
     */
    model_pos_val1.encodeSymbol(&rc, (pos >> 0 ) & 0xff);
    model_pos_val2.encodeSymbol(&rc, (pos >> 8 ) & 0xff);
    model_pos_val3.encodeSymbol(&rc, (pos >> 16) & 0xff);
    model_pos_val4.encodeSymbol(&rc, (pos >> 24) & 0xff);
}

static uint32_t decode_pos() {
    int p1 = model_pos_val1.decodeSymbol(&rc);
    int p2 = model_pos_val2.decodeSymbol(&rc);
    int p3 = model_pos_val3.decodeSymbol(&rc);
    int p4 = model_pos_val4.decodeSymbol(&rc);

    return p1 | (p2 << 8) | (p3 << 16) | (p4 << 24);
}

static void encode_mapQ(int mq) {
    model_mapQ.encodeSymbol(&rc, mq);
}

static int decode_mapQ() {
    return model_mapQ.decodeSymbol(&rc);
}

static void encode_flags(uint32_t flag) {
    model_flag1.encodeSymbol(&rc, flag & 0xff);
    model_flag2.encodeSymbol(&rc, flag >> 8);
}

static uint32_t decode_flags(void) {
    uint32_t f1 = model_flag1.decodeSymbol(&rc);
    uint32_t f2 = model_flag2.decodeSymbol(&rc);

    return f1 + (f2<<8);
}

enum name_type {N_UNK = 0, N_ALPHA, N_CHAR,
		N_ZERO, N_DIGITS, N_DDELTA, N_MATCH, N_END};
static void encode_name2(char *name, int len) {
    int i, j, k;

    static int last_token_type[1024];
    static int last_token_int[1024];
    static int last_token_str[1024];

    static char last_name[MAX_NAME];

    //fprintf(stderr, "NAME: %.*s\n", len, name);

    int ntok = 0;
    for (i = j = 0, k = 0; i < len; i++, j++, k++) {
	/* Determine data type of this segment */
	if (isalpha(name[i])) {
	    int s = i+1;
	    while (s < len && isalpha(name[s]))
		s++;

	    if (last_token_type[ntok] == N_ALPHA) {
		if (s-i == last_token_int[ntok] &&
		    memcmp(&name[i], 
			   &last_name[last_token_str[ntok]],
			   s-i) == 0) {
		    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
		    model_name_type[ntok].encodeSymbol(&rc, N_MATCH);
		} else {
		    //fprintf(stderr, "Tok %d (alpha)\n", N_ALPHA);
		    model_name_type[ntok].encodeSymbol(&rc, N_ALPHA);
		    model_name_alpha_len[ntok].encodeSymbol(&rc, s-i);
		    for (int x = 0; x < s-i; x++) {
			model_name_alpha[ntok].encodeSymbol(&rc, name[i+x]);
		    }
		}
	    } else {
		//fprintf(stderr, "Tok %d (alpha)\n", N_ALPHA);
		model_name_type[ntok].encodeSymbol(&rc, N_ALPHA);
		model_name_alpha_len[ntok].encodeSymbol(&rc, s-i);
		for (int x = 0; x < s-i; x++) {
		    model_name_alpha[ntok].encodeSymbol(&rc, name[i+x]);
		}
	    }

	    last_token_int[ntok] = s-i;
	    last_token_str[ntok] = i;
	    last_token_type[ntok] = N_ALPHA;

	    i = s-1;
	} else if (name[i] == '0') {
	    int s = i, v;
	    while (s < len && name[s] == '0')
		s++;
	    v = s-i;

	    if (last_token_type[ntok] == N_ZERO) {
		if (last_token_int[ntok] == v) {
		    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
		    model_name_type[ntok].encodeSymbol(&rc, N_MATCH);
		} else {
		    //fprintf(stderr, "Tok %d (0)\n", N_ZERO);
		    model_name_type[ntok].encodeSymbol(&rc, N_ZERO);
		    model_name_zero[ntok].encodeSymbol(&rc, v);
		}
	    } else {
		//fprintf(stderr, "Tok %d (0)\n", N_ZERO);
		model_name_type[ntok].encodeSymbol(&rc, N_ZERO);
		model_name_zero[ntok].encodeSymbol(&rc, v);
	    }

	    last_token_int[ntok] = v;
	    last_token_type[ntok] = N_ZERO;

	    i = s-1;
	} else if (isdigit(name[i])) {
	    int s = i;
	    int v = 0;
	    int d = 0;
	    while (s < len && isdigit(name[s]) && v < (1<<27)) {
		v = v*10 + name[s] - '0';
		s++;
	    }

	    if (last_token_type[ntok] == N_DIGITS) {
		if ((d = v - last_token_int[ntok]) == 0) {
		    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
		    model_name_type[ntok].encodeSymbol(&rc, N_MATCH);
		} else if (d < 256 && d > 0) {
		    //fprintf(stderr, "Tok %d (delta)\n", N_DDELTA);
		    model_name_type[ntok].encodeSymbol(&rc, N_DDELTA);
		    model_name_ddelta[ntok].encodeSymbol(&rc, d);
		} else {
		    //fprintf(stderr, "Tok %d (dig)\n", N_DIGITS);
		    model_name_type[ntok].encodeSymbol(&rc, N_DIGITS);
		    model_name_digit0[ntok].encodeSymbol(&rc, (v>> 0) & 0xff);
		    model_name_digit1[ntok].encodeSymbol(&rc, (v>> 8) & 0xff);
		    model_name_digit2[ntok].encodeSymbol(&rc, (v>>16) & 0xff);
		    model_name_digit3[ntok].encodeSymbol(&rc, (v>>24) & 0xff);
		}
	    } else {
		//fprintf(stderr, "Tok %d (dig)\n", N_DIGITS);
		model_name_type[ntok].encodeSymbol(&rc, N_DIGITS);
		model_name_digit0[ntok].encodeSymbol(&rc, (v>> 0) & 0xff);
		model_name_digit1[ntok].encodeSymbol(&rc, (v>> 8) & 0xff);
		model_name_digit2[ntok].encodeSymbol(&rc, (v>>16) & 0xff);
		model_name_digit3[ntok].encodeSymbol(&rc, (v>>24) & 0xff);
	    }

	    last_token_int[ntok] = v;
	    last_token_type[ntok] = N_DIGITS;

	    i = s-1;
	} else {
	    if (last_token_type[ntok] == N_CHAR) {
		if (name[i] == last_token_int[ntok]) {
		    //fprintf(stderr, "Tok %d (mat)\n", N_MATCH);
		    model_name_type[ntok].encodeSymbol(&rc, N_MATCH);
		} else {
		    //fprintf(stderr, "Tok %d (chr)\n", N_CHAR);
		    model_name_type[ntok].encodeSymbol(&rc, N_CHAR);
		    model_name_char[ntok].encodeSymbol(&rc, name[i]);
		}
	    } else {
		//fprintf(stderr, "Tok %d (chr)\n", N_CHAR);
		model_name_type[ntok].encodeSymbol(&rc, N_CHAR);
		model_name_char[ntok].encodeSymbol(&rc, name[i]);
	    }

	    last_token_int[ntok] = name[i];
	    last_token_type[ntok] = N_CHAR;
	}

	ntok++;
    }
    //fprintf(stderr, "Tok %d (end)\n", N_END);
    model_name_type[ntok].encodeSymbol(&rc, N_END);
    
    memcpy(last_name, name, len);
}

static void decode_name2(char *name) {
    enum name_type tok;
    int ntok = 0, i = 0, v;

    static int last_token_type[1024];
    static int last_token_int[1024];
    static int last_token_str[1024];

    static char last_name[MAX_NAME];

    for (;;) {
	tok = (enum name_type)model_name_type[ntok].decodeSymbol(&rc);
	//fprintf(stderr, "tok=%d, last type=%d int=%d str=%d\n",
	//	tok, last_token_type[ntok], last_token_int[ntok],
	//	last_token_str[ntok]);
	if (tok == N_END)
	    break;

	switch (tok) {
	    /* Str delta too? */
	case N_ALPHA:
	    v = model_name_alpha_len[ntok].decodeSymbol(&rc);
	    last_token_int[ntok] = v; // len
	    last_token_str[ntok] = i;
	    for (int x = 0; x < v; x++)
		// also per 'x'; per pos in name? */
		name[i++] = model_name_alpha[ntok].decodeSymbol(&rc);
	    last_token_type[ntok] = N_ALPHA;
	    break;

	case N_CHAR:
	    v = model_name_char[ntok].decodeSymbol(&rc);
	    name[i++] = v;
	    last_token_int[ntok] = v;
	    last_token_type[ntok] = N_CHAR;
	    break;

	case N_ZERO:
	    v = model_name_zero[ntok].decodeSymbol(&rc);
	    last_token_int[ntok] = v;
	    for (int x = 0; x < v; x++)
		name[i++] = '0';
	    last_token_type[ntok] = N_ZERO;
	    break;

	case N_DIGITS: {
	    char rev[100];
	    int ri = 0, v0, v1, v2, v3;

	    v0 = model_name_digit0[ntok].decodeSymbol(&rc);
	    v1 = model_name_digit1[ntok].decodeSymbol(&rc);
	    v2 = model_name_digit2[ntok].decodeSymbol(&rc);
	    v3 = model_name_digit3[ntok].decodeSymbol(&rc);
	    v = v0 + (v1<<8) + (v2<<16) + (v3<<24);
	    last_token_int[ntok] = v;
	    while (v > 0) {
		rev[ri++] = '0' + (v%10);
		v /= 10;
	    }
	    while (ri > 0)
		name[i++] = rev[--ri];
	    last_token_type[ntok] = N_DIGITS;
	    break;
	}

	case N_DDELTA: {
	    char rev[100];
	    int ri = 0;

	    v = model_name_ddelta[ntok].decodeSymbol(&rc);
	    v += last_token_int[ntok];
	    last_token_int[ntok] = v;
	    while (v > 0) {
		rev[ri++] = '0' + (v%10);
		v /= 10;
	    }
	    while (ri > 0)
		name[i++] = rev[--ri];
	    last_token_type[ntok] = N_DIGITS;
	    break;
	}

	case N_MATCH:
	    switch (last_token_type[ntok]) {
	    case N_CHAR:
		name[i++] = last_token_int[ntok];
		break;

	    case N_ALPHA:
		v = last_token_int[ntok];
		for (int x = 0; x < v; x++)
		    name[i++] = last_name[last_token_str[ntok]+x];
		last_token_str[ntok] = i-v;
		break;

	    case N_ZERO:
		v = last_token_int[ntok];
		for (int x = 0; x < v; x++)
		    name[i++] = '0';
		break;

	    case N_DIGITS: {
		char rev[100];
		int ri = 0;
		v = last_token_int[ntok];

		while (v > 0) {
		    rev[ri++] = '0' + (v%10);
		    v /= 10;
		}
		while (ri > 0)
		    name[i++] = rev[--ri];
		break;
	    }
	    }
	    break;

	default:
	    fprintf(stderr, "Unexpected name token %d\n", tok);
	    return;
	}

	ntok++;
    }

    name[i] = '\0';
    memcpy(last_name, name, i);

    //fprintf(stderr, "%s\n", name);

    return;
}

static void encode_qual(char *qual, int len, int rev) {
    unsigned int last = 0;
    int i, j, len2 = len;
    int delta = 5;
    unsigned char q1 = 0, q2 = 0;

    if (rev) {
	for (i = 0, j = len-1; i < j; i++, j--) {
	    char t = qual[i];
	    qual[i] = qual[j];
	    qual[j] = t;
	}
    }

    for (i = 0; i < len2; i++) {
	unsigned char q = qual[i] & 0x7f;
	model_qual[last].encodeSymbol(&rc, q);
	last   = ((MAX(q1,q2)<<6) + q) & ((1<<QBITS)-1);
	last  += (q1==q2)<<QBITS;
	delta += (q1>q)*(q1-q);
	last  += MIN(7,delta>>3) << (QBITS+1);
	last  += (MIN(i+15,127)&(15<<3))<<(QBITS+1);

	q2 = q1; q1 = q;
    }

    if (len != len2)
	model_qual[last].encodeSymbol(&rc, 0x7f); /* terminator */
}

static void decode_qual(char *qual, /*char *seq,*/ int len, int rev) {
    unsigned int last = 0;
    int i, j;
    int delta = 5;
    unsigned char q1 = 0, q2 = 0;

    for (i = 0; i < len; i++) {
	unsigned char q = model_qual[last].decodeSymbol(&rc);

	qual[i] = q + '!';
	last   = ((MAX(q1,q2)<<6) + q) & ((1<<QBITS)-1);
	last  += (q1==q2) << QBITS;
	delta += (q1>q)*(q1-q);
	last  += MIN(7,delta>>3) << (QBITS+1);
	last  += (MIN(i+15,127)&(15<<3))<<(QBITS+1);

	q2 = q1; q1 = q;

    }

    if (rev) {
	for (i = 0, j = len-1; i < j; i++, j--) {
	    char t = qual[i];
	    qual[i] = qual[j];
	    qual[j] = t;
	}
    }
}


static unsigned char bam_base4(bam_seq_t *b, int pos) {
    // NACM GRSV TWYH KDBN
    static int L[] = {0, 0, 1, 0,
		      2, 0, 0, 0,
		      3, 0, 0, 0,
		      0, 0, 0, 0};
    return L[bam_seqi((unsigned char *)bam_seq(b), pos)];
}

static void encode_seq8(bam_seq_t *seq, int len) {
    /*
     * Replace N with the most likely base call (ie the one with 
     * the largest range). It's very minimal benefit though so may not
     * be worth the time.
     *
     * Or... just remove N completely? We still need to push something into
     * 'last' though to preserve hash-word alignment.
     *
     * Needs reordering of seq/qual code.
     */

    /* Corresponds to a 12-mer word that doesn't occur in human genome. */
    int last  = 0x7616c7 & ((1<<(2*NS))-1);
//    int last2 = 0x1c9791 & ((1<<(2*NS))-1);
    for (int i = 0; i < len; i++) {
	//if (seq[i] == 'N') continue;

	unsigned char  b = bam_base4(seq, i);
	model_seq8[last].encodeSymbol(&rc, b);

	last = (last*4 + b) & ((1<<(2*NS))-1);

	/*
	 * On deep data hashing both strands works well. On shallow data
	 * it adds a significant CPU hit for very minimal gains (at best
	 * 0.5% smaller seq portion).
	 * Eg: -6=>513049628, -6b=>510382143, -8b=501978520 (Seq portion only)
	 *
	 * Pre-seeding the hash table with double-stranded human genome
	 * hashes seems like a faster starting point and will help more
	 * for shallow data too. However even this isn't as significant
	 * as it sounds.
	 * -5=>516624591, -5h=>514002730
	 */
//	if (both_strands) {
//	    int b2 = last2 & 3;
//	    last2 = last2/4 + ((3-b) << (2*NS-2));
//	    model_seq8[last2].updateSymbol(b2);
//	}
    }
}

void decode_seq8(char *qual, char *seq, int len) {
    int last  = 0x7616c7 & ((1<<(2*NS))-1);
//   int last2 = 0x1c9791 & ((1<<(2*NS))-1);
    for (int i = 0; i < len; i++) {
	unsigned char b = model_seq8[last].decodeSymbol(&rc);
	seq[i] = "ACGTN"[b];
	last = (last*4 + b) & ((1<<(2*NS))-1);

	if (qual[i] == '!') seq[i] = 'N';

//	if (both_strands) {
//	    int b2 = last2 & 3;
//	    last2 = last2/4 + ((3-b) << (2*NS-2));
//	    model_seq8[last2].updateSymbol(b2);
//	}
    }
}

#if 0
static void print_cigar(FILE *fp, uint32_t *cigar, int ncigar) {
    fprintf(fp, "cigar ");
    for (int i = 0; i < ncigar; i++)
	fprintf(fp, "%d%c",
	       cigar[i] >> BAM_CIGAR_SHIFT,
	       "MIDNSHP=X"[cigar[i] & BAM_CIGAR_MASK]);
    fprintf(fp, "\n");
}

static void format_cigar(char *cig, uint32_t *cigar, int ncigar) {
    for (int i = 0; i < ncigar; i++)
	cig += sprintf(cig, "%d%c",
	       cigar[i] >> BAM_CIGAR_SHIFT,
	       "MIDNSHP=X"[cigar[i] & BAM_CIGAR_MASK]);
}
#endif

static void encode_cigar(bam_seq_t *b, int seq_len) {
    int i, nc = bam_cigar_len(b);
    uint32_t *cig = bam_cigar(b);
    int last_op = 0;

    /* Special op = remainder M */
    /* Special len = remainder len */

    if (nc == 1 && ((cigar_op)(*cig & BAM_CIGAR_MASK)) == BAM_CMATCH) {
	model_cigar_allM.encodeSymbol(&rc, 1);
	return;
    } else {
	model_cigar_allM.encodeSymbol(&rc, 0);
    }

    for (i = 0; i < nc; i++) {
	cigar_op op = (cigar_op)(cig[i] & BAM_CIGAR_MASK);
	int len = cig[i] >> BAM_CIGAR_SHIFT;

	switch (op) {
	case BAM_CMATCH:
	    seq_len -= len;
	    //if (seq_len == 0) op=BAM_CMATCH_TO_END;
	    break;

	case BAM_CBASE_MATCH:
	case BAM_CBASE_MISMATCH:
	case BAM_CINS:
	case BAM_CSOFT_CLIP:
	    seq_len -= len;
	    break;

	default:
	    break;
	}

	if (seq_len == 0 && (op != BAM_CHARD_CLIP || op != BAM_CDEL))
	    len = 0; /* Special case for remainder */

	model_cigar_op[last_op].encodeSymbol(&rc, op);
	if (len >= 0x80) {
	    do {
		model_cigar_len[op].encodeSymbol(&rc, (len & 0x7f) | ((len >= 0x80)?0x80:0));
		len >>= 7;
	    } while (len);
	} else {
	    model_cigar_len[op].encodeSymbol(&rc, len);
	}
	last_op = op;
    }

    model_cigar_op[last_op].encodeSymbol(&rc, 15);
}

/*
 * Decompresses the next cigar string and stores in p.
 * Seq_len is context required for the decompression.
 *
 * Debug: also return a static string form of cigar
 */
static void decode_cigar(int seq_len, uint32_t *cigar, int *nc) {
    int last_op = 0;

    *nc = 0;

    if (model_cigar_allM.decodeSymbol(&rc)) {
	cigar[(*nc)++] = (seq_len << BAM_CIGAR_SHIFT) | BAM_CMATCH;
	return;
    }

    //while (seq_len > 0) {
    for (;;) {
	int op  = model_cigar_op[last_op].decodeSymbol(&rc);
	int len;
	int len_tmp, len_shift;

	if (op == 15)
	    break;

	len = 0, len_shift = 0;
	do {
	    len_tmp = model_cigar_len[op].decodeSymbol(&rc);
	    len |= (len_tmp & 0x7f) << len_shift;
	    len_shift += 7;
	} while (len_tmp & 0x80);

	last_op = op;

	if (!len)
	    len = seq_len;

	cigar[(*nc)++] = (len << BAM_CIGAR_SHIFT) | op;

	switch(op) {
	case BAM_CMATCH:
	    seq_len -= len;
	    break;

	case BAM_CBASE_MATCH:
	case BAM_CBASE_MISMATCH:
	case BAM_CINS:
	case BAM_CSOFT_CLIP:
	    seq_len -= len;
	    break;
	}
    }

    return;
}

/* Encodes a new consensus base b vs expected ref. r */
static void encode_consensus(int spos, char b, char r) {
    static int last = 0, last_ct = 0;

    //model_base[r].encodeSymbol(&rc, b); return; //680701

//    /* Whether the consensus matches ref */
//    model_cons_type[last_ct].encodeSymbol(&rc, b == r); //470784

    //453495
    assert(last_ct*128 + MIN(spos,127) <= 256*128);
    //fprintf(stderr, "M[%5d] ", last_ct*128 + MIN(spos,127));
    model_cons_type[last_ct*128 + MIN(spos,127)].encodeSymbol(&rc, b == r);
    last_ct = (last_ct*2 + (b == r)) % 256;
    if (b==r) return;

    /* And if not, how it differs using last base and cons as ctx */
    if ((b^r) >= 8) {
	fprintf(stderr, "b=%d r=%d\n", b, r);
	abort();
    }

    model_base[last*6 + r].encodeSymbol(&rc, b ^ r); //104786
    last = (last*6+b)%36;
}

/* Decode a consensus base when given reference r */
static unsigned char decode_consensus(int spos, unsigned char r) {
    static int last = 0, last_ct = 0;
    int ct, b;

    //return model_base[r].decodeSymbol(&rc);
    
    /* Does it match */
    //ct = model_cons_type[last_ct].decodeSymbol(&rc);
    assert(last_ct*128 + MIN(spos,127) < 256*128);
    //fprintf(stderr, "M[%5d] ", last_ct*128 + MIN(spos,127));
    ct = model_cons_type[last_ct*128 + MIN(spos,127)].decodeSymbol(&rc);
    last_ct = (last_ct*2 + ct) % 256;
    if (ct)
	return r;

    /* No, so find the new base */
    b = model_base[last*6 + r].decodeSymbol(&rc) ^ r;
    last = (last*6+b)%36;

    return b;
}


/* --------------------------------------------------------------------------
 * Fast printf replacements.
 */

char *append_int(char *cp, int i) {
    int j, k = 0;

    if (i < 0) {
	if (i == INT_MIN)
	    return cp + sprintf(cp, "%d", INT_MIN);
	*cp++ = '-';
	i = -i;
    } else if (i == 0) {
	*cp++ = '0';
	return cp;
    }

    if (i < 1000)
	goto b1;
    if (i < 100000)
	goto b2;
    if (i < 100000000)
	goto b3;

    j = i / 1000000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000000000;

    j = i / 100000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 100000000;
    
 b3:
    j = i / 10000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 10000000;
    
    j = i / 1000000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000000;
    
    j = i / 100000;
    if (j || k) *cp++ = j + '0', k=1, i %= 100000;
    
 b2:
    j = i / 10000;
    if (j || k) *cp++ = j + '0', k=1, i %= 10000;

    j = i / 1000;
    if (j || k) *cp++ = j + '0', k=1, i %= 1000;

 b1:
    j = i / 100;
    if (j || k) *cp++ = j + '0', k=1, i %= 100;

    j = i / 10;
    if (j || k) *cp++ = j + '0', k=1, i %= 10;

    if (i || k) *cp++ = i + '0';

    return cp;
}

//#define EVEN_SIMPLER
int fast_fprintf(FILE *fp, char *fmt, ...) {
    char *cp;
    long l;
    int i;
    double d; 
    char tmp[128], *ct;
    va_list ap;

    va_start(ap, fmt);

    for (cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    cp++;
#ifndef EVEN_SIMPLER
	    {
		char c;

		/* Firstly, strip the modifier flags (+-#0 and [space]) */
		for(; (c=*cp);) {
		    if ('#' == c)
			;
		    else if ('-' == c || '+' == c || ' ' == c)
			;
		    else
			break;
		}
	    }

	    /* Width specifier */
	    if (*cp >= '0' && *cp <= '9') {
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }
#endif

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len2 = l;
		if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

#ifndef EVEN_SIMPLER
	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }
#else
	    arg_size = 0;
#endif

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		putc('%', fp);
		break;

	    case 'd':
	    case 'i':
	    case 'u':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else 
		    l = (long)va_arg(ap, int);

		ct = append_int(tmp, (int)l);
		*ct++ = 0;
		fputs(tmp, fp);
		break;

	    case 'c':
		i = va_arg(ap, int);
		putc(i, fp);
		break;

	    case 'f':
		d = va_arg(ap, double);
		fprintf(fp, "%f", d);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		fprintf(fp, "%g", d);
		break;

	    case 'p':
	    case 'x':
	    case 'X':
		l = (long)va_arg(ap, void *);
		puts("TODO");
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		if (conv_len2) {
		    fwrite(s, conv_len2, 1, fp);
		} else
		    fputs(s, fp);
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		puts("TODO");
	    }
	    
	}

	case '\0':
	    break;

	default:
	    putc(*cp, fp);
	}
    }

    va_end(ap);
    
    return 0;
}

/* Max 8k of output */
int faster_fprintf(FILE *fp, const char *fmt, ...) {
    const char *cp;
    long l;
    int i;
    double d; 
    va_list ap;
    char max_line[8192], *out = max_line;

    va_start(ap, fmt);

    for (cp = fmt; *cp; cp++) {
	switch(*cp) {

	/* A format specifier */
	case '%': {
	    long conv_len1=0, conv_len2=0, conv_len=0;
	    signed int arg_size;

	    cp++;

#ifndef EVEN_SIMPLER
	    {
		char c;

		/* Firstly, strip the modifier flags (+-#0 and [space]) */
		for(; (c=*cp);) {
		    if ('#' == c)
			;
		    else if ('-' == c || '+' == c || ' ' == c)
			;
		    else
			break;
		}
	    }

	    /* Width specifier */
	    if (*cp >= '0' && *cp <= '9') {
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len = conv_len1 = l;
	    } else if (*cp == '*') {
		conv_len = conv_len1 = (int)va_arg(ap, int);
		cp++;
	    }
#endif

	    /* Precision specifier */
	    if ('.' == *cp) {
		cp++;
		for (l = 0; *cp >= '0' && *cp <= '9'; cp++)
		    l = l*10 + *cp - '0';
		conv_len2 = l;
		if (*cp == '*') {
		    conv_len2 = (int)va_arg(ap, int);
		    cp++;
		}
		conv_len = MAX(conv_len1, conv_len2);
	    }

#ifndef EVEN_SIMPLER
	    /* Short/long identifier */
	    if ('h' == *cp) {
		arg_size = -1; /* short */
		cp++;
	    } else if ('l' == *cp) {
		arg_size = 1; /* long */
		cp++;
	    } else {
		arg_size = 0; /* int */
	    }
#else
	    arg_size = 0;
#endif

	    /* The actual type */
	    switch (*cp) {
	    case '%':
		/*
		 * Not real ANSI I suspect, but we'll allow for the
		 * completely daft "%10%" example.
		 */
		*out++ = '%';
		break;

	    case 'd':
	    case 'i':
	    case 'u':
		/* Remember: char and short are sent as int on the stack */
		if (arg_size == -1)
		    l = (long)va_arg(ap, int);
		else if (arg_size == 1)
		    l = va_arg(ap, long); 
		else 
		    l = (long)va_arg(ap, int);

		out = append_int(out, (int)l);
		break;

	    case 'c':
		i = va_arg(ap, int);
		*out++ = i;
		break;

	    case 'f':
		d = va_arg(ap, double);
		out += sprintf(out, "%f", d);
		break;

	    case 'e':
	    case 'E':
	    case 'g':
	    case 'G':
		d = va_arg(ap, double);
		out += sprintf(out, "%g", d);
		break;

	    case 'p':
	    case 'x':
	    case 'X':
		l = (long)va_arg(ap, void *);
		puts("TODO");
		break;

	    case 'n':
		/* produces no output */
		break;

	    case 's': {
		char *s = (char *)va_arg(ap, char *);
		if (conv_len2) {
		    //memcpy(out, s, conv_len2);
		    //out += conv_len2;

		    //strncpy(out, s, conv_len2);
		    //out += MIN(conv_len2, strnlen(s));

		    while(conv_len2 && *s) {
		        *out++ = *s++;
		        conv_len2--;
		    }
		} else {
		    //strcpy(out, s);
		    //out += strlen(s);

		    while(*s) *out++ = *s++;
		}
		//		if (conv_len2) {
		//		    fwrite(s, conv_len2, 1, fp);
		//		} else
		//		    fputs(s, fp);
		break;
	    }

	    default:
		/* wchar_t types of 'C' and 'S' aren't supported */
		puts("TODO");
	    }
	    
	}

	case '\0':
	    break;

	default:
	    *out++ = *cp;
	}
    }

    *out = 0;
    fwrite(max_line, 1, out-max_line, fp);

    va_end(ap);
    
    return 0;
}

#define fprintf faster_fprintf

/* --------------------------------------------------------------------------
 */

/*
 * Fast conversion from encoded SAM base nibble to a printable character
 */
static char tab[256][2];
static void init_tab(void) {
    int i, j;
    unsigned char b2;
    static int done = 0;

    if (done)
	return;

    for (i = 0; i < 16; i++) {
	for (j = 0; j < 16; j++) {
	    b2 = (i<<4) | j;
	    tab[b2][0] = "NACMGRSVTWYHKDBN"[i];
	    tab[b2][1] = "NACMGRSVTWYHKDBN"[j];
	}
    }

    done = 1;
}

static char *ref_dir = NULL;

#ifdef FASTA_REFS
/* Obtains a reference, via fasta format files */
static unsigned char get_ref_base(const char *fname, int rnum, int pos) {
    static int last_rnum = -1;
    static char *seq = NULL;

    if (rnum != last_rnum) {
	char path[PATH_MAX];
	struct stat sb;
	FILE *fa_fp;
	size_t j = 0;

	last_rnum = rnum;

	/* alloc & load it */
	sprintf(path, "%s/%s.fa", ref_dir, fname);
	fprintf(stderr, "Loading %s...", path);
	fflush(stderr);

	if (-1 == stat(path, &sb)) {
	    perror(path);
	    return 0;
	}

	if (seq)
	    free(seq);

	if (NULL == (seq = (char *)malloc(sb.st_size)))
	    return 0;

	fa_fp = fopen(path, "r");

	while (fgets(&seq[j], 1024, fa_fp)) {
	    if (seq[j] == '>')
		continue;

	    while (seq[j] != '\n')
		j++;
	}

	fclose(fa_fp);
	fprintf(stderr, "done\n");

	//assert(j == fp->ref[rnum].len);
    }

    return seq ? lookup[(uc)seq[pos-1]] : 0 /* A */;
}

#else

/* Obtains a reference, via 2-bit encoded files */
static unsigned char get_ref_base(const char *fname, int rnum, int pos) {
    static int last_rnum = -1;
    static char *seq = NULL;

    if (rnum != last_rnum) {
	char path[PATH_MAX];
	struct stat sb;
	FILE *fa_fp;

	last_rnum = rnum;

	/* alloc & load it */
	sprintf(path, "%s/%s.fb", ref_dir, fname);
	fprintf(stderr, "Loading %s...", path);
	fflush(stderr);

	if (-1 == stat(path, &sb)) {
	    perror(path);
	    return 0;
	}

	if (seq)
	    free(seq);

	if (NULL == (seq = (char *)malloc(sb.st_size)))
	    return 0;

	fa_fp = fopen(path, "r");
	if (sb.st_size != (off_t)fread(seq, 1, sb.st_size, fa_fp)) {
	    perror(path);
	    free(seq);
	    seq = NULL;
	}

	fclose(fa_fp);
	fprintf(stderr, "done\n");
    }

    return seq ? ((uc)seq[(pos-1)>>2] >> (2*((pos-1)&3))) & 3 : 0;
}
#endif

static unsigned char bam_ref_base(bam_file_t *fp, int rnum, int pos) {
    if (!ref_dir || rnum < 0 || rnum >= (int)fp->nref)
	return 0; /* A, guess */

    return get_ref_base(fp->ref[rnum].name, rnum, pos);
}

static char bam_base(bam_seq_t *b, int pos) {
    unsigned char *seq = (unsigned char *)bam_seq(b);
    return tab[seq[pos/2]][pos&1];
}

//#define bam_base(b,p) (tab[(unsigned char)(bam_seq((b)))[(p)/2]][(p)&1])

#define RLEN 300000000

static int sam_compress(bam_file_t *fp, int no_MQ) {
    bam_seq_t *b = NULL;
    int r;
    int last_pos = 1, last_cons = 0, last_ref = -1;
    int max_rpos = 0;
    bit_vec done(RLEN);

    init_tab();

    do {
	r = bam_next_seq(fp, &b);
	if (r == -1) {
	    fprintf(stderr, "bam_next_seq() failure on line %d\n", fp->line);
	    return -1;
	}
	if (!r)
	    break;

	if (b->ref != -1 && b->ref != last_ref) {
	    last_ref = b->ref;
	    //memset(done, 0, max_rpos);
	    done.clear_range(max_rpos);
	    max_rpos = 0;
	}

	encode_len(b->len);
	encode_flags(bam_flag(b));
	encode_pos(b->pos - last_pos); last_pos = b->pos;

	if (!(bam_flag(b) & 4)) {
	    encode_rname(fp, b->ref);
	    encode_cigar(b, b->len);
	}

	encode_mapQ(no_MQ ? 255 : bam_map_qual(b));
	//encode_name(bam_name(b));
	encode_name2(bam_name(b), strlen(bam_name(b)));
	encode_qual(bam_qual(b), b->len, bam_flag(b) & BAM_FREVERSE);

#if 1
	if (bam_flag(b) & 4) {
	    /* No alignment, so use fqzcomp seq method */
	    encode_seq8(b, b->len);
	    
	} else {
	    /* Output bases aligned to ref */
	    int rpos = b->pos, spos = 0;
#ifndef NDEBUG
	    int nc = bam_cigar_len(b);
#endif
	    uint32_t *cig = bam_cigar(b);
	    int cig_op = BAM_CUNKNOWN, cig_len = 0, cig_ind = 0;

	    while (spos < b->len) {
		int apos = bam_flag(b) & 0x10 ? b->len-1 - spos : spos;

		if (cig_len == 0) {
		    assert(cig_ind < nc);
		    cig_op = cig[cig_ind] & BAM_CIGAR_MASK;
		    cig_len = cig[cig_ind] >> BAM_CIGAR_SHIFT;
		    cig_ind++;
		    //printf("%d%c\n", cig_len, "MIDNSHP=X"[cig_op]);
		}

		switch (cig_op) {
		case BAM_CMATCH:
		case BAM_CBASE_MATCH:
		case BAM_CBASE_MISMATCH: {
		    //printf("M %d/%d => %c\n", rpos, spos, bam_base(b, spos));
		    char c = lookup[(uc)bam_base(b, spos)];

		    //if (!done[rpos]) {
		    if (!done.get(rpos)) {
			int s[6];
			char r;

			//done[rpos] = 1;
			done.set(rpos);

			if (ref_dir) {
			    r = bam_ref_base(fp, b->ref, rpos+1);
			    if (r > 3) r = 0;
			    //fprintf(stderr, "%d %d %d/%d\n", rpos, apos, r, c);
			    encode_consensus(apos, c, r);
			} else {
			    r = c;
			    model_cons_base[last_cons].encodeSymbol(&rc, c);
			    last_cons = (last_cons * 5 + c) % MODEL_CONS_SZ;
			}

			/*
			 * These are parameters we need to learn.
			 * See above for example methods.
			 */
			if (sizeof(base_t) == 1) {
			    s[0] = 1;
			    s[1] = 1;
			    s[2] = 1;
			    s[3] = 1;
			    s[4] = 1;
			    s[5] = 0;
			    s[(uc)c] += 5;
			    s[(uc)r] += 1;
			} else {
			    s[0] = 5;
			    s[1] = 5;
			    s[2] = 5;
			    s[3] = 5;
			    s[4] = 1; // 0 if we remove N
			    s[5] = 0;
			    s[(uc)c] += 200;
			    s[(uc)r] += 1800;
			}
			//fprintf(stderr, "%d -> %d %d %d %d %d\n", rpos, s[0], s[1], s[2], s[3], s[4]);
			model_ref_pos[rpos & (MODSZ-1)].reset(s);
			s[0] = 5;
		    } else {
			//fprintf(stderr, "%d %d \n", rpos, c);
			model_ref_pos[rpos & (MODSZ-1)].encodeSymbol(&rc, c);
		    }

		    if (max_rpos < rpos)
			max_rpos = rpos;

		    /*n/a model_ref_pos[rpos].encodeSymbol(&rc, c);*/
		    rpos++;
		    spos++;
		    cig_len--;
		    break;
		}

		case BAM_CDEL:
		case BAM_CREF_SKIP:
		    //printf("D %d/%d\n", rpos, spos);
		    rpos += cig_len;
		    cig_len = 0;
		    break;

		case BAM_CINS: {
		    char c = lookup[(uc)bam_base(b, spos)];
		    //printf("I %d/%d => %c\n",  rpos, spos, bam_base(b, spos));
		    model_ins_pos[(rpos + cig_len) & (MODSZ-1)].encodeSymbol(&rc, c);
		    spos++;
		    cig_len--;
		    break;
		}

		case BAM_CSOFT_CLIP:
		    //fprintf(stderr, "Error: soft clips not yet supported\n");
		    for (int i = 0; i < cig_len; i++) {
			char c = lookup[(uc)bam_base(b, spos+i)];
			model_sclip_pos[(rpos + i + MODSZ) & (MODSZ-1)].encodeSymbol(&rc, c);
			//fprintf(stderr, "Soft clip %d=%d\n", i, c);
		    }
		    //printf("S %d/%d => ...\n", rpos, spos);
		    spos += cig_len;
		    cig_len = 0;
		    break;

		case BAM_CHARD_CLIP:
		    cig_len = 0;
		    break;

		default:
		    fprintf(stderr, "Unhandled cigar_op %d\n", cig_op);
		    return -1;
		}
	    }
	}
#endif
    } while (r > 0);

    if (b)
	free(b);

    encode_len(0);

    return 0;
}

/* Inline reversal of str */
char *rev(char *str, int len) {
    for (int i = 0, j = len-1; i < j; i++, j--) {
	char c = str[i];
	str[i] = str[j];
	str[j] = c;
    }

    return str;
}

/* Inline reverse complementing of str */
char *revcomp(char *str, int len) {
    int i, j;
    static char M[256];

    if (M[0] == 0) {
	memset(M, 'N', 256);
	M[(uc)'A'] = 'T';
	M[(uc)'C'] = 'G';
	M[(uc)'G'] = 'C';
	M[(uc)'T'] = 'A';
    }

    for (i = 0, j = len-1; i < j; i++, j--) {
	char c = str[i];
	str[i] = M[(uc)str[j]];
	str[j] = M[(uc)c];
    }
    if (i == j)
	str[i] = M[(uc)str[i]];

    return str;
}

static int sam_decompress(format_t format) {
    const char *ref_seq;

    int seq_len = 100;
    char seq[MAX_SEQ], qual[MAX_SEQ];
    char cigar[MAX_CIG];
    char name[MAX_SEQ];
    int last_pos = 1, last_cons = 0;
    int ref_num = 0;
    int max_rpos = 0;
    bit_vec done(RLEN);

    init_tab();

    for (;;) {
	uint32_t cig[MAX_CIG];
	int nc = 0, pos, flags, mq;
	char *cp = cigar;

	if ((seq_len = decode_len()) == 0) break;
	flags = decode_flags();
	pos = decode_pos() + last_pos; last_pos = pos;
	pos++;

	if (flags & 4) {
	    ref_seq = "*";
	    cigar[0] = '*';
	    cigar[1] = '\0';
	} else {
	    int diff;
	    ref_seq = decode_rname(&diff);
	    decode_cigar(seq_len, cig, &nc);

	    if (diff) {
		//fprintf(stderr, "Decoding for %s\n", ref_seq);
		ref_num++;
		//memset(done, 0, max_rpos);
		done.clear_range(max_rpos);
		max_rpos = 0;
	    }
	}
	//print_cigar(stdout, cig, nc);
	mq = decode_mapQ();
	//decode_name(name);
	decode_name2(name);
	decode_qual(qual, seq_len, flags & BAM_FREVERSE);

#if 1
	if (flags & 4) {
	    decode_seq8(qual, seq, seq_len);
	} else {
	    /* Align bases to ref */
	    int rpos = pos-1, spos = 0;
	    int cig_op = BAM_CUNKNOWN, cig_len = 0, cig_ind = 0;

	    while (cig_ind < nc || spos < seq_len) {
		int apos = flags & 0x10 ? seq_len-1 - spos : spos;

		if (cig_len == 0) {
		    assert(cig_ind < nc);
		    cig_op = cig[cig_ind] & BAM_CIGAR_MASK;
		    cig_len = cig[cig_ind] >> BAM_CIGAR_SHIFT;
		    cig_ind++;
		    cp += sprintf(cp, "%d%c", cig_len, "MIDNSHP=X"[cig_op]);
		}

		switch (cig_op) {
		case BAM_CMATCH:
		case BAM_CBASE_MATCH:
		case BAM_CBASE_MISMATCH: {
		    //printf("M %d/%d => %c\n", rpos, spos, bam_base(b, spos));
		    char c;

		    //if (!done[rpos]) {
		    if (!done.get(rpos)) {
			int s[6];
			char r;

			//done[rpos] = 1;
			done.set(rpos);

			if (ref_dir) {
			    r = get_ref_base(ref_seq, ref_num, rpos+1);
			    if (r > 3) r = 0;
			    c = decode_consensus(apos, r);
			    //fprintf(stderr, "%d %d %d/%d\n", rpos, apos, r, c);
			} else {
			    c = model_cons_base[last_cons].decodeSymbol(&rc);
			    last_cons = (last_cons * 5 + c) % MODEL_CONS_SZ;
			    r = c;
			}

			if (sizeof(base_t) == 1) {
			    s[0] = 1;
			    s[1] = 1;
			    s[2] = 1;
			    s[3] = 1;
			    s[4] = 1;
			    s[5] = 0;
			    s[(uc)c] += 5;
			    s[(uc)r] += 1;
			} else {
			    s[0] = 5;
			    s[1] = 5;
			    s[2] = 5;
			    s[3] = 5;
			    s[4] = 1;
			    s[5] = 0;
			    s[(uc)c] += 200;
			    s[(uc)r] += 1800;
			}
			//fprintf(stderr, "%d -> %d %d %d %d %d\n", rpos, s[0], s[1], s[2], s[3], s[4]);
			model_ref_pos[rpos & (MODSZ-1)].reset(s);
		    } else {
			c = model_ref_pos[rpos & (MODSZ-1)].decodeSymbol(&rc);
			//fprintf(stderr, "%d %d\n", rpos, c);
		    }

		    if (max_rpos < rpos)
			max_rpos = rpos;

		    rpos++;

		    if (spos >= MAX_SEQ)
			abort();

		    seq[spos++] = "ACGTN*"[(uc)c];
		    cig_len--;
		    break;
		}

		case BAM_CDEL:
		case BAM_CREF_SKIP:
		    //printf("D %d/%d\n", rpos, spos);
		    rpos += cig_len;
		    cig_len = 0;
		    break;

		case BAM_CINS: {
		    //printf("I %d/%d => %c\n",  rpos, spos, bam_base(b, spos));
		    char c = model_ins_pos[(rpos+cig_len) & (MODSZ-1)].decodeSymbol(&rc);
		    seq[spos++] = "ACGTN*"[(uc)c];
		    cig_len--;
		    break;
		}

		case BAM_CSOFT_CLIP:
		    //printf("S %d/%d => ...\n", rpos, spos);
		    for (int i = 0; i < cig_len; i++) {
			char c = model_sclip_pos[(rpos+i+MODSZ) & (MODSZ-1)].decodeSymbol(&rc);
			seq[spos++] = "ACGTN*"[(uc)c];
		    }
		    cig_len = 0;
		    break;

		case BAM_CHARD_CLIP:
		    cig_len = 0;
		    break;

		default:
		    fprintf(stderr, "Unhandled cigar_op %d\n", cig_op);
		    return -1;
		}
	    }
	}
#else
	memset(seq, 'A', seq_len);
#endif

	switch (format) {
	case SAM:
	case BAM:
	    printf("%s\t%d\t%s\t%d\t%d\t%s\t*\t0\t0\t%.*s\t%.*s\n",
		   name, flags, ref_seq, pos, mq, cigar, seq_len, seq,
		   seq_len, qual);
	    break;

	case FASTQ:
	    if (flags & 16) {
		printf("@%s\n%.*s\n+\n%.*s\n",
		       name,
		       seq_len, revcomp(seq, seq_len),
		       seq_len, rev(qual, seq_len));
	    } else {
		printf("@%s\n%.*s\n+\n%.*s\n",
		       name, seq_len, seq, seq_len, qual);
	    }
	    break;

	case FQ1:
	    /* FQ1 is designed to be piped into tr '\011' '\012', but
	     * being on one line allows for greps, sorting, etc before
	     * the final transformation to fastq.
	     */
	    if (flags & 16) {
		printf("@%s\t%.*s\t+\t%.*s\n",
		       name,
		       seq_len, revcomp(seq, seq_len),
		       seq_len, rev(qual, seq_len));
	    } else {
		printf("@%s\t%.*s\t+\t%.*s\n",
		       name, seq_len, seq, seq_len, qual);
	    }
	    break;
	}
    }

    return 0;
}

void usage(int err) {
    FILE *fp = err ? stderr : stdout;

    fprintf(fp, "sam_comp SAM/BAM compression program.\n");
    fprintf(fp, "James Bonfield, 2012, 2014.\n\n");

    fprintf(fp, "sam_comp [options] < fn_in > fn_out\n\n");

    fprintf(fp, "To compress:\n");
    fprintf(fp, "    sam_comp [-r ref_dir] < foo.sam > foo.zam\n\n");

    fprintf(fp, "To uncompress:\n");
    fprintf(fp, "    sam_comp -d [-r ref_dir] [-f format] foo.zam > foo.sam\n\n");
    fprintf(fp, "Format defaults to SAM, but can be specified as one of\n");
    fprintf(fp, "SAM, fastq or fq1. The latter being a 1-line fastq.\n");


    exit(err);
}

int main(int argc, char **argv) {
    int decode = 0;
    format_t format = SAM;
    int opt;
    int no_MQ = 0;

    while ((opt = getopt(argc, argv, "hdr:f:M")) != -1) {
	switch (opt) {
	case 'h':
	    usage(0);

	case 'd':
	    decode = 1;
	    break;

	case 'r':
	    ref_dir = optarg;
	    break;
	    
	case 'f':
	    if (strcmp(optarg, "sam") == 0)
		format = SAM;
	    else if (strcmp(optarg, "SAM") == 0)
		format = SAM;
	    else if (strcmp(optarg, "bam") == 0)
		format = BAM;
	    else if (strcmp(optarg, "BAM") == 0)
		format = BAM;
	    else if (strcmp(optarg, "fastq") == 0)
		format = FASTQ;
	    else if (strcmp(optarg, "fq1") == 0)
		format = FQ1;
	    else
		usage(1);
	    break;

	case 'M':
	    no_MQ = 1;
	    break;

	default:
	    usage(1);
	}
    }

    //if (optind == argc)
    //    usage(1);

    lookup_init();

    if (decode) {
	rc.StartDecode();
	sam_decompress(format);
	rc.FinishDecode();
    } else {
	const char *fn;

	/* FIXME, write bam_fdopen() too */
	if (format != SAM && format != BAM) {
	    fprintf(stderr, "Input format must be SAM or BAM\n");
	    exit(1);
	}

	if (optind == argc)
	    fn = "/dev/stdin";
	else
	    fn = argv[optind];

	bam_file_t *fp = bam_open(fn, format == SAM ? "r" : "rb");
	if (!fp) {
	    perror(fn);
	    return 1;
	}
    
	rc.StartEncode();
	sam_compress(fp, no_MQ);
	rc.FinishEncode();

	bam_close(fp);
    }

    f_flush();

    return 0;
}
