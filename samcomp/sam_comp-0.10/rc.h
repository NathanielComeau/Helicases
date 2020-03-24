/*
 * Range Coder.
 *
 * This implements a basic carryless range coder by mixing 64-bit and
 * 32-bit code together (an idea I took from Eugene Shelwien's many
 * range coder implementations).
 *
 * To encode, call EncodeStart() to initialise and then multiple calls
 * to EncodeRange() for new symbol, ending on EncodeFlush() to flush
 * out any bits still left in the "low" value.
 *
 * Decoding needs DecodeStart() to prime the variables and then one
 * DecodeFreq() and DecodeRange() call per symbol. The first of these
 * returns the current encoded frequency to the model so it can look
 * up in its tables and identify the associated symbol. The
 * DecodeRange() call then updates the range coder variables in the
 * same manner than EncodeRange does.
 *
 * The input frequencies for {Encode,Decode}Range() should all be less
 * than 1<<16. Failure to do this will silently generate incorrect
 * results, or crash.
 */

#define RC RangeCoder

#define StartEncode EncodeStart
#define Encode EncodeRange
#define FinishEncode EncodeFlush

#define StartDecode DecodeStart
#define GetFreq DecodeFreq
#define Decode DecodeRange

typedef unsigned char uc;

#define TOP   0x00ffffff                   // lowest 24-bit set
#define TALL  0xffffffff                   // 32-bit all set
#define LMASK (uint64_t)0xff00000000000000 // bits 56-63 set

class RC {
 public:
    char *in_buf;
    char *out_buf;

    /* Buffer handling */
    void input (char *in)  { out_buf = in_buf = in;  }
    void output(char *out) { in_buf = out_buf = out; }
    char *input(void)  {return in_buf;}
    char *output(void) {return out_buf;}
    int size_out(void) {return out_buf - in_buf;}
    int size_in(void)  {return in_buf - out_buf;}

    /* Encoding */
    void EncodeStart(void);
    void EncodeRange(uint32_t cumFreq, uint32_t symFreq, uint32_t totFreq);
    void EncodeFlush(void);

    /* Decoding */
    void DecodeStart(void);
    uint32_t DecodeFreq(uint32_t totFreq);
    void DecodeRange(uint32_t cumFreq, uint32_t symFreq, uint32_t totFreq);

protected:
    /*
     * low being larger than range means that after low<<=8 range<<=8
     * the sum low+range still doesn't wrap around and become lower
     * than the previous range. This avoids needing to handle the
     * overflow explicitly.
     */
    uint64_t low;
    uint32_t range;
    //uint64_t inFreq;
    uint32_t inFreq; // significant bytes of encoded stream, for decoding
};

/* ----------------------------------------------------------------------------
 * Encoder
 */

void RC::EncodeStart(void) {
    low   = 0;
    range = TALL;
}

/*
 * Given a starting range from low to high (range = high-low), we
 * subdivide it into a smaller range in the ratios cumFreq/totFreq to
 * (cumFreq+symFreq)/totFreq. ie:
 *
 *       low                high    (range == high-low)
 *       |                  |
 *   ....///////////////////...
 * +
 *   .............XXXXXXXXX....
 *   <cumFreq----><symFreq>
 *   <totFreq----------------->
 * =>
 *   ....,,,,,,,,,,/////,,,,...
 *                 |    |
 *                 low  high
 *
 * We track when the range gets sufficiently small that we can output
 * the top bits of low and then shift these bits out from low & range,
 * with special care taken when low & high cross a byte boundary.
 */
void RC::EncodeRange(uint32_t cumFreq, uint32_t symFreq, uint32_t totFreq) {
    range /= totFreq;
    low  += range * cumFreq;
    range = range * symFreq;

    //fprintf(stderr, "low=%08x%08x\n", (uint32_t)(low>>32), (uint32_t)(low&0xffffffff));

    /*
     * If the top-bits of range are zero, then emit the top bits of
     * low. There is an assumption here that low and high share the
     * same top bits; ((low^high)&0xff000000) == 0.
     * If this is false then we have to shrink range until it is true,
     * so that decoding can be achieved.
     *
     * Eg low/high may be matching top-8
     *       88332211/88342211 (range 00010000)
     *
     * or mismatching top-8
     *       88ff2211/89002211 (range 00010000)
     *    => 88ff2211/88ffffff (adjusted high => range 0000ddee)
     * 
     * the latter is simply computing new high from low|ffffff (88ffffff)
     * and then subtracting from low again to get the new range.
     * (A consequence of this is often that we'll need to shift out
     * more low bits again.)
     */
    while (range <= TOP) {
	if ((low^(low+range)) & LMASK) {
	    uint32_t l32 = low;
	    range = (l32 | TOP) - l32;
	}
	OutTgtByte(uc(low>>56));
	low   <<= 8;
	range <<= 8;
    }
}

void RC::EncodeFlush( void ) {
    for (int i = 0; i < 8; i++) {
	OutTgtByte(uc(low>>56));
	low<<=8;
    }
}


/* ----------------------------------------------------------------------------
 * Decoder
 *
 * We keep track of the current encoded frequence (inFreq). This
 * frequency is used by the model to identify the appropriate symbol
 * X, as it will be somewhere between sum(freq(0) to freq(X-1)) and
 * sum(freq(0) to freq(X)).
 *
 * Having found the symbol and the corresponding cumulative frequency
 * and symbol frequency, we then pass them into DecodeRange function
 * to update the frequencies in an identical manner to the EncodeRange
 * function. (The only difference is editing inFreq too.)
 */
void RC::DecodeStart(void) { 
    inFreq = low = 0;
    range = TALL;

    inFreq  = ((uint64_t)InpSrcByte()) << 56;
    inFreq |= ((uint64_t)InpSrcByte()) << 48;
    inFreq |= ((uint64_t)InpSrcByte()) << 40;
    inFreq |= ((uint64_t)InpSrcByte()) << 32;
    inFreq |= ((uint64_t)InpSrcByte()) << 24;
    inFreq |= ((uint64_t)InpSrcByte()) << 16;
    inFreq |= ((uint64_t)InpSrcByte()) <<  8;
    inFreq |= ((uint64_t)InpSrcByte());
}

uint32_t RC::DecodeFreq(uint32_t totFreq) {
    return inFreq / (range/=totFreq);
}

void RC::DecodeRange(uint32_t cumFreq, uint32_t symFreq, uint32_t totFreq) {
    /*
     * See EncodeRange() too.
     *
     * inFreq starts off as the top-bits of the 'low' value we
     * previously had prior to a call to EncodeRange, but we're
     * modelling expanding our range instead of shrinking it. We also
     * know it won't suffer overflow (unlike low + range), so it's
     * sufficient to be 32-bit still.
     */
    uint32_t delta = cumFreq*range;
    //fprintf(stderr, "inFreq=%08x low=%08x%08x\n", inFreq,  (uint32_t)(low>>32), (uint32_t)(low&0xffffffff));
    inFreq -= delta;
    low    += delta;
    range  *= symFreq;
    

    while (range <= TOP) {
	/*
	 * FIXME: we only need low here to spot the overflow case of
	 * low and high top-bits being different.
	 * Can we force this to never be true in the encoder by proper
	 * carry handling, thus making ther decode have no need for low?
	 *
	 * See: http://encode.ru/threads/?p=9987&viewfull=1#post9987
	 */
	if ((low^(low+range)) & LMASK) {
	    uint32_t l32 = low;
	    range = (l32 | TOP) - l32;
	}
	
	/* Fetch more input */
	inFreq = (inFreq << 8) + InpSrcByte();
	low   <<= 8;
	range <<= 8;
    }
}


