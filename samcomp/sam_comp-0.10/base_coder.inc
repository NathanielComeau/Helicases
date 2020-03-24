/*
 * A fixed alphabet encoder for symbols 0 to 5 inclusive.
 * We have no escape symbol or resorting of data. It simply
 * accumulates stats and encodes in proportion.
 *
 * The intention is to use this per column for encoding
 * bases that differ to the consensus.
 */

/*
 * Fast mode uses 8-bit counters. It's generally 40% faster or so, but less
 * accurate and so larger file size (maybe only 1% though unless deep data).
 */

// Small enough to not overflow uint16_t even after +STEP
#ifdef WSIZ
#  undef WSIZ
#endif

template <typename st_t>
struct BASE_CODER {
    enum { STEP = sizeof(st_t) == 1 ? 1 : 300 };
    enum { WSIZ = (1<<8*sizeof(st_t))-3*STEP };

    BASE_CODER();
    BASE_CODER(int *start);
    void reset();
    void reset(int *start);
    inline void encodeSymbol(RangeCoder *rc, uint sym);
    inline void updateSymbol(uint sym);
    inline uint decodeSymbol(RangeCoder *rc);

//protected:
    void   rescaleRare();

    st_t SummFreq;
    st_t Stats[5];
};

template <typename st_t>
BASE_CODER<st_t>::BASE_CODER()
{
    reset();
}

template <typename st_t>
BASE_CODER<st_t>::BASE_CODER(int *start) {
    SummFreq = 0;
    for (int i = 0; i < 5; i++) {
	Stats[i] =  start[i];
	SummFreq += start[i];
    }
    assert(SummFreq > 0);
}

template <typename st_t>
void BASE_CODER<st_t>::reset() {
    SummFreq = 0;

    for ( int i=0; i<5; i++ )
	SummFreq += (Stats[i] = 3*STEP);
    assert(SummFreq > 0);
}

template <typename st_t>
void BASE_CODER<st_t>::reset(int *start) {
    SummFreq = 0;
    for (int i = 0; i < 5; i++) {
	Stats[i] =  start[i];
	SummFreq += start[i];
    }
    assert(SummFreq > 0);
}

template <typename st_t>
void BASE_CODER<st_t>::rescaleRare()
{
    Stats[0] -= Stats[0]>>1;
    Stats[1] -= Stats[1]>>1;
    Stats[2] -= Stats[2]>>1;
    Stats[3] -= Stats[3]>>1;
    Stats[4] -= Stats[4]>>1;
    SummFreq = Stats[0] + Stats[1] + Stats[2] + Stats[3] + Stats[4];
    assert(SummFreq > 0 && SummFreq < WSIZ);
}

template <typename st_t>
inline void BASE_CODER<st_t>::encodeSymbol(RangeCoder *rc, uint sym) {
    if ( SummFreq>WSIZ ) rescaleRare();

    switch(sym) {
    case 0:
	rc->Encode(0,
		   Stats[0], SummFreq);
	Stats[0] += STEP; 
	SummFreq += STEP;
	assert(SummFreq > 0);
	return;
    case 1:
	rc->Encode(Stats[0],
		   Stats[1], SummFreq);
	Stats[1] += STEP;
	SummFreq += STEP;
	assert(SummFreq > 0);
	return;
    case 2:
	rc->Encode(Stats[0] + Stats[1],
		   Stats[2], SummFreq);
	Stats[2] += STEP;
	SummFreq += STEP;
	assert(SummFreq > 0);
	return;
    case 3:
	rc->Encode(Stats[0] + Stats[1] + Stats[2],
		   Stats[3], SummFreq);
	Stats[3] += STEP;
	SummFreq += STEP;
	assert(SummFreq > 0);
	return;
    case 4:
	rc->Encode(Stats[0] + Stats[1] + Stats[2] + Stats[3],
		   Stats[4], SummFreq);
	Stats[4] += STEP;
	SummFreq += STEP;
	assert(SummFreq > 0);
	return;
    }

    return;
}

template <typename st_t>
inline void BASE_CODER<st_t>::updateSymbol(uint sym) {
  if ( SummFreq>WSIZ ) rescaleRare();

  if (sym == 4)
      return;

  /* known symbol */
  Stats[sym] += STEP;            
  SummFreq   += STEP;
  assert(SummFreq > 0);
}

template <typename st_t>
inline uint BASE_CODER<st_t>::decodeSymbol(RangeCoder *rc) {
    if ( SummFreq>WSIZ) rescaleRare();

    uint count=rc->GetFreq(SummFreq);
    uint HiCount=0;             

    st_t* p=Stats;
    if ((HiCount += *p) > count) {
	rc->Decode(0, *p, SummFreq);
	Stats[0] += STEP;
	SummFreq += STEP;
	return 0;
    }

    if ((HiCount += *++p) > count) {
	rc->Decode(HiCount-*p, *p, SummFreq);
	Stats[1] += STEP;
	SummFreq += STEP;
	return 1;
    }

    if ((HiCount += *++p) > count) {
	rc->Decode(HiCount-*p, *p, SummFreq);
	Stats[2] += STEP;
	SummFreq += STEP;
	return 2;
    }

    if ((HiCount += *++p) > count) {
	rc->Decode(HiCount-*p, *p, SummFreq);
	Stats[3] += STEP;
	SummFreq += STEP;
	return 3;
    }

    rc->Decode(HiCount, Stats[4], SummFreq);
    Stats[4] += STEP;
    SummFreq += STEP;
    return 4;
}
