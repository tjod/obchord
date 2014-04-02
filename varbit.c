/* Copyright TJ O'Donnell 2009, 2010 */

/* useful bit string functions */

#include "postgres.h"			/* general Postgres declarations */

#include "fmgr.h"			/* for argument/result macros */
#include "executor/executor.h"		/* for GetAttributeByName() */

#include "utils/varbit.h"

PG_MODULE_MAGIC;

/* These prototypes just prevent possible warnings from gcc. */

Datum	nbits_set(PG_FUNCTION_ARGS);
Datum	isbit_set(PG_FUNCTION_ARGS);
Datum	bit_set(PG_FUNCTION_ARGS);
Datum	bit_contains(PG_FUNCTION_ARGS);

PG_FUNCTION_INFO_V1(nbits_set);
PG_FUNCTION_INFO_V1(isbit_set);
PG_FUNCTION_INFO_V1(bit_set);
PG_FUNCTION_INFO_V1(bit_contains);

#ifndef SET_VARSIZE
#define SET_VARSIZE(v,l) (VARATT_SIZEP(v) = (l))
#endif 

Datum
nbits_set(PG_FUNCTION_ARGS)
{
/* how many bits are set in a bitstring? */

	VarBit	   *a = PG_GETARG_VARBIT_P(0);
	int n=0;

	/*
	 * VARBITLENTOTAL is the total size of the struct in bytes.
	 * VARBITLEN is the size of the string in bits.
	 * VARBITBYTES is the size of the string in bytes.
	 * VARBITHDRSZ is the total size of the header in bytes.
	 * VARBITS is a pointer to the data region of the struct.
	 */
	int i;
	unsigned char *ap = VARBITS(a);
	unsigned char aval;
	for (i=0; i < VARBITBYTES(a); ++i) {
		aval = *ap; ++ap;
		if (aval == 0) continue;
		if (aval & 1) ++n;
		if (aval & 2) ++n;
		if (aval & 4) ++n;
		if (aval & 8) ++n;
		if (aval & 16) ++n;
		if (aval & 32) ++n;
		if (aval & 64) ++n;
		if (aval & 128) ++n;
	}
	PG_RETURN_INT32(n);
}

Datum
bit_set(PG_FUNCTION_ARGS)
{
/* set the n'th (arg 2) bit in a bitstring (arg 1) */

	VarBit  *a = PG_GETARG_VARBIT_P(0);
	int      n = PG_GETARG_INT32(1);
	VarBit  *result;

	/*
	 * VARBITLENTOTAL is the total size of the struct in bytes.
	 * VARBITLEN is the size of the string in bits.
	 * VARBITBYTES is the size of the string in bytes.
	 * VARBITHDRSZ is the total size of the header in bytes.
	 * VARBITS is a pointer to the data region of the struct.
	 */

	int alen = VARBITLEN(a);
	if (n < 1 || n > alen) {
		elog(ERROR, "no bit #%d in bitstring",  n);
		PG_RETURN_NULL();
	}

/* create result bitstring and copy input bitstring to result */
        int rlen = VARBITTOTALLEN(alen);
	result = (VarBit *) palloc0(rlen);
        SET_VARSIZE(result, rlen);
        VARBITLEN(result) = alen;
	unsigned char *ap = VARBITS(a);
	unsigned char *rp = VARBITS(result);
        memcpy(rp, ap, rlen);

/* set appropriate bit in appropriate byte of result */
	int i = (n-1)/8;        /* i'th byte */
	int j = 8 - (n - i*8) ; /* j'th bit */
	rp[i] |= 1 << j;

        PG_RETURN_VARBIT_P(result);
}

Datum
isbit_set(PG_FUNCTION_ARGS)
{
/* is the n'th (arg 2) bit set in a bitstring (arg 1) */

	VarBit  *a = PG_GETARG_VARBIT_P(0);
	int      n = PG_GETARG_INT32(1);

	/*
	 * VARBITLENTOTAL is the total size of the struct in bytes.
	 * VARBITLEN is the size of the string in bits.
	 * VARBITBYTES is the size of the string in bytes.
	 * VARBITHDRSZ is the total size of the header in bytes.
	 * VARBITS is a pointer to the data region of the struct.
	 */

	int alen = VARBITLEN(a);
	if (n < 1 || n > alen) {
		elog(ERROR, "no bit #%d in bitstring",  n);
		PG_RETURN_NULL();
	}

	unsigned char *ap = VARBITS(a);

/* check appropriate bit in appropriate byte of result */
	int i = (n-1)/8;        /* i'th byte */
	unsigned char j = 8 - (n - i*8) ; /* j'th bit */
	if ( ap[i] & (1 << j) ) {
	        PG_RETURN_BOOL(true);
	} else {
		PG_RETURN_BOOL(false);
	}
}

Datum
bit_contains(PG_FUNCTION_ARGS)
{
/* Is bitstring (arg2) contained within bitstring (arg 1) */

	VarBit  *a = PG_GETARG_VARBIT_P(0);
	VarBit  *b = PG_GETARG_VARBIT_P(1);

	/*
	 * VARBITLENTOTAL is the total size of the struct in bytes.
	 * VARBITLEN is the size of the string in bits.
	 * VARBITBYTES is the size of the string in bytes.
	 * VARBITHDRSZ is the total size of the header in bytes.
	 * VARBITS is a pointer to the data region of the struct.
	 */

	int alen = VARBITLEN(a);
	int blen = VARBITLEN(b);
	//if (alen == 0 || blen == 0) PG_RETURN_NULL();
	if (alen != blen) {
/*
		elog(ERROR, "bitstrings must be of equal length");
		PG_RETURN_NULL();
*/
		int ai,bi;
		for (ai=0; ai<32; ++ai) if (alen == 1<<ai) break;
		for (bi=0; bi<32; ++bi) if (blen == 1<<bi) break;
		if (ai == 32 || bi == 32) {
			//fprintf(stderr, "%d %d %d %d\n",alen,blen,ai,bi);
			elog(ERROR, "bitstrings of different length and not both power of 2 in length");
			PG_RETURN_NULL();
		}
		register unsigned char *ap = VARBITS(a);
		register unsigned char *bp = VARBITS(b);
		register int i,j;
		alen = VARBITBYTES(a);
		blen = VARBITBYTES(b);
		unsigned char abyte;
		unsigned char bbyte;
		for (i=0; i<alen && i<blen; ++i) {
			if (alen<blen) {
				abyte = *ap;
/* fold the longer bitstring */
				for (bbyte=0,j=0;j<blen;j+=alen) {
				//fprintf(stderr, "%x(%d) ",bbyte,j);
				 bbyte |= *(bp+j);
				}
			} else {
				bbyte = *bp;
/* fold the longer bitstring */
				for (abyte=0,j=0;j<alen;j+=blen) {
				 abyte |= *(ap+j);
				//fprintf(stderr, "%x(%d) ",abyte,j);
				}
			}
			//fprintf(stderr, " : %x %x %d %d\n",abyte,bbyte,alen,blen);
			if ( (abyte & bbyte) != bbyte ) PG_RETURN_BOOL(false);
			++ap; ++bp;
		}

/* equal length bitstrings */
	} else {
		register unsigned char *ap = VARBITS(a);
		register unsigned char *bp = VARBITS(b);
		alen = VARBITBYTES(a);
		register int i;
		for (i=0; i<alen; ++i) {
			if ( (*ap & *bp) != *bp ) PG_RETURN_BOOL(false);
			++ap; ++bp;
		}

	        PG_RETURN_BOOL(true);
	}
}
