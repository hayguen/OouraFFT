/* test of fft*g_h.c */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <ooura/fftxg_h.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))

/* random number generator, 0 <= RND < 1 */
#define RND(p) ((*(p) = (*(p) * 7141 + 54773) % 259200) * (1.0 / 259200.0))

#ifdef OOURA_SINGLE_PREC
#define ERR_LIMIT 2.0e-6    /* must be sufficient for 512 and 65536 */
#else
#define ERR_LIMIT 3.0e-15   /* must be sufficient for 512 and 65536 */
#endif

void putdata(int nini, int nend, OouraReal *a);
OouraReal errorcheck(int nini, int nend, OouraReal scale, OouraReal *a);


int main(int argc, char *argv[])
{
    int n, ret, retCode;
    OouraReal *a, err;

    printf("Ooura %s-precision FFT-library and sizeof(OouraReal) = %u\n", ooura_prec(), (unsigned)sizeof(OouraReal));

    if (1 < argc)
    {
        n = atoi(argv[1]);
        printf("running %s with N = %d\n", argv[0], n);
    }
    else
    {
        fprintf(stderr, "usage: %s <N>\n", argv[0]);
        fprintf(stderr, "  N = data length: must be 2^m\n");
        return 1;
    }

    a = (OouraReal*)malloc( (n) * sizeof(OouraReal));
    retCode = 0;

    /* check of CDFT */
    putdata(0, n - 1, a);
    cdft(n, 1, a);
    cdft(n, -1, a);
    err = errorcheck(0, n - 1, 2.0 / n, a);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("cdft err= %g\n", err);

    /* check of RDFT */
    putdata(0, n - 1, a);
    rdft(n, 1, a);
    rdft(n, -1, a);
    err = errorcheck(0, n - 1, 2.0 / n, a);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("rdft err= %g\n", err);

    /* check of DDCT */
    putdata(0, n - 1, a);
    ddct(n, 1, a);
    ddct(n, -1, a);
    a[0] *= 0.5;
    err = errorcheck(0, n - 1, 2.0 / n, a);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("ddct err= %g\n", err);

    /* check of DDST */
    putdata(0, n - 1, a);
    ddst(n, 1, a);
    ddst(n, -1, a);
    a[0] *= 0.5;
    err = errorcheck(0, n - 1, 2.0 / n, a);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("ddst err= %g\n", err);

    /* check of DFCT */
    putdata(0, n, a);
    a[0] *= 0.5;
    a[n] *= 0.5;
    dfct(n, a);
    a[0] *= 0.5;
    a[n] *= 0.5;
    dfct(n, a);
    err = errorcheck(0, n, 2.0 / n, a);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("dfct err= %g\n", err);

    /* check of DFST */
    putdata(1, n - 1, a);
    dfst(n, a);
    dfst(n, a);
    err = errorcheck(1, n - 1, 2.0 / n, a);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("dfst err= %g\n", err);

    free(a);

    return retCode;
}


void putdata(int nini, int nend, OouraReal *a)
{
    int j, seed = 0;

    for (j = nini; j <= nend; j++) {
        a[j] = RND(&seed);
    }
}


OouraReal errorcheck(int nini, int nend, OouraReal scale, OouraReal *a)
{
    int j, seed = 0;
    OouraReal err = 0, e;

    for (j = nini; j <= nend; j++) {
        e = RND(&seed) - a[j] * scale;
        err = MAX(err, fabs(e));
    }
    return err;
}

