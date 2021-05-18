/*  Miscelaneous test functions.
 *
 *  Rather than deal with parallel I/O, we just list the tests here:
 *     1.  Griewank
 *     2.  Rastrigin
 *     3.  Schwefel
 *     4.  Ackley
 */
#include <pgapack.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265354
#endif

double griewank(PGAContext *, int, int);
double rastrigin(PGAContext *, int, int);
double schwefel(PGAContext *, int, int);
double ackley(PGAContext *, int, int);

void   printResultInterpretation(PGAContext *, int);
int    GetIntegerParameter(char *query);

int    NumCoords[4]  = { 10, 20, 10, 20};
double Lower[4]      = { -512.0, -5.12, -512.0, -32.768};
double Upper[4]      = { 511.0, 5.11, 511.0, 32.767};

const int maxN = 100;

int ackleyA = 100;
double ackleyB = 0.3;
double ackleyC = 1.85 * M_PI;

/*******************************************************************
 *                   user main program                              *
 *******************************************************************/
int main( int argc, char **argv ) {
    PGAContext *ctx;     /* the context variable */
    int testnum;         /* the DeJong test to run */
    int maxiter;         /* the maximum number of iterations */
    int n;
    double l[maxN], u[maxN]; /* for initializing lu ranges */
    int i;

    clock_t start = clock();

    MPI_Init(&argc, &argv);

    // testnum = GetIntegerParameter("Which test? (1-Griewank, 2-Rastrigin, 3-Schwefel, 4-Ackley)\n") - 1;
    // maxiter = GetIntegerParameter("How many iterations?\n");

    testnum = 3;
    maxiter = 100000;

    n = NumCoords[testnum];

    // if(testnum == 3)
    // {
    //     int overrideParams = GetIntegerParameter("Want to override parameters(0 - false, 1 - true)?\nDefault Parameters are (n = 20, a = 70, b = 0.2, c = 1.0 * PI)\n");
    //     if(overrideParams == 1)
    //     {
    //         n = GetIntegerParameter("N (5-100)\n");
    //         ackleyA = GetIntegerParameter("A(int)\n");
    //         ackleyB = ((double)GetIntegerParameter("B / 100\n")) / 100.0;
    //         ackleyC = ((double)GetIntegerParameter("C / 100 * PI\n")) / 100.0 * M_PI;
    //     }
    // }

    for (i = 0; i < n; i++)
    {
        l[i] = Lower[testnum];
        u[i] = Upper[testnum];
    }

    ctx = PGACreate(&argc, argv, PGA_DATATYPE_REAL, 
		    n, PGA_MINIMIZE);

    PGASetRandomSeed(ctx, 1);

    PGASetRealInitRange(ctx, l, u);
    PGASetMaxGAIterValue(ctx, maxiter);
    
    PGASetUp(ctx);

    if (testnum == 0)    PGARun(ctx, griewank);
    if (testnum == 1)    PGARun(ctx, rastrigin);
    if (testnum == 2)    PGARun(ctx, schwefel);
    if (testnum == 3)    PGARun(ctx, ackley);

    PGADestroy(ctx);
    
    MPI_Finalize();

    clock_t difference = clock() - start;
    int msec = difference * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %ds:%dms\n", msec / 1000, msec % 1000);
    return 0;
}


double griewank(PGAContext *ctx, int p, int pop) {
    int i, len;
    double term, sum, product;

    sum = 0;
    product = 1;
    len = PGAGetStringLength(ctx);
    for (i = 0; i < len; i++) {
        term = PGAGetRealAllele(ctx, p, pop, i);
        sum = sum + term * term / 4000.0;
        product = product * cos(term / sqrt(((double)i + 1)));
    }

    return (1 + sum - product);
}


double rastrigin(PGAContext *ctx, int p, int pop)
{
    int i, len;
    double term, sum;

    sum = 0;
    len = PGAGetStringLength(ctx);
    for (i = 0; i < len; i++) {
        term = PGAGetRealAllele(ctx, p, pop, i);
        sum = sum + term * term - 10 * cos(2 * M_PI * term);
    }
    return (len * 10 + sum);
}


double schwefel(PGAContext *ctx, int p, int pop) {
    int i, len;
    double term, sum;

    sum = 0;
    len = PGAGetStringLength(ctx);
    for (i = 0; i < len; i++) {
        term = PGAGetRealAllele(ctx, p, pop, i);
        sum = sum - term * sin(sqrt(fabs(term)));
    }
    return (sum);
}


double ackley(PGAContext *ctx, int p, int pop) {
    int i, len;
    double term, sum, s1, s2;
    double a, b, c;

    sum = 0;
    len = PGAGetStringLength(ctx);
    s1 = 0;
    s2 = 0;
    for (i = 0; i < len; i++)
    {
        term = PGAGetRealAllele(ctx, p, pop, i);
        s1 = s1 + term * term;
        s2 = s2 + cos(ackleyC * term);
    }
    sum = -ackleyA * exp(-ackleyB * sqrt(1.0 / len * s1)) - exp(1.0 / len * s2) + ackleyA + exp(1.0);
    return (sum);
}

/*  Get an integer parameter from the user.  Since this is
 *  typically a parallel program, we must only do I/O on the
 *  "master" process -- process 0.  Once we read the parameter,
 *  we broadcast it to all the other processes, then every 
 *  process returns the correct value.
 */
int GetIntegerParameter(char *query) {
    int  rank, tmp;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf(query);
        scanf("%d", &tmp);
    }
    MPI_Bcast(&tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return(tmp);
}



