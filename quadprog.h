#ifndef QUADPROG_H
#define QUADPROG_H
extern void init_genrand(unsigned long s);
extern unsigned long genrand_int32(void);
extern double genrand_real2(void);
extern double genrand_res53(void);
extern int quadprog_ineq(const double *H, const double *f, const double *A, const double *b, const double *x0, int objectivesLen, int constraintsLen, double *ans, double *fval);
#endif