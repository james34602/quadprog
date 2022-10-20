#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
//#include <vld.h>
#include "solvopt.h"
#include "quadprog.h"
void matrix_matrix_product(const double *a, const double *b, double *product, int rows1, int cols1, int rows2, int cols2, double weight)
{
	int i, j, k;
	for (i = 0; i < rows1; ++i)
		for (j = 0; j < cols2; ++j)
		{
			double res = 0.0;
			for (k = 0; k < cols1; ++k)
				res += a[i * cols1 + k] * b[k * cols2 + j] * weight;
			product[i * cols2 + j] = res;
		}
}
double unroll_dot_product(const double *x, const double *y, int n, double weight)
{
	double res = 0.0;
	int i = 0;
	for (; i <= n - 4; i += 4)
		res += (x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] + x[i + 3] * y[i + 3]) * weight;
	for (; i < n; i++)
		res += (x[i] * y[i]) * weight;
	return res;
}
void matrix_vector_mult(const double *mat, const double *vec, double *c, int rows, int cols, double weight)
{
	for (int i = 0; i < rows; i++)
	{
		double res = 0.0;
		for (int j = 0; j < cols; j++)
			res += (mat[i * cols + j] * vec[j] * weight);
		c[i] = res;
	}
}
// Quadprog
double obj(double x[], SharedFunctionVariables *pass)
{
	double ans1 = unroll_dot_product(pass->f, x, pass->dimY, 1.0);
	matrix_vector_mult(pass->H, x, pass->tmp1, pass->dimY, pass->dimY, 1.0);
	double ans2 = unroll_dot_product(x, pass->tmp1, pass->dimY, 0.5);
	return ans1 + ans2;
}
void gobj(double x[], double g[], SharedFunctionVariables *pass)
{
	matrix_vector_mult(pass->H, x, pass->tmp1, pass->dimY, pass->dimY, 1.0);
	for (int i = 0; i < pass->dimY; i++)
		g[i] = pass->tmp1[i] + pass->f[i];
}
double pen(double x[], SharedFunctionVariables *pass)
{
	matrix_vector_mult(pass->A, x, pass->tmp2, pass->dimX, pass->dimY, 1.0);
	double val, ans = 0.0;
	for (int i = 0; i < pass->dimX; i++)
	{
		val = pass->tmp2[i] - pass->b[i];
		if (val > 0.0)
			ans += val;
	}
	return ans;
}
void gpen(double x[], double g[], SharedFunctionVariables *pass)
{
	matrix_vector_mult(pass->A, x, pass->tmp2, pass->dimX, pass->dimY, 1.0);
	double val;
	for (int i = 0; i < pass->dimX; i++)
	{
		val = pass->tmp2[i] - pass->b[i];
		if (val < 0.0)
			val = 0.0;
		pass->tmp2[i] = val;
	}
	double ans = sqrt(unroll_dot_product(pass->tmp2, pass->tmp2, pass->dimX, 1.0));
	matrix_matrix_product(pass->tmp2, pass->A, g, 1, pass->dimX, pass->dimX, pass->dimY, ans > 0.0 ? 1.0 / ans : 1.0);
}
// Matlab inv(), warm start, b_x_data is nxn double, perm1 and perm2 is n int
void myinvWarm(const double *A, double *Y, double *b_x_data, int *perm1, int *perm2, int n)
{
	if (n == 1)
	{
		Y[0] = 1.0 / A[0];
		return;
	}
	int i, j, k, yk, i2, jA, ldap1, mmj_tmp, jp1j, jj, ix;
	double smax, s;
	memset(Y, 0, n * n * sizeof(double));
	memcpy(b_x_data, A, n * n * sizeof(double));
	perm1[0] = 1;
	yk = 1;
	for (k = 2; k <= n; k++)
		perm1[k - 1] = ++yk;
	ldap1 = n + 1;
	for (j = 0; j < n - 1; j++)
	{
		mmj_tmp = n - j;
		jA = j * ldap1;
		jj = j * ldap1;
		jp1j = jA + 2;
		yk = 0;
		ix = jA;
		smax = fabs(b_x_data[jA]);
		for (k = 2; k <= mmj_tmp; k++)
		{
			ix++;
			s = fabs(b_x_data[ix]);
			if (s > smax)
			{
				yk = k - 1;
				smax = s;
			}
		}
		if (b_x_data[jj + yk] != 0.0)
		{
			if (yk != 0)
			{
				yk += j;
				perm1[j] = yk + 1;
				ix = j;
				for (k = 0; k < n; k++)
				{
					smax = b_x_data[ix];
					b_x_data[ix] = b_x_data[yk];
					b_x_data[yk] = smax;
					ix += n;
					yk += n;
				}
			}
			for (i = jp1j; i <= jj + mmj_tmp; i++)
				b_x_data[i - 1] /= b_x_data[jj];
		}
		yk = jA + n;
		jA = jj + ldap1;
		for (jp1j = 0; jp1j <= mmj_tmp - 2; jp1j++)
		{
			smax = b_x_data[yk];
			if (b_x_data[yk] != 0.0)
			{
				ix = jj + 1;
				i2 = jA + 1;
				k = mmj_tmp + jA;
				for (i = i2; i < k; i++)
				{
					b_x_data[i - 1] += b_x_data[ix] * -smax;
					ix++;
				}
			}
			yk += n;
			jA += n;
		}
	}
	perm2[0] = 1;
	yk = 1;
	for (k = 2; k <= n; k++)
		perm2[k - 1] = ++yk;
	for (k = 0; k < n; k++)
	{
		if (perm1[k] > 1 + k)
		{
			yk = perm2[perm1[k] - 1];
			perm2[perm1[k] - 1] = perm2[k];
			perm2[k] = yk;
		}
	}
	for (k = 0; k < n; k++)
	{
		jp1j = n * (perm2[k] - 1);
		Y[k + jp1j] = 1.0;
		for (j = k + 1; j <= n; j++)
			if (Y[(j + jp1j) - 1] != 0.0)
				for (i = j + 1; i <= n; i++)
					Y[(i + jp1j) - 1] -= Y[(j + n * (perm2[k] - 1)) - 1] * b_x_data[(i + n * (j - 1)) - 1];
	}
	for (j = 0; j < n; j++)
	{
		yk = n * j - 1;
		for (k = n; k >= 1; k--)
		{
			jA = n * (k - 1) - 1;
			i2 = k + yk;
			if (Y[i2] != 0.0)
			{
				Y[i2] /= b_x_data[k + jA];
				for (i = 0; i <= k - 2; i++)
					Y[i + yk + 1] -= Y[i2] * b_x_data[i + jA + 1];
			}
		}
	}
}
void DenmanBeaversSqrtm(const double *A, int n, int nSteps, double *solution, char *workingBuffer)
{
	int i, j;
	memcpy(solution, A, n * n * sizeof(double));
	double *Y = (double*)workingBuffer;
	double *Xk = Y + n * n;
	double *Yk = Xk + n * n;
	double *tmp = Yk + n * n;
	double *bPtr = tmp + n * n;
	int *permPtr = (int*)(bPtr + n * n);
	int *perPtr2 = permPtr + n * n;
	memset(Y, 0, n * n * sizeof(double));
	for (i = 0; i < n; i++)
		Y[i * n + i] = 1.0;
	for (j = 0; j < nSteps; j++)
	{
		memcpy(Xk, solution, n * n * sizeof(double));
		memcpy(Yk, Y, n * n * sizeof(double));
		myinvWarm(Yk, tmp, bPtr, permPtr, perPtr2, n); // inv(Yk)
		for (i = 0; i < n * n; i++)
			solution[i] = (Xk[i] + tmp[i]) * 0.5;
		myinvWarm(Xk, tmp, bPtr, permPtr, perPtr2, n); // inv(Xk)
		for (i = 0; i < n * n; i++)
			Y[i] = (Yk[i] + tmp[i]) * 0.5;
	}
}
void SharedFunctionVariablesInit(SharedFunctionVariables *pass, const double *H, const double *f, const double *A, const double *b, int objectivesLen, int constraintsLen)
{
	pass->H = H;
	pass->f = f;
	pass->A = A;
	pass->b = b;
	pass->dimX = constraintsLen;
	pass->dimY = objectivesLen;
	pass->tmp1 = (double*)malloc((pass->dimX + pass->dimY) * sizeof(double));
	pass->tmp2 = pass->tmp1 + pass->dimY;
}
void SharedFunctionVariablesFree(SharedFunctionVariables *pass)
{
	if (pass->tmp1)
		free(pass->tmp1);
}
double defaultOptions[13] = { -1.0, DBL_EPSILON, DBL_EPSILON, 20e5, -1.0, DBL_EPSILON, 2.5, DBL_EPSILON }; // very high precision
// Return error code
int quadprog_ineq(const double *H, const double *f, const double *A, const double *b, const double *x0, int objectivesLen, int constraintsLen, double *ans, double *fval)
{
	int i, j;
	double solvopt_options[13];
	memcpy(solvopt_options, defaultOptions, sizeof(defaultOptions));
	//if any(eig(H) < 0)
	//error('Non-convex. (Quadratic form is not positive definite.)')
	//end
	for (i = 0; i < objectivesLen; i++)
		for (j = 0; j < objectivesLen; j++)
			if (fabs(H[i * objectivesLen + j] - H[j * objectivesLen + i]) > DBL_EPSILON)
				return -1;
	size_t workingBufferLength = (objectivesLen * objectivesLen
		+ objectivesLen
		+ constraintsLen
		+ constraintsLen
		+ objectivesLen
		+ objectivesLen * objectivesLen
		+ constraintsLen * objectivesLen
		+ constraintsLen
		+ objectivesLen * objectivesLen) * sizeof(double)
		+ (objectivesLen * objectivesLen * 5 * sizeof(double)) // DenmanBeaversSqrtm
		+ (objectivesLen * objectivesLen * 2 * sizeof(int)); // DenmanBeaversSqrtm
	char *workingMatrix = (char*)malloc(workingBufferLength * sizeof(char));
	double *HI = (double*)workingMatrix;
	// try a quick check on unconstrained solution
	double *pos = (double*)((workingMatrix + workingBufferLength) - (objectivesLen * objectivesLen * 2 * sizeof(int)) - (objectivesLen * objectivesLen * sizeof(double)));
	int *permPtr = (int*)(pos + objectivesLen * objectivesLen);
	int *perPtr2 = permPtr + objectivesLen * objectivesLen;
	// Check Hessian matrix is all zero, all zero -> reduce to linear programming problem
	int counter = 0;
	for (i = 0; i < objectivesLen * objectivesLen; i++)
		if (fabs(H[i]) > 0.0)
			counter++;
	if (counter == objectivesLen * objectivesLen)
		myinvWarm(H, HI, pos, permPtr, perPtr2, objectivesLen);
	else
		memcpy(HI, H, objectivesLen * objectivesLen * sizeof(double));
	double *mqf = HI + objectivesLen * objectivesLen; // inverse exists since H is positive definite
	matrix_vector_mult(HI, f, mqf, objectivesLen, objectivesLen, -1.0);
	double *quickCheckOfConstraintSatisfied = mqf + objectivesLen;
	matrix_vector_mult(A, mqf, quickCheckOfConstraintSatisfied, constraintsLen, objectivesLen, 1.0);
	counter = 0;
	for (i = 0; i < constraintsLen; i++)
		if (quickCheckOfConstraintSatisfied[i] <= b[i])
			counter++;
	SharedFunctionVariables pass = { 0 };
	if (counter == constraintsLen)
	{
		memcpy(ans, mqf, objectivesLen * sizeof(double));
		SharedFunctionVariablesInit(&pass, H, f, A, b, objectivesLen, constraintsLen);
		*fval = obj(mqf, &pass);
		SharedFunctionVariablesFree(&pass);
		free(workingMatrix);
		return 0;
	}
	if (!x0)
		for (i = 0; i < objectivesLen; i++)
			ans[i] = mqf[i] + genrand_res53() * 0.5;
	else
		memcpy(ans, x0, objectivesLen * sizeof(double));
	double *bn = quickCheckOfConstraintSatisfied + constraintsLen;
	double *fn = bn + constraintsLen;
	double *T = fn + objectivesLen;
	double *An = T + objectivesLen * objectivesLen;
	double *nA = An + constraintsLen * objectivesLen;
	double *Hn = nA + constraintsLen;
	char *sqrtmWorkingBuffer = (char*)(Hn + objectivesLen * objectivesLen);
	double sum;
	if (objectivesLen >= 10)
	{
		// geometric pre-conditioning:
		for (i = 0; i < objectivesLen; i++)
			bn[i] = b[i] - quickCheckOfConstraintSatisfied[i];
		memset(fn, 0, objectivesLen * sizeof(double));
		DenmanBeaversSqrtm(HI, objectivesLen, 50, T, sqrtmWorkingBuffer);
		matrix_matrix_product(A, T, An, constraintsLen, objectivesLen, objectivesLen, objectivesLen, 1.0);
		for (i = 0; i < constraintsLen; i++)
		{
			sum = 0.0;
			for (j = 0; j < objectivesLen; j++)
				sum += (An[i * objectivesLen + j] * An[i * objectivesLen + j]);
			nA[i] = sqrt(sum);
		}
		for (i = 0; i < constraintsLen; i++)
		{
			for (j = 0; j < objectivesLen; j++)
				An[i * objectivesLen + j] = An[i * objectivesLen + j] / nA[i];
			bn[i] = bn[i] / nA[i];
		}
		memset(Hn, 0, objectivesLen * objectivesLen * sizeof(double));
		for (i = 0; i < objectivesLen; i++)
			Hn[i * objectivesLen + i] = 1.0;
		SharedFunctionVariablesInit(&pass, Hn, fn, An, bn, objectivesLen, constraintsLen);
		solvopt_options[5] *= 0.1;
		solvopt(objectivesLen, ans, obj, gobj, solvopt_options, pen, gpen, &pass);
		matrix_vector_mult(T, ans, fn, objectivesLen, objectivesLen, 1.0);
		for (i = 0; i < objectivesLen; i++)
			ans[i] = fn[i] + mqf[i];
		pass.H = H;
		pass.f = f;
		*fval = obj(ans, &pass);
		matrix_vector_mult(A, ans, bn, constraintsLen, objectivesLen, 1.0);
		for (i = 0; i < constraintsLen; i++)
			bn[i] -= b[i];
		double max = bn[0];
		for (i = 1; i < constraintsLen; i++)
			if (bn[i] > max)
				max = bn[i];
		if (max < DBL_EPSILON)
		{
			SharedFunctionVariablesFree(&pass);
			free(workingMatrix);
			return 0;
		}
		else
		{
//			Cond 1
			pass.A = A;
			pass.b = b;
			solvopt_options[5] *= 0.1;
			*fval = solvopt(objectivesLen, ans, obj, gobj, solvopt_options, pen, gpen, &pass);
			matrix_vector_mult(A, ans, bn, constraintsLen, objectivesLen, 1.0);
			for (i = 0; i < constraintsLen; i++)
				bn[i] -= b[i];
			max = bn[0];
			for (i = 1; i < constraintsLen; i++)
				if (bn[i] > max)
					max = bn[i];
			while (max > solvopt_options[5])
			{
//				Cond 2
				for (i = 0; i < objectivesLen; i++)
					ans[i] = mqf[i] + genrand_res53() * 0.5;
				*fval = solvopt(objectivesLen, ans, obj, gobj, solvopt_options, pen, gpen, &pass);
				matrix_vector_mult(A, ans, bn, constraintsLen, objectivesLen, 1.0);
				for (i = 0; i < constraintsLen; i++)
					bn[i] -= b[i];
				max = bn[0];
				for (i = 1; i < constraintsLen; i++)
					if (bn[i] > max)
						max = bn[i];
			}
			SharedFunctionVariablesFree(&pass);
			free(workingMatrix);
			return 0;
		}
	}
	else
	{
		SharedFunctionVariablesInit(&pass, H, f, A, b, objectivesLen, constraintsLen);
		*fval = solvopt(objectivesLen, ans, obj, gobj, solvopt_options, pen, gpen, &pass);
		matrix_vector_mult(A, ans, bn, constraintsLen, objectivesLen, 1.0);
		for (i = 0; i < constraintsLen; i++)
			bn[i] -= b[i];
		double max = bn[0];
		for (i = 1; i < constraintsLen; i++)
			if (bn[i] > max)
				max = bn[i];
		while (max > solvopt_options[5])
		{
//			Cond 2
			for (i = 0; i < objectivesLen; i++)
				ans[i] = mqf[i] + genrand_res53() * 0.5;
			*fval = solvopt(objectivesLen, ans, obj, gobj, solvopt_options, pen, gpen, &pass);
			matrix_vector_mult(A, ans, bn, constraintsLen, objectivesLen, 1.0);
			for (i = 0; i < constraintsLen; i++)
				bn[i] -= b[i];
			max = bn[0];
			for (i = 1; i < constraintsLen; i++)
				if (bn[i] > max)
					max = bn[i];
		}
		SharedFunctionVariablesFree(&pass);
		free(workingMatrix);
		return 0;
	}
	free(workingMatrix);
	return 0;
}
void quadprog(int problemLen, double *H, double *f, int inequalityLen, double *A, double *b, int equalityLen, double *Aeq, double *beq, double *lb, double *ub)
{
	int i;
	int lowerBoundLen = 0;
	int upperBoundLen = 0;
	if (lb)
		lowerBoundLen = problemLen;
	if (ub)
		upperBoundLen = problemLen;
	if (!H || !f)
	{
		printf("Fatal error\n");
		return;
	}
	if ((!A || !b) && inequalityLen)
	{
		printf("wtf with inequalityLen greater than 0 while the inequality matrix is empty\nCheck input argument, now assuming no inequality constraints\n");
		inequalityLen = 0;
	}
	if ((!Aeq || !beq) && equalityLen)
	{
		printf("wtf with equalityLen greater than 0 while the equality matrix is empty\nCheck input argument, now assuming no equality constraints\n");
		equalityLen = 0;
	}
	int constraintLength = inequalityLen + (equalityLen << 1) + lowerBoundLen + upperBoundLen;
	double *synthesizedA = (double*)malloc(constraintLength * problemLen * sizeof(double));
	double *synthesizedb = (double*)malloc(constraintLength * sizeof(double));
	if (inequalityLen)
	{
		memcpy(synthesizedA, A, inequalityLen * problemLen * sizeof(double));
		memcpy(synthesizedb, b, inequalityLen * sizeof(double));
	}
	if (equalityLen)
	{
		for (i = 0; i < equalityLen; i++)
		{
			synthesizedb[inequalityLen + i] = beq[i];
			synthesizedb[inequalityLen + equalityLen + i] = -beq[i];
			for (int j = 0; j < problemLen; j++)
			{
				synthesizedA[inequalityLen * problemLen + i * problemLen + j] = Aeq[i * problemLen + j];
				synthesizedA[inequalityLen * problemLen + equalityLen * problemLen + i * problemLen + j] = -Aeq[i * problemLen + j];
			}
		}
	}
	if (lb)
	{
		memset(synthesizedA + inequalityLen * problemLen + (equalityLen << 1) * problemLen, 0, lowerBoundLen * lowerBoundLen * sizeof(double));
		for (i = 0; i < lowerBoundLen; i++)
		{
			synthesizedA[inequalityLen * problemLen + (equalityLen << 1) * problemLen + i * lowerBoundLen + i] = -1.0;
			synthesizedb[inequalityLen + (equalityLen << 1) + i] = -lb[i];
		}
	}
	if (ub)
	{
		memset(synthesizedA + inequalityLen * problemLen + (equalityLen << 1) * problemLen + lowerBoundLen * lowerBoundLen, 0, upperBoundLen * upperBoundLen * sizeof(double));
		for (i = 0; i < upperBoundLen; i++)
		{
			synthesizedA[inequalityLen * problemLen + (equalityLen << 1) * problemLen + lowerBoundLen * lowerBoundLen + i * upperBoundLen + i] = 1.0;
			synthesizedb[inequalityLen + (equalityLen << 1) + lowerBoundLen + i] = ub[i];
		}
	}
	double *ans = (double*)malloc(problemLen * sizeof(double));
	double fval;
	int err = quadprog_ineq(H, f, synthesizedA, synthesizedb, 0, problemLen, constraintLength, ans, &fval);
	for (i = 0; i < problemLen; i++)
		printf("%d %1.15lf\n", i + 1, ans[i]);
	free(ans);
	if (err == -1)
		printf("Quadratic form must be symmetric.\n");
	if (err == 0)
		printf("SolvOpt find the solution\n");
	free(synthesizedA);
	free(synthesizedb);
}
#include "testingConstant.h"
int main()
{
	int i;
	init_genrand(1337ul);
	int problemLength = 11;
	int equalityLen = 3;
	int inequalityLen = 102;
	//quadprog(problemLength, Hqp, fqp, inequalityLen, Aqp, bqp, equalityLen, Aeqqp, beqqp, lbqp, ubqp);
	quadprog(problemLength, Hqp, fqp, inequalityLen, Aqp, 0, equalityLen, Aeqqp, beqqp, 0, 0);
	/*
	system("pause");*/
	return 0;
}