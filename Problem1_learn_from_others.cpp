#include<bits/stdc++.h>

using namespace std;
const int N = 510;
const double eps = 1e-12;
const int n = 501;
double a[5][N], e[N], cond_a2, E_minnorm;

double Power(double antip) {
	static double ret[5][N], u[N], y[N];
	memcpy(ret, a, sizeof a);
	for(int j = 0; j < n; ++j)
		ret[2][j] = a[2][j] - antip;
	for(int i = 0; i < n; ++i) u[i] = 0.5;
	memset(y, 0, sizeof y);
	double beta = 0. ,prev = 0.;
	for(int k = 1; ; k++) {
		double s = 0;
		for(int i = 0; i < n; ++i) s += u[i] * u[i];
		s = sqrt(s);
		for(int i = 0; i < n; ++i) y[i] = u[i] / s;
		for(int i = 0; i < n; ++i) {
			double sum = 0.;
			for(int j = max(0, i - 2); j <= min(i + 2, 500); ++j)
				sum += ret[i - j + 2][j] * y[j];
			u[i] = sum;
		}
		beta = 0.;
		for(int i = 0; i < n; ++i) beta += y[i] * u[i];
		if(fabs((beta - prev) / beta) > eps) prev = beta;
		else break;
	}
	return beta + antip;
}

void Elimination(double mat[][N]) {
	for(int k = 0; k < n; ++k) {
		for(int j = k; j <= min(k + 2, n - 1); ++j) {
			double s = 0.;
			for(int t = max(0, max(j - 2, k - 2)); t < k; ++t)
				s += mat[k - t + 2][t] * mat[t - j + 2][j];
			mat[k - j + 2][j] -= s;
		}
		if(k < n - 1) {
			for(int i = k + 1; i <= min(k + 2, n - 1); ++i) {
				double s = 0.;
				for(int t = max(0, max(i - 2, k - 2)); t < k; ++t)
					s += mat[i - t + 2][t] * mat[t - k + 2][k];
				mat[i - k + 2][k] = (mat[i - k + 2][k] - s) / mat[2][k];
			}
		}
	}
}

double iPower(double ip) {
	static double ret[5][N], u[N], y[N], b[N];
	double beta, prev = 0.;
	memcpy(ret, a, sizeof a);
	for(int j = 0; j < n; ++j)
		ret[2][j] = a[2][j] - ip;
	for(int i = 0; i < n; ++i) u[i] = 0.5;
	memset(y, 0, sizeof y);
	
	Elimination(ret);
	
	for(int k = 1; ; k++){
		double s = 0;
		for(int i = 0; i < n; ++i) s += u[i] * u[i];
		s = sqrt(s);
		for(int i = 0; i < n; ++i) b[i] = y[i] = u[i] / s;
		for(int i = 1; i < n; ++i) {
			double sum = 0.;
			for(int j = max(0, i - 2); j < i; ++j)
				sum += ret[i - j + 2][j] * b[j];
			b[i] -= sum;
		}
		u[n - 1] = b[n - 1] / ret[2][n - 1];
		for(int i = n - 2; i >= 0; --i) {
			double sum = 0.;
			for(int j = i + 1; j <= min(i + 2, n - 1); ++j)
				sum += ret[i - j + 2][j] * u[j];
			u[i] = (b[i] - sum) / ret[2][i];
		}
		beta = 0.;
		for(int i = 0; i < n; ++i) beta += y[i] * u[i];
		beta = 1 / beta + ip;
		if(fabs((beta - prev) / beta) > eps) prev = beta;
		else break;
	}
	return beta;
}

double det(double mat[][N]) {
	static double a[5][N];
	memcpy(a, mat, sizeof a);
	Elimination(a);
	double ret = 1.0;
	for(int i = 0; i < n; ++i) ret *= a[2][i];
	return ret;
}


int main() {
	freopen("Problem1.out", "w", stdout);
	for(int i = 2; i < n; ++i) a[0][i] = -0.064;
	for(int i = 1; i < n; ++i) a[1][i] = 0.16;
	for(int i = 1; i <= n; ++i) {
		a[2][i - 1] = (1.64 - 0.024 * i) * sin(0.2 * i) - 
			0.64 * exp(0.1 / (double)i);
	}
	for(int i = 0; i < n - 1; ++i) a[3][i] = 0.16;
	for(int i = 0; i < n - 2; ++i) a[4][i] = -0.064;
	e[0] = Power(0);
	e[n - 1] = Power(e[0]);
	double tmp_e = e[0];
	E_minnorm = e[0] * e[n - 1] > 0. ? e[n - 1] :
		iPower(0);
	if(e[0] > 0. && e[n - 1] != 0.) swap(e[0], e[500]);
	printf("最小特征值: %1.11E\n\n", e[0]);
	printf("最大特征值: %1.11E\n\n", e[n - 1]);
	printf("按模最小特征值: %1.11E\n\n",E_minnorm);
	for(int k = 1; k < 40; ++k) {
		double p = e[0] + k * (e[n - 1] - e[0]) / 40.;
		e[k] = iPower(p);
		printf("λ[i, %d]: %17.11E\n",k, e[k]);
	}
	puts("");
	cond_a2 = fabs(tmp_e / E_minnorm);
	printf("cond(A)_2: %1.11E\n", cond_a2);
	printf("det(A): %1.11E\n", det(a));
	return 0;
}

程序输出：
最小特征值: -1.07001136150E+001

最大特征值: 9.72463409878E+000

按模最小特征值: -5.55791079423E-003

λ[i, 1]: -1.01829340331E+001
λ[i, 2]: -9.58570742498E+000
λ[i, 3]: -9.17267242393E+000
λ[i, 4]: -8.65228400790E+000
λ[i, 5]: -8.09348380865E+000
λ[i, 6]: -7.65940540769E+000
λ[i, 7]: -7.11968464869E+000
λ[i, 8]: -6.61176433940E+000
λ[i, 9]: -6.06610322659E+000
λ[i, 10]: -5.58510105263E+000
λ[i, 11]: -5.11408352981E+000
λ[i, 12]: -4.57887217687E+000
λ[i, 13]: -4.09647092631E+000
λ[i, 14]: -3.55421121575E+000
λ[i, 15]: -3.04109001813E+000
λ[i, 16]: -2.53397031130E+000
λ[i, 17]: -2.00323076956E+000
λ[i, 18]: -1.50355761123E+000
λ[i, 19]: -9.93558606007E-001
λ[i, 20]: -4.87042673885E-001
λ[i, 21]: 2.23173624957E-002
λ[i, 22]: 5.32417474207E-001
λ[i, 23]: 1.05289896269E+000
λ[i, 24]: 1.58944588188E+000
λ[i, 25]: 2.06033046027E+000
λ[i, 26]: 2.55807559707E+000
λ[i, 27]: 3.08024050931E+000
λ[i, 28]: 3.61362086769E+000
λ[i, 29]: 4.09137851045E+000
λ[i, 30]: 4.60303537828E+000
λ[i, 31]: 5.13292428390E+000
λ[i, 32]: 5.59490634808E+000
λ[i, 33]: 6.08093385702E+000
λ[i, 34]: 6.68035409211E+000
λ[i, 35]: 7.29387744814E+000
λ[i, 36]: 7.71711171424E+000
λ[i, 37]: 8.22522001405E+000
λ[i, 38]: 8.64866606518E+000
λ[i, 39]: 9.25420034458E+000

cond(A)_2: 1.92520427390E+003
det(A): 2.77278614175E+118
 
