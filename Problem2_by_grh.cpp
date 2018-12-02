#include<bits/stdc++.h>

const int N = 30;
const int M = 1000010;
const double eps = 1e-12;
double a[N][N];
const int n = 10;

struct Matrix {
	double a[N][N];
	Matrix() {
		memset(a, 0, sizeof a);
	}
	void out() {
		for(int i = 1; i <= n; ++i)
		for(int j = 1; j <= n; ++j)
			printf("%1.11E%c", a[i][j], "\t\n"[j == n]);
		puts("");
	}
} A0, A0_inv;

Matrix operator *(const Matrix &a, const Matrix &b) {
	Matrix ret;
	for(int i = 1; i <= n; ++i)
	for(int j = 1; j <= n; ++j)
	for(int k = 1; k <= n; ++k)
		ret.a[i][j] += a.a[i][k] * b.a[k][j];
	return ret;
}

struct Vector {
	double a[N];
	Vector() {
		memset(a, 0, sizeof a);
	}
	void out() {
		for(int i = 1; i <= n; ++i)
			printf("%1.11E%c", a[i], ",\n"[i == n]);
		putchar('\n');
	}
} ;

Vector operator /(const Vector &a, double k) {
	Vector ret;
	for(int i = 1; i <= n; ++i)
		ret.a[i] = a.a[i] / k;
	return ret;
}

double dot(const Vector &a, const Vector &b) {
	double ret = 0.;
	for(int i = 1; i <= n; ++i)
		ret += a.a[i] * b.a[i];
	return ret;
}

Vector operator *(const Matrix &a, const Vector &b) {
	Vector ret;
	for(int i = 1; i <= n; ++i) {
		double t = 0.;
		for(int j = 1; j <= n; ++j)
			t += a.a[i][j] * b.a[j];
		ret.a[i] = t;
	}
	return ret;
}

void Power() {
	static Vector u[M], y[M];
	static Matrix A = A0;
	static double yita[M], beita[M];
	for(int i = 1; i <= n; ++i)
		u[0].a[i] = rand() * 1.0;
	int k = 0;
	while(!k || fabs((beita[k] - beita[k - 1]) / beita[k]) > eps)  {
		yita[k] = sqrt(dot(u[k], u[k]));
		y[k] = u[k] / yita[k];
		u[k + 1] = A * y[k];
		beita[k + 1] = dot(y[k], u[k + 1]);
		++k;
	}
	printf("幂法 迭代步数: %d\n特征值：%1.11E\n", k, beita[k - 1]);
	printf("特征向量: "); 
	y[k - 1].out();
}

Matrix inverse(Matrix A) {
	Matrix ret;
	static double a[N][N << 1];
	memset(a, 0, sizeof a);
	for(int i = 1; i <= n; ++i) {
		for(int j = 1; j <= n; ++j)
			a[i][j] = A.a[i][j];
		a[i][i + n] = 1.0;
	}
	int m = n + n;
	for(int t = 1; t <= n; ++t) {
		for(int i = t + 1; i <= n; ++i) {
			if(a[i][t] == 0.0) continue;
			double k = a[i][t] / a[t][t];
			for(int j = t; j <= m; ++j)
				a[i][j] -= a[t][j] * k;
		}
	}
	for(int t = n; t > 1; --t)
	for(int i = 1; i < t; ++i) {
		double k = a[i][t] / a[t][t];
		for(int j = t; j <= m; ++j)
			a[i][j] -= a[t][j] * k;
	}
	for(int i = 1; i <= n; ++i) {
		double k = 1 / a[i][i];
		for(int j = n + 1; j <= m; ++j)
			a[i][j] *= k;
	}
	for(int i = 1; i <= n; ++i)
	for(int j = 1; j <= n; ++j)
		ret.a[i][j] = a[i][j + n];
	return ret;
}

void iPower() {
	static Vector u[M], y[M];
	static Matrix A_inv = A0_inv;
	static double yita[M], beita[M];
	for(int i = 1; i <= n; ++i)
		u[0].a[i] = rand() * 1.0;
	int k = 0;
	while(!k || fabs((beita[k] - beita[k - 1]) / beita[k]) > eps)  {
		yita[k] = sqrt(dot(u[k], u[k]));
		y[k] = u[k] / yita[k];
		u[k + 1] = A_inv * y[k];
		beita[k + 1] = dot(y[k], u[k + 1]);
		++k;
	}
	printf("反幂法 迭代步数: %d\n", k);
	printf("特征值 %1.11E\n", 1 / beita[k - 1]);
	printf("特征向量: ");
	y[k - 1].out();
}

Matrix Quasi_upper_tri(Matrix A) {
	for(int j = 1; j <= 8; ++j) {
		int m = 0;
		for(int i = j + 2; i <= 10; ++i) {
			m += A.a[i][j] != 0.;
		}
		if(!m) continue;
		Vector U, P, Q, W;
		double s = 0.;
		for(int i = j + 1; i <= n; ++i)
			s += A.a[i][j] * A.a[i][j];
		s = sqrt(s);
		double c = A.a[j + 1][j] <= 0 ? s : -s;
		double h = c * (c - A.a[j + 1][j]);
		U.a[j + 1] = A.a[j + 1][j] - c;
		for(int i = j + 2; i <= n; ++i)
			U.a[i] = A.a[i][j];
		for(int i = 1; i <= n; ++i) {
			for(int k = 1; k <= n; ++k)
				P.a[i] += U.a[k] * A.a[k][i];
			for(int k = 1; k <= n; ++k)
				Q.a[i] += U.a[k] * A.a[i][k];
			P.a[i] /= h;
			Q.a[i] /= h;
		}
		double t = dot(P, U) / h;
		
		for(int i = 1; i <= n; ++i) {
			W.a[i] = Q.a[i] - t * U.a[i];
			for(int k = 1; k <= n; ++k) {
				A.a[i][k] -= W.a[i] * U.a[k] + U.a[i] * P.a[k];
				if(fabs(A.a[i][k]) < eps) A.a[i][k] = 0.;
			}
		}
	}
	printf("拟上三角化 A(n-1): \n");
	A.out();
	return A;
}

bool Continue(const Matrix &A, const Matrix &B) {
	for(int i = 1; i <= n; ++i)
		if(fabs((A.a[i][i] - B.a[i][i])) > eps)
			return true;
	return false;
}

void QR(Matrix A, Matrix &Q, Matrix &R) {
	for(int i = 1; i <= n; ++i) Q.a[i][i] = 1.;
	for(int r = 1; r < n; ++r) {
		int m = 0;
		for(int i = r + 1; i <= n; ++i)
			m += A.a[i][r] != 0.;
		if(!m) continue;
		double d = 0.;
		for(int i = r; i <= n; ++i)
			d += A.a[i][r] * A.a[i][r];
		d = sqrt(d);
		double c = A.a[r][r] <= 0 ? d : -d;
		double h = c * (c - A.a[r][r]);
		Vector U, W, P;
		U.a[r] = A.a[r][r] - c;
		for(int i = r + 1; i <= n; ++i)
			U.a[i] = A.a[i][r];
		
		for(int i = 1; i <= n; ++i) {
			for(int j = 1; j <= n; ++j)
				W.a[i] += Q.a[i][j] * U.a[j];
		}
		for(int i = 1; i <= n; ++i) {
			for(int j = 1; j <= n; ++j)
				Q.a[i][j] -= W.a[i] * U.a[j] / h;
		}
		for(int i = 1; i <= n; ++i) {
			for(int j = 1; j <= n; ++j)
				P.a[i] += A.a[j][i] * U.a[j];
			P.a[i] /= h;
		}
		for(int i = 1; i <= n; ++i) {
			for(int j = 1; j <= n; ++j) {
				A.a[i][j] -= U.a[i] * P.a[j];
				if(fabs(A.a[i][j]) < eps) A.a[i][j] = 0.;
			}
		}
	}
	R = A;
}

void Calc(double a, double b, double c, double d) {
	double A = 1, B = -a - d, C = a * d - b * c;
	double imag = sqrt(-(B * B - 4 * A * C)) / 2.;
	double real = -B / 2;
	printf("%1.11E + %1.11Ei\n", real, imag); 
	printf("%1.11E - %1.11Ei\n", real, imag); 
}

void Gauss(Matrix A, double r) {
	static double a[N][N];
	memset(a, 0, sizeof a);
	for(int i = 1; i <= n; ++i) {
		for(int j = 1; j <= n; ++j)
			a[i][j] = A.a[i][j];
		a[i][i] -= r;
	}
	for(int t = 1; t <= n; ++t) {
		if(fabs(a[t][t]) <= eps) {
			for(int i = t + 1; i <= n; ++i)
			if(fabs(a[i][t]) > eps) {
				for(int j = t; j <= n; ++j)
					std::swap(a[t][j], a[i][j]);
				break;
			}
		}
		for(int i = t + 1; i <= n; ++i) {
			if(a[i][t] == 0.0) continue;
			double k = a[i][t] / a[t][t];
			for(int j = t; j <= n; ++j) {
				a[i][j] -= a[t][j] * k;
				if(fabs(a[i][j]) <= eps) a[i][j] = 0.;
			}
		}
	}
	int m = n + 1;
	a[n][n] = a[n][m] = 1.;
	for(int t = n; t > 1; --t) {
		double k = 1. / a[t][t];
		a[t][t] *= k;
		a[t][m] *= k;
		for(int i = 1; i < t; ++i) {
			double k = a[i][t] / a[t][t];
			for(int j = t; j <= m; ++j)
				a[i][j] -= a[t][j] * k;
		}
	}
	for(int i = 1; i <= n; ++i)
		printf("%1.11E%c", a[i][m], ",\n"[i == n]);
	puts("");
}

void LoopQR(Matrix &A, double *r) {
	Matrix prev = A, mat = A;
	int cnt = 0;
	do {
		Matrix Q, R;
		QR(A, Q, R);
		A = R * Q;
		++cnt;
	} while(cnt < 100000 && Continue(prev, A));
	cnt = 0;
	puts("虚特征值有:"); 
	for(int i = 1; i <= n; ++i) {
		if(fabs(A.a[i + 1][i] > eps)) {
			Calc(A.a[i][i], A.a[i][i + 1], A.a[i + 1][i], A.a[i + 1][i + 1]);
			++i;
		} else r[++cnt] = A.a[i][i];
	}
	puts("实特征值有:");
	for(int i = 1; i <= cnt; ++i) {
		printf("%1.11E\n", r[i]);
	}
	for(int i = 1; i <= cnt; ++i) {
		printf("实特征值 %1.11E 对应的特征向量:\n", r[i]);
		Gauss(mat, r[i]);
	}
}

int main() {
	freopen("2.output", "w", stdout); 
	for(int i = 1; i <= n; ++i)
	for(int j = 1; j <= n; ++j) {
		A0.a[i][j] = i != j ? sin(0.5 * i + 0.2 * j) : 1.52 * cos(i + 1.2 * j);
	}
	A0_inv = inverse(A0);
	Power();
	iPower();
	Matrix A = Quasi_upper_tri(A0);
	double r[N];
	LoopQR(A, r);
	return 0;
}

程序输出如下： 
幂法 迭代步数: 96
特征值：3.38961343881E+000
特征向量: 1.04871999320E-001,2.17676976319E-001,4.74694012242E-001,2.59383624651E-001,3.04665248521E-001,2.59451746662E-001,-8.68664182726E-002,-4.05258126692E-001,-5.09628289643E-001,-2.39514692166E-001

反幂法 迭代步数: 12
特征值 4.95499092363E-002
特征向量: 2.13767977959E-001,2.06773621699E-001,-3.86828983510E-001,3.11123946362E-002,3.80938960236E-001,1.25173726812E-001,-6.44715735839E-001,3.08201272967E-001,2.95976727013E-001,-4.37229510135E-002

拟上三角化 A(n-1): 
-8.94521698228E-001	-9.93313649183E-002	-1.09983175888E+000	-7.66503870908E-001	1.70760114146E-001	-1.93488255889E+000	-8.39020870525E-002	9.13256511314E-001	-6.40797700919E-001	1.94673367868E-001
-2.34787836242E+000	2.37205792160E+000	1.82799855232E+000	3.26655688471E-001	2.08236058364E-001	2.08898700994E+000	1.84786191029E-001	-1.26301526608E+000	6.79069466850E-001	-4.67215088650E-001
0.00000000000E+000	1.73595446995E+000	-1.16502336748E+000	-1.24674444352E+000	-6.29822548908E-001	-1.98482018099E+000	2.97575006080E-001	6.33930059659E-001	-1.30851892877E-001	3.04030103610E-001
0.00000000000E+000	0.00000000000E+000	-1.29293756392E+000	-1.12623922590E+000	1.19078291192E+000	-1.30877298390E+000	1.86015166267E-001	4.23673393688E-001	-1.01960082655E-001	1.94366091451E-001
0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	1.57771115303E+000	8.16935832816E-001	4.46153172383E-001	-4.36509254161E-002	-4.66597916719E-001	2.94123156618E-001	-1.03442111366E-001
0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	-7.72897513499E-001	-1.60102824405E+000	-2.91268547483E-001	-2.43433785832E-001	6.73628608451E-001	2.62477290494E-001
0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	-7.29677394636E-001	-7.96545627982E-003	9.71073910201E-001	-1.29896736857E-001	2.78024208124E-002
0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	7.94553961298E-001	-4.52514345461E-001	5.04890152758E-001	-1.21121019351E-001
0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	7.03991137351E-001	1.26753552350E-001	-3.71469673551E-001
0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	0.00000000000E+000	-4.91958687221E-001	4.08150976640E-001

虚特征值有:
-2.33686593224E+000 + 8.93437921022E-001i
-2.33686593224E+000 - 8.93437921022E-001i
-9.89114346473E-001 + 1.08475863151E-001i
-9.89114346473E-001 - 1.08475863151E-001i
实特征值有:
3.38961343882E+000
1.59031345881E+000
-1.49314708091E+000
9.43287957277E-001
6.48948820211E-001
4.95499092362E-002

实特征值 3.38961343882E+000 对应的特征向量:
1.99983436817E+004,4.01689189193E+004,1.72052065162E+004,-6.04838571264E+003,-3.60453538921E+003,5.66355647144E+002,-1.29270201738E+002,-2.75610573222E+001,-6.06039194218E+000,1.00000000000E+000

实特征值 1.59031345881E+000 对应的特征向量:
1.24327530147E+001,-2.78128607069E+001,-1.25447464607E-001,-3.63537417195E+001,-6.27290100308E+001,1.60030740169E+001,-9.80797852888E+000,-4.46798580484E+000,-2.40297104782E+000,1.00000000000E+000

实特征值 -1.49314708091E+000 对应的特征向量:
7.54938217963E+000,4.98489836385E+000,-7.69218518915E+000,-7.85912225608E+000,2.26834045153E+000,5.82893570608E+000,8.65261832518E+000,-8.36522416402E+000,3.86475146588E+000,1.00000000000E+000

实特征值 9.43287957277E-001 对应的特征向量:
-2.64885719357E-001,3.70597383671E-001,2.22791818920E-001,1.11758282132E-001,3.30381335193E-001,-1.63935175750E-001,-4.45780050176E-001,-7.34001325819E-001,-1.08776812878E+000,1.00000000000E+000

实特征值 6.48948820211E-001 对应的特征向量:
-2.69762334705E-001,3.01764324536E-001,6.46188426614E-001,2.59596558072E-001,4.10916310578E-001,-2.78383891527E-001,6.91635224435E-001,1.64592985517E-001,-4.89467611462E-001,1.00000000000E+000

实特征值 4.95499092362E-002 对应的特征向量:
3.15983398129E-001,-1.29762490172E-001,-5.25704519239E-001,3.34078340085E-002,-2.37024609949E-001,5.06377292754E-001,-2.78391560757E-002,4.47724379626E-001,7.28925165300E-001,1.00000000000E+000


特别地， 
从网上云计算得出的矩阵特征值：
	特征值1：	3.3896
	特征值2：	1.5903
	特征值3：	0.9433
	特征值4：	0.6489
	特征值5：	0.0495
	特征值6：	-0.9891 + 0.1085i
	特征值7：	-0.9891 - 0.1085i
	特征值8：	-1.4931
	特征值9：	-2.3369 + 0.8934i
	特征值10：	-2.3369 - 0.8934i
因而可以初步验证：特征值计算基本正确。 

