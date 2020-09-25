#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>

using namespace std;

typedef double** Matrix;
typedef double* Vector;

bool WithGauss(Matrix M, Vector A, int n,double& Lmin, double& Lmax){
	int l;
	int bi,r,c,i;
	Vector Temp;
	double a1,b1,sum;
	double el;
	bool Result;

	Temp=0;
	Result=false;
	if (n==0) return Result;
	l=n+1;

	Lmin=1.0E300;
	Lmax=-1.0E300;
	for (i=0;i<n;i++){
	  bi=i;
	  for (r=i;r<n;r++)
		if (fabsl(M[r][i])>fabsl(M[bi][i])) bi=r;
		  Temp=M[bi];
		  M[bi]=M[i];
		  M[i]=Temp;
		  if (M[i][i]==0.0l) return Result;
		  el=M[i][i];
		  if (Lmin>el) {
			Lmin=el;
		  }
		  if (Lmax<el) {
			Lmax=el;
		  }
		  M[i][i]=1;
		  for(c=i+1;c<l;c++)
			if (M[i][c]==0) M[i][c]=0;
			else M[i][c]=M[i][c]/el;
		  for(r=i+1;r<n;r++){
			if (M[r][i]==0) M[r][i]=0;
			else {
			  el=M[r][i];
			  M[r][i]=0;
			  for(c=i+1;c<l;c++){
				if (M[r][c]==0) a1=0;
				else a1=M[r][c];
				if (M[i][c]==0) b1=0;
				else b1=M[i][c]*el;
				M[r][c]=a1-b1;
			  }
			}
		  }
	}

	for(i=n-1;i>=0;i--){
	  a1=0;
	  for(c=i+1;c<n;c++)
		a1=a1+M[i][c]*A[c];
	  if (a1==0) a1=0;
	  if (M[i][n]==0) A[i]=-a1;
	  else A[i]=M[i][n]-a1;
	}

	Result=true;
	return Result;
}

void Inverse(Matrix M, int n){
	int l;
	int bi,r,c,i;
	double a1,b1,sum;
	double el;
	l=2*n;
	for (i=0;i<n;i++){
		bi=i;
		el=M[i][i];
		M[i][i]=1;
		for(c=i+1;c<l;c++)
			if (M[i][c]==0) M[i][c]=0;
			else M[i][c]=M[i][c]/el;
		for(r=0;r<n;r++){
			if (r==i) continue;
			if (M[r][i]==0) M[r][i]=0;
			else {
				el=M[r][i];
				M[r][i]=0;
				for(c=i+1;c<l;c++){
					if (M[r][c]==0) a1=0;
					else a1=M[r][c];
					if (M[i][c]==0) b1=0;
					else b1=M[i][c]*el;
					M[r][c]=a1-b1;
				}
			}
		}
	}
}

int GaussSeidel(Matrix a, Vector x, int n){
	int i, j, m = 0;
	int method;
	double eps=1.0E-12;
	double* p=new double[n];
	for (int i = 0; i < n; i++) {
		x[i] = 0;
		p[i] = 0;
	}
	int iteration=0;
	while (true)
	{
		iteration++;
		for (int i = 0; i < n; i++)
		{
			x[i] = a[i][n];
			for (int j = 0; j < n; j++)
			{
				if (j < i)
				{
					x[i] -= a[i][j] * x[j];
				}
				if (j > i)
				{
					x[i] -= a[i][j] * p[j];
				}
			}
			x[i] /= a[i][i];
		}
		double error = 0.0;
		for (int i = 0; i < n; i++)
		{
			error += fabs (x[i] - p[i]);
		}
		if (error < eps)
		{
			break;
		}
		for (int i = 0; i < n; i++){
			p[i] = x[i];
		}
	}

	delete [] p;

	return iteration;
}

//double EMatrNorm(Matrix a, int n)
//{
//	double sum = 0;
//	for(int i=0;i<n;i++)
//	{
//		for(int j=0;j<n;j++)
//		{
//			sum += a[i][j]*a[i][j];
//		}
//	}
//	return sqrt(sum);
//}

int main (int argc, char **argv) {

	int n=99;
	vector<double> a(n);
	vector<double> b(n);
	vector<double> c(n);
	vector<double> f(n+1);
	vector<double> p(n+1);
	for (int i = 0; i < n; i++) {
		a[i]=1;
		b[i]=10+(i+1);
		c[i]=1;
		p[i]=2;
	}
	for (int i = 0; i < n+1; i++) {
		f[i]=(i+1)/(double)n;
	}
	p[0]=p[n]=1;
	Matrix m=new double*[n+1];
	for (int i = 0; i < n+1; i++) {
		m[i]=new double[n+2];
		for (int j = 0; j < n+2; j++) {
			m[i][j]=0;
		}
	}
	for (int i = 0; i < n+1; i++) {
		if (i<n) {
			if (i>0) {
				m[i][i-1]=a[i];
			}
			m[i][i]=b[i];
			m[i][i+1]=c[i];
		} else {
			for (int j = 0; j < n+1; j++) {
				m[i][j]=p[j];
			}
		}
		m[i][n+1]=f[i];
	}

	Matrix m2=new double*[n+1];
	for (int i = 0; i < n+1; i++) {
		m2[i]=new double[n+2];
		for (int j = 0; j < n+2; j++) {
			m2[i][j]=m[i][j];
		}
	}

	Matrix m3=new double*[n+1];
	for (int i = 0; i < n+1; i++) {
		m3[i]=new double[n+2];
		for (int j = 0; j < n+2; j++) {
			m3[i][j]=m[i][j];
		}
	}

//	for (int i = 0; i < n+1; i++) {
//		for (int j = 0; j < n+2; j++) {
//			cout<<m3[i][j]<<" ";
//		}
//		cout<<endl;
//	}

	Matrix m4=new double*[n+1];
	for (int i = 0; i < n+1; i++) {
		m4[i]=new double[2*(n+1)];
		for (int j = 0; j < n+1; j++) {
			m4[i][j]=m[i][j];
		}
		for (int j = n+1; j < 2*(n+1); j++) {
			if (j-n-1==i) {
				m4[i][j]=1;
			} else {
				m4[i][j]=0;
			}
		}
	}

	Inverse(m4,n+1);

//	cout << "Inverse matrix" << endl;
//	for (int i = 0; i < n+1; i++) {
//		for (int j = n+1; j < 2*n+1; j++) {
//			cout<<fixed<<setprecision(3)<<m4[i][j]<<" ";
//		}
//		cout<<endl;
//	}


	/*
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; j < n+1; j++) {
			double sum=0;
			for (int k = 0; k < n+1; k++) {
				sum+=m3[i][k]*m4[k][j+n+1];
			}
			cout<<sum<<" ";
		}
		cout<<endl;
	}*/

//	double sum1=0;
//	for (int i = 0; i < n+1; i++) {
//		double max=0;
//		for (int j = 0; j<n; j++) {
//			if (max<m4[i][j]) {
//				max=m4[i][j];
//			}
//		}
//		sum1+=max;
//	}
//
//	double sum2=0;
//	for (int i = 0; i < n+1; i++) {
//		double max=0;
//		for (int j = 0; j<n; j++) {
//			if (max<m4[i][j+n+1]) {
//				max=m4[i][j+n+1];
//			}
//		}
//		sum2+=max;
//	}

	double sum1=0;
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; j<n+1; j++) {
			sum1+=m3[i][j]*m3[i][j];
		}
	}
	sum1 = sqrt(sum1);

	double sum2=0;
	for (int i = 0; i < n+1; i++) {
		for (int j = 0; j<n+1; j++) {
			sum2+=m4[i][j+n+1]*m4[i][j+n+1];
		}
	}
	sum2 = sqrt(sum2);

	Vector res=new double[n+1];

	double Lmin,Lmax;
	WithGauss(m,res,n+1,Lmin,Lmax);

	Vector res2=new double[n+1];
	int iteration=GaussSeidel(m2,res2,n+1);

	Vector error2=new double[n+1];
	Vector error=new double[n+1];

	// nevazka Gauss
	double err1_norm = 0;
	for (int i = 0; i < n+1; i++) {
		double sum=0;
		for (int j = 0; j < n+1; j++) {
			sum+=m3[i][j]*res[j];
		}
		error[i]=sum-m3[i][n+1];
		err1_norm += error[i]*error[i];
	}
	err1_norm = sqrtf(err1_norm);

	// nevazka Seidel
	double err2_norm = 0;
	for (int i = 0; i < n+1; i++) {
		double sum=0;
		for (int j = 0; j < n+1; j++) {
			sum+=m3[i][j]*res2[j];
		}
		error2[i]=sum-m3[i][n+1];
		err2_norm += error2[i]*error2[i];
	}
	err2_norm = sqrtf(err2_norm);

	cout<<"Gauss,              Gauss-Seidel ("
	<<iteration<<" iterations)"<<endl;
	for (int i = 0; i < n+1; i++) {
		cout.precision(16);
		cout<<fixed;

		cout<<res[i]<<", " << res2[i] << endl;

//		if (error[i]>=0)
//			cout<<res[i]<<"+"<<error[i];
//		else
//			cout<<res[i]<<error[i];

//		cout<<", ";

//		if (error2[i]>=0)
//			cout<<res2[i]<<"+"<<error2[i];
//		else
//			cout<<res2[i]<<error2[i];
//		cout<<endl;
	}
	cout<<endl;

	cout << "Norma nevazki Gauss=" << err1_norm << endl;
	cout << "Norma nevazki Gauss-Seidel=" << err2_norm << endl;
	cout << endl;
	
	cout<<"Lambda min="<<Lmin<<endl;
	cout<<"Lambda max="<<Lmax<<endl;
	cout << endl;

	cout << "||A||= " << sum1 << endl;
	cout << "||A-1|| = " << sum2 << endl;
	cout<<"||A||*||A-1||="<<sum1*sum2<<endl;

	for (int i = 0; i < n+1; i++) {
		delete [] m[i];
		delete [] m2[i];
		delete [] m3[i];
	}

	delete [] m;
	delete [] m2;
	delete [] m3;
	delete [] res;
	delete [] res2;

	system("pause");

	return 0;
}
