#include<stdio.h>
#include<math.h>

	double p_method(double pm_mat[][501],double sub_epsilon,double antip);
	double ip_method(double ipm_mat[][501],double sub_epsilon,double ip);
	void lua(double lua_mat[][501]);
	void det(double d_mat[][501]);
	int max(int x1,int x2,int x3);
	int min(int x1,int x2,int x3);

int main() 
{
	int k,t;
	int r=2,s=2;
	double a[5][501]={0};
	int m,n;
	double p;
	double absemin; 

	double epsilon=1e-12;

	m=0;
	for(n=2;n<501;n++) a[m][n]=-0.064;
	m++;
	for(n=1;n<501;n++) a[m][n]=0.16;
	m++;
	for(n=0;n<501;n++) 
	{
		t=n+1;
		a[m][n]=(1.64-0.024*t)*sin(0.2*t)-0.64*exp(0.1/t);
	}
	m++;
	for(n=0;n<500;n++) a[m][n]=0.16;
	m++;
	for(n=0;n<499;n++) a[m][n]=-0.064;

	double e[501]={0};
	
	det(a);

	double cond_a2;

	e[0]=p_method(a,epsilon,0);
	double e1=e[0];

	e[500]=p_method(a,epsilon,e[0]);

	if(e[0]*e[500]>0)
		absemin=e[500];
    else
		absemin=ip_method(a,epsilon,0);

	if((e[0]>0&&e[500]>0)|| (e[0]>0&&e[500]<0))
	{
		double etemp;
		etemp=e[0];
		e[0]=e[500];
		e[500]=etemp;
	}
	printf("最小特征值λ1=%17.11e\n\n",e[0]);
	printf("最大特征值λ501=%17.11e\n\n",e[500]);
	printf("按模最小特征值λs=%17.11e\n",absemin);

    for(k=1;k<=39;k++)
	{
		p=e[0]+k*(e[500]-e[0])/40;
		e[k]=ip_method(a,epsilon,p);
		printf("λi%d的近似值=%lf;其附近特征值=%17.11e\n",k,p,e[k]);
	}

	cond_a2=fabs(e1/absemin);
	printf("\nA的谱范数条件数cond(A)=%17.11e\n\n",cond_a2);

}

double p_method(double pm_mat[][501],double sub_epsilon,double antip)//带位移幂法
{
	int i,j,t,k;
	int r=2,s=2;
	int n=501;
	double beta,eta,betak_1=0;
	double pm_mt[5][501];

	for(i=0;i<5;i++)
	{
		for(j=0;j<501;j++)
			pm_mt[i][j]=pm_mat[i][j];
		if(i==2)
			for(j=0;j<n;j++) pm_mt[i][j]=pm_mt[i][j]-antip;
	}

	double u[501];
	for(t=0;t<501;t++) u[t]=0.5;

	double y[501]={0};

	for(k=1;;k++)
	{
		double s_sum=0;
		for(t=0;t<501;t++) s_sum=s_sum+u[t]*u[t];
		eta=sqrt(s_sum);
		for(t=0;t<501;t++) y[t]=u[t]/eta;

		for(i=0;i<501;i++)
		{
			double u_sum=0;
			if(i<=r)
			{
				for(j=0;j<=i+s;j++) u_sum=u_sum+pm_mt[i-j+s][j]*y[j];
				u[i]=u_sum;
			}
			else if((i>r)&&(i<(500-s)))
			{
				for(j=i-2;j<=i+2;j++) u_sum=u_sum+pm_mt[i-j+s][j]*y[j];
				u[i]=u_sum;
			}
			else if(i>=(500-s))
			{
				for(j=i-2;j<=500;j++) u_sum=u_sum+pm_mt[i-j+s][j]*y[j];
				u[i]=u_sum;
			}
		}

		double b_sum=0;
		for(t=0;t<501;t++) b_sum=b_sum+y[t]*u[t];
		beta=b_sum;

		if((fabs(beta-betak_1)/fabs(beta))<=sub_epsilon) break;
		else betak_1=beta;
	}

	beta=beta+antip;
	return beta;
}

double ip_method(double ipm_mat[][501],double sub_epsilon,double ip)
{
	int i,j,t,k;
	int r=2,s=2;
	int n=501;
	double beta,eta,betak_1=0;
	double ipm_mt[5][501]={0};
	double b[501]={0};

	for(i=0;i<5;i++)
	{
		for(j=0;j<501;j++)
			ipm_mt[i][j]=ipm_mat[i][j];
		if(i==2)
			for(j=0;j<n;j++) ipm_mt[i][j]=ipm_mt[i][j]-ip;
	}

	double u[501];
	for(t=0;t<501;t++) u[t]=0.5;

	double y[501]={0};

    lua(ipm_mt);

	for(k=1;;k++)
	{
		double s_sum=0;
		for(t=0;t<501;t++) s_sum=s_sum+u[t]*u[t];
		eta=sqrt(s_sum);
		for(t=0;t<501;t++) y[t]=u[t]/eta;

		for(i=0;i<n;i++) b[i]=y[i];      
		for(i=1;i<=n-1;i++)              
		{
		double sum=0;
		for(t=max(0,i-r,0);t<=i-1;t++)
			sum=sum+ipm_mt[i-t+s][t]*b[t];
		b[i]=b[i]-sum;
		}

	    u[n-1]=b[n-1]/ipm_mt[s][n-1];   
	    for(i=n-2;i>=0;i--)
		{
		double sum=0;
		for(t=i+1;t<=min(i+s,n-1,n-1);t++)
			sum=sum+ipm_mt[i-t+s][t]*u[t];
		u[i]=(b[i]-sum)/ipm_mt[s][i];
		}

		double b_sum=0;
		for(t=0;t<501;t++) b_sum=b_sum+y[t]*u[t];

		beta=1/b_sum+ip;

		if((fabs((beta-betak_1)/beta))<=sub_epsilon) break;
		else betak_1=beta;
	}

	return beta;
}
	
void lua(double lua_mat[][501]) 
{
	int i,j,t,k;
	int r=2,s=2;
	int n=501;

	for(k=0;k<=n-1;k++)
	{
		for(j=k;j<=min(k+s,n-1,n-1);j++)  
		{
			double sum=0;
			for(t=max(0,k-r,j-s);t<=k-1;t++)
				sum=sum+lua_mat[k-t+s][t]*lua_mat[t-j+s][j];
			lua_mat[k-j+s][j]=lua_mat[k-j+s][j]-sum;

		} 

		if(k<(n-1)) 
		{
			for(i=k+1;i<=min(k+r,n-1,n-1);i++)
			{
				double sum=0;
			    for(t=max(0,i-r,k-s);t<=k-1;t++)
				    sum=sum+lua_mat[i-t+s][t]*lua_mat[t-k+s][k];
			    lua_mat[i-k+s][k]=(lua_mat[i-k+s][k]-sum)/lua_mat[s][k];
			}
		}
	}
}

void det(double d_mat[][501]) 
{
	int i,j;
	double d_mt[5][501]={0};
	int n=501;
	for(i=0;i<5;i++)
		for(j=0;j<n;j++)
			d_mt[i][j]=d_mat[i][j];

	lua(d_mt);

	double det_prd=1;
	for(i=0;i<501;i++)
		det_prd=det_prd*d_mt[2][i];
	printf("A的行列式值detA=%17.11e\n",det_prd);
}

int max(int x1,int x2,int x3) 
{
	int tx;
	tx=x1;
	if(tx<x2) tx=x2;
	if(tx<x3) tx=x3;
	return tx;
}

int min(int x1,int x2,int x3)
{
	int tx;
	tx=x1;
	if(tx>x2) tx=x2;
	if(tx>x3) tx=x3;
	return tx;
}

