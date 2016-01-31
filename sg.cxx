#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
const double pi=M_PI;
const cmplx im= cmplx(0,1.0);
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);


void step(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          const double h, const double m, const double k,const double x, const int N);
//-----------------------------------
int main(){

	const double h=1, m=1;
	const int Nx = 300;
	const double xmin = -40;
  const double xmax = 40;
	const double Tend = 10*pi;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = dx/10.0  ;
  double t = 0;
	const double omega=0.2;
	const double k=omega*omega*m;
	const double alpha=pow(m*k/h,1.0/4.0);
	const int Na = 10;
  int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
        cmplx* psi1 = new cmplx[Nx];
        cmplx* temp=NULL;

	init(psi0,alpha,lambda,dx,dt,Nx,xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
            step(psi1,psi0,dt,dx,h,m,k,xmin,Nx);
            temp=psi1; psi1=psi0; psi0=temp;
         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

        delete[] psi0;
        delete[] psi1;
	return 0;
}
//-----------------------------------
void step(cmplx* const f1, cmplx* const f0,
          const double dt, const double dx,
          const double h, const double m, const double k,const double x, const int N)
{

  cmplx* d=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* f0_prime=new cmplx[N];
    //initial conditions
  for(int i=0;i<N;i++) d[i] = 1.0 + im*(h*dt/2.0/m/dx/dx+dt/4.0/h*k*pow(x+i*dx,2));
  for(int i=0;i<N;i++) u[i] = - im*h*dt/4.0/m/dx/dx;
    //New vector conj(A)*f0
  f0_prime[0] = conj(d[0])*f0[0]+conj(u[0])*f0[1];
  for(int i=1;i<N-1;i++) f0_prime[i] = conj(u[i])*f0[i-1]+conj(d[i])*f0[i]+conj(u[i])*f0[i+1];
  f0_prime[N-1]=conj(u[N-1])*f0[N-2]+conj(d[N-1])*f0[N-1];
    //solve via backwards substitution
 for(int i=1; i<N; i++){
    d[i]-=u[i]/d[i-1]*u[i-1];
    f0_prime[i]-=u[i]/d[i-1]*f0_prime[i-1];
  }
  f1[N-1]=f0_prime[N-1]/d[N-1];
  for(int i=N-2; i>0; i--) f1[i]=(f0_prime[i]-u[i]*f1[i+1])/d[i];

  delete[] d;
  delete[] u;
  delete[] f0_prime;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
