#include<iostream>
#include<cmath>

using fptr = double(double);

double f(double x);
double gaussian(double x);
double deriv(double x);
double bisection( double xl, double xu, double eps, fptr fun );
double newton(double x0, double eps, fptr fun, fptr fderiv);

int main(int argc, char **argv)
{
  double xl=-1.0;
  double xu= 1.0;
  double eps= 1.0e-9;
  double inversa_bisection= bisection(xl,xu,eps,f);
  double x0= 0.0;
  double inversa_newton= newton(x0,eps,f,deriv);
  std::cout<< inversa_bisection<<"  "<< f(inversa_bisection)<<"\n";
  std::cout<< inversa_newton<<"  "<< f(inversa_newton)<<"\n";
  
  
  return 0;
  
}

double gaussian( double x)
{
  return  (2/M_PI)*std::exp(-(x*x)/2);
  
}

double deriv(double x)
{
  return 2*std::exp(-x*x)/std::sqrt(M_PI);
  // double h = 0.0001;
  //return (f(x+h/2) - f(x-h/2))/h;
}   

double bisection( double xl, double xu, double eps, fptr fun )
{
  const int NITERMAX = 10000;
  double xr = 0.0;
  int niter = 0;
  while (niter <= NITERMAX) {
        xr = 0.5*(xl + xu);
        if (std::fabs(fun(xr)) <= eps) {
            break;
        } else if (fun(xr)*fun(xu) > 0) {
            xu = xr;
        } else {
            xl = xr;
        }
        niter++;
    }
    std::cout << "biseccion Info -> Niter: " << niter << "\n";
    return xr;
}

double newton(double x0, double eps, fptr fun, fptr fderiv)
{
    const int NITERMAX = 1000;
    double xr = x0;
    int niter = 0;
    while (niter <= NITERMAX) {
        xr = xr - fun(xr)/fderiv(xr);
        if (std::fabs(fun(xr)) <= eps) {
            break;
        } 
        niter++;
    }
    std::cout << "Newton Info -> Niter: " << niter << "\n";
    return xr;
}

double f(double x)
{

  double F= (std::erf(x))- gaussian(x);
  return F;
  
}
