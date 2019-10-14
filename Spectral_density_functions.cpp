#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <complex>
#include <fftw3.h>

void drude (double * spectr, double * freqarr, int numbofpoints, double a, double b);
void ohmic (double * spectr,double * freqarr, int numbofpoints, double a, double b);
void superohmic (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c);
void lognormal (double * spectr, double * freqarr, int numbofpoints, double a, double b, double d);
void fractionalfunction (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c, double d);

void type_of_sp_density (double * spectr, std::string type, double * freqarr, int numbofpoints, double a, double b, double c, double d)
{
	if (!type.compare("drude"))
	{
	std::cout << " drude spectral density " <<std::endl;
	std::cout << " parameters: " <<std::endl;
	std::cout << a << " " << b << std::endl;
	drude (spectr,freqarr, numbofpoints, a,b);		
	} 
	else if (!type.compare("ohmic"))
	{
	std::cout << " ohmic spectral density " <<std::endl;
	std::cout << " parameters: " <<std::endl;
	std::cout << a << " " << b << " " << c <<std::endl;
	ohmic (spectr, freqarr, numbofpoints, a, b);	
	}
	else if(!type.compare("superohmic"))
	{
	std::cout << " superohmic spectral density " <<std::endl;
	std::cout << " parameters: " <<std::endl;
	std::cout << a << " " << b << " " << c <<std::endl;
	superohmic (spectr, freqarr, numbofpoints, a, b, c);		
	}
	else if (!type.compare("lognormal"))
	{
	std::cout << " lognormal spectral density " <<std::endl;
	std::cout << " parameters: " <<std::endl;
	std::cout << a << " " << b << " " << c <<std::endl;
	lognormal (spectr, freqarr, numbofpoints, a, b, c);
    }
	else if (!type.compare("fractionalfunction"))
	{
	std::cout << " fractional function spectral density " <<std::endl;
	std::cout << " parameters: " <<std::endl;
	std::cout << a << " " << b << " " << c << " " << d <<std::endl;
	fractionalfunction(spectr, freqarr, numbofpoints, a, b, c, d);
	}
	else std::cout << "no spectral density with this name" <<std::endl;
}

void Fouriertransformfunction (double * spectr, double * x0, int numbofpoints, double * freqarr, double b)
{
 int n = 2*numbofpoints-1;
 //input
 fftw_complex * x = new fftw_complex[n];
 // output
 fftw_complex * y = new fftw_complex[n];
 for (int p = 0; p < n; p++)
 {
  x[p][0]= x0[p];
  x[p][1] = 0;
 }
 //plan for FFT
 fftw_plan plan =fftw_plan_dft_1d(n, x, y, FFTW_BACKWARD, FFTW_ESTIMATE);
 fftw_execute(plan);
 fftw_destroy_plan(plan);
 fftw_cleanup();
 for (int i =0; i<n; i++)
 {
  y[i][0]=y[i][0]*(((2*M_PI)/((freqarr[1]-freqarr[0])*numbofpoints)));
 }
 for (int p =0; p<numbofpoints; p++)
 {
  spectr[p]=(freqarr[p]/(2*b))*y[p][0];
 }
 delete[] x;
 delete[] y;
}

// reorgnization energy calc
void ReorganizationEnergy (double * freqarr, int numbofpoints, double * spectr)
{
 double reorgenergy = 0;
 for (int reorgn = 0; reorgn<numbofpoints-1; reorgn++)
 {
 	if (freqarr[reorgn] > 0)
 	{
 	 reorgenergy += (freqarr[reorgn+1]-freqarr[reorgn])*(((spectr[reorgn]/freqarr[reorgn])+(spectr[reorgn+1]/freqarr[reorgn+1]))/2);
	}
	else
	{
	 reorgenergy += (freqarr[reorgn+1]-freqarr[reorgn])*((spectr[reorgn+1]/freqarr[reorgn+1])/2);	
	}
 }
 reorgenergy =reorgenergy/M_PI;
 std::cout << "Reorgnization Energy: " << reorgenergy <<std::endl;
}

// drude
void drude (double * spectr, double * freqarr, int numbofpoints, double a, double b)
{
for (int l=0; l<numbofpoints; l++)
{
spectr[l]=(2*b*freqarr[l]*a)/(pow(freqarr[l],2)+pow(a,2));
}
ReorganizationEnergy(freqarr,numbofpoints, spectr);		
}

// ohmic function
void ohmic (double * spectr,double * freqarr, int numbofpoints, double a, double b)
{
for (int l=0; l<numbofpoints; l++)
{
spectr[l]=M_PI*b*freqarr[l]*exp(-freqarr[l]/a);	
}
ReorganizationEnergy(freqarr,numbofpoints, spectr);		
}

// super ohmic function
void superohmic (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c)
{
if (b==2)
{
for (int l=0; l<numbofpoints; l++)
{
spectr[l]=M_PI*c*(pow(freqarr[l],2)/pow(a,1))*exp(-freqarr[l]/a);	
}		
} 
else if (b==3)
{
for (int l=0; l<numbofpoints; l++)
{
spectr[l]=M_PI*c*(pow(freqarr[l],3)/pow(a,2))*exp(-freqarr[l]/a);	
}	
}
else 
{
std::cout << "wrong alpha parameter. alpha parameter should be 2 or 3" <<std::endl;	
}
ReorganizationEnergy(freqarr,numbofpoints, spectr);		
}

// lognormal function
void lognormal (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c)
{
for (int l=0; l<numbofpoints; l++)
{
spectr[l]=((c*M_PI*freqarr[l])/(b*sqrt(2*M_PI)))*exp(-(pow(log(freqarr[l]/a),2)/(2*pow(b,2))));	
}
ReorganizationEnergy(freqarr,numbofpoints, spectr);		
}

// fractional function
void fractionalf(double * x0, int numbofpoints, double * freqarr, double a, double b, double c, double d)
{
 int n = 2*numbofpoints-1;
 double * T=new double[n];
 for (int p = 0; p <numbofpoints; p++)
 {
  T[p]=(p*(2*M_PI))/((freqarr[1]-freqarr[0])*numbofpoints);		
 }
 for (int p = 1; p<numbofpoints; p++)
 {
 x0[p]=a*b*pow((tanh(T[p]/d)/(T[p]/d)),1-c);
 }
 x0[0] = a*b;
 for (int p = 0; p<numbofpoints; p++)
 {
  x0[n-1-p]=x0[p+1];	
 }
 delete[] T;
}

void fractionalfunction (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c, double d)
{
double * x0 = new double[2*numbofpoints-1];
fractionalf(x0, numbofpoints, freqarr, a, b, c, d);
std::cout << "fractional function" <<std::endl;
std::cout << a << " " << b <<" " << c << " " << d << std::endl;
Fouriertransformfunction(spectr, x0, numbofpoints, freqarr, b);
ReorganizationEnergy(freqarr,numbofpoints, spectr);
}

std::string Stringtonumber(int i)
{
std::string number;
std::ostringstream convert;
convert << i-1;
number = convert.str();	
return number;	
}

void spectral_density_file (std::string filename, double * spectr, int numbofpoints, int initfreq, double freqsteps)
{	
std::ofstream spfile(filename.c_str());
spfile << "# Spectral density file" <<std::endl;
spfile << "# Number of points" <<std::endl;
spfile << numbofpoints <<std::endl;
spfile << "# Initial frequency" <<std::endl;
spfile << initfreq <<std::endl;
spfile << "# Frequency step" <<std::endl;
spfile << freqsteps <<std::endl;
spfile.precision(15);
for (int r=0; r<numbofpoints; r++)
{
spfile << spectr[r] << std::endl;
}
spfile.close();				
}
