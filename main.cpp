#include <iostream>
#include <fstream>
#include <string>

#include "spectral_density_functions.h"

int main ()
{	
std::ifstream file("complexsp.txt");
std::string line, type, comments, number;
int spnumb, numbofpoints, initfreq;
double freqsteps, a, b, c, d;
int i =1;

getline(file, line);
file >> spnumb;
getline(file, line);
getline(file, line);
file >> numbofpoints;
getline(file, line);
getline(file, line);
file >> initfreq;
getline(file, line);
getline(file, line);
file >> freqsteps;
getline(file, line);

std::cout << "number of spectral densities :" <<std::endl;
std::cout << spnumb <<std::endl;
std::cout << "number of points :" <<std::endl;
std::cout << numbofpoints <<std::endl;
std::cout << "initial frequency :" <<std::endl;
std::cout << initfreq <<std::endl;
std::cout << "frequency step :" <<std::endl;
std::cout << freqsteps <<std::endl;

double * freqarr = new double[numbofpoints];
freqarr[0]=initfreq;
for (int k=1; k<numbofpoints; k++)
{
 freqarr[k]=freqarr[k-1]+freqsteps;
}

while(i++ < spnumb+1)
{
getline(file, line);
type = line;
file >> a >> b >> c >> d;
getline(file, line);

double * spectr = new double[numbofpoints];

std::cout << "calculating spectral density..." <<std::endl;
type_of_sp_density (spectr, type, freqarr, numbofpoints, a, b, c, d);

number = Stringtonumber(i);
std::string filename = type+"_spectral_density_"+number+".txt";

std::cout << "making spectral density file:"<<std::endl;
std::cout << filename <<std::endl;
spectral_density_file (filename, spectr, numbofpoints, initfreq, freqsteps);

delete[] spectr;
}
std::cout << "spectral density files has been made" <<std::endl;
file.close();
return 0;
}
