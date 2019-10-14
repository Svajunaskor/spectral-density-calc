#include <string>

// spectral density functions
void type_of_sp_density (double * spectr, std::string type, double * freqarr, int numbofpoints, double a, double b, double c, double d);
void drude (double * spectr, double * freqarr, int numbofpoints, double a, double b);
void ohmic (double * spectr,double * freqarr, int numbofpoints, double a, double b);
void superohmic (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c);
void lognormal (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c);
void fractionalfunction (double * spectr, double * freqarr, int numbofpoints, double a, double b, double c, double d);
void Fouriertransformfunction (double * spectr, double * x0, int numbofpoints, double * freqarr, double b);
void ReorganizationEnergy (double * freqarr, int numbofpoints, double * spectr);
std::string Stringtonumber(int i);
void spectral_density_file (std::string filename, double * spectr, int numbofpoints, int initfreq, double freqsteps);


