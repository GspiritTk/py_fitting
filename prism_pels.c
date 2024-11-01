
#include <math.h>

static double
form_volume(double length, double kuhn_length, double radius)
{
    return 1.0;
}

static double cq_pel(double q, double kuhn_length, double radius) 
{
    double sigma = (exp(radius / (2.0 * kuhn_length)) - 1.0) / kuhn_length;
    double qRc = q * radius;
    return (sin(qRc) / qRc) * (exp(-pow(q, 2.0) * pow(sigma, 2.0)) - 1.0);
}


static double
Iq(double q,
   double length,
   double kuhn_length,
   double radius,
   double sld,
   double solvent_sld,
   double parav,
   double conc,
   double wamw)
{
    const double nump = 6.009956e-3 * conc / wamw;
    const double contrast = sld - solvent_sld;
    const double cross_section = sas_2J1x_x(q*radius);
    const double volume = M_PI*radius*radius*length;
    const double pwlc = Sk_WR(q, length, kuhn_length);
    const double cq = cq_pel(q, kuhn_length, radius);
    return nump * pow(volume, 2) * pow(contrast, 2) * pwlc * cross_section / (1 + parav * cq * pwlc);
}