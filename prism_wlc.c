
#include <math.h>

static double
form_volume(double length, double kuhn_length, double radius)
{
    return 1.0;
}

static double cq_pel(double q, double kuhn_length, double radius, double sigma, double rc) 
{
    double qRc = q * rc;
    return (sin(qRc) / qRc) * (exp(-pow(q, 2.0) * pow(sigma, 2.0)));
}

#define N 10000
double Si(double x) {
    if (x == 0.0) {
        return 0.0;
    }

    double a = 0.0;
    double b = x;
    double h = (b - a) / N;
    double sum = 0.0;
    for (int i = 0; i <= N; i++) {
        double t = a + i * h;
        double coef;

        if (i == 0 || i == N) {
            coef = 1;
        } else if (i % 2 == 0) {
            coef = 2;
        } else {
            coef = 4;
        }

        double f_t = (t == 0) ? 1 : sin(t) / t;
        sum += coef * f_t;
    }

    sum = (h / 3.0) * sum;
    return sum;
}

static double cq_rod(double q, double lr)
{
    double term1 = 2 * Si(q * lr) / (q * lr);
    double term2 = 4 * pow(sin(q * lr / 2), 2) / (pow(q, 2) * pow(lr, 2));
    return term1 - term2;

}

static double
Iq(double q,
   double length,
   double kuhn_length,
   double radius,
   double sld,
   double solvent_sld,
   double lr,
   double parav,
   double conc,
   double wamw)
{
    const double nump = 6.009956e-3 * conc / wamw;
    const double contrast = sld - solvent_sld;
    const double pcs = pow(sas_2J1x_x(q*radius), 2);
    const double volume = M_PI*radius*radius*length;
    const double pwlc = Sk_WR(q, length, kuhn_length);
    const double cq = cq_rod(q, lr);
    return 1e8 * nump * pow(contrast, 2) * pwlc * pcs / (1 + parav * cq * pwlc);
}