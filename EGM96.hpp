#include "EGM96_data.h"
#include <math.h>

#define _coeffs (65341) //!< Size of correction and harmonic coefficients arrays (361*181)
#define _nmax   (360)   //!< Maximum degree and orders of harmonic coefficients.
#define _361    (361)

double hundu(double p[_coeffs+1],
             double sinml[_361+1], double cosml[_361+1],
             double gr, double re)
{
    // WGS 84 gravitational constant in m³/s² (mass of Earth’s atmosphere included)
    const double GM = 0.3986004418e15;
    // WGS 84 datum surface equatorial radius
    const double ae = 6378137.0;

    double ar = ae/re;
    double arn = ar;
    double ac = 0;
    double a = 0;

    unsigned k = 3;
    for (unsigned n = 2; n <= _nmax; n++)
    {
        arn *= ar;
        k++;
        double sum = p[k]*egm96_data[k][2];
        double sumc = p[k]*egm96_data[k][0];

        for (unsigned m = 1; m <= n; m++)
        {
            k++;
            double tempc = egm96_data[k][0]*cosml[m] + egm96_data[k][1]*sinml[m];
            double temp  = egm96_data[k][2]*cosml[m] + egm96_data[k][3]*sinml[m];
            sumc += p[k]*tempc;
            sum  += p[k]*temp;
        }
        ac += sumc;
        a += sum*arn;
    }
    ac += egm96_data[1][0] + (p[2]*egm96_data[2][0]) + (p[3] * (egm96_data[3][0]*cosml[1] + egm96_data[3][1]*sinml[1]));

    // Add haco = ac/100 to convert height anomaly on the ellipsoid to the undulation
    // Add -0.53m to make undulation refer to the WGS84 ellipsoid

    return ((a * GM) / (gr * re)) + (ac / 100.0) - 0.53;
}

void dscml(double rlon, double sinml[_361+1], double cosml[_361+1])
{
    double a = sin(rlon);
    double b = cos(rlon);

    sinml[1] = a;
    cosml[1] = b;
    sinml[2] = 2*b*a;
    cosml[2] = 2*b*b - 1;

    for (unsigned m = 3; m <= _nmax; m++)
    {
        sinml[m] = 2*b*sinml[m-1] - sinml[m-2];
        cosml[m] = 2*b*cosml[m-1] - cosml[m-2];
    }
}

/*!
 * \param m: order.
 * \param theta: Colatitude (radians).
 * \param rleg: Normalized legendre function.
 *
 * This subroutine computes all normalized legendre function in 'rleg'.
 * The dimensions of array 'rleg' must be at least equal to nmax+1.
 * All calculations are in double precision.
 *
 * Original programmer: Oscar L. Colombo, Dept. of Geodetic Science the Ohio State University, August 1980.
 * ineiev: I removed the derivatives, for they are never computed here.
 */
void legfdn(unsigned m, double theta, double rleg[_361+1])
{
    static double drts[1301], dirt[1301], cothet, sithet, rlnn[_361+1];
    static int ir; // TODO 'ir' must be set to zero before the first call to this sub.

    unsigned nmax1 = _nmax + 1;
    unsigned nmax2p = (2 * _nmax) + 1;
    unsigned m1 = m + 1;
    unsigned m2 = m + 2;
    unsigned m3 = m + 3;
    unsigned n, n1, n2;

    if (!ir)
    {
        ir = 1;
        for (n = 1; n <= nmax2p; n++)
        {
            drts[n] = sqrt(n);
            dirt[n] = 1 / drts[n];
        }
    }

    cothet = cos(theta);
    sithet = sin(theta);

    // compute the legendre functions
    rlnn[1] = 1;
    rlnn[2] = sithet * drts[3];
    for (n1 = 3; n1 <= m1; n1++)
    {
        n = n1 - 1;
        n2 = 2 * n;
        rlnn[n1] = drts[n2 + 1] * dirt[n2] * sithet * rlnn[n];
    }

    switch (m)
    {
        case 1:
            rleg[2] = rlnn[2];
            rleg[3] = drts[5] * cothet * rleg[2];
            break;
        case 0:
            rleg[1] = 1;
            rleg[2] = cothet * drts[3];
            break;
    }
    rleg[m1] = rlnn[m1];

    if (m2 <= nmax1)
    {
        rleg[m2] = drts[m1*2 + 1] * cothet * rleg[m1];
        if (m3 <= nmax1)
        {
            for (n1 = m3; n1 <= nmax1; n1++)
            {
                n = n1 - 1;
                if ((!m && n < 2) || (m == 1 && n < 3)) continue;
                n2 = 2 * n;
                rleg[n1] = drts[n2+1] * dirt[n+m] * dirt[n-m] * (drts[n2-1] * cothet * rleg[n1-1] - drts[n+m-1] * drts[n-m-1] * dirt[n2-3] * rleg[n1-2]);
            }
        }
    }
}

/*!
 * \param lat: Latitude in radians.
 * \param lon: Longitude in radians.
 * \param re: Geocentric radius.
 * \param rlat: Geocentric latitude.
 * \param gr: Normal gravity (m/sec²).
 *
 * This subroutine computes geocentric distance to the point, the geocentric
 * latitude, and an approximate value of normal gravity at the point based the
 * constants of the WGS84(g873) system are used.
 */
void radgra(double lat, double lon, double *rlat, double *gr, double *re)
{
    const double a = 6378137.0;
    const double e2 = 0.00669437999013;
    const double geqt = 9.7803253359;
    const double k = 0.00193185265246;
    double t1 = sin(lat) * sin(lat);
    double n = a / sqrt(1.0 - (e2 * t1));
    double t2 = n * cos(lat);
    double x = t2 * cos(lon);
    double y = t2 * sin(lon);
    double z = (n * (1 - e2)) * sin(lat);

    *re = sqrt((x * x) + (y * y) + (z * z));            // compute the geocentric radius
    *rlat = atan(z / sqrt((x * x) + (y * y)));          // compute the geocentric latitude
    *gr = geqt * (1 + (k * t1)) / sqrt(1 - (e2 * t1));  // compute normal gravity (m/sec²)
}

/*!
 * \brief Compute the geoid undulation from the EGM96 potential coefficient model, for a given latitude and longitude.
 * \param lat: Latitude in radians.
 * \param lon: Longitude in radians.
 * \return The geoid undulation / altitude offset (in meters).
 */
double undulation(double lat, double lon)
{
    double p[_coeffs+1], sinml[_361+1], cosml[_361+1], rleg[_361+1];

    double rlat, gr, re;
    unsigned nmax1 = _nmax + 1;

    // compute the geocentric latitude, geocentric radius, normal gravity
    radgra(lat, lon, &rlat, &gr, &re);
    rlat = (M_PI / 2) - rlat;

    for (unsigned j = 1; j <= nmax1; j++)
    {
        unsigned m = j - 1;
        legfdn(m, rlat, rleg);
        for (unsigned i = j ; i <= nmax1; i++)
        {
            p[(((i - 1) * i) / 2) + m + 1] = rleg[i];
        }
     }
     dscml(lon, sinml, cosml);

     return hundu(p, sinml, cosml, gr, re);
}

/*!
 * \brief Compute the geoid undulation from the EGM96 potential coefficient model, for a given latitude and longitude.
 * \param latitude: Latitude (in degrees).
 * \param longitude: Longitude (in degrees).
 * \return The geoid undulation / altitude offset (in meters).
 */
double egm96_compute_altitude_offset(double lat, double lon)
{
    const double rad = (180.0 / M_PI);
    return undulation(lat/rad, lon/rad);
}
