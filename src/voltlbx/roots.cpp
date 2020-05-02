#include "roots.h"

#include <algorithm> 

namespace voltlbx
{
    std::pair<double, double> bracket_root(const std::function<double(double)>& f, double x1, double x2)
    {
        const int NTRY = 50;
        const double FACTOR = 1.6;
            
        if(std::abs(x2-x1) == 0.0)
            throw std::invalid_argument("bracket_root : invalid initial bracket");
        
        double xa = x1;
        double xb = x2;
        double fa = f(xa);
        double fb = f(xb);

        double last = std::numeric_limits<double>::quiet_NaN();
        double lastf = std::numeric_limits<double>::quiet_NaN();

        for (int j = 0; j < NTRY; j++)
        {
            if (fa * fb < 0.0)
            {
                if (lastf * fa < 0.0)
                {
                    xb = last;
                    fb = lastf;
                }
                if (lastf * fb < 0.0)
                {
                    xa = last;
                    fa = lastf;
                }
                return { xa, xb };
            }
            if (std::abs(fa) < std::abs(fb))
            {
                last = xa;
                lastf = fa;
                xa += FACTOR * (xa - xb);
                fa = f(xa);
            }
            else
            {
                last = xb;
                lastf = fb;
                xb += FACTOR * (xb - xa);
                fb = f(xb);
            }
        }
        throw std::runtime_error("bracket_root : too many iterations");
    }
       

    double brenth(const std::function<double(double)>& f, 
                  double xa, double xb,
                  double xtol, double rtol, std::size_t max_iter)
    {
        double xpre = xa;
        double xcur = xb;
        double fpre = f(xa);
        double fcur = f(xb);
        std::size_t funcalls = 0;
        std::size_t iterations = 0;

        if (fpre * fcur > 0)
            throw std::invalid_argument("Brent : root must be bracketed");

        if (std::abs(fpre)==0.0) return xpre;
        if (std::abs(fcur)==0.0) return xcur;

        double xblk = 0.0, fblk = 0.0, spre = 0.0, scur = 0.0;
        for (int i = 0; i < max_iter; i++)
        {
            iterations++;
            if (fpre * fcur < 0)
            {
                xblk = xpre;
                fblk = fpre;
                spre = scur = xcur - xpre;
            }
            if (std::abs(fblk) < std::abs(fcur))
            {
                xpre = xcur;
                xcur = xblk;
                xblk = xpre;
                fpre = fcur;
                fcur = fblk;
                fblk = fpre;
            }

            double tol = xtol + rtol * std::abs(xcur);
            double sbis = (xblk - xcur) / 2.0;
            if (std::abs(fcur)==0.0 || std::abs(sbis) < tol)
                return xcur;

            if (std::abs(spre) > tol && std::abs(fcur) < std::abs(fpre))
            {
                double stry;
                if (xpre == xblk)
                {
                    // interpolate
                    stry = -fcur * (xcur - xpre) / (fcur - fpre);
                }
                else
                {
                    // extrapolate
                    double dpre = (fpre - fcur) / (xpre - xcur);
                    double dblk = (fblk - fcur) / (xblk - xcur);

                    stry = -fcur * (fblk - fpre) / (fblk * dpre - fpre * dblk);
                }

                if (2.0 * std::abs(stry) < std::min(std::abs(spre), 3.0 * std::abs(sbis) - tol))
                {
                    // accept step
                    spre = scur;
                    scur = stry;
                }
                else
                {
                    // bisect 
                    spre = sbis;
                    scur = sbis;
                }
            }
            else
            {
                // bisect 
                spre = sbis;
                scur = sbis;
            }

            xpre = xcur;
            fpre = fcur;
            if (std::abs(scur) > tol)
                xcur += scur;
            else
                xcur += (sbis > 0.0 ? tol : -tol);

            fcur = f(xcur);
            funcalls++;
        }
        throw std::runtime_error("Brent : max iteration excedeed");
    }

}