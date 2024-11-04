#include "openmc/tpms.h"

namespace openmc {

TPMS::TPMS(double _cst, double _pitch, double _x0, double _y0, double _z0, double _a, double _b, double _c, double _d, double _e, double _f, double _g, double _h, double _i)
{
    cst = _cst; pitch = _pitch;
    x0 = _x0; y0 = _y0; z0 = _z0;
    a = _a; b = _b; c = _c; d = _d; e = _e; f = _f; g = _g; h = _h; i = _i;
}

double TPMS::fk2(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position xyz = Position(x, y, z);
    return this->fxyz(xyz);
}

double TPMS::fpk2(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position xyz = Position(x, y, z);
    std::vector<double> grad = this->grad(xyz);
    return u.x*grad[0] + u.y*grad[1] +u.x*grad[2];
}

double TPMS::fppk2(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position xyz = Position(x, y, z);
    std::vector<std::vector<double>> hess = this->hess(xyz);
    double fpx = u.x*( u.x*hess[0][0] + u.y*hess[0][1] + u.z*hess[0][2] );
    double fpy = u.y*( u.x*hess[1][0] + u.y*hess[1][1] + u.z*hess[1][2] );
    double fpz = u.z*( u.x*hess[2][0] + u.y*hess[2][1] + u.z*hess[2][2] );
    return fpx + fpy + fpz;
}

TPMS::rootFinding TPMS::root_in_interval(double L0, double L1, Position r, Direction u)
{
    bool solFound = false;
    double xa = L0;
    double xb = L1;
    rootFinding solution;
    while (solFound==false)
    {
        solution.xa = xa;
        solution.xb = xb;
        double fa   = this->fk(xa, r, u);
        double fpa  = this->fpk(xa, r, u);
        double fppa = this->fppk(xa, r, u);
        double fb   = this->fk(xb, r, u);
        double fpb  = this->fpk(xb, r, u);
        double fppb = this->fppk(xb, r, u);
        //std::cout << "[" << L0 << "-" << L1 <<"] XA: " << xa <<" XB: " << xb << " - FA " << fa << " FPA " << fpa << " FPPA " << fppa << " FB " << fb << " FPB " << fpb << " FPPB " << fppb << " " << solFound << std::endl;
        if (std::signbit(fa) != std::signbit(fb)) {solFound=true; solution.isRoot = true;solution.status=1;} // if the two sampled point are of different sign, there is a root in the interval
        else if (std::signbit(fpa) == std::signbit(fpb)) // if the two root are of same sign and the two derivative are of same sign, there is no root in the interval
        {
            solFound=true; 
            solution.isRoot = false;
        }
        else if (std::signbit(fppa) != std::signbit(fppb)) // If inflexion change, split the interval in half to study root existance. (rare case)
        {
            TPMS::rootFinding firstInterval = this->root_in_interval(L0, 0.5*(L0+L1), r, u);
            if (firstInterval.isRoot == true) 
            {
                solFound=true; 
                solution.isRoot = true; 
                solution.xa = firstInterval.xa; 
                solution.xb = firstInterval.xb;
                solution.status = firstInterval.status;
            }
            else 
            {
                solFound=true; 
                solution = this->root_in_interval(0.5*(L0+L1), L1, r, u);
            }
        }
        else if (std::signbit(fppa) != std::signbit(fa)) // If no inflexion change and the sign of the second derivative is not the one of the sampled points, there is no root
        {
            solFound=true; 
            solution.isRoot = false;
        }
        else
        {
            double xn = -(fa - fpa*xa - fb + fpb*xb)/(fpa - fpb); // x de l'intersection des droites
            double fn = fpa * xn + fa - fpa*xa; // y de l'intersection des droites
            double fxn = this->fk(xn, r, u); // evaluation de la tpms en x
            if (std::signbit(fa) == std::signbit(fn)) // if the intersection ordinate is of same sign as the sampled points, there is no root
            {
                solFound=true; 
                solution.isRoot = false;
            } 
            else if (std::signbit(fa) != std::signbit(fxn)) // if the evaluate function in x is of different sign of the sampled points, there is a root
            {
                solFound=true; 
                solution.isRoot = true; 
                solution.xb = xn;
                solution.status = 2;
            }
            else // Otherwise, we don't know and need to refine.
            {
                double fpn = this->fpk(xn, r, u); 
                if (std::signbit(fpn) != std::signbit(fpa)) {xb = xn;} 
                else {xa = xn;}
            }
        }
    }
    return solution;
}

double TPMS::ray_tracing(Position r, Direction u)
{
    std::uintmax_t max_iter = 1000000;
    const double w0 = this->sampling_frequency(u);
    double L0 = 0.;
    double L1 = L0 + w0;
    double root;
    bool rootFound = false;
    while (L0 < TPMS::XLIM && rootFound==false)
    {
        TPMS::rootFinding solution = this->root_in_interval(L0, L1, r, u);
        if (solution.isRoot) 
        {
            // std::pair<double, double> sol = boost::math::tools::bisect([this,r,u](double k){return this->fk(k, r, u);}, solution.xa, solution.xb, [](double l, double r){return abs(l-r) < 1e-8;});
            std::pair<double, double> sol = bisect([this,r,u](double k){return this->fk(k, r, u);}, solution.xa, solution.xb, [](double l, double r){return abs(l-r) < 1e-8;}, max_iter);
            root = sol.second;
            rootFound = true;
        }
        else {
            L0 += w0;
            L1 += w0;
        }
    }
    return root;
}

double SchwarzP::fxyz(Position xyz) const
{
    const double x = xyz.x;
    const double y = xyz.y;
    const double z = xyz.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    return cos(l*xx) + cos(l*yy) + cos(l*zz) - cst;
}

std::vector<double> SchwarzP::grad(Position xyz) const
{
    const double x = xyz.x;
    const double y = xyz.y;
    const double z = xyz.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double dfdx = -a*l*sin(l*xx) -d*l*sin(l*yy) -g*l*sin(l*zz);
    const double dfdy = -b*l*sin(l*xx) -e*l*sin(l*yy) -h*l*sin(l*zz);
    const double dfdz = -c*l*sin(l*xx) -f*l*sin(l*yy) -i*l*sin(l*zz);
    return {dfdx, dfdy, dfdz};
}

std::vector<std::vector<double>> SchwarzP::hess(Position xyz) const
{
    const double x = xyz.x;
    const double y = xyz.y;
    const double z = xyz.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double dfdxx = -a*a*l*l*cos(l*xx) -d*d*l*l*cos(l*yy) -g*g*l*l*cos(l*zz);
    const double dfdxy = -a*b*l*l*cos(l*xx) -d*e*l*l*cos(l*yy) -g*h*l*l*cos(l*zz);
    const double dfdxz = -a*c*l*l*cos(l*xx) -d*f*l*l*cos(l*yy) -g*i*l*l*cos(l*zz);
    const double dfdyy = -b*b*l*l*cos(l*xx) -e*e*l*l*cos(l*yy) -h*h*l*l*cos(l*zz);
    const double dfdyz = -b*c*l*l*cos(l*xx) -e*f*l*l*cos(l*yy) -h*i*l*l*cos(l*zz);
    const double dfdzz = -c*c*l*l*cos(l*xx) -f*f*l*l*cos(l*yy) -i*i*l*l*cos(l*zz);
    return {
        {dfdxx, dfdxy, dfdxz},
        {dfdxy, dfdyy, dfdyz},
        {dfdxz, dfdyz, dfdzz},
    };
}

double SchwarzP::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double L = pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    return cos(2*M_PI*xx/L) + cos(2*M_PI*yy/L) + cos(2*M_PI*zz/L) - cst;
}

double SchwarzP::fpk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double L = pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return -2*M_PI*xxp*sin(2*M_PI*xx/L)/L -2*M_PI*yyp*sin(2*M_PI*yy/L)/L -2*M_PI*zzp*sin(2*M_PI*zz/L)/L;
}

double SchwarzP::fppk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double L = pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    const double cX = pow(2*M_PI*xxp/L, 2);
    const double cY = pow(2*M_PI*yyp/L, 2);
    const double cZ = pow(2*M_PI*zzp/L, 2);
    return -cX*cos(2*M_PI*xx/L) -cY*cos(2*M_PI*yy/L) -cZ*cos(2*M_PI*zz/L);
}

double SchwarzP::sampling_frequency(Direction u) const
{
    const double xxp = abs(a*u.x + b*u.y + c*u.z);
    const double yyp = abs(d*u.x + e*u.y + f*u.z);
    const double zzp = abs(g*u.x + h*u.y + i*u.z);
    const double pmax = std::max(std::max(xxp, yyp), zzp);
    return 0.125*pitch/pmax;
}

}