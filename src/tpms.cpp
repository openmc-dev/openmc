#include "openmc/tpms.h"

#include "openmc/error.h"

#include "openmc/constants.h"

namespace openmc {

// *****************************************************************************
// *   TPMS GENERAL DEFINITION AND RAY TRACING SOLVER
// *****************************************************************************

TPMS::TPMS(double _x0, double _y0, double _z0, double _a, double _b, double _c, double _d, double _e, double _f, double _g, double _h, double _i)
{
    x0 = _x0; y0 = _y0; z0 = _z0;
    a = _a; b = _b; c = _c; d = _d; e = _e; f = _f; g = _g; h = _h; i = _i;
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

double TPMS::ray_tracing(Position r, Direction u, double max_range)
{
    std::uintmax_t max_iter = 1000000;
    const double w0 = this->sampling_frequency(u);
    double L0 = 1.e-7; // A tolerance is set to discard eventual zero solutions
    double L1 = L0 + w0;
    double root;
    bool rootFound = false;
    while (L0 < max_range && rootFound==false)
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
    if (L0 >= max_range) {root = INFTY;}
    return root;
}

// *****************************************************************************
// *   TPMS CLASSIC DEFINITION WITH CONSTANT PITCH AND THICKNESS
// *****************************************************************************

TPMSClassic::TPMSClassic(double _cst, double _pitch, double _x0, double _y0, double _z0, double _a, double _b, double _c, double _d, double _e, double _f, double _g, double _h, double _i)
: TPMS( _x0, _y0, _z0, _a, _b, _c, _d, _e, _f, _g, _h, _i), cst(_cst), pitch(_pitch)
{}

// *****************************************************************************
// *   TPMS FUNCTION DEFINITION WITH NON-CONSTANT PITCH AND THICKNESS
// *****************************************************************************

TPMSFunction::TPMSFunction(const FunctionForTPMS& _fThickness, const FunctionForTPMS& _fPitch, double _x0, double _y0, double _z0, double _a, double _b, double _c, double _d, double _e, double _f, double _g, double _h, double _i)
: fThickness(_fThickness), fPitch(_fPitch), TPMS(_x0, _y0, _z0, _a, _b, _c, _d, _e, _f, _g, _h, _i)
{}

double TPMSFunction::get_pitch(Position r) const
{
    double pitch = fPitch.fxyz(r.x, r.y, r.z);
    return pitch;
}

double TPMSFunction::get_thickness(Position r) const
{
    double thickness = fThickness.fxyz(r.x, r.y, r.z);
    return thickness;
}

std::vector<double> TPMSFunction::get_pitch_first_partial_derivatives(Position r) const
{
    std::vector<double> derivatives = fPitch.first_partial_derivatives(r.x, r.y, r.z);
    return derivatives;
}

std::vector<double> TPMSFunction::get_thickness_first_partial_derivatives(Position r) const
{
    std::vector<double> derivatives = fThickness.first_partial_derivatives(r.x, r.y, r.z);
    return derivatives;
}

std::vector<std::vector<double>> TPMSFunction::get_pitch_second_partial_derivatives(Position r) const
{
    std::vector<std::vector<double>> derivatives = fPitch.second_partial_derivatives(r.x, r.y, r.z);
    return derivatives;
}

std::vector<std::vector<double>> TPMSFunction::get_thickness_second_partial_derivatives(Position r) const
{
    std::vector<std::vector<double>> derivatives = fThickness.second_partial_derivatives(r.x, r.y, r.z);
    return derivatives;
}

// *****************************************************************************
// *   SCHWARZ P DEFINITION
// *****************************************************************************

double SchwarzP::evaluate(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    return cos(l*xx) + cos(l*yy) + cos(l*zz) - cst;
}

Direction SchwarzP::normal(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double dx = l*(-a*sin(l*xx) -d*sin(l*yy) -g*sin(l*zz));
    const double dy = l*(-b*sin(l*xx) -e*sin(l*yy) -h*sin(l*zz));
    const double dz = l*(-c*sin(l*xx) -f*sin(l*yy) -i*sin(l*zz));
    const double norm = pow(pow(dx,2)+pow(dy,2)+pow(dz,2),0.5);
    Direction grad = Direction(dx/norm, dy/norm, dz/norm);
    return grad;
}

double SchwarzP::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position rk = Position(x,y,z);
    return this->evaluate(rk);
}

double SchwarzP::fpk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return -xxp*l*sin(l*xx) -yyp*l*sin(l*yy) -zzp*l*sin(l*zz);
}

double SchwarzP::fppk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return -xxp*xxp*l*l*cos(l*xx) -yyp*yyp*l*l*cos(l*yy) -zzp*zzp*l*l*cos(l*zz);
}

double SchwarzP::sampling_frequency(Direction u) const
{
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    std::vector<double> pulses = {abs(xxp), abs(yyp), abs(zzp)};
    std::vector<double>::iterator pmax;
    pmax = std::max_element(pulses.begin(), pulses.end());
    return 0.125*pitch / *pmax;
}

// *****************************************************************************
// *   GYROID DEFINITION
// *****************************************************************************

// Linearized expression : -0.5*sin(x - y) + 0.5*sin(x + y) + 0.5*sin(x - z) + 0.5*sin(x + z) - 0.5*sin(y - z) + 0.5*sin(y + z)

double Gyroid::evaluate(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    return sin(l*xx)*cos(l*zz) + sin(l*yy)*cos(l*xx) + sin(l*zz)*cos(l*yy) - cst;
}

Direction Gyroid::normal(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double cosxmy = cos(l*(xx - yy));
    const double cosxpy = cos(l*(xx + yy));
    const double cosxmz = cos(l*(xx - zz));
    const double cosxpz = cos(l*(xx + zz));
    const double cosymz = cos(l*(yy - zz));
    const double cosypz = cos(l*(yy + zz));
    const double dx = 0.5*l*((-a + d)*cosxmy + (a + d)*cosxpy + (a - g)*cosxmz + (a + g)*cosxpz + (-d + g)*cosymz + (d + g)*cosypz);
    const double dy = 0.5*l*((-b + e)*cosxmy + (b + e)*cosxpy + (b - h)*cosxmz + (b + h)*cosxpz + (-e + h)*cosymz + (e + h)*cosypz);
    const double dz = 0.5*l*((-c + f)*cosxmy + (c + f)*cosxpy + (c - i)*cosxmz + (c + i)*cosxpz + (-f + i)*cosymz + (f + i)*cosypz);
    const double norm = pow(pow(dx,2)+pow(dy,2)+pow(dz,2),0.5);
    Direction grad = Direction(dx/norm, dy/norm, dz/norm);
    return grad;
}

double Gyroid::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position rk = Position(x,y,z);
    return this->evaluate(rk);
}

double Gyroid::fpk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return 0.5*l*((-xxp+yyp)*cos(l*(xx-yy)) + (xxp+yyp)*cos(l*(xx+yy)) + (xxp-zzp)*cos(l*(xx-zz)) + (xxp+zzp)*cos(l*(xx+zz)) + (-yyp+zzp)*cos(l*(yy-zz)) + (yyp+zzp)*cos(l*(yy+zz)));
}

double Gyroid::fppk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return -0.5*l*l*(-(xxp - yyp)*(xxp - yyp)*sin(l*(xx - yy)) + (xxp + yyp)*(xxp + yyp)*sin(l*(xx + yy)) + (xxp - zzp)*(xxp - zzp)*sin(l*(xx - zz)) + (xxp + zzp)*(xxp + zzp)*sin(l*(xx + zz)) - (yyp - zzp)*(yyp - zzp)*sin(l*(yy - zz)) + (yyp + zzp)*(yyp + zzp)*sin(l*(yy + zz)));
}

double Gyroid::sampling_frequency(Direction u) const
{
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    std::vector<double> pulses = {abs(xxp+yyp), abs(xxp-yyp), abs(xxp+zzp), abs(xxp-zzp), abs(yyp+zzp), abs(yyp-zzp)};
    std::vector<double>::iterator pmax;
    pmax = std::max_element(pulses.begin(), pulses.end());
    return 0.125*pitch / *pmax;
}

// *****************************************************************************
// *   DIAMOND DEFINITION
// *****************************************************************************

// Linearized expression : 0.5*sin(-x + y + z) + 0.5*sin(x - y + z) + 0.5*sin(x + y - z) + 0.5*sin(x + y + z)

double Diamond::evaluate(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    return 0.5*sin(l*(-xx+yy+zz)) + 0.5*sin(l*(+xx-yy+zz)) + 0.5*sin(l*(+xx+yy-zz)) + 0.5*sin(l*(+xx+yy+zz)) - cst;
}

Direction Diamond::normal(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double cosypzmx = cos(l*(-xx+yy+zz));
    const double coszpxmy = cos(l*(+xx-yy+zz));
    const double cosxpymz = cos(l*(+xx+yy-zz));
    const double cosxpypz = cos(l*(+xx+yy+zz));
    const double dx = 0.5*l*((-a + d + g)*cosypzmx + (a - d + g)*coszpxmy + (a + d - g)*cosxpymz + (a + d + g)*cosxpypz);
    const double dy = 0.5*l*((-b + e + h)*cosypzmx + (b - e + h)*coszpxmy + (b + e - h)*cosxpymz + (b + e + h)*cosxpypz);
    const double dz = 0.5*l*((-c + f + i)*cosypzmx + (c - f + i)*coszpxmy + (c + f - i)*cosxpymz + (c + f + i)*cosxpypz);
    const double norm = pow(pow(dx,2)+pow(dy,2)+pow(dz,2),0.5);
    Direction grad = Direction(dx/norm, dy/norm, dz/norm);
    return grad;
}

double Diamond::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position rk = Position(x,y,z);
    return this->evaluate(rk);
}

double Diamond::fpk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return 0.5*l*((-xxp + yyp + zzp)*cos(l*(-xx + yy + zz)) + (xxp - yyp + zzp)*cos(l*(xx - yy + zz)) + (xxp + yyp - zzp)*cos(l*(xx + yy - zz)) + (xxp + yyp + zzp)*cos(l*(xx + yy + zz)));
}

double Diamond::fppk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    const double l = 2*M_PI/pitch;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    return -0.5*l*l*(-(xxp - yyp)*(xxp - yyp)*sin(l*(xx - yy)) + (xxp + yyp)*(xxp + yyp)*sin(l*(xx + yy)) + (xxp - zzp)*(xxp - zzp)*sin(l*(xx - zz)) + (xxp + zzp)*(xxp + zzp)*sin(l*(xx + zz)) - (yyp - zzp)*(yyp - zzp)*sin(l*(yy - zz)) + (yyp + zzp)*(yyp + zzp)*sin(l*(yy + zz)));
}

double Diamond::sampling_frequency(Direction u) const
{
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    std::vector<double> pulses = {abs(+xxp+yyp+zzp), abs(+xxp+yyp-zzp), abs(+xxp-yyp+zzp), abs(-xxp+yyp+zzp)};
    std::vector<double>::iterator pmax;
    pmax = std::max_element(pulses.begin(), pulses.end());
    return 0.125*pitch / *pmax;
}

// *****************************************************************************
// *   FUNCTION PITCH & THICKNESS SCHWARZ P DEFINITION
// *****************************************************************************

double FunctionSchwarzP::evaluate(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    Position rr = Position(xx,yy,zz);
    const double local_pitch = this->get_pitch(rr);
    const double local_thickness = this->get_thickness(rr);
    const double l = 2*M_PI/local_pitch;
    const double local_cst = l*local_thickness;
    return cos(l*xx) + cos(l*yy) + cos(l*zz) - local_cst;
}

Direction FunctionSchwarzP::normal(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    Position rr = Position(xx,yy,zz);
    const double local_pitch = this->get_pitch(rr);
    const double local_thickness = this->get_thickness(rr);
    const std::vector<double> der_pitch = this->get_pitch_first_partial_derivatives(rr);
    const std::vector<double> der_thick = this->get_thickness_first_partial_derivatives(rr);
    const double l = 2*M_PI/local_pitch;
    const double dxx_l = -l*der_pitch[0]/local_pitch;
    const double dyy_l = -l*der_pitch[1]/local_pitch;
    const double dzz_l = -l*der_pitch[2]/local_pitch;
    const double dx_l = a*dxx_l + d*dyy_l + g*dzz_l;
    const double dy_l = b*dxx_l + e*dyy_l + h*dzz_l;
    const double dz_l = c*dxx_l + f*dyy_l + i*dzz_l;
    const double dxx_cst = dxx_l*local_thickness + l * der_thick[0];
    const double dyy_cst = dyy_l*local_thickness + l * der_thick[1];
    const double dzz_cst = dzz_l*local_thickness + l * der_thick[2];
    const double dx_cst = a*dxx_cst + d*dyy_cst + g*dzz_cst;
    const double dy_cst = b*dxx_cst + e*dyy_cst + h*dzz_cst;
    const double dz_cst = c*dxx_cst + f*dyy_cst + i*dzz_cst;    
    const double dx = -(xx*dx_l + a*l)*sin(l*xx) -(yy*dx_l + d*l)*sin(l*yy) -(zz*dx_l + g*l)*sin(l*zz) - dx_cst;
    const double dy = -(xx*dy_l + b*l)*sin(l*xx) -(yy*dy_l + e*l)*sin(l*yy) -(zz*dy_l + h*l)*sin(l*zz) - dy_cst;
    const double dz = -(xx*dz_l + c*l)*sin(l*xx) -(yy*dz_l + f*l)*sin(l*yy) -(zz*dz_l + i*l)*sin(l*zz) - dz_cst;
    const double norm = pow(pow(dx,2)+pow(dy,2)+pow(dz,2),0.5);
    Direction grad = Direction(dx/norm, dy/norm, dz/norm);
    return grad;
}

double FunctionSchwarzP::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position rk = Position(x,y,z);
    return this->evaluate(rk);
}

double FunctionSchwarzP::fpk(double k, Position r, Direction u) const
{
    if (fPitch.useFirstDerivatives == false)
    {
        double dk = 1.e-9;
        double f1 = this->fk(k+dk, r, u);
        double f0 = this->fk(k-dk, r, u);
        double brut_value = (f1-f0)/(2*dk);
        return (f1-f0)/(2*dk);
    }
    else
    {
        // to code explicit derivatives the day functions not using approximate derivatives are implemented
        fatal_error("Analitical first derivate is not implemented for FunctionSchwarzP.");
    }
}

double FunctionSchwarzP::fppk(double k, Position r, Direction u) const
{
    if (fPitch.useSecondDerivatives == false)
    {
        double dk = 1.e-9;
        double f2 = this->fk(k+dk, r, u);
        double f1 = this->fk(k, r, u);
        double f0 = this->fk(k-dk, r, u);
        return (f2-2*f1+f0)/(pow(dk,2));
    }
    else
    {
        // to code explicit derivatives the day functions not using approximate derivatives are implemented
        fatal_error("Analitical second derivate is not implemented for FunctionSchwarzP.");
    }
}

double FunctionSchwarzP::sampling_frequency(Direction u) const
{
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    std::vector<double> pulses = {abs(xxp), abs(yyp), abs(zzp)};
    std::vector<double>::iterator pmax;
    pmax = std::max_element(pulses.begin(), pulses.end());
    const double min_pitch = fPitch.minimalValue;
    if (min_pitch <= 0.) {fatal_error("Pitch matrix for interpolation contains zero or negative values.");}
    return 0.125*min_pitch / *pmax;
}

// *****************************************************************************
// *   FUNCTION PITCH & THICKNESS GYROID DEFINITION
// *****************************************************************************

double FunctionGyroid::evaluate(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    Position rr = Position(xx,yy,zz);
    const double local_pitch = this->get_pitch(rr);
    const double local_thickness = this->get_thickness(rr);
    const double l = 2*M_PI/local_pitch;
    const double local_cst = l*local_thickness;
    return sin(l*xx)*cos(l*zz) + sin(l*yy)*cos(l*xx) + sin(l*zz)*cos(l*yy) - local_cst;
}

Direction FunctionGyroid::normal(Position r) const
{
    // TODO: a coder quand j'aurais un peu moins la flemme.
    return Direction();
}

double FunctionGyroid::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position rk = Position(x,y,z);
    return this->evaluate(rk);
}

double FunctionGyroid::fpk(double k, Position r, Direction u) const
{
    if (fPitch.useFirstDerivatives == false)
    {
        double dk = 1.e-9;
        double f1 = this->fk(k+dk, r, u);
        double f0 = this->fk(k-dk, r, u);
        double brut_value = (f1-f0)/(2*dk);
        return (f1-f0)/(2*dk);
    }
    else
    {
        // to code explicit derivatives the day functions not using approximate derivatives are implemented
        fatal_error("Analitical first derivate is not implemented for FunctionSchwarzP.");
    }
}

double FunctionGyroid::fppk(double k, Position r, Direction u) const
{
    if (fPitch.useSecondDerivatives == false)
    {
        double dk = 1.e-9;
        double f2 = this->fk(k+dk, r, u);
        double f1 = this->fk(k, r, u);
        double f0 = this->fk(k-dk, r, u);
        return (f2-2*f1+f0)/(pow(dk,2));
    }
    else
    {
        // to code explicit derivatives the day functions not using approximate derivatives are implemented
        fatal_error("Analitical second derivate is not implemented for FunctionSchwarzP.");
    }
}

double FunctionGyroid::sampling_frequency(Direction u) const
{
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    std::vector<double> pulses = {abs(xxp+yyp), abs(xxp-yyp), abs(xxp+zzp), abs(xxp-zzp), abs(yyp+zzp), abs(yyp-zzp)};
    std::vector<double>::iterator pmax;
    pmax = std::max_element(pulses.begin(), pulses.end());
    const double min_pitch = fPitch.minimalValue;
    if (min_pitch <= 0.) {fatal_error("Pitch matrix for interpolation contains zero or negative values.");}
    return 0.125*min_pitch / *pmax;
}

// *****************************************************************************
// *   FUNCTION PITCH & THICKNESS DIAMOND DEFINITION
// *****************************************************************************

double FunctionDiamond::evaluate(Position r) const
{
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double xx = a*(x-x0) + b*(y-y0) + c*(z-z0);
    const double yy = d*(x-x0) + e*(y-y0) + f*(z-z0);
    const double zz = g*(x-x0) + h*(y-y0) + i*(z-z0);
    Position rr = Position(xx,yy,zz);
    const double local_pitch = this->get_pitch(rr);
    const double local_thickness = this->get_thickness(rr);
    const double l = 2*M_PI/local_pitch;
    const double local_cst = l*local_thickness;
    return 0.5*sin(l*(-xx+yy+zz)) + 0.5*sin(l*(+xx-yy+zz)) + 0.5*sin(l*(+xx+yy-zz)) + 0.5*sin(l*(+xx+yy+zz)) - local_cst;
}

Direction FunctionDiamond::normal(Position r) const
{
    // TODO: a coder quand j'aurais un peu moins la flemme.
    return Direction();
}

double FunctionDiamond::fk(double k, Position r, Direction u) const
{
    const double x = r.x + k*u.x;
    const double y = r.y + k*u.y;
    const double z = r.z + k*u.z;
    Position rk = Position(x,y,z);
    return this->evaluate(rk);
}

double FunctionDiamond::fpk(double k, Position r, Direction u) const
{
    if (fPitch.useFirstDerivatives == false)
    {
        double dk = 1.e-9;
        double f1 = this->fk(k+dk, r, u);
        double f0 = this->fk(k-dk, r, u);
        double brut_value = (f1-f0)/(2*dk);
        return (f1-f0)/(2*dk);
    }
    else
    {
        // to code explicit derivatives the day functions not using approximate derivatives are implemented
        fatal_error("Analitical first derivate is not implemented for FunctionSchwarzP.");
    }
}

double FunctionDiamond::fppk(double k, Position r, Direction u) const
{
    if (fPitch.useSecondDerivatives == false)
    {
        double dk = 1.e-9;
        double f2 = this->fk(k+dk, r, u);
        double f1 = this->fk(k, r, u);
        double f0 = this->fk(k-dk, r, u);
        return (f2-2*f1+f0)/(pow(dk,2));
    }
    else
    {
        // to code explicit derivatives the day functions not using approximate derivatives are implemented
        fatal_error("Analitical second derivate is not implemented for FunctionSchwarzP.");
    }
}

double FunctionDiamond::sampling_frequency(Direction u) const
{
    const double xxp = a*u.x + b*u.y + c*u.z;
    const double yyp = d*u.x + e*u.y + f*u.z;
    const double zzp = g*u.x + h*u.y + i*u.z;
    std::vector<double> pulses = {abs(+xxp+yyp+zzp), abs(+xxp+yyp-zzp), abs(+xxp-yyp+zzp), abs(-xxp+yyp+zzp)};
    std::vector<double>::iterator pmax;
    pmax = std::max_element(pulses.begin(), pulses.end());
    const double min_pitch = fPitch.minimalValue;
    if (min_pitch <= 0.) {fatal_error("Pitch matrix for interpolation contains zero or negative values.");}
    return 0.125*min_pitch / *pmax;
}


}