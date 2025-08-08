import numpy as np
import scipy
import scipy.optimize

cos, sin, pi = np.cos, np.sin, np.pi

class TPMS:

    TYPES:list = ['Schwarz_P', 'Gyroid', 'Diamond']

    def __init__(self, name:  str, pitch:float, position:np.ndarray, direction:np.ndarray) -> None:
        assert name in TPMS.TYPES
        assert abs(direction[0]**2 + direction[1]**2 + direction[2]**2 - 1.) < 1e-15
        
        self.name = name
        self.pitch = pitch
        self.position = position
        self.direction = direction

    def get_ABCDEF(self) -> tuple[float]:
        return self.direction[0],self.position[0],self.direction[1],self.position[1],self.direction[2],self.position[2]

    def f(self, k:float) -> float:
        A,B,C,D,E,F = self.get_ABCDEF()
        x = 2*pi*(A*k+B)/self.pitch
        y = 2*pi*(C*k+D)/self.pitch
        z = 2*pi*(E*k+F)/self.pitch
        if self.name == 'Schwarz_P':
            return cos(x) + cos(y) + cos(z)
        elif self.name == 'Gyroid':
            return cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x)
        elif self.name == 'Diamond':
            return sin(x)*sin(y)*sin(z) + sin(x)*cos(y)*cos(z) + cos(x)*sin(y)*cos(z) + cos(x)*cos(y)*sin(z)
        else:
            raise ValueError(f"TPMS: {self.name} not implemented")

    def fp(self, k:float):
        A,B,C,D,E,F = self.get_ABCDEF()
        x = 2*pi*(A*k+B)/self.pitch
        y = 2*pi*(C*k+D)/self.pitch
        z = 2*pi*(E*k+F)/self.pitch
        xp = 2*pi*A/self.pitch
        yp = 2*pi*C/self.pitch
        zp = 2*pi*E/self.pitch
        if self.name == 'Schwarz_P':
            return -xp*sin(x) -yp*sin(y) -zp*sin(z)
        elif self.name == 'Gyroid':
            return xp*(-sin(x)*sin(y) + cos(x)*cos(z)) + yp*(-sin(y)*sin(z) + cos(x)*cos(y)) + zp*(-sin(x)*sin(z) + cos(y)*cos(z))
        elif self.name == 'Diamond':
            part1 = xp*(-sin(x)*sin(y)*cos(z) - sin(x)*sin(z)*cos(y) + sin(y)*sin(z)*cos(x) + cos(x)*cos(y)*cos(z))
            part2 = yp*(-sin(x)*sin(y)*cos(z) + sin(x)*sin(z)*cos(y) - sin(y)*sin(z)*cos(x) + cos(x)*cos(y)*cos(z))
            part3 = zp*(+sin(x)*sin(y)*cos(z) - sin(x)*sin(z)*cos(y) - sin(y)*sin(z)*cos(x) + cos(x)*cos(y)*cos(z))
            return part1+part2+part3
        else:
            raise ValueError(f"TPMS: {self.name} not implemented")

    def fpp(self, k:float):
        A,B,C,D,E,F = self.get_ABCDEF()
        x = 2*pi*(A*k+B)/self.pitch
        y = 2*pi*(C*k+D)/self.pitch
        z = 2*pi*(E*k+F)/self.pitch
        xp = 2*pi*A/self.pitch
        yp = 2*pi*C/self.pitch
        zp = 2*pi*E/self.pitch
        if self.name == 'Schwarz_P':
            return -(xp)**2*cos(x) -(yp)**2*cos(y) -(zp)**2*cos(z)
        elif self.name == 'Gyroid':
            part1 = xp*(xp*(-sin(x)*cos(z) - sin(y)*cos(x)) + yp*(-sin(x)*cos(y))                 + zp*(-sin(z)*cos(x)))
            part2 = yp*(xp*(-sin(x)*cos(y))                 + yp*(-sin(y)*cos(x) - sin(z)*cos(y)) + zp*(-sin(y)*cos(z)))
            part3 = zp*(xp*(-sin(z)*cos(x))                 + yp*(-sin(y)*cos(z))                 + zp*(-sin(x)*cos(z) - sin(z)*cos(y)))
            return part1 + part2 + part3
        elif self.name == 'Diamond':
            part1 = -sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            part2 = +sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) + sin(z)*cos(x)*cos(y)
            part3 = +sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) + sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            part4 = +sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) + sin(z)*cos(x)*cos(y)
            part5 = -sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            part6 = +sin(x)*sin(y)*sin(z) + sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            part7 = +sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) + sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            part8 = +sin(x)*sin(y)*sin(z) + sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            part9 = -sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)
            parta = xp*(xp*part1 + yp*part2 + zp*part3)
            partb = xp*(xp*part4 + yp*part5 + zp*part6)
            partc = xp*(xp*part7 + yp*part8 + zp*part9)
            return xp*(xp*(-sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)) + yp*(sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) + sin(z)*cos(x)*cos(y)) + zp*(sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) + sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y))) + yp*(xp*(sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) + sin(z)*cos(x)*cos(y)) + yp*(-sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)) + zp*(sin(x)*sin(y)*sin(z) + sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y))) + zp*(xp*(sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) + sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)) + yp*(sin(x)*sin(y)*sin(z) + sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)) + zp*(-sin(x)*sin(y)*sin(z) - sin(x)*cos(y)*cos(z) - sin(y)*cos(x)*cos(z) - sin(z)*cos(x)*cos(y)))
        else:
            raise ValueError(f"TPMS: {self.name} not implemented")

    def get_pol(self, k0):
        f0,fp0,fpp0 = self.f(k0), self.fp(k0), self.fpp(k0)
        a = fpp0/2
        b = fp0 - 2*a*k0
        c = f0 - a*k0**2 -b*k0
        return a,b,c

    def brut_root(self, Nb = 100000, rang = 10.):
        t = np.linspace(0.,rang,Nb)
        truth = self.f(t)>0.
        roots = np.array([x for x in range(0,Nb-1)])[truth[1:]!=truth[:-1]]
        if len(roots) > 0:
            return rang*(roots[0]/Nb)
        else:
            return None

    def get_sample_step_length(self):
        u,v,w = self.direction
        if self.name == 'Schwarz_P':
            return 0.5* self.pitch * 0.25/np.max(np.abs([u,v,w]))
        if self.name == 'Gyroid':
            return 0.5* self.pitch * 0.25/np.max(np.abs([u+v,u-v,v+w,v-w,w+u,w-u]))
        if self.name == 'Diamond':
            return 0.5* self.pitch * 0.25/np.max(np.abs([u+v+w,u+v-w,u-v+w,u-v-w]))
        else:
            return 0.5* self.pitch * 0.25/np.max(np.abs([u,v,w,u+v,u-v,v+w,v-w,w+u,w-u,u+v+w,u+v-w,u-v+w,u-v-w]))
        # else:
        #     raise ValueError(f"TPMS: {self.name} not implemented")

    @classmethod
    def from_seed(cls, name:str, pitch:float=1., seed=None) -> 'TPMS':
        if seed is not None:
            np.random.seed(seed)
        angle0 = np.random.uniform(-pi,+pi, 1)[0]
        angle1 = np.random.uniform(-pi/2,+pi/2, 1)[0]
        position  = [float(x) for x in np.random.uniform(-1.,+1., 3)]
        direction = [float(x) for x in [cos(angle1)*cos(angle0),cos(angle1)*sin(angle0),sin(angle1)]]
        tpms= cls(name, pitch, position, direction)
        return tpms

class Solver:

    XLIM = 1e6

    def __init__(self, tpms:TPMS) -> None:
        self.tpms = tpms
        self.w0 = tpms.get_sample_step_length()

    def root_in_interval(self, L0:float, L1: float, itermax = 10):
        f, fp, fpp = self.tpms.f,self.tpms.fp,self.tpms.fpp
        assert L0 <= L1
        xa = L0
        xb = L1
        isRoot = None
        ITER = 0
        while isRoot is None and ITER < itermax:
            fa = f(xa)
            fpa = fp(xa)
            fppa = fpp(xa)
            fb  = f(xb)
            fpb = fp(xb)
            fppb = fpp(xb)
            # print(f"[{L0:1f}-{L1:1f}] XA: {xa:3f} XB: {xb:3f} - FA {fa:+.2E} FPA {fpa:+.2E} FPPA {fppa:+.2E} FB {fb:+.2E} FPB {fpb:+.2E} FPPB {fppb:+.2E}") # uncomment to debug pathological cases
            if fa*fb < 0.: # if the two sampled point are of different sign, there is a root in the interval
                isRoot = True
            elif fpa*fpb > 0.: # if the two root are of same sign and the two derivative are of same sign, there is no root in the interval
                isRoot = False
            elif fppa*fppb < 0.: # If inflexion change, split the interval in half to study root existance. (rare case)
                inInterval1, xa1, xb1 = self.root_in_interval(L0, 0.5*(L0+L1), itermax)
                if inInterval1:
                    return inInterval1, xa1, xb1
                else:
                    return self.root_in_interval(0.5*(L0+L1), L1, itermax)
            elif fppa*fa < 0.: # If no inflexion change and the sign of the second derivative is not the one of the sampled points, there is no root
                isRoot = False
            else: # Otherwise, we don't know if there is a root
                polX = np.linspace(L0,L1,100)
                polY1 = fpa * polX + fa - fpa*xa # tangeante au point1
                polY2 = fpb * polX + fb - fpb*xb # tangeante au point2
                xn = -(fa - fpa*xa - fb + fpb*xb)/(fpa - fpb) # x de l'intersection des droites
                fn = fpa * xn + fa - fpa*xa # y de l'intersection des droites
                fxn = f(xn) # evaluation de la tpms en x
                if fa*fn > 0.: # if the intersection ordinate is of same sign as the sampled points, there is no root
                    isRoot = False
                elif fa*fxn < 0.:
                    isRoot = True # if the evaluate function in x is of different sign of the sampled points, there is a root
                    xb = xn
                else: # Otherwise, we don't know and need to refine.
                    fpn = fp(xn)
                    if fpn*fpa <0.:
                        xb = xn
                    else:
                        xa = xn
        if isRoot is None:
            isRoot = False
        return isRoot, xa, xb

    def get_interval(self, L0=None):
        w0 = self.w0
        if L0 is None:
            L0, L1 = 0., w0
        else:
            L1 = L0 + w0
        while L0 < Solver.XLIM:
            rootInInterval, xa, xb = self.root_in_interval(L0, L1)
            if rootInInterval:
                break
            else:
                L0+=w0
                L1+=w0
        return xa, xb

    def solve(self, L0=None):
        f = self.tpms.f
        x0, x1 = self.get_interval(L0)
        root = scipy.optimize.root_scalar(f, x0=x0, x1=x1, method='bisect',bracket=(x0,x1)).root
        return root


# INPUT ************************************************************************
seeds = [0+i for i in range(0,10000)] 
# ******************************************************************************

skips, oks, fails = [], [], []

for i, seed in enumerate(seeds):
    tpms = TPMS.from_seed('Diamond', 10.5, seed)
    brut = tpms.brut_root(Nb = 100000, rang = 10.)
    sol = Solver(tpms)
    root = sol.solve()
    # PRINT
    if brut is None:
        skips.append(seed)
        print(f"{i:5d} - {seed}: ROOT: >100.        vs. FOUND: {root:.6E} difference: ----------",'\033[96m',"SKIP",'\033[0m')
    elif abs(root-brut)<2e-4:
        oks.append(seed)
        print(f"{i:5d} - {seed}: ROOT: {brut:.6E} vs. FOUND: {root:.6E} difference: {root-brut:+.3E}",'\033[92m',"OK",'\033[0m')
    else:
        fails.append(seed)
        print(f"{i:5d} - {seed}: ROOT: {brut:.6E} vs. FOUND: {root:.6E} difference: {root-brut:+.3E}",'\033[91m',"FAIL",'\033[0m')
        
print(fails, skips)