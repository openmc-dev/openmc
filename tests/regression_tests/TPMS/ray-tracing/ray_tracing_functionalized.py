import numpy as np
import canvaslib2 as can
import scipy
import scipy.optimize
import scipy.interpolate
import time
import matplotlib
import matplotlib.pyplot as plt
import warnings
try:
    matplotlib.use('TkAgg')
except:
    warnings.warn("failed to import 'TkAgg' as matplotlib backend. Falling back to default backend. Interactive display may not work.",ImportWarning)

cos, sin, pi = np.cos, np.sin, np.pi

def fPitch(x,y,z):
    return 0.08*z**2 + 1.

def fThick(x,y,z):
    return 0.3*(x**2 + y**2)**0.5

# def interp_3d(x,y,z, points, values):
#     def _interp_1d(v, V):
#         vlen = np.shape(V)[0]
#         lower = (v <= V[0])
#         higher = (v > V[-1])
#         iV = np.arange(vlen)
#         iv = np.interp(v, V, iV).astype(int)
#         iv0 = iv
#         iv1 = np.minimum(iv+1, np.full_like(iv, vlen-1))
#         v0, v1 = V[iv0], V[iv1]
#         rv0, rv1 = np.zeros_like(v), np.zeros_like(v)
#         rv0[higher], rv1[higher] = 0.5, 0.5
#         rv0[lower], rv1[lower] = 0.5, 0.5
#         middle = (higher==False)*(lower==False)
#         rv0[middle] = (v1[middle]-v[middle])/(v1[middle]-v0[middle])
#         rv1[middle] = (v[middle]-v0[middle])/(v1[middle]-v0[middle])
#         return iv0, iv1, rv0, rv1
#     X,Y,Z = points
#     if isinstance(x, float): x = np.array([x])
#     if isinstance(y, float): y = np.array([y])
#     if isinstance(z, float): z = np.array([z])
#     ix0, ix1, rx0, rx1 = _interp_1d(x, X)
#     iy0, iy1, ry0, ry1 = _interp_1d(y, Y)
#     iz0, iz1, rz0, rz1 = _interp_1d(z, Z)
#     v000 = rx0*ry0*rz0*values[ix0,iy0,iz0]
#     v001 = rx0*ry0*rz1*values[ix0,iy0,iz1]
#     v010 = rx0*ry1*rz0*values[ix0,iy1,iz0]
#     v011 = rx0*ry1*rz1*values[ix0,iy1,iz1]
#     v100 = rx1*ry0*rz0*values[ix1,iy0,iz0]
#     v101 = rx1*ry0*rz1*values[ix1,iy0,iz1]
#     v110 = rx1*ry1*rz0*values[ix1,iy1,iz0]
#     v111 = rx1*ry1*rz1*values[ix1,iy1,iz1]
#     return v000+v001+v010+v011+v100+v101+v110+v111

def interp_3d(x, y, z, points, values):
    def _interp_1d(v, V):
        vlen = np.shape(V)[0]
        iV = np.arange(vlen)
        iv = np.interp(v, V, iV).astype(int)
        # Combine iv0 and iv1 assignment into one step to avoid multiple array creations
        iv0 = iv
        iv1 = np.clip(iv + 1, 0, vlen - 1)  # Use np.clip to avoid np.minimum
        v0, v1 = V[iv0], V[iv1]
        # Pre-calculate conditions
        lower = (v <= V[0])
        higher = (v > V[-1])
        middle = ~lower & ~higher  # Use bitwise not and and for conditions
        # Initialize rv0 and rv1 to 0.5, then modify only the middle cases
        rv0 = np.full_like(v, 0.5)
        rv1 = np.full_like(v, 0.5)
        rv0[middle] = (v1[middle] - v[middle]) / (v1[middle] - v0[middle])
        rv1[middle] = (v[middle] - v0[middle]) / (v1[middle] - v0[middle])
        return iv0, iv1, rv0, rv1
    X, Y, Z = points
    # Convert x, y, z to arrays if they are not already
    if isinstance(x, float): x = np.array([x])
    if isinstance(y, float): y = np.array([y])
    if isinstance(z, float): z = np.array([z])
    # Interpolate each dimension
    ix0, ix1, rx0, rx1 = _interp_1d(x, X)
    iy0, iy1, ry0, ry1 = _interp_1d(y, Y)
    iz0, iz1, rz0, rz1 = _interp_1d(z, Z)
    # Use broadcasting to avoid repetitive multiplications
    rxy00 = rx0*ry0
    rxy01 = rx0*ry1
    rxy10 = rx1*ry0
    rxy11 = rx1*ry1
    v000 = rxy00*rz0*values[ix0,iy0,iz0]
    v001 = rxy00*rz1*values[ix0,iy0,iz1]
    v010 = rxy01*rz0*values[ix0,iy1,iz0]
    v011 = rxy01*rz1*values[ix0,iy1,iz1]
    v100 = rxy10*rz0*values[ix1,iy0,iz0]
    v101 = rxy10*rz1*values[ix1,iy0,iz1]
    v110 = rxy11*rz0*values[ix1,iy1,iz0]
    v111 = rxy11*rz1*values[ix1,iy1,iz1]
    # Sum up all contributions
    return v000 + v001 + v010 + v011 + v100 + v101 + v110 + v111

# ['Gyroid','Double_Gyroid1','Double_Gyroid2','G_prime1','G_prime2','Lidinoid',
#          'Diamond','Double_Diamond1','Double_Diamond2','D_prime','Schwarz_P','Double_Schwarz_P','IWP1',
#          'Neovius','Octo1','Octo2','PN','KP','FRD1','FRD2','Split_P','IWP2','L-Type','Skeletal_1',
#          'Skeletal_2','Tubular_G','Tubular_P','I2-Y','G','Double_Diamond3','Fischer-Koch_S','IWP3']

class TPMS:

    TYPES:list = ['fSchwarz_P']

    def __init__(self, name:  str, position:np.ndarray, direction:np.ndarray, xref, yref, zref, mPitch, mThick) -> None:
        assert name in TPMS.TYPES
        assert abs(direction[0]**2 + direction[1]**2 + direction[2]**2 - 1.) < 1e-15
        assert zpitch0 > 0.
        assert zpitch1 > 0.
        self.name = name
        self.position = position
        self.direction = direction
        self._xref, self._yref, self._zref = xref, yref, zref
        self._mPitch = mPitch
        self._mThick = mThick
        self.interPitch = scipy.interpolate.RegularGridInterpolator((xref, yref, zref), mPitch)
        self.interThick = scipy.interpolate.RegularGridInterpolator((xref, yref, zref), mThick)

    def get_ABCDEF(self) -> tuple[float]:
        return self.direction[0],self.position[0],self.direction[1],self.position[1],self.direction[2],self.position[2]

    def get_pitch(self, x, y, z):
        #pitches = interp_3d(x,y,z, (self._xref, self._yref, self._zref), self._mPitch)
        if isinstance(x, float): x = np.array([x])
        if isinstance(y, float): y = np.array([y])
        if isinstance(z, float): z = np.array([z])
        x = np.clip(x, self._xref[0], self._xref[-1])
        y = np.clip(y, self._yref[0], self._yref[-1])
        z = np.clip(z, self._zref[0], self._zref[-1])
        pts = np.array([[x[i], y[i], z[i]] for i in range(0,np.shape(x)[0])])
        pitches = self.interPitch(pts)
        return pitches

    def get_thickness(self, x, y, z):
        #thicknesses = interp_3d(x,y,z, (self._xref, self._yref, self._zref), self._mThick)
        if isinstance(x, float): x = np.array([x])
        if isinstance(y, float): y = np.array([y])
        if isinstance(z, float): z = np.array([z])
        x = np.clip(x, self._xref[0], self._xref[-1])
        y = np.clip(y, self._yref[0], self._yref[-1])
        z = np.clip(z, self._zref[0], self._zref[-1])
        pts = np.array([[x[i], y[i], z[i]] for i in range(0,np.shape(x)[0])])
        thicknesses = self.interThick(pts)
        return thicknesses

    def f(self, k:float) -> float:
        A,B,C,D,E,F = self.get_ABCDEF()
        x = A*k+B
        y = C*k+D
        z = E*k+F
        pitch = self.get_pitch(x, y, z)
        thickness = self.get_thickness(x, y, z)
        l = 2*pi/pitch
        cst = thickness*2*pi/pitch
        if self.name == 'fSchwarz_P':
            return cos(l*x) + cos(l*y) + cos(l*z) - cst
        else:
            raise ValueError(f"TPMS: {self.name} not implemented")

    def fp(self, k:float):
        dk = 1.e-9
        if self.name == 'fSchwarz_P':
            return (self.f(k+dk) - self.f(k-dk)) / (2*dk)
        else:
            raise ValueError(f"TPMS: {self.name} not implemented")

    def fpp(self, k:float):
        dk = 1.e-5
        if self.name == 'fSchwarz_P':
            return (self.fp(k+dk) - self.fp(k-dk)) / (2*dk)
        else:
            raise ValueError(f"TPMS: {self.name} not implemented")

    def brut_root(self, Nb = 100000, rang = 10.):
        t = np.linspace(0.,rang,Nb)
        truth = self.f(t)>0.
        roots = np.array([x for x in range(0,Nb-1)])[truth[1:]!=truth[:-1]]
        if len(roots) > 0:
            return rang*(roots[0]/Nb)
        else:
            return None

    def get_pitch_min(self):
        return np.min(self._mPitch)

    def get_sample_step_length(self):
        u,v,w = self.direction
        if self.name == 'fSchwarz_P':
            return 0.5* self.get_pitch_min() * 0.25/np.max(np.abs([u,v,w]))
        # else:
        #     raise ValueError(f"TPMS: {self.name} not implemented")

    @classmethod
    def from_seed(cls, name:str, seed=None, xref=np.array([0.]), yref=np.array([0.]), zref=np.array([0.]), mPitch=np.array([[[2.]]]), mThick=np.array([[[0.]]])) -> 'TPMS':
        if seed is not None:
            np.random.seed(seed)
        angle0 = np.random.uniform(-pi,+pi, 1)[0]
        angle1 = np.random.uniform(-pi/2,+pi/2, 1)[0]
        position  = [float(x) for x in np.random.uniform(-1.,+1., 3)]
        direction = [float(x) for x in [cos(angle1)*cos(angle0),cos(angle1)*sin(angle0),sin(angle1)]]
        tpms= cls(name, position, direction, xref, yref, zref, mPitch, mThick)
        return tpms

class Solver:

    XLIM = 1e2

    def __init__(self, tpms:TPMS) -> None:
        self.tpms = tpms
        self.w0 = tpms.get_sample_step_length()

    def root_in_interval(self, L0:float, L1: float, interactive, itermax = 10):
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
            if interactive:
                print(f"[{L0:1f}-{L1:1f}] XA: {xa:3f} XB: {xb:3f} - FA {float(fa):+.2E} FPA {float(fpa):+.2E} FPPA {float(fppa):+.2E} FB {float(fb):+.2E} FPB {float(fpb):+.2E} FPPB {float(fppb):+.2E}")
            if fa*fb < 0.: # if the two sampled point are of different sign, there is a root in the interval
                isRoot = True
            elif fpa*fpb > 0.: # if the two root are of same sign and the two derivative are of same sign, there is no root in the interval
                isRoot = False
            elif fppa*fppb < 0.: # If inflexion change, split the interval in half to study root existance. (rare case)
                inInterval1, xa1, xb1 = self.root_in_interval(L0, 0.5*(L0+L1), interactive, itermax)
                if inInterval1:
                    return inInterval1, xa1, xb1
                else:
                    return self.root_in_interval(0.5*(L0+L1), L1, interactive, itermax)
            elif fppa*fa < 0.: # If no inflexion change and the sign of the second derivative is not the one of the sampled points, there is no root
                isRoot = False
            else: # Otherwise, we don't know if there is a root
                polX = np.linspace(L0,L1,100)
                polY1 = fpa * polX + fa - fpa*xa # tangeante au point1
                polY2 = fpb * polX + fb - fpb*xb # tangeante au point2
                xn = -(fa - fpa*xa - fb + fpb*xb)/(fpa - fpb) # x de l'intersection des droites
                fn = fpa * xn + fa - fpa*xa # y de l'intersection des droites
                fxn = f(xn) # evaluation de la tpms en x
                self._add_plot(polX, polY1, interactive, color = 3, lw=0.5)
                self._add_plot(polX, polY2, interactive, color = 3, lw=0.5)
                self._add_plot([xn, xn], [fn,fxn], interactive, color = 4, marker="+", markersize=3, ls= ' ')
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

    def get_interval(self, interactive, L0=None):
        w0 = self.w0
        if L0 is None:
            L0, L1 = 0., w0
        else:
            L1 = L0 + w0
        self._add_plot([L0, L1],[self.tpms.f(x) for x in [L0, L1]],interactive, marker='x', color = 2, markersize = 3, ls=' ')
        while L0 < Solver.XLIM:
            rootInInterval, xa, xb = self.root_in_interval(L0, L1, interactive)
            if rootInInterval:
                return xa, xb
            else:
                L0+=w0
                L1+=w0
                self._add_plot([L1],[self.tpms.f(L1)],interactive, marker='x', color = 2, markersize = 3, ls=' ')
        return None, None

    def solve(self, interactive=True, L0=None):
        if interactive:
            self._baseplot()
        f = self.tpms.f
        x0, x1 = self.get_interval(interactive, L0)
        if x0 is not None:
            root = scipy.optimize.root_scalar(f, x0=x0, x1=x1, method='bisect',bracket=(x0,x1)).root
        else:
            root = None
        self._add_plot([root],[0.],interactive, marker='o', color = 1, markersize = 5, ls=' ')
        return root

    def _baseplot(self):
        xmin, xmax = self.IPLOT.ax.get_xlim()
        X = np.linspace(xmin, xmax, 500)
        f, w0 = self.tpms.f, self.w0
        self.IPLOT.clear()
        self.IPLOT.interactive_plot(X, f(X))
        nSample = int(xmax/w0)+1
        xSample = np.arange(start=0.,stop=nSample*w0,step=w0)
        self.IPLOT.interactive_plot(xSample, f(xSample),marker='x',markersize=3,ls=' ',color=5)
        self.IPLOT.display()

    def _add_plot(self, X, Y, interactive, **kwargs):
        if interactive:
            self.IPLOT.interactive_plot(X, Y, **kwargs)
            self.IPLOT.display()

# INPUT ************************************************************************
interactive = False
seeds = [1000+i for i in range(0,10000)] # 1 Million de cas
zpitch0=3.
zpitch1=1.
zthick0=-0.5
zthick1=+0.5
# ******************************************************************************

skips, oks, fails = [], [], []


t0 = time.time()
xg, yg, zg = np.linspace(-5.,+5.,11, endpoint=True),np.linspace(-5.,+5.,11, endpoint=True),np.linspace(-5.,+5.,11, endpoint=True)
X,Y,Z = np.meshgrid(xg,yg,zg)
mpitch = fPitch(X,Y,Z)
mthick = fThick(X,Y,Z)
for i, seed in enumerate(seeds):
    tpms = TPMS.from_seed('fSchwarz_P', seed, xg, yg, zg, mpitch, mthick)
    brut = tpms.brut_root(Nb = 100000, rang = 10.)
    sol = Solver(tpms)
    root = sol.solve(interactive)
    # PRINT
    if root is not None:
        if brut is None:
            skips.append(seed)
            print(f"{i:5d} - {seed}: ROOT: >100.        vs. FOUND: {root:.6E} difference: ----------",'\033[96m',"SKIP",'\033[0m')
        elif abs(root-brut)<2e-4:
            oks.append(seed)
            print(f"{i:5d} - {seed}: ROOT: {brut:.6E} vs. FOUND: {root:.6E} difference: {root-brut:+.3E}",'\033[92m',"OK",'\033[0m')
        else:
            fails.append(seed)
            print(f"{i:5d} - {seed}: ROOT: {brut:.6E} vs. FOUND: {root:.6E} difference: {root-brut:+.3E}",'\033[91m',"FAIL",'\033[0m')
        if interactive:
            input()
    else:
        print(f"{i:5d} - {seed}: ROOT: >100.        vs. FOUND: >100.        difference: ----------",'\033[96m',"SKIP",'\033[0m')
t1 = time.time()
print(f"TIME: {t1-t0}")
print(fails, skips)