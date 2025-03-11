# distutils: language = c++
import numpy as np
from scipy.special import spence
from libcpp.vector cimport vector
import gepard as g
from gepard.fits import th_KM15

cdef double M = 0.938272081  # Proton mass (GeV)
cdef double me = 0.5109989461 * 0.001  # Electron mass in GeV
cdef double ycolcut = 0.000001  # P1 cut
cdef double alpha = 1.0/137.036

cpdef double printKM(double xB, double Q2, double t, double phi,
                     int pol=0, double E=10.604, str model='km15'):
    """
    Computes cross section for KM15 or KM15 Bethe-Heitler model
    """
    phi = np.pi - phi
    pt1 = g.DataPoint(
        xB=xB, t=-t, Q2=Q2, phi=phi,
        process='ep2epgamma', exptype='fixed target', frame='trento',
        in1energy=E, in1charge=-1, in1polarization=pol
    )
    pt1.prepare()
    if model == 'km15':
        return th_KM15.XS(pt1)
    elif model == 'km15_bh':
        return th_KM15.PreFacSigma(pt1)*th_KM15.TBH2unp(pt1)


cpdef double nu(double xB, double Q2):
    return Q2/(2.0*M*xB)

cpdef double y(double xB, double Q2):
    return nu(xB, Q2)/10.604

cpdef double W2(double xB, double Q2):
    return M*M + 2.0*M*nu(xB, Q2) - Q2

cpdef double tmin2(double xB, double Q2):
    return -0.5 * (
        (Q2/xB - Q2)*(Q2/xB - np.sqrt((Q2/xB)**2 + 4.0*M*M*Q2))
        + 2.*M*M*Q2
    ) / W2(xB, Q2)

cpdef double tmax2(double xB, double Q2):
    return -0.5 * (
        (Q2/xB - Q2)*(Q2/xB + np.sqrt((Q2/xB)**2 + 4.*M*M*Q2))
        + 2.*M*M*Q2
    ) / W2(xB, Q2)

cpdef double del2(double t):
    return -t

cpdef double del2q2(double t, double Q2):
    return del2(t)/Q2

cpdef double eps(double xB, double Q2):
    return 2.0*xB*M/np.sqrt(Q2)

cpdef double eps2(double xB, double Q2):
    return eps(xB, Q2)**2

cpdef double qeps2(double xB, double Q2):
    return 1 + eps2(xB, Q2)

cpdef double sqeps2(double xB, double Q2):
    return np.sqrt(qeps2(xB, Q2))

cpdef double y1eps(double xB, double Q2):
    cdef double yd = y(xB, Q2)
    return 1 - yd - yd*yd*eps2(xB, Q2)/4.0

cpdef double Kfac2(double xB, double Q2, double t):
    cdef double tmind = -tmin2(xB, Q2)
    cdef double eps2d = eps2(xB, Q2)
    return (
        -del2q2(t, Q2)
        * (1 - xB)
        * y1eps(xB, Q2)
        * (1 - (tmind)/t)
        * (
            np.sqrt(1 + eps2d)
            + ((4*xB*(1 - xB) + eps2d)/(4*(1 - xB)))
              * (-(t - (tmind))/Q2)
        )
    )

cpdef double Kfac(double xB, double Q2, double t):
    return np.sqrt(Kfac2(xB, Q2, t))

cpdef double Jfac(double xB, double Q2, double t):
    cdef double yd = y(xB, Q2)
    cdef double eps2d = eps2(xB, Q2)
    cdef double del2q2d = del2q2(t, Q2)
    return (
        (1 - yd - yd*eps2d/2.0)*(1 + del2q2d)
        - (1 - xB)*(2 - yd)*del2q2d
    )

cpdef double P1(double xB, double Q2, double t, double phi):
    cdef double yd = y(xB, Q2)
    cdef double eps2d = eps2(xB, Q2)
    return -(
        Jfac(xB, Q2, t)
        + 2*Kfac(xB, Q2, t)*np.cos(np.pi - np.radians(phi))
    ) / (yd*(1 + eps2d))

cpdef double getScale(double xBmin, double xBmax,
                      double Q2min, double Q2max,
                      double tmin, double tmax,
                      double ymin, double ymax,
                      double w2min, int rad=0, double Ed=10.604):
    """
    Scans the phase space to find maximum XS for rejection sampling (unused by default).
    """
    cdef int nx=40
    cdef int nq=20
    cdef int nt=40

    cdef double dx = (xBmax - xBmin)/nx
    cdef double dq = (Q2max - Q2min)/nq
    cdef double dt

    cdef double xBd
    cdef double Q2d
    cdef double yd
    cdef double w2d
    cdef double td

    cdef int elPold
    cdef double phigd
    cdef double dstot = 0.0
    cdef double xs

    for elPold in [-1, 1]:
        for ix in range(1, nx+1):
            xBd = xBmin + ix*dx
            for iq in range(1, nq+1):
                Q2d = Q2min + iq*dq
                yd = y(xBd, Q2d)
                w2d = W2(xBd, Q2d)

                if (yd < ymin) or (yd > ymax):
                    continue
                if (w2d < w2min):
                    continue

                for it in range(1, nt+1):
                    tmax = min(-tmax2(xBd, Q2d), tmax)
                    tmin = max(-tmin2(xBd, Q2d), tmin)
                    if (tmax < tmin):
                        continue
                    dt = (tmax - tmin)/nt
                    td = tmin + it*dt
                    if (abs(P1(xBd, Q2d, td, 0)) < ycolcut):
                        continue

                    xs = printKM(xBd, Q2d, td, 0, pol=elPold, E=Ed)
                    if xs > dstot:
                        dstot = xs
    return dstot


cpdef str genOneEvent(double xBmin, double xBmax,
                      double Q2min, double Q2max,
                      double tmin, double tmax,
                      double ymin, double ymax,
                      double w2min, double xsmax,
                      int rad=0, double Ed=10.604,
                      str filename="km15gen", str model="km15"):
    """
    Generates one event in Lund format (KM15 model),
    returning its string representation.
    """
    cdef vector[double] kine
    cdef double cl_be, costel, afac, dE1, nud, Esc, dE2
    cdef double Eprime_el_e, E_el_e, eta, deltaEs
    cdef double delta_vertex, delta_vac, delta_R, delta_vvr
    cdef double deld, rho, aks, delta_1, Eprime_p
    cdef double pprime_p, delta_2, xs_born, xs_part, xs
    cdef double xBd_tr, Q2d_tr, nud_tr
    cdef double V3l1, V3l2, V3l3, V3gam1, V3gam2, V3gam3
    cdef double V3p1, V3p2, V3p3
    cdef double sintel, cosphe, sinphe
    cdef double vx, vy, vz
    cdef double radQ2, radxB, radEd
    cdef str result

    result = ""

    # small vertex smearing
    vx = 0.025*(np.random.rand() - 0.5)
    vy = 0.025*(np.random.rand() - 0.5)
    vz = np.random.rand()*5.0  # dummy offset for radiator simulation
    cl_be = Ed

    cdef int elPold = 2*np.random.randint(2) - 1

    cdef double xBd = xBmin + (xBmax - xBmin)*np.random.rand()
    cdef double Q2d = Q2min + (Q2max - Q2min)*np.random.rand()
    cdef double yd  = y(xBd, Q2d)
    cdef double w2d = W2(xBd, Q2d)
    cdef double td  = tmin + (tmax - tmin)*np.random.rand()
    cdef double phigd  = np.random.rand()*2.0*np.pi
    cdef double phield = np.random.rand()*2.0*np.pi

    radQ2d = Q2d
    radxBd = xBd

    if (yd < ymin) or (yd > ymax):
        return result
    if (w2d < w2min):
        return result

    nud = Q2d/(2.0*M*xBd)
    Esc = Ed - nud
    costel = 1.0 - Q2d/(2.0*Ed*Esc)

    # check t-range
    if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
        return result
    if (-P1(xBd, Q2d, td, phigd) < ycolcut):
        return result

    xs_born = printKM(xBd, Q2d, td, phigd, pol=elPold, E=Ed, model=model)
    if np.isnan(xs_born) or np.isinf(xs_born) or (xs_born == 0):
        return result

    # -------------- handle radiative corrections if enabled --------------
    if rad:
        Ed = Ed - Ed*np.random.rand()**(3.0/4.0/(vz/929.0 + 0.003/8.897))
        # subtract some energy for external radiator thickness

        nud = nud - (cl_be - Ed)
        if nud <= 0:
            return result

        Q2d = Q2d*(Ed/cl_be)
        xBd = Q2d/(2.0*M*nud)

        # re-check t
        if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
            return result
        if (-P1(xBd, Q2d, td, phigd) < ycolcut):
            return result

        # pick electron internal rad
        afac = alpha/np.pi*(np.log(Q2d/me**2) - 1.0)
        dE1 = np.random.rand()**(1.0/afac)*Ed
        E_el_e = Ed - dE1

        if dE1 >= nud:
            return result
        nud_tr = nud - dE1

        afac = alpha/np.pi*(np.log(Q2d/me**2) - 1.0)
        dE2 = np.random.rand()**(1.0/afac)*Esc
        Eprime_el_e = Esc + dE2

        if (dE1 + dE2 >= nud):
            return result
        nud_tr = nud_tr - dE2

        Q2d_tr = Q2d*Eprime_el_e/Esc*E_el_e/Ed
        xBd_tr = Q2d_tr/(2.0*M*nud_tr)

        if (xBd_tr > 1.0) or (xBd_tr < 0.0):
            return result

        if (td < -tmin2(xBd_tr, Q2d_tr)) or (td > -tmax2(xBd_tr, Q2d_tr)):
            return result
        if (-P1(xBd_tr, Q2d_tr, td, phigd) < ycolcut):
            return result

        xs_part = printKM(xBd_tr, Q2d_tr, td, phigd, pol=elPold,
                          E=(Ed - dE1), model=model)
        if np.isnan(xs_part) or np.isinf(xs_part) or (xs_part == 0):
            return result

        # vertex + vac. + ...
        delta_vertex = alpha/np.pi*(1.5*np.log(Q2d/me**2) - 2
                                    - 0.5*np.log(Q2d/me**2)**2
                                    + np.pi**2/6.0)
        delta_vac = alpha/np.pi*(2.0/3.0)*(-5.0/3.0 + 1.0*np.log(Q2d/me**2))
        delta_R = alpha/np.pi*(
            -0.5*np.log(Ed/Esc)**2
            + 0.5*np.log(Q2d/me**2)**2
            - np.pi**2/3.0
            + spence(1-((1+costel)/2.0))
        )

        delta_vvr = alpha/np.pi*(
            (3.0/2.0 + 2.0/3.0)*np.log(Q2d/me**2)
            - 28./9.0 - 0.5*np.log(Ed/Esc)**2 - np.pi**2/6.0
            + spence(1-((1+costel)/2.0))
        )

        xs = xs_part*np.exp(delta_vertex + delta_R)/(1 - delta_vac*0.5)**2
        if np.isnan(xs) or np.isinf(xs) or (xs == 0):
            return result

    # -------------- if radiative is off or after rad, proceed --------------
    else:
        if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
            return result
        if abs(P1(xBd, Q2d, td, phigd)) < ycolcut:
            return result
        xBd_tr = xBd
        Q2d_tr = Q2d
        xs = xs_born

    vz = vz - 5.5  # shift the vertex a bit more
    # ========== Now we actually compute final 4-vectors and write them ==========

    # retrieve final 3-vectors
    if rad:
        # in rad block we call getphoton with Ed - dE1, then tweak V3l*, etc.
        kine = getphoton(xBd_tr, Q2d_tr, td, phigd, phield, Ed=(Ed - dE1))
        V3l1, V3l2, V3l3, V3gam1, V3gam2, V3gam3, V3p1, V3p2, V3p3, costgg = kine

        # subtract dE2 in direction of scattered electron
        sintel = np.sqrt(1 - costel**2)
        cosphe = np.cos(phield)
        sinphe = np.sin(phield)
        V3l1 = V3l1 - dE2*sintel*cosphe
        V3l2 = V3l2 - dE2*sintel*sinphe
        V3l3 = V3l3 - dE2*costel
    else:
        # no rad
        kine = getphoton(xBd_tr, Q2d_tr, td, phigd, phield, Ed=Ed)
        V3l1, V3l2, V3l3, V3gam1, V3gam2, V3gam3, V3p1, V3p2, V3p3, costgg = kine

    # Now compute energies + masses
    cdef double El  = np.sqrt(V3l1**2 + V3l2**2 + V3l3**2 + me**2)
    cdef double Ep  = np.sqrt(V3p1**2 + V3p2**2 + V3p3**2 + M**2)
    cdef double Egam= np.sqrt(V3gam1**2 + V3gam2**2 + V3gam3**2 + 0.0)  # photon mass=0

    # ========== Write out event to LUND ==========

    # We'll have different # of particles depending on s/p peak, etc.
    # The code already uses "if (dE1>=0.1) and (dE2>=0.1): #both s and p" etc.
    # We replace only the final 7–14 columns with px,py,pz, E, mass, vx, vy, vz

    with open("{}.dat".format(filename), "a") as file_out:
        #  -------- Radiative block with 5 or 4 or 3 final states -----------
        if rad:
            # El, Ep, Egam computed above
            if (dE1 >= 0.1) and (dE2 >= 0.1):
                # 5-particle final
                # header line (#particles=5, etc.)
                result += f"5   1       1  0.0{elPold:>4}   11   {cl_be:.3f}   1       1      {xs:6f}\n"

                # 1) electron
                # index=1, lifetime=0, type=1, PID=11, parent=0, firstDaughter=4
                # px, py, pz, E, mass, vx, vy, vz
                result += (
                    f"1  0.0000  1   11   0    4"
                    f"  {V3l1: .4f}  {V3l2: .4f}  {V3l3: .4f}  {El: .4f}  {me: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # 2) proton
                result += (
                    f"2  0.0000  1  2212  0    0"
                    f"  {V3p1: .4f}  {V3p2: .4f}  {V3p3: .4f}  {Ep: .4f}  {M: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # 3) photon
                result += (
                    f"3  0.0000  1   22   0    0"
                    f"  {V3gam1: .4f}  {V3gam2: .4f}  {V3gam3: .4f}  {Egam: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # 4) 'photon' for s peak
                result += (
                    f"4  0.0000  1   22   0    0"
                    f"  0.0000  0.0000  {dE1: .4f}  {dE1: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # 5) 'photon' for p peak
                result += (
                    f"5  0.0000  1   22   0    0"
                    f"  {dE2*sintel*cosphe: .4f}  {dE2*sintel*sinphe: .4f}  {dE2*costel: .4f}"
                    f"  {dE2: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )

            elif (dE1 >= 0.1) and (dE2 < 0.1):
                # 4-particle final
                result += f"4   1       1  0.0{elPold:>4}   11   {cl_be:.3f}   1       1      {xs:6f}\n"
                # electron
                result += (
                    f"1  0.0000  1   11   0    2"
                    f"  {V3l1: .4f}  {V3l2: .4f}  {V3l3: .4f}  {El: .4f}  {me: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # proton
                result += (
                    f"2  0.0000  1  2212  0    0"
                    f"  {V3p1: .4f}  {V3p2: .4f}  {V3p3: .4f}  {Ep: .4f}  {M: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # photon
                result += (
                    f"3  0.0000  1   22   0    0"
                    f"  {V3gam1: .4f}  {V3gam2: .4f}  {V3gam3: .4f}  {Egam: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # extra 'photon'
                result += (
                    f"4  0.0000  1   22   0    0"
                    f"  0.0000  0.0000  {dE1: .4f}  {dE1: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )

            elif (dE1 < 0.1) and (dE2 >= 0.1):
                # 4-particle final
                result += f"4   1       1  0.0{elPold:>4}   11   {cl_be:.3f}   1       1      {xs:6f}\n"
                # electron
                result += (
                    f"1  0.0000  1   11   0    3"
                    f"  {V3l1: .4f}  {V3l2: .4f}  {V3l3: .4f}  {El: .4f}  {me: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # proton
                result += (
                    f"2  0.0000  1  2212  0    0"
                    f"  {V3p1: .4f}  {V3p2: .4f}  {V3p3: .4f}  {Ep: .4f}  {M: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # photon
                result += (
                    f"3  0.0000  1   22   0    0"
                    f"  {V3gam1: .4f}  {V3gam2: .4f}  {V3gam3: .4f}  {Egam: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # extra 'photon'
                result += (
                    f"4  0.0000  1   22   0    0"
                    f"  {dE2*sintel*cosphe: .4f}  {dE2*sintel*sinphe: .4f}  {dE2*costel: .4f}"
                    f"  {dE2: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )

            else:
                # (dE1 < 0.1) and (dE2 < 0.1) => 3-particle final
                result += f"3   1       1  0.0{elPold:>4}   11   {cl_be:.3f}   1       1      {xs:6f}\n"
                # electron
                result += (
                    f"1  0.0000  1   11   0    1"
                    f"  {V3l1: .4f}  {V3l2: .4f}  {V3l3: .4f}  {El: .4f}  {me: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # proton
                result += (
                    f"2  0.0000  1  2212  0    0"
                    f"  {V3p1: .4f}  {V3p2: .4f}  {V3p3: .4f}  {Ep: .4f}  {M: .4f}"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )
                # photon
                result += (
                    f"3  0.0000  1   22   0    0"
                    f"  {V3gam1: .4f}  {V3gam2: .4f}  {V3gam3: .4f}  {Egam: .4f}  0.0000"
                    f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
                )

        else:
            # ----------------------- Non-radiative block -----------------------
            # 3 final-state particles: electron, proton, photon
            result += f"3   1       1  0.0{elPold:>4}   11   {cl_be:.3f}   1       1      {xs:6f}\n"
            # electron
            result += (
                f"1  0.0000  1   11   0    1"
                f"  {V3l1: .4f}  {V3l2: .4f}  {V3l3: .4f}  {El: .4f}  {me: .4f}"
                f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
            )
            # proton
            result += (
                f"2  0.0000  1  2212  0    0"
                f"  {V3p1: .4f}  {V3p2: .4f}  {V3p3: .4f}  {Ep: .4f}  {M: .4f}"
                f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
            )
            # photon
            result += (
                f"3  0.0000  1   22   0    0"
                f"  {V3gam1: .4f}  {V3gam2: .4f}  {V3gam3: .4f}  {Egam: .4f}  0.0000"
                f"  {vx: .4f}  {vy: .4f}  {vz: .4f}\n"
            )

    return result


cpdef vector[double] bhdvcs(double xBd, double Q2d, double td,
                            double phigd, int elPold, int rad=0,
                            double Ed=10.604):
    """
    BH + DVCS cross section code (unused in main?).
    """
    cdef vector[double] vec
    cdef double costel, afac, dE1, nud, Esc, dE2, Eprime_e
    cdef double E_e, eta, deltaEs, delta_vertex, delta_vac, delta_R
    cdef double delta_vvr, deld, rho, aks, delta_1, Eprime_p
    cdef double pprime_p, delta_2, xs_born, xs

    if rad:
        afac = alpha/np.pi*(np.log(Q2d/me**2) - 1.0)
        dE1 = np.random.rand()**(1.0/afac)*Ed
        E_e = Ed - dE1

        nud = Q2d/(2.0*M*xBd)
        Esc = Ed - nud
        costel = 1.0 - Q2d/(2.0*E_e*Esc)

        Q2d = Q2d*(E_e/Ed)
        nud = nud - dE1
        xBd = Q2d/(2.0*M*nud)

        xs_born = printKM(xBd, Q2d, td, phigd, pol=elPold, E=(Ed - dE1))

        afac = alpha/np.pi*(np.log(Q2d/me**2) - 1.0)
        dE2 = np.random.rand()**(1.0/afac)*Esc
        Eprime_e = Esc - dE2

        delta_vertex = alpha/np.pi*(
            1.5*np.log(Q2d/me**2)
            - 2.0
            - 0.5*np.log(Q2d/me**2)**2
            + np.pi**2/6.0
        )
        delta_vac = alpha/np.pi*(2./3.)*(-5./3. + np.log(Q2d/me**2))
        delta_R = alpha/np.pi*(
            -0.5*np.log(E_e/Eprime_e)**2
            + 0.5*np.log(Q2d/me**2)**2
            - np.pi**2/3.0
            + spence(1 - ((1+costel)/2.0))
        )

        delta_vvr = alpha/np.pi*(
            (3./2. + 2./3.)*np.log(Q2d/me**2)
            - 28./9.
            - 0.5*np.log(E_e/Eprime_e)**2
            - np.pi**2/6.
            + spence(1 - ((1+costel)/2.0))
        )

        xs = xs_born*np.exp(delta_vertex + delta_R)/(1 - delta_vac*0.5)**2

        vec.push_back(xs)
        vec.push_back(xs_born)
        vec.push_back(xBd)
        vec.push_back(Q2d)
        vec.push_back(dE1)
        vec.push_back(dE2)
        return vec

    else:
        xs = printKM(xBd, Q2d, td, phigd, pol=elPold, E=Ed)
        vec.push_back(xs)
        return vec


cpdef vector[double] getphoton(double xBd, double Q2d, double td,
                               double phigd, double phield,
                               double Ed=10.604):
    """
    Returns the 3-vectors for electron, photon, and proton in the final state
    from the chosen kinematics. V3gam*, V3p*, etc.
    """
    cdef double nud = Q2d/(2.0*M*xBd)
    cdef double Esc = Ed - nud
    cdef double yb = nud/Ed
    cdef double costel = 1.0 - Q2d/(2.0*Ed*Esc)
    cdef double sintel = np.sqrt(1.0 - costel**2)
    cdef double cosphe = np.cos(phield)
    cdef double sinphe = np.sin(phield)

    # incoming beam 4-vector in z-direction
    cdef double V3k1 = 0.0
    cdef double V3k2 = 0.0
    cdef double V3k3 = Ed

    # scattered electron
    cdef double V3l1 = Esc*sintel*cosphe
    cdef double V3l2 = Esc*sintel*sinphe
    cdef double V3l3 = Esc*costel

    # q-vector
    cdef double V3q1 = V3k1 - V3l1
    cdef double V3q2 = V3k2 - V3l2
    cdef double V3q3 = V3k3 - V3l3

    cdef double Ep = M + td/(2.0*M)
    cdef double Egam = nud - td/(2.0*M)

    cdef double qmod = np.sqrt(V3q1**2 + V3q2**2 + V3q3**2)
    cdef double costVq = (Ed - Esc*costel)/qmod
    cdef double sintVq = np.sqrt(1.0 - costVq**2)

    cdef double costgg = (
        2.0*Egam*(M + nud)
        + Q2d
        - 2.0*M*nud
    ) / (2.0*Egam*qmod)
    cdef double sintgg = np.sqrt(1.0 - costgg**2)

    cdef double Vgx = Egam*sintgg*np.cos(phigd)
    cdef double Vgy = Egam*sintgg*np.sin(phigd)
    cdef double Vgz = Egam*costgg

    cdef double V3gam1 = (
        Vgx*costVq*cosphe
        - Vgz*sintVq*cosphe
        - Vgy*sinphe
    )
    cdef double V3gam2 = (
        Vgx*costVq*sinphe
        - Vgz*sintVq*sinphe
        + Vgy*cosphe
    )
    cdef double V3gam3 = Vgx*sintel + Vgz*costVq

    cdef double V3p1 = V3q1 - V3gam1
    cdef double V3p2 = V3q2 - V3gam2
    cdef double V3p3 = V3q3 - V3gam3

    cdef vector[double] vec
    vec.push_back(V3l1)
    vec.push_back(V3l2)
    vec.push_back(V3l3)
    vec.push_back(V3gam1)
    vec.push_back(V3gam2)
    vec.push_back(V3gam3)
    vec.push_back(V3p1)
    vec.push_back(V3p2)
    vec.push_back(V3p3)
    vec.push_back(costgg)

    return vec