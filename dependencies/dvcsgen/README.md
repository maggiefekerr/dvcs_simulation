# dvcsgen
dvcs/pi0/eta  generator using pdfs and gpds. 

```
git clone https://github.com/JeffersonLab/dvcsgen.git
cd dvcsgen
make
```

Please cite the following articles for this generator.
```
V. A. Korotkov and W. D. Nowak, Eur. Phys. J. C 23, 455–
461 (2002).
I. Akushevich and A. Ilyichev, Phys. Rev. D 98, 013005 (2018).
```

The supplementary material of the following article contains the global fitting of the EM form factors.
```
Z. Ye, J. R. Arrington, R. J.Hill and G. Lee, Phys. Lett. B 777 (2018) 8–15
```

To get command line options `./dvcsgen --help`

Code requires GPD grid files. Set the `CLASDVCS_PDF` variable to point to directory where the
gpd.dat file is located (ex. CLASDVCS_PDF=/scratch/username/dvcsgen )

Example
`./dvcsgen --trig 10000000 --nmax 10000`
will write `gemc lund type` data files with 10K DVCS events (--pi0 for pi0, --eta for eta) in the current directory  with total of 10M events.

`./dvcsgen --trig 10000 --docker`
will write `gemc lund type` single data file dvcs.dat with 10K events 



```
 dvcspi0gen [options]
       option  value      default    comment
      --pi0                    exclusive pi-0 on
      --eta                    exclusive eta on
      --nodvcs                 DVCS off
      --v    verbos_level    0  additional printout
      --trig nevents  10      number of triggers
      --seed random_seed 0    use the local time
      --nmax nevents   40k     # of triggers per file
      --beam beam_energy   5.754 e- momentum in GeV
      --ycol P_1 cut        0.005      P_1>ycol_cut
      --x  x-min   x-max   0.1 0.65   min max x_Bj
      --q2 Q2-min Q2-max 1.0  10. min max Q2 in GeV^2
      --y y-min y-max 0.05 0.9    min max for y=nu/E
      --w w^2-min         4.0 min for w^2=M^2+2MyE-Q^2
      --t tmin tmax  0 1.0      t  min/max in GeV^2
      --th thmin thmax  0.2 1 theta min/max \theta rad
      --xpos x-position      0 beam x position in cm
      --ypos y-position      0 beam y position in cm
      --zpos z-position      0 target z position in cm
      --zwidth z-width 0  width z in cm (zpos+/-zwidth/2)
       --raster diameter 0.75   raster diameter in cm
      --weight   flat distributions with weight(part12)
      --printgpd               print gpds and exit
      --nont               do not write out the ntuple
      --file              dvcspi0gen   filename
      --gpd  Igpd 3  GPD model(1-A,2-B,3-C,4-D) 101 VGG+Dterm
      --scale  scale      1.0   scale the sinphi-mom
      --targ target       proton   deut/neut possible
      --lpol                    Long.pol.target
      --tpol                    Trans.pol.target
      --writef    format      0-lund12, 1-lundgsim, 2-radiative correction
      --mod     write-mod      0-all, 1-cut on events
      --mom                include moments in ntuple
      --proloss                  add proton loss
      --ktcor          FALSE   turn on k_t cor for A_LU
      --globalfit      FALSE   use the global fitting results of EM form factors
      --radgen                   include radgen
      --radext                   include external radiative corrections
      --radstable                use born cross sections for rejection sampling
      --nodat               do not write a data file
      --acce16            include e16 acceptance for e-
      --acceg1           include eg1 acceptance for e-
      --acc12         include clas++ acceptance for e-
      --smear                   smear moments
      --A    value         0.006   smear param-A
      --B    value         0.001  smear param-B
      --C    value         0.0008  smear param-C
      --D    value         0.001  smear param-D
      --nmax   value     2000  maximum events per file
      --print nprint     1000   print ev nprint event
      --bh  value      3 BH status:3-All, 1-only BH
      --delta  value      0.01 Minimum rad photon energy (GeV)'
      --vv2cut value      0.1 cuts on missing mass ep squared (GeV^2)'

```

For the writef 2 option, the file format is still lund, [https://gemc.jlab.org/gemc/html/documentation/generator/lund.html](https://gemc.jlab.org/gemc/html/documentation/generator/lund.html).

But the contents are changed for the radiative corrections.
So, ```--writef 2``` is only useful when ```--radgen``` is on.

The header's event weight is still the radiative cross section (weight of MC::Header).

The lund particles have three user-defined values that are not used by the dvcsgen and the gemc.
These are (2) lifetime, (10) energy, (11) mass, and will be saved in MC::Lund.

The electron: (2) xB, (6) radiation mode (1: nonrad, 2:s-peak, 3:p-peak), (10) Q2, (11) -t.

The proton: (2) phi (radians), (10) shifted xB of the virtual photon, (11) shifted Q2 of the virtual photon.

The photon: (2) actual beam electron energy (See ```External Radiative Corrections``` below), (11) born cross section.

## External Radiative Corrections

For the realistic radiative generators, the external radiative corrections are required for the incoming electron because the interactions between the beam electron and the target materials happen prior to the vertex. So far, the program is reducing the energy only, without considering the angle smearing from the multiple scattering. The energy loss formula is from the approximation by Mo & Tsai (1969) [https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.41.205](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.41.205). The materials depends on the target configuration. So far, this option assumes the RG-A target geometry only. One can turn on the external RC by adding ```--radext``` option. If ```--writef 2``` is turned on, the photon's (2) lifetime column will record the actual incoming electron energy.