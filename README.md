## Fortran codes for spherical harmonic synthesis of all-element gravity field from global geopotential coefficient model
https://www.zcyphygeodesy.com/en/h-nd-120.html
## [Algorithm purpose]
    From global geopotential coefficient model, calculate the model value of the height anomaly (m), gravity anomaly (mGal), gravity disturbance (mGal), vertical deflection vector (ʺ, south, west), disturbing gravity gradient (E, radial), tangential gravity gradient vector (E, north, west) or Laplace operator (E).
    When the minimum and maximum degree n to be set is equal, the program calculates the contribution of the degree n geopotential coefficients to the anomalous gravity field element, which can be employed to analyze and evaluate the spectral and space properties of the geopotential coefficient model.
## [Main program for test entrance]
    GGMallelementgrfd.f90
    The record format of the input calculation point file: ID (point no / point name), longitude (decimal degrees), latitude (decimal degrees), ellipsoidal height (m)......
    The record format of the output file reslt.txt: Behind the record of the calculation point file, appends 9 columns of model values of anomalous field elements which inculde the height anomaly (m), gravity anomaly (mGal), gravity disturbance (mGal), vertical deflection vector (ʺ, south, west), disturbing gravity gradient (E, radial), tangential gravity gradient vector (E, north, west) and Laplace operator (E).
## (1) Algorithm module for spherical harmonic synthesis of anomalous gravity field elements
    RntGravFdpnm(nmin,maxn,rln,cnm,snm,gvm,GRS,pnm,dpt1,dpt2,gr)
    Input parameters: nmin,maxn - the minimum and maximum calculation degree
    Input parameters: cnm, snm - geopotential coefficients
    Input parameters: rln(3) - the spherical coordinates of calculation point in IERS.
    Input parameters: GRS(6) - gm,ae,j2,omega,1/f, default value.
    Returnparameters: gvm(8) - height anomaly (m), gravity anomaly (mGal), gravity disturbance (mGal), vertical deflection vector (ʺ, south, west), disturbing gravity gradient (E, radial), tangential gravity gradient vector (E, north, west) 
## (2) Calculation module for normal geopotential coefficients
    normdjn(GRS,djn)   !GRS(6) - gm,ae,j2,omega,1/f, default value
    cstonorm (maxn,djn,GRS,HC,HS)
## (3) Calculation module for the normal gravity field elements
    GNormalfd(BLH,NFD,GRS)
    Input：BLH(3) - Latitude, longitude (degree decimal) and ellipsoidal height (m) of the calculation point
    Return parameters: NFD(5) - the normal geopotential (m2/s2), normal gravity (mGal), normal gravity gradient (E), normal gravity line direction (', expressed by its north declination relative to the center of the Earth center of mass) or normal gravity gradient direction (', expressed by its north declination relative to the Earth center of mass)..
## (4) Algorithm module for normalized associative Legendre functions and their derivatives
    BelPnmdt(pnm,dpt1,dpt2,maxn,t)
    The normalized associated Legendre functions and thier derivatives. Improved Belikov recursion algorithm for mnm and non-singular recursive algorithm for derivative of pnm.
## (5) Calculation module for Legendre functions and their derivatives to ψ
    LegPn_dt2(pn,dp1,dp2,n,t) ! t=cos ψ
## (6) Calculation module for the basic parameters of the Earth ellipsoid
    ELLIPSOIDPARA(GRS)
    Input parameters: GRS(1), GRS(2), GRS(4) and one of GRS(3), GRS(5) and GRS(6).
    GRS(1) - Geocentric gravitational constant GM. 
    GRS(2) - Major semi axis of the Earth. 
    GRS(3) - Dynamic form factor.
    GRS(4) - Mean angular velocity omega (e-5/s) of the Earth.
    GRS(5) - ellipsoid flattening f.
    GRS(6) - Normal ellipsoid geopotential.
    Return parameters: GRS(1:6)
## (7) Algorithm module for transforming ellipsoid geodetic coordinates into spherical coordinates
    BLH_RLAT(GRS,BLH,RLAT)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] PAGravf4.5 User Reference https://www.zcyphygeodesy.com/en/
    7.2 Calculation formulas of Earth gravity field from geopotential coefficient model
    7.3 Algorithms of normalized associative Legendre function and its derivative
    It is suggested that at the poles of the earth, the vertical deflection westward and the horizontal gravity gradient westward are not defined, and the relative calculation results are meaningless.
    DOS executable test file and all input and output data.
