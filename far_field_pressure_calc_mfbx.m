function [pmfb,pmfbij,pmfbgij] = far_field_pressure_calc_mfbx(Z,m,Stidx,md,absx,chi,rs,xs,usewinx,opts)

    % DESCRIPTION ---------------------------------------------------------
    % Far field spectral pressure calculation p\hat(x,f) for one azimuthal mode (m),
    % frequency (f) and block (b) of the lighthill tensor, and one far-field position (x).
    %
    % Currently, this function ignores azimuthal variation of the observer
    % location. The observer azimuthal coordinate is fixed at phi=0. 
    %
    % This function is to be used in the folder 
    % ~/AWC/SPOD-sound/lighthill-modes-pure
    %
    % INPUTS --------------------------------------------------------------
    % Z:       Complex array of size (nx,nr,3), where nx is the number of axial
    %          grid points and nr is the number of radial grid points. Contains 
    %          the components of the Lighthill source for the chosen azimuthal 
    %          mode m, Strouhal index Stidx and block. The final 
    %          dimension has size 3 to accommodate T_xx, T_xy and T_yy.
    % m:       Azimuthal mode of the source. Nonnegative integer.
    % Stidx:   Strouhal index. Integer between 1 and length of the Strouhal
    %          vector.
    % md:      String. If md is 'r', calculate for the round jet. 'c' is
    %          for the chevron jet
    % absx:    The distance |x| between the nozzle and the observer in 
    %          meters.
    % chi:     Polar angle of the observer in degrees
    % rs:      Vector of dimensionless radial points of size (1,nr).
    % xs:      Vector of dimensionless axial points of size (nx,1).
    % usewinx: Boolean. 1 = use one-sided axial taper. 0 = don't use. 
    % opts:    See INPUT section of function 'dft_parameters' below.
    %
    % OUTPUTS -------------------------------------------------------------
    % pmfb:    Contribution of a chosen azimuthal and frequency mode and
    %          block of the Lighthill tensor to the spectral far-field pressure 
    %          p\hat(x,\omega). Complex scalar. 
    % pmfbij:  Contribution of each of the components T_xx, T_xy and T_yy
    %          to p\hat(x,\omega). Complex vector of size (1,3)
    % pmfbgij: The Lighthill source weighted with the Green's function.
    %          Complex array of size (nx,nr,3). 
    %
    % FUNCTION BODY -------------------------------------------------------

    % 1) Get frequency parameter and axial and radial wavenumbers
    U        = opts.U;               % Jet exit speed (m/s)
    D        = opts.D;               % Jet diameter (m)
    c0       = opts.c0;              % Ambient speed of sound (m/s)
    Sts      = dft_parameters(opts); % Sts = Strouhal vector.
    St       = Sts(Stidx);           % St corresponding to Stidx
    M        = U/c0;                 % Mach number of the jet
    k        = 2*pi*St*M;            % Dimensionless wavenumber (i.e k*D).
    omega    = (2*pi*St*U)/D;        % Dimensional angular frequency
    chi      = chi*pi/180;           % Convert polar angle to radians
    kx       = k*cos(chi);           % Axial wavenumber
    kr       = k*sin(chi);           % Radial wavenumber
    nx       = length(xs);           % Number of axial grid points
    nr       = length(rs);           % Numer of radial grid points

    % 2) Calculate Green's function frequency factor
    ffac = -(omega^2)/(2*(c0^2));
    if m>0
        ffac = ffac*2;               % Accounts for -m. *2 assumes plane symmetry
    end

    % 3) Calculate Green's function directivity factor
    Dir    = zeros([1,3]);           % Initialise
    x1     = absx*cos(chi); 
    x2     = absx*sin(chi);
    Dir(1) = x1*x1;
    Dir(2) = 2*x1*x2;                % 2* because T_xy = T_yx. 
    Dir(3) = x2*x2;
    Dir    = Dir/(absx^3);           % This is the term xi*xj/|x^3|

    % 4) Calculate Green's function phase factor
    pfac   = exp(-1i*k*absx)*((1i)^m);

    % 5) Calculate overall Green's function prefactor
    prefac = ffac*Dir*pfac;                 % (1,3)
    prefac = reshape(prefac,[1,1,3]);       % (1,1,3)

    % 6) Get trapezoidal quadrature weights and axial taper
    weight  = trapzWeightsPolar(rs,xs,md);  % (nx,nr)
    winx    = onesided_hann(nx);            % (nx,1)

    % 7) Compute the spectral far-field pressure
    Jmkr    = besselj(m,kr*rs);                                     % Bessel function. (1,nr)
    eikx    = exp(1i*kx*xs);                                        % Axial phase. (nx,1) 
    G       = (1/(2*pi))*prefac.*Jmkr.*eikx.*weight;                % Green's function. (nx,nr,3)
    if usewinx == 1
        G  = winx.*G;                                               % Optional use of axial taper. 
    end
    pmfbgij              = Z.*G;                                    % (nx,nr,3)
    pmfbij               = transpose(squeeze(sum(pmfbgij,[1,2])));  % Integrating over x and r. (1,3) 
    pmfb                 = sum(pmfbij);                             % Integrating over x and r and summing over i,j. 

end



function [weight] = trapzWeightsPolar(r,z,md)

    % DESCRIPTION ---------------------------------------------------------
    % Based on:
    % TRAPZWEIGHTSPOLAR Integration weight matrix for cylindical 
    % coordinates using trapazoidal rule OTS, 2015
    % 
    % INPUTS --------------------------------------------------------------
    % r:      Vector of radial coordinates.
    % x:      Vector of axial coordinates.
    % md:     String which is either 'r' or 'c', which controls whether to 
    %         apply the round blank or the chevron blank for the nozzle.
    %
    % OUTPUTS -------------------------------------------------------------
    % weight: Trapezoidal quadrature weight matrix of size (nx,nr), where
    %         nx is the length of x and nr is the length of r.
    %
    % FUNCTION BODY -------------------------------------------------------
    %
    % 1) Create the 1D weight in the radial direction
    nothetar = length(r);
    weight_thetar = zeros(nothetar,1);
    weight_thetar(1) = pi*( r(1) + (r(2)-r(1))/2)^2;
    for i=2:nothetar-1
        weight_thetar(i) = pi*( r(i) + (r(i+1)-r(i))/2 )^2 - pi*( r(i) - (r(i)-r(i-1))/2 )^2;
    end
    weight_thetar(nothetar) = pi*r(end)^2 - pi*( r(end) - (r(end)-r(end-1))/2 )^2;
    
    % 2) Create the 1D weight in the axial direction
    noz = length(z);
    weight_z = zeros(noz,1);
    weight_z(1) = (z(2)-z(1))/2;
    for i=2:noz-1
        weight_z(i) = (z(i)-z(i-1))/2 + (z(i+1)-z(i))/2;
    end
    weight_z(noz) = (z(noz)-z(noz-1))/2;
    
    % 3) Combine radial and axial weights
    weight    = weight_z*weight_thetar'; % (nx,nr)
    
    % 4) Blank out nozzle
    if strcmp(md,'r')
        idx=importdata('~/JER_round/JER_round/SPOD/blank_round_nozzle.dat');
    elseif strcmp(md,'c')
        idx=importdata('~/JER_chev/JER_chev/SPOD_modes/blank_chev_nozzle.dat');
    end
    
    for i=1:size(idx,1)
        weight(idx(i,1),idx(i,2))=0; 
    end

end




function w = onesided_hann(N, frac)

    % DESCRIPTION -------------------------------------------------------------
    % One-sided Hanning taper of length N. This is optionally when integrating
    % the Green's function with source if the source has significant amplitude
    % towards the edge of the domain. The taper smoothly rolls off the source
    % amplitude so that there is no sharp amplitude cutoff at the edge of the
    % domain. In this way, the aim of the taper is to reduce spectral leakage. 
    %
    % INPUTS ------------------------------------------------------------------
    % N:     Number of points in the taper
    % frac:  fraction of the record to taper at the end. A default value of
    %        0.95 is set for this variable.
    %
    % OUTPUTS -----------------------------------------------------------------
    % w:     Vector of size (N,1) containing the taper weights.
    %
    % FUNCTION BODY -----------------------------------------------------------
    if nargin<2, frac = 0.95; end
    M = max(1, round(frac*N));
    w = ones(N,1);
    i0 = N-M+1:N;
    t = linspace(0,1,numel(i0)).';
    w(i0) = 0.5*(1 + cos(pi*t));   % 1 → 0 over last M points

end