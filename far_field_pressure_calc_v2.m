function [pmfb,pmfbij,pmfbgij] = far_field_pressure_calc_v2(Z,m,Stidx,md,absx,chis,rs,xs,usewinx,opts)

    % DESCRIPTION ---------------------------------------------------------
    % Version 2 of the far field pressure calculation of the full (no SPOD)
    % lighthill modes. 
    %
    % Calculates the contribution of a chosen azimuthal and frequency mode and 
    % block of the Lighthill tensor to the spectral far-field pressure 
    % p\hat(x,\omega). The observer location x is a vector of positions,
    % all at a chosen fixed distance from the nozzle, and spanning a range
    % of chosen polar angles. The minimum polar angle is 20°, while the 
    % maximum is 160°. 
    %
    % To obtain the total far-field pressure for this frequency, one should
    % sum over all available modes.
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
    % chis:    Vector of observer polar angles in degrees. The minimum and 
    %          maximum angles that should be used are 20° and 160°.
    % rs:      Vector of dimensionless radial points of size (1,nr).
    % xs:      Vector of dimensionless axial points of size (nx,1).
    % usewinx: Boolean. 1 = use one-sided axial taper. 0 = don't use. 
    % opts:    See INPUT section of function 'dft_parameters' below.

    % OUTPUTS -------------------------------------------------------------
    % pmfb:    Contribution of a chosen azimuthal and frequency mode and
    %          block of the Lighthill tensor to the spectral far-field pressure 
    %          p\hat(x,\omega). Complex vector of size (n,1), where n is the 
    %          length of the vector chis. 
    % pmfbij:  Contribution of each of the components T_xx, T_xy and T_yy
    %          to p\hat(x,\omega). Complex matrix of size (n,3)
    % pmfbgij: The Lighthill source weighted with the Green's function.
    %          Complex array of size (n,nx,nr,3). 
    %
    % FUNCTION BODY -------------------------------------------------------
    
    % 1) Get frequency parameter and axial and radial wavenumbers
    n        = length(chis);         % Number of observer polar angles
    nx       = length(xs);           % Number of axial grid points
    nr       = length(rs);           % Numer of radial grid points

    % 7) Loop over the polar angles
    pmfb       = zeros([n,1]);       % Initialise
    pmfbij     = zeros([n,3]);       % Initialise
    pmfbgij    = zeros([n,nx,nr,3]); % Initialise
    for chidx  = 1:n
        chi                   = chis(chidx);
        [Pmfb,Pmfbij,Pmfbgij] = far_field_pressure_calc_mfbx(Z,m,Stidx,md,absx, ...
                                                            chi,rs,xs,usewinx,opts);
        pmfbgij(chidx,:,:,:)  = Pmfbgij;       % Assign
        pmfbij(chidx,:)       = Pmfbij;        % Assign
        pmfb(chidx)           = Pmfb;          % Assign
    end 

end