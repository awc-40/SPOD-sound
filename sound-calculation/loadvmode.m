function [Q_blk_hat_fi] = loadvmode(m,fidx,bidx,md)

    % DESCRIPTION ---------------------------------------------------------
    % Loads velocity mode  at desired azimuthal mode m, frequency index
    % fidx and block index bidx.
    %
    % INPUTS --------------------------------------------------------------
    % m:    Azimuthal mode
    % fidx: Frequency index. Integer between 1 and 52. fidx=52 corresponds
    %       to St = 1.98. Only 52 saved modes. 
    % bidx: Block index. Integer between 1 and 22. Based on block
    %       parameters nDFT=256, novlp=128, nT=3000. 
    %
    % OUTPUTS -------------------------------------------------------------
    % Z:    Lighthill mode. Complex array of shape (nx,nr,2), where nx=576
    %       is the number of axial grid points and nr=81 is the number of 
    %       radial grid points. 
    %
    % FUNCTION BODY -------------------------------------------------------
    
    if strcmp(md,'r')
        base = ['~/AWC/SPOD-sound/velocity-modes/' ...
        'azi-dft-realisations/m' sprintf('%02d',m) '/nfft256_novlp128_nblks22/'];
    elseif strcmp(md,'c')
        base = ['~/AWC/SPOD-sound/velocity-modes/' ...
            'azi-dft-realisations-chvmodes/m' sprintf('%02d',m) '/nfft256_novlp128_nblks22/']; 
    end
    file = ['/fft_block' ...
        sprintf('%04d',bidx) '_freq' sprintf('%04d',fidx) '.mat'];
    filepath = [base,file];
    load(filepath,'Q_blk_hat_fi');

end