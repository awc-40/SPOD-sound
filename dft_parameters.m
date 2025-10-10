function Sts = dft_parameters(opts)

    % DESCRIPTION ---------------------------------------------------------
    % Calculate the Strouhal vector given block parameters

    % INPUTS --------------------------------------------------------------
    % opts.dt:    Timestep between snaphots.
    % opts.nDFT:  Number of snapshots per block.
    % opts.D:     Nozzle diameter (m).
    % opts.U:     Jet exit speed (m/s).

    % OUTPUTS -------------------------------------------------------------
    % Sts:        Strouhal vector, defined such that St=f*D/U. 

    % FUNCTION BODY -------------------------------------------------------
    
    % 1) Define variables
    dt = opts.dt;
    nDFT = opts.nDFT;
    D = opts.D;
    U = opts.U;

    % 2) Calculate frequency interval and the Nyquist frequency
    df = 1/(nDFT*dt);
    fn = 1/(2*dt);

    % 3) Construct Strouhal vector
    f = 0:df:2*fn-df;
    f(nDFT/2+1:end) = f(nDFT/2+1:end) - 2*fn;
    Sts = f*D/U;

end