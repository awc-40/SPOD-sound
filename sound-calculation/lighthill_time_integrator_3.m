function p = lighthill_time_integrator_3(X,theta,mup,opts,varargin)

    % DESCRIPTION -------------------------------------------------------------
    % See the document 'lighthill-time-domain' for more details. This function
    % calculates the timeseries the double derivative of which gives the
    % far-field pressure timeseries using Lighthill's acoustic analogy in a
    % cylindrical domain.
    
    % INPUTS ------------------------------------------------------------------
    % X:       Observer distance in meters
    % theta:   Observer polar angle in degrees
    % mup:     Highest azimuthal mode of Lighthill tensor to include
    % opts.xs: 1D vector of dimensional axial grid points of shape (nx,1)
    % opts.rs: 1D vector of dimensional radial gridpoints of shape (1,nr)
    % opts.as: 1D vector of azimuthal gridpoints of shape (1,na) in radians
    % opts.ts: 1D vector of timepoints of shape (nt,1) in seconds
    % opts.D:  Nozzle diameter in meters
    % opts.pr: Progress display control. If opts.pr == 0, use a waitbar - this is for GUI. 
    %          If opts.pr == 1, display progress bar in terminal. 
    % meth:    varargin{1}. Temporal interpolation methd for the retarded time. Full
    %          range of options detailed in MATLAB function interp1.
    %          Default unspecified is linear, but can specify e.g
    %          'spline','cubic' etc. 
    
    
    % OUTPUTS -----------------------------------------------------------------
    % p:     1D vector, the far-field pressure timeseries
    
    % FUNCTION BODY -----------------------------------------------------------

    nvarargin = length(varargin);
    
    % Initial parameters
    theta_r = deg2rad(theta);
    c0      = 340.015;
    nx      = length(opts.xs);
    nr      = length(opts.rs);
    ncols   = nx*nr;
    na      = length(opts.as);
    nt      = length(opts.ts);
    dt      = opts.ts(2)-opts.ts(1);
    T       = opts.ts(end);
    weight  = trapzWeightsPolarnoblank(opts.rs,opts.xs);
    winx    = hamming(nx); 
    ww      = 1/mean(winx);
    winx    = ww*winx;
    weight  = winx.*weight;
    bs      = '~/AWC/SPOD-sound/lighthill-modes-pure/azimodes-cyl-ss/';

    % Calculate the separation vector Rv, its magnitude R and all observer
    % times t_o. Calculate also the cylindrical components Rx and Rr
    Rv = zeros([3,nx,nr,na]);
    Xv = [X*cos(theta_r),X*sin(theta_r),0];
    
    for i = 1:nx
        for j = 1:nr
            for l = 1:na
                xp = opts.xs(i);
                rp = opts.rs(j);
                ap = opts.as(l);
    
                y  = [xp,rp*cos(ap),rp*sin(ap)];
                Rv(:,i,j,l) = Xv - y;
            end
        end
    end

    Rx  = Rv(1,:,:,:); % (1,nx,nr,na)
    R   = sqrt(sum(Rv.^2,1)); % (1,nx,nr,na)
    Rr  = sqrt(R.^2 - Rx.^2); % (1,nx,nr,na)
    t_o = reshape(squeeze(R/c0),[ncols,na]); % (ncols,na)

    % Calculate indices of the limits of the pressure timeseries
    t0_max  = max(t_o,[],'all');
    t0_min  = min(t_o,[],'all');
    t_aug   = opts.ts;
    bl      = t_aug(end) < T + t0_min;
    while bl
        t_aug(end+1) = t_aug(end) + dt;
        bl = t_aug(end) < T + t0_min;
    end
    n       = find(opts.ts > t0_max,1,'first');
    s       = length(t_aug) - 1; 
    k       = s-n+1; % length of pressure timeseries
    cq      = repmat(1:ncols,[k,1]); % Index grid for 2D interpolant (k,ncols)
    tseg    = t_aug(n:s); % (k,1)

    % Calculate prefactor
    g = 1/(4*pi*(c0^2)*X); % Don't need opts.D^3 because opts.xs and opts.rs are dimensional

    % Initialise S_i
    S = zeros([k,1]);

    % -------------------- PROGRESS DISPLAY SETUP -----------------------------
    useGUI = ~isempty(opts.pr) && opts.pr == 0;
    useTTY = ~isempty(opts.pr) && opts.pr == 1;

    totalIts = (mup+1) * na;   % total (m,l) combinations
    it       = 0;              % completed iterations counter

    if useGUI
        updateEvery = 1;  % increase if needed
        hwb = waitbar(0,'Starting...','Name','Lighthill time integrator');
        cwb = onCleanup(@() safeCloseWaitbar(hwb)); %#ok<NASGU>
    elseif useTTY
        pb = termProgressBar(totalIts, ...
            'Prefix', 'Lighthill', ...
            'BarWidth', 34, ...
            'UseUnicode', true, ...
            'MinUpdatePeriod', 0.10); % seconds; prevents excessive printing
        pb.start();
    end
    % ------------------------------------------------------------------------

    % Loop over azimuthal wavenumbers
    for m = 0:mup
        mfl = ['mode' sprintf('%02d',m) 'tp.mat'];
        load([bs,mfl],'q'); disp(['Loaded m = ' num2str(m)])
        mn  = mean(q);
        q   = q - mn;
        q   = reshape(q,[nt,ncols,6]);

        % Form the gridded interpolant
        F = cell(1,6);
        for c = 1:6
            if nvarargin == 0
                F{c} = griddedInterpolant({opts.ts,1:ncols},q(:,:,c),'linear','none');
            elseif nvarargin == 1
                F{c} = griddedInterpolant({opts.ts,1:ncols},q(:,:,c),varargin{1},'none');
            end
        end

        for l = 1:na

            % -------------------- PROGRESS UPDATE ---------------------------
            it = it + 1;

            if useGUI
                if mod(it, updateEvery) == 0 || it == 1 || it == totalIts
                    frac = it/totalIts;
                    waitbar(frac, hwb, sprintf('m = %d/%d, l = %d/%d (%.1f%%)', ...
                        m, mup, l, na, 100*frac));
                end
            elseif useTTY
                pb.update(it, sprintf('m %d/%d | l %d/%d', m, mup, l, na));
            end
            % ---------------------------------------------------------------

            a_l = opts.as(l);
            t_ob = transpose(t_o(:,l)); % (1,ncol)
            % Query times ti - tobs
            tq = tseg - t_ob; % (k,ncol)
            J = zeros([k,nx,nr,6],'like',q);

            % Interpolate the 6 components of T
            for c = 1:6

                % Evaluate at (tq,cq)
                F_c = F{c};
                J(:,:,:,c) = reshape(F_c(tq,cq),[k,nx,nr]);

            end

            % Calculate the double time derivative
            G = zeros(size(J),'like',J); % (k,nx,nr,6)
            G(2:k-1,:,:,:) = (J(3:k,:,:,:) - 2*J(2:k-1,:,:,:) + J(1:k-2,:,:,:))/(dt^2); % Central for interior points
            G(1,:,:,:)     = (J(3,:,:,:) - 2*J(2,:,:,:) + J(1,:,:,:))/(dt^2);           % Forward for first point
            G(k,:,:,:)     = (J(k,:,:,:) - 2*J(k-1,:,:,:) + J(k-2,:,:,:))/(dt^2);       % Backward for last point

            % Project onto observer direction
            R_x = Rx(:,:,:,l); % (1,nx,nr)
            R_r = Rr(:,:,:,l); % (1,nx,nr)
            R_  = R(:,:,:,l);  % (1,nx,nr)

            G(:,:,:,1)  = (R_x.^2).*G(:,:,:,1);
            G(:,:,:,2)  = 2*(R_x.*R_r*cos(a_l)).*G(:,:,:,2);
            G(:,:,:,3)  = -2*(R_x.*R_r*sin(a_l)).*G(:,:,:,3);
            G(:,:,:,4)  = ((R_r.^2)*((cos(a_l))^2)).*G(:,:,:,4);
            G(:,:,:,5)  = -((R_r.^2)*sin(2*a_l)).*G(:,:,:,5);
            G(:,:,:,6)  = ((R_r.^2)*((sin(a_l))^2)).*G(:,:,:,6);

            G = G./(R_.^3);
            % Account for azimuthal variation of nonzero azimuthal modes
            if m > 0
                G = 2*real(exp(1i*m*a_l)*G); % (k,nx,nr,6)
            end

            S = S + squeeze(sum(permute(weight,[3,1,2]).*G,[2,3,4]));

        end      
    end

    % Finish terminal progress bar (if used)
    if exist('useTTY','var') && useTTY
        pb.finish('Done');
    end

    % Multiply by prefactor
    p = g*S;

end

function safeCloseWaitbar(hwb)
    if ~isempty(hwb) && isvalid(hwb)
        close(hwb);
    end
end

function pb = termProgressBar(total, varargin)
%TERMPROGRESSBAR Professional single-line terminal progress bar.
% Usage:
%   pb = termProgressBar(N,'Prefix','Task','BarWidth',30,'UseUnicode',true);
%   pb.start(); pb.update(i,'status'); pb.finish('Done');

    p = inputParser;
    p.addRequired('total', @(x)isnumeric(x) && isscalar(x) && x>=1);
    p.addParameter('Prefix', '', @(s)ischar(s) || isstring(s));
    p.addParameter('BarWidth', 30, @(x)isnumeric(x) && isscalar(x) && x>=10);
    p.addParameter('UseUnicode', true, @(x)islogical(x) && isscalar(x));
    p.addParameter('MinUpdatePeriod', 0.1, @(x)isnumeric(x) && isscalar(x) && x>=0);
    p.parse(total, varargin{:});
    cfg = p.Results;

    % Choose glyphs
    if cfg.UseUnicode
        fullChar = char(9608);  % '█'
        emptyChar = char(9617); % '░'
    else
        fullChar = '#';
        emptyChar = '-';
    end

    state.started = false;
    state.t0 = [];
    state.lastPrint = -inf;
    state.lastLen = 0;

    pb.start  = @start;
    pb.update = @update;
    pb.finish = @finish;

    function start()
        if state.started, return; end
        state.started = true;
        state.t0 = tic;
        state.lastPrint = -inf;
        state.lastLen = 0;
        printLine(0, 0, '', 0, NaN);
    end

    function update(i, status)
        if ~state.started
            start();
        end
        if nargin < 2, status = ''; end
        i = max(0, min(total, i));
        tNow = toc(state.t0);

        % Throttle printing
        if (tNow - state.lastPrint) < cfg.MinUpdatePeriod && i < total
            return
        end

        frac = i / total;
        eta = NaN;
        if i > 0
            rate = tNow / i;           % sec/it
            eta = rate * (total - i);  % sec
        end

        printLine(frac, i, status, tNow, eta);
        state.lastPrint = tNow;
    end

    function finish(msg)
        if nargin < 1, msg = 'Done'; end
        update(total, msg);
        fprintf('\n');
        state.started = false;
    end

    function printLine(frac, i, status, elapsed, eta)
        filled = round(cfg.BarWidth * frac);
        filled = min(cfg.BarWidth, max(0, filled));
        barStr = [repmat(fullChar, 1, filled), repmat(emptyChar, 1, cfg.BarWidth - filled)];

        pct = 100 * frac;

        elStr = fmtTime(elapsed);
        if isnan(eta)
            etaStr = '--:--';
        else
            etaStr = fmtTime(eta);
        end

        prefix = string(cfg.Prefix);
        if strlength(prefix) > 0
            prefix = prefix + " ";
        end

        status = string(status);
        if strlength(status) > 0
            status = " | " + status;
        end

        line = sprintf('%s[%s] %6.2f%%  %d/%d  elapsed %s  ETA %s%s', ...
            char(prefix), barStr, pct, i, total, elStr, etaStr, char(status));

        % Clear remnants of previous longer line
        pad = max(0, state.lastLen - strlength(line));
        fprintf(1, '\r%s%s', line, repmat(' ', 1, pad));
        state.lastLen = strlength(line);
    end

    function s = fmtTime(tsec)
        if ~isfinite(tsec) || tsec < 0
            s = '--:--';
            return
        end
        tsec = round(tsec);
        hh = floor(tsec/3600);
        mm = floor((tsec - 3600*hh)/60);
        ss = tsec - 3600*hh - 60*mm;
        if hh > 0
            s = sprintf('%d:%02d:%02d', hh, mm, ss);
        else
            s = sprintf('%02d:%02d', mm, ss);
        end
    end
end

function [weight] = trapzWeightsPolarnoblank(r,z)
    %TRAPZWEIGHTSPOLAR Integration weight matrix for cylindical coordinates using trapazoidal rule
    % OTS, 2015
    
    nothetar = length(r);
    weight_thetar = zeros(nothetar,1);
    weight_thetar(1) = pi*( r(1) + (r(2)-r(1))/2)^2;
    for i=2:nothetar-1
        weight_thetar(i) = pi*( r(i) + (r(i+1)-r(i))/2 )^2 - pi*( r(i) - (r(i)-r(i-1))/2 )^2;
    end
    weight_thetar(nothetar) = pi*r(end)^2 - pi*( r(end) - (r(end)-r(end-1))/2 )^2;
    
    % dz
    noz = length(z);
    weight_z = zeros(noz,1);
    weight_z(1) = (z(2)-z(1))/2;
    for i=2:noz-1
        weight_z(i) = (z(i)-z(i-1))/2 + (z(i+1)-z(i))/2;
    end
    weight_z(noz) = (z(noz)-z(noz-1))/2;
    
    weight    = weight_z*weight_thetar'; % (nx,nr)
end