%% EXAMPLE 4: Calculate the realisations of DFT of large data and save results on hard drive.

clc, clear variables functions
clear all; clear classes
dt=0.0001;

% Add path to function
addpath('~/AWC/SPOD-sound/')

for m = 0:8

    opts.savefft    = true;            
    opts.deletefft  = false;            
    opts.savedir    = ['./azi-dft-realisations-chv/m' sprintf('%02d', m) '/'];        
    opts.savefreqs  = 1:52;          % Save frequencies from St=0 to St=1.98
    opts.nt         = 3000;                        
    %opts.mean       = p_mean;           % provide a long-time mean
    % opts.nsave      = 3;                % save the 5 most energetic modes
    
    file = matfile(['./azimodes-chv/lhmode' sprintf('%02d', m) 'tp.mat']);
    disp(['Loaded lighthill mode ' sprintf('%02d', m)])
    
    input_function = @(i) squeeze(file.q(i,:,:,:));
    
    [f] = create_dft_blocks(input_function,256,128,dt,opts);
    disp(['Saved in ' opts.savedir])

end

return;

