%% EXAMPLE 4: Calculate the realisations of DFT of large data and save results on hard drive.

clc, clear variables functions
clear all; clear classes
dt=1.5e-5;

% Add path to function
addpath('~/AWC/SPOD-sound/')

% This current version of the file is for m = -1 and m = -2

dirnames = {['~/AWC/SPOD-sound/lighthill-modes-pure' ...
        '/azi-dft-realisations-cyl-ss/m-1/'],['~/AWC/SPOD-sound/lighthill-modes-pure' ...
        '/azi-dft-realisations-cyl-ss/m-2/']};

for m = 1:2

    opts.savefft    = true;            
    opts.deletefft  = false;            
%    opts.savedir    = ['~/AWC/SPOD-sound/lighthill-modes-pure' ...
%        '/azi-dft-realisations-cyl-ss/m' sprintf('%02d', m) '/'];        
    opts.savedir = dirnames{m};  
    opts.savefreqs  = 1:256;          % Save all positive frequencies
    opts.nt         = 3600;                        
    
    file = matfile(['~/AWC/SPOD-sound/lighthill-modes-pure/azimodes-cyl-ss/'...
        'mode' sprintf('%02d', m) 'tp.mat']);
%    opts.mean = squeeze(mean(file.q,1));
    opts.mean = squeeze(mean(conj(file.q),1));
    disp(['Loaded lighthill mode ' sprintf('%02d', m)])
    
%    input_function = @(i) squeeze(file.q(i,:,:,:));
    input_function = @(i) squeeze(conj(file.q(i,:,:,:)));
    nDFT = 512;
    ovlp = 0.75*nDFT;
    
    [f] = create_dft_blocks_no_double_ifft(input_function,nDFT,ovlp,dt,opts);

end

return;

