clear all;

% %  initialisation
%--------------------------------------------------------------------------
fprintf('%s \n','initialising variables...')
% intialise variables
ist_file=1;
ied_file=3;
nfiles=ied_file-ist_file+1;
itt=0;

nt=128;
nt1=nt+1;
nr=81;
nx=576;
nmode=6;
rgas=286.8973;
rho0 = 1.226; % Density of air

 % read hdf data dimension, dims = [nnodes, nsmpales, npdes]
fid=H5F.open('/rds/projects/w/wangzt-hifi-cfd/data/JER_chev/ldt/ldt-01.h5');
dset_name='ldt';
dset_id=H5D.open(fid,dset_name);
type_id=H5D.get_type(dset_id);    
dspace_id = H5D.get_space(dset_id);
[~,dims]=H5S.get_simple_extent_dims(dspace_id); 

% declare variable size
nsamples=nfiles*dims(2);
mid=round((nt+1)/2);
if(nmode > nt) then
   fprintf('%s%03d%s \n','the number of output mode is larger than the maxium mode number (', nt, ').'); 
   return;
end
tmode=nan(nx,nr,nt,3);

% load mean flow file
fprintf('%s \n','loading mean flows...');
load('/rds/projects/w/wangzt-hifi-cfd/data/JER_chev/meanflows_tp.mat','flows');

% initialise output file for azimuthal modes

modeM = 0:8;

for im=1:length(modeM)
    matfilename=['~/AWC/SPOD-sound/lighthill-modes-pure/azimodes-chv/',...
        sprintf('%s%02d%s','lhmode',modeM(im),'tp.mat')];
    mode{im}=matfile(matfilename);
end

%% azimuthal decomposition for instantaneous flow field
%--------------------------------------------------------------------------
fprintf('%s \n','begin to perform azimuthal decomposition');
% loop over files
for ifile=ist_file:ied_file

    % set hdf file name
    filename=['/rds/projects/w/wangzt-hifi-cfd/data/JER_chev/ldt/' sprintf('%s%02d%s','ldt-',ifile,'.h5')];
    fprintf('%s%s \n','begin to process ',filename);
    
    % open hdf file
    fid=H5F.open(filename);
    
    % open dataset and get data type
    dset_name='ldt';    
    dset_id=H5D.open(fid,dset_name);
    type_id=H5D.get_type(dset_id);
    
    % get dataset dataspace and its dimension
    dspace_id = H5D.get_space(dset_id);
    [~,dims]=H5S.get_simple_extent_dims(dspace_id);
    
    % create memory space
    dim=fliplr([5 1 dims(1)]); % [pdes ntime nprobes]
    maxdim = dim;
    mspace_id = H5S.create_simple(3, dim, maxdim);
    
    % loop over time samples
   for it=1:dims(2)
        fprintf('%s%3d%s%s \n','processing the ',it,' samples of ',filename);
        % select elements to be copied for the dataspace
        start = fliplr([0 it-1 0]);
        stride= fliplr([1 1 1]);
        count = fliplr([1 1 1]);
        block = fliplr([5 1 dims(1)]);
        H5S.select_hyperslab(dspace_id, 'H5S_SELECT_SET', start, stride, count, block);
    
        % read data
        data = H5D.read(dset_id, type_id, mspace_id, dspace_id, 'H5P_DEFAULT');
        data = data-flows(1:5,1,:);
        
        % restructrue the data
        idx=0;
        for k=1:nt1
            theta=(k-1)*2*pi/nt;
            yt=cos(theta);
            zt=sin(theta);   
            for j=1:nr
                for i=1:nx
                    idx=idx+1; 

                    vxdash = data(2,1,idx);
                    vydash = data(3,1,idx);

                    pp(k,1,i,j) = rho0*vxdash*vxdash; % Txx
                    pp(k,2,i,j) = rho0*vxdash*vydash; % Txy
                    pp(k,3,i,j) = rho0*vydash*vydash; % Tyy
  
                end
            end
        end        
        
        % azimuthal decomposition
        fprintf('%s%3d%s%s \n','saving the ',it,' samples of ',filename);
        itt=itt+1;
        i=58; j=18;
        for j=1:nr
            for i=1:nx                
                tmode(i,j,1:nt,1:3)=fft(pp(1:nt,1:3,i,j),nt,1)/nt;
            end
        end
        
        % save modes to files        
         for im=1:length(modeM)
             mode{im}.q(itt,1:nx,1:nr,1:3)=reshape(tmode(1:nx,1:nr,modeM(im)+1,1:3),[1,nx,nr,3]);
         end
        
    end % loop over samples
end % loop over files
return;
clear data pp flows filename tmode
return;
