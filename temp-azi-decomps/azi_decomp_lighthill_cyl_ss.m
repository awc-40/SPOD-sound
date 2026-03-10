clear all;

% %  initialisation
%--------------------------------------------------------------------------
fprintf('%s \n','initialising variables...')
% intialise variables
ist_file=4;
ied_file=6;
nfiles=ied_file-ist_file+1;
itt=0;

nt=128;
nt1=nt+1;
nr=100;
nx=512;
nmode=9;
rgas=286.8973;
rho0 = 1.226; % Density of air

 % read hdf data dimension, dims = [nnodes, nsmpales, npdes]
fid=H5F.open('/rds/projects/w/wangzt-hifi-cfd/data/isolated_jets/ldt-04.h5');
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
tmode=nan(nx,nr,nt,6);
%ampmode=nan(nsamples,nx,nr,nmode);
%angmode=nan(nsamples,nx,nr,nmode);
%pmode=ones(nsamples,nx,nr,5,nmode)+i*ones(nsamples,nx,nr,5,nmode);

% load mean flow file
fprintf('%s \n','loading mean flows...');
load('/rds/projects/w/wangzt-hifi-cfd/data/isolated_jets/meanflows_tp.mat','flows');

% initialise output file for azimuthal mode
for im=1:nmode
    matfilename=sprintf('%s%02d%s','mode',im-1,'tp.mat');
    mode{im}=matfile(matfilename,'Writable',true);
end

%% azimuthal decomposition for instantaneous flow field
%--------------------------------------------------------------------------
fprintf('%s \n','begin to perform azimuthal decomposition');
% loop over files
for ifile=ist_file:ied_file

    % set hdf file name
    filename=['/rds/projects/w/wangzt-hifi-cfd/data/isolated_jets/' sprintf('%s%02d%s','ldt-',ifile,'.h5')];
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
        pres(1,1,:) = data(5,1,:) - flows(6,1,:);
        data(5,1,:) = data(5,1,:)./(rgas*data(1,1,:));
        % data = data-flows(1:5,1,:);

        % restructrue the data
        idx=0;
        for k=1:nt1
            alpha=(k-1)*2*pi/nt;
            ca=cos(alpha);
            sa=sin(alpha);   
            for j=1:nr
                for i=1:nx
                    idx = idx+1;  

                    ux = data(2,1,idx); 
                    uy = data(3,1,idx);
                    uz = data(4,1,idx);
                    ur = uy*ca + uz*sa;
                    ua = -uy*sa + uz*ca;

                    Txx = rho0*ux*ux;
                    Txr = rho0*ux*ur;
                    Txa = rho0*ux*ua;
                    Trr = rho0*ur*ur;
                    Tra = rho0*ur*ua;
                    Taa = rho0*ua*ua;

                    pp(k,1,i,j) = Txx;
                    pp(k,2,i,j) = Txr;
                    pp(k,3,i,j) = Txa;
                    pp(k,4,i,j) = Trr;
                    pp(k,5,i,j) = Tra;
                    pp(k,6,i,j) = Taa;

                end
            end
        end        
        
        % azimuthal decomposition
        fprintf('%s%3d%s%s \n','saving the ',it,' samples of ',filename);
        itt=itt+1;
        for j=1:nr
            for i=1:nx                
                tmode(i,j,1:nt,1:6)=fft(pp(1:nt,1:6,i,j),nt,1)/nt;
            end
        end
        
        % save modes to files        
        for im=1:nmode
            mode{im}.q(itt,1:nx,1:nr,1:6)=reshape(tmode(1:nx,1:nr,im,1:6),[1,nx,nr,6]);
        end
        
    end % loop over samples
end % loop over files
clear data pp flows filename tmode
return;
