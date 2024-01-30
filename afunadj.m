function [fmodel, amodel] = afunadj(u)
%
% Purpose: afun is a function handle in FGGK to run the forward and 
% adjoint of the atmospheric model. In this script, I've created a 
% function to launch GEOS-Chem and wait for the outputs.
% The first iteration of FGGK requires an adjoint model run but no forward model run.
% This script is designed to execute that run.
% **** For these runs, make sure that only the vector u gets passed into the adjoint within
% the GEOS-Chem observation operator. I.e., ignore the forward model outputs.
%
% Inputs:
%	u: The initial guess for u in the first iteration. This vector becomes the observations
%		that will get run through the adjoint model.
%
% Outputs:
% 	amodel: Outputs of the adjoint model. This object is typically labeled "gdt.01" in
%		 	the GEOS-Chem outputs.


%--------------------------------------------%
% Write the vector "u" to observation files  %
%--------------------------------------------%
    % confirm the dimension of u
    % TODO: do a asset with dimensions from driver script 
    % assert 
    % 

	% convert u to obs nc files 
    % TODO: VERY IMPORTANT!!!!!!!!!
    %       check the order of u, 
    %obsdir_template = '/Users/leyang/Documents/JHU-projects/ATD/gcadj-support/tasks/tropomi-processing/out/US-2022/C/';
    obs_fns = dir(obsdir_template);
    obs_fns = {obs_fns.name};
    obs_fns = obs_fns(3:end);
    
    u = flip(u);

    currentdir = pwd;
	cd(obsdir);
    
    u_ith_start = 1;
    u_ith_end   = 0;
    for i = 1:length(obs_fns)
        fn = strcat(obsdir, obs_fns{1});

        dat_NOBS_MAX                        = ncread(fn, 'NOBS_MAX'                       );  
        dat_londim_one                      = ncread(fn, 'londim_one'                     );  
        dat_level                           = ncread(fn, 'level'                          );  
        dat_index                           = ncread(fn, 'index'                          );  
        dat_NOBS                            = ncread(fn, 'NOBS'                           );   
        dat_lat                             = ncread(fn, 'lat'                            );  
        dat_lon                             = ncread(fn, 'lon'                            );  
        dat_ch4_true_state_no_error         = ncread(fn, 'ch4_true_state_no_error'        );   
        dat_ch4_true_state_obs_error        = ncread(fn, 'ch4_true_state_obs_error'       );   
        %dat_ch4_from_tropomi_bias_corrected = ncread(fn, 'ch4_from_tropomi_bias_corrected');   
        dat_QFLAG                           = ncread(fn, 'QFLAG'                          );   
        dat_surf_pres                       = ncread(fn, 'surf_pres'                      );   
        dat_pres_interval                   = ncread(fn, 'pres_interval'                  );   
        dat_AK                              = ncread(fn, 'AK'                             );  
        dat_profile_priori                  = ncread(fn, 'profile_priori'                 );   
        dat_dry_air_column                  = ncread(fn, 'dry_air_column'                 );   
        dat_GCFRAC                          = ncread(fn, 'GCFRAC'                         );   
        dat_GCII                            = ncread(fn, 'GCII'                           );   
        dat_GCJJ                            = ncread(fn, 'GCJJ'                           );   
        dat_YYYYMMDD                        = ncread(fn, 'YYYYMMDD'                       );   
        dat_HHMMSS                          = ncread(fn, 'HHMMSS'                         );

        ncid = netcdf.create(fn,'NETCDF4');

        dim_londim_one = netcdf.defDim(ncid,'londim_one', length(dat_londim_one));
        dim_level      = netcdf.defDim(ncid,'level',      length(dat_level)     );     
        dim_index      = netcdf.defDim(ncid,'index',      length(dat_index)     );

        var_londim_one                      = netcdf.defVar(ncid,'londim_one'                      , 'double', dim_londim_one        );
        var_level                           = netcdf.defVar(ncid,'level'                           , 'double', dim_level             );
        var_index                           = netcdf.defVar(ncid,'index'                           , 'double', dim_index             );
        var_NOBS_MAX                        = netcdf.defVar(ncid, 'NOBS_MAX'                       , 'double', dim_londim_one        );     
        var_NOBS                            = netcdf.defVar(ncid, 'NOBS'                           , 'double', dim_index             );
        var_lat                             = netcdf.defVar(ncid, 'lat'                            , 'double', dim_index             );
        var_lon                             = netcdf.defVar(ncid, 'lon'                            , 'double', dim_index             );
        var_ch4_true_state_no_error         = netcdf.defVar(ncid, 'ch4_true_state_no_error'        , 'double', dim_index             );
        var_ch4_true_state_obs_error        = netcdf.defVar(ncid, 'ch4_true_state_obs_error'       , 'double', dim_index             );
        var_ch4_from_tropomi_bias_corrected = netcdf.defVar(ncid, 'ch4_from_tropomi_bias_corrected', 'double', dim_index             );
        var_QFLAG                           = netcdf.defVar(ncid, 'QFLAG'                          , 'double', dim_index             );
        var_surf_pres                       = netcdf.defVar(ncid, 'surf_pres'                      , 'double', dim_index             );
        var_pres_interval                   = netcdf.defVar(ncid, 'pres_interval'                  , 'double', dim_index             );
        var_AK                              = netcdf.defVar(ncid, 'AK'                             , 'double', dim_index             );
        var_profile_priori                  = netcdf.defVar(ncid, 'profile_priori'                 , 'double', dim_index             );
        var_dry_air_column                  = netcdf.defVar(ncid, 'dry_air_column'                 , 'double', [dim_level dim_index] );            
        var_GCFRAC                          = netcdf.defVar(ncid, 'GCFRAC'                         , 'double', [dim_level dim_index] );            
        var_GCII                            = netcdf.defVar(ncid, 'GCII'                           , 'double', [dim_level dim_index] );            
        var_GCJJ                            = netcdf.defVar(ncid, 'GCJJ'                           , 'double', dim_index             );
        var_YYYYMMDD                        = netcdf.defVar(ncid, 'YYYYMMDD'                       , 'double', dim_index             );
        var_HHMMSS                          = netcdf.defVar(ncid, 'HHMMSS'                         , 'double', dim_index             );
        
        netcdf.endDef(ncid);

        netcdf.putVar(var_londim_one                     ,  dat_londim_one                      );
        netcdf.putVar(var_level                          ,  dat_level                           );
        netcdf.putVar(var_index                          ,  dat_index                           );
        netcdf.putVar(var_NOBS_MAX                       ,  dat_NOBS_MAX                        );
        netcdf.putVar(var_NOBS                           ,  dat_NOBS                            );
        netcdf.putVar(var_lat                            ,  dat_lat                             );
        netcdf.putVar(var_lon                            ,  dat_lon                             );
        netcdf.putVar(var_ch4_true_state_no_error        ,  dat_ch4_true_state_no_error         );
        netcdf.putVar(var_ch4_true_state_obs_error       ,  dat_ch4_true_state_obs_error        );
        netcdf.putVar(var_QFLAG                          ,  dat_QFLAG                           );
        netcdf.putVar(var_surf_pres                      ,  dat_surf_pres                       );
        netcdf.putVar(var_pres_interval                  ,  dat_pres_interval                   );
        netcdf.putVar(var_AK                             ,  dat_AK                              );
        netcdf.putVar(var_profile_priori                 ,  dat_profile_priori                  );
        netcdf.putVar(var_dry_air_column                 ,  dat_dry_air_column                  );
        netcdf.putVar(var_GCFRAC                         ,  dat_GCFRAC                          );
        netcdf.putVar(var_GCII                           ,  dat_GCII                            );
        netcdf.putVar(var_GCJJ                           ,  dat_GCJJ                            );
        netcdf.putVar(var_YYYYMMDD                       ,  dat_YYYYMMDD                        );
        netcdf.putVar(var_HHMMSS                         ,  dat_HHMMSS                          );        
	    
        u_ith_end = u_ith_end + dat_NOBS_MAX;
        netcdf.putVar(var_ch4_from_tropomi_bias_corrected,  u(u_ith_start:u_ith_end));
        
        netcdf.close(ncid);	

        u_ith_start = u_ith_end+1;
    end
    
    cd(currentdir);   

%------------------------------------------------%
% Launch GEOS-Chem and wait for it to complete   %
%------------------------------------------------%

    % `bash run_simple` takes care of the cleaning of the dir temp files. 
    disp(['GEOS-Chem launched. Current iteration: ',num2str(count-1)]);
	disp(datetime)
	unix(['bash run_simple &> ADJOINT.out']); % LY: in obs_operator: theres no substraction!
	disp('GEOS-Chem model finished running. Continue calculating cost function.')
	disp(datetime)

%------------------------------%
% Read in GEOS-Chem outputs    %
%------------------------------%

	% Read in the gradient file and create object 'amodel'
    % TODO: check dimension of adj maybe using assert
    adj = ncread(strcat(geosdir,'OptData/gctm.gdt.01'), 'dJ_dESF');    
	adj=reshape(adj,[],1);

	% Only keep land cells
	amodel = adj(landall);
	
