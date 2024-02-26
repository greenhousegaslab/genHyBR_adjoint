function [fmodel, amodel] = afun(emiss)
%
% Purpose: afun is a function handle in FGGK to run the forward and 
% adjoint of the atmospheric model. In this script, I've created a 
% function to launch GEOS-Chem and wait for the outputs.
%
% Inputs:
%	emiss: a vector of numbers that will be fed into GEOS-Chem in place of the emissions.
%	  	   This vector can be emissions, or it can be anything that will be fed into 
%		   GEOS-Chem in place of emissions.
%
% Outputs:
%	fmodel: Outputs of the forward model that have been interpolated to the observations.
% 	amodel: Outputs of the adjoint model. This object is typically labeled "gdt.01" in
%		 	the GEOS-Chem outputs. For the GEOS-Chem simulations here, we want to run 
%		    R^-1*emiss through the adjoint to compute A'*R^-1*emiss. I.e., don't
%		    add or subtract the value of the observations in the GESO-Chem observation
%		    operator.


%--------------------------------------------%
% Write the emissions to file for GEOS-Chem  %
%--------------------------------------------%

	% code for methane.

	disp('Write fluxes to netcdf file for GEOS-Chem model');

	% One additional, important note: the function "netcdf.create" only works if you are currently in the folder where you intend to write
	% the netcdf file. Hence, I've included a command to "cd" into that folder.
	currentdir = pwd;
	cd(fluxpath);

    %%% Keep this block for now
    % For flexibility, lat and lon data are read-in in the for loop blow, 
    % but if the performance is really bad, consider using the hard-code
    % defined lat lon here
	% Define latitudes and longitudes for the netcdf file
    % lons = -126.875:0.625:-65;
    % lats = 23:0.5:51;
    %%% End - keep this block for now

    % Check the incoming emiss dimension, and error out if need to 
    % TODO: add dimension info in the driver script so to avoid reading of nc file now
    % assert([100, 57, 265] == size(emiss), 'the dimension does not look right!')

	% Loop over each time period in the inverse model
    for j=1:ntimes

	% Convert day of year to month and day of month 
	[yy, month, day, ~, ~, ~] = datevec(datetime(year,1,j));

    year_str = num2str(yy);
	if month<10; month_str = strcat('0',num2str(month)); else; month_str =num2str(month); end
	if day<10;   day_str   = strcat('0',num2str(day));   else; day_str   =num2str(day);   end
	
    % Construct nc file name
    nc_fn = strcat('HEMCO_sa_diagnostics.', year_str, month_str, day_str, '0000.nc');
	
    % Read lat/lon/AREA/hyam/hybm/lev/P0/time/... from template emissions
    dat_EmisCH4_SoilAbsorb = ncread(strcat(fluxpath_template, '/', nc_fn), 'EmisCH4_SoilAbsorb');
    dat_lat = ncread(strcat(fluxpath_template, '/', nc_fn), 'lat');
    dat_lon = ncread(strcat(fluxpath_template, '/', nc_fn), 'lon');

    % Open the new nc file
    ncid = netcdf.create(nc_fn,'NETCDF4');

	% Define the dimensions of the netcdf file
    dim_lon = netcdf.defDim(ncid,'lon', length(dat_lon));
    dim_lat = netcdf.defDim(ncid,'lon', length(dat_lat));

	% Define new variables for the netcdf file
    var_EmisCH4_Oil         = netcdf.defVar(ncid, 'EmisCH4_Oil',         'double',[dim_lon dim_lat]);
    var_EmisCH4_Gas         = netcdf.defVar(ncid, 'EmisCH4_Gas',         'double',[dim_lon dim_lat]);
    var_EmisCH4_Coal        = netcdf.defVar(ncid, 'EmisCH4_Coal',        'double',[dim_lon dim_lat]);
    var_EmisCH4_Livestock   = netcdf.defVar(ncid, 'EmisCH4_Livestock',   'double',[dim_lon dim_lat]);
    var_EmisCH4_Wastewater  = netcdf.defVar(ncid, 'EmisCH4_Wastewater',  'double',[dim_lon dim_lat]);
    var_EmisCH4_Landfills   = netcdf.defVar(ncid, 'EmisCH4_Landfills',   'double',[dim_lon dim_lat]);
    var_EmisCH4_Rice        = netcdf.defVar(ncid, 'EmisCH4_Rice',        'double',[dim_lon dim_lat]);
    var_EmisCH4_OtherAnth   = netcdf.defVar(ncid, 'EmisCH4_OtherAnth',   'double',[dim_lon dim_lat]);
    var_EmisCH4_BiomassBurn = netcdf.defVar(ncid, 'EmisCH4_BiomassBurn', 'double',[dim_lon dim_lat]);
    var_EmisCH4_Wetlands    = netcdf.defVar(ncid, 'EmisCH4_Wetlands',    'double',[dim_lon dim_lat]);
    var_EmisCH4_Termites    = netcdf.defVar(ncid, 'EmisCH4_Termites',    'double',[dim_lon dim_lat]);
    var_EmisCH4_Lakes       = netcdf.defVar(ncid, 'EmisCH4_Lakes',       'double',[dim_lon dim_lat]);
    var_EmisCH4_Seeps       = netcdf.defVar(ncid, 'EmisCH4_Seeps',       'double',[dim_lon dim_lat]);
    var_EmisCH4_SoilAbsorb  = netcdf.defVar(ncid, 'EmisCH4_SoilAbsorb',  'double',[dim_lon dim_lat]);
    var_lon = netcdf.defVar(ncid,'lon','double',dim_lon);
    var_lat = netcdf.defVar(ncid,'lat','double',dim_lat);
	
	netcdf.endDef(ncid);

	% Write data to the netcdf variable
    netcdf.putVar(ncid,var_lon, dat_lon);
    netcdf.putVar(ncid,var_lat, dat_lat);
    netcdf.putVar(ncid, var_EmisCH4_SoilAbsorb, dat_EmisCH4_SoilAbsorb);
    % useless note: picking oil to be the lucky one to store all values 
    netcdf.putVar(ncid, var_EmisCH4_Oil, emiss(:,:,j));
	dat_dummy = zeros(length(dat_lon), length(dat_lat));
    netcdf.putVar(ncid, var_EmisCH4_Gas         , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Coal        , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Livestock   , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Wastewater  , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Landfills   , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Rice        , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_OtherAnth   , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_BiomassBurn , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Wetlands    , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Termites    , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Lakes       , dat_dummy);
    netcdf.putVar(ncid, var_EmisCH4_Seeps       , dat_dummy);
	
	% Close the netcdf file
	netcdf.close(ncid);	

    end % End of j loop
	clear emiss;

	% Move from the flux directory back to the previous directory
	cd(currentdir);


%------------------------------------------------%
% Launch GEOS-Chem and wait for it to complete   %
%------------------------------------------------%
    
    % `bash run_simple` takes care of the cleaning of the dir temp files. 
	
	disp(['GEOS-Chem launched. Current iteration: ',num2str(count-1)]);
	disp(datetime)
	unix(['bash run_simple &> ADJOINT.out']);
	disp('GEOS-Chem model finished running. Continue calculating cost function.')
	disp(datetime)


%------------------------------%
% Read in GEOS-Chem outputs    %
%------------------------------%

	% Read in the forward model output and assign to object 'fmodel'
    % TODO: double check with Scot to confirm fmodel and amodel 
    % TODO: think about the performance of reallines 
    % TODO: check the dimensions of fmodel and amodel in fggk 
    % TODO: VERY IMPORTANT!!!!!!
    %       think and check the order of fmodel in here (aka the order of
    %       diff.01.m), what I have right now could be wrong!!!!!
    fn = strcat(geosdir,'diagadj/diff.CH4.01.m');
    %fn = '/Users/leyang/Documents/JHU-projects/ATD/gcadj-support/tasks/tropomi-gen-syn/diff-files/NA-2022/diff_20180101-20180301/diff.CH4.01.m'
    
    dat_lines = readlines(fn);
    [n_lines, ~] = size(dat_lines);

    fmodel = zeors(n_lines, 1);
    for i = 1:n_lines
        line_parts = strsplit(dat_lines(i));
        tango = line_parts(2);
        tango = str2double(tango);
        fmodel(i, 1) = tango;
    end 

	% Read in the gradient file and create object 'amodel'
    % TODO: check dimension of adj maybe using assert
    adj = ncread(strcat(geosdir,'OptData/gctm.gdt.01'), 'dJ_dESF');    
	adj=reshape(adj,[],1);

	% Only keep land cells
	amodel = adj(landall);
