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

	% Note from Scot: The code below is what I use for CO2. 

	disp('Write fluxes to netcdf file for GEOS-Chem model');

	% Add ocean fluxes back into the flux vector
	landall = repmat(landmap,ntimes,1);
	landall = landall==1;
	temp = ocean;
	temp(landall) = temp(landall) + shat_backtransform;

	% Convert the fluxes to netcdf. The GEOS-Chem adjoint will read in the fluxes in netcdf format.
	temp = reshape(temp,72,46,ntimes); % 72 lon and 46 lat

	% One additional, important note: the function "netcdf.create" only works if you are currently in the folder where you intend to write
	% the netcdf file. Hence, I've included a command to "cd" into that folder.
	currentdir = pwd;
	cd(fluxpath);

	% Define latitudes and longitudes for the netcdf file
        lons = -180:5:175;
        lats = [-89 -86:4:86 89];

	% Loop over each time period in the inverse model
    for j=1:ntimes,

	% Convert day of year to month and day of month 
	[yy month day HH MM] = datevec(datenum(year,1,j)); end;
 
	if month<10; month1=strcat('0',num2str(month)); else; month1=num2str(month); end;
	if day<10; day1=strcat('0',num2str(day)); else; day1=num2str(day); end;

	% Open the netcdf file
	ncid = netcdf.create(strcat('CO2.daily.geos.4x5.',num2str(yy),".",month1,".",day1,'.nc'),'NETCDF4');

	% Define the dimensions of the netcdf file
	dimid1 = netcdf.defDim(ncid,'lon', length(lons));
	dimid2= netcdf.defDim(ncid,'lat', length(lats));

	% Define new variables for the netcdf file
	my_varID = netcdf.defVar(ncid,'CO2_flux','double',[dimid1 dimid2]);
	my_varID1= netcdf.defVar(ncid,'lon','double',[dimid1]);
	my_varID2= netcdf.defVar(ncid,'lat','double',[dimid2]);
	netcdf.endDef(ncid);

	% Write data to the netcdf variable
	netcdf.putVar(ncid,my_varID,temp(:,:,j));
	netcdf.putVar(ncid,my_varID1,lons);
	netcdf.putVar(ncid,my_varID2,lats);
	
	% Close the netcdf file
	netcdf.close(ncid);	

	end; % End of j loop
	clear temp;

	% Move from the flux directory back to the previous directory
	cd(currentdir);


%------------------------------------------------%
% Launch GEOS-Chem and wait for it to complete   %
%------------------------------------------------%

	% First, it's possible you'll need to delete some files from the run directory before
	% you start GEOS-Chem (e.g., I'm not sure if the observation operator writes new
	% model-data comparisons each time it runs, or if it appends to any existin model-data
	% comparison file).
	
	disp(['GEOS-Chem launched. Current iteration: ',num2str(count-1)]);
	disp(fix(clock))
	unix(['./run &> ADJOINT.out']);
	disp('GEOS-Chem model finished running. Continue calculating cost function.')
	disp(fix(clock))


%------------------------------%
% Read in GEOS-Chem outputs    %
%------------------------------%

	% Read in the forward model output and assign to object 'fmodel'
	% [ FILL IN HERE!!!!! ]

	% Read in the gradient file and create object 'amodel'
	% Here, I've assumed the gradient file is binary punch, and I have the function
	% readBPCHSingle() if needed.
	[ adj ] = readBPCHSingle(strcat(geosdir,'OptData/gctm.gdt.01'),'C_IJ_GDE','CO2bal', ...
	strcat(geosdir,'tracerinfo.dat'),strcat(geosdir,'diaginfo.dat'),true,true,false);

	adj=reshape(adj,[],1);

	% Only keep land cells
	amodel = adj(landall);
	
