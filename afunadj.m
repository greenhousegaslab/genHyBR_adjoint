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

	% [ FILL IN CODE HERE ]

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

	% Read in the gradient file and create object 'amodel'
	% Here, I've assumed the gradient file is binary punch, and I have the function
	% readBPCHSingle() if needed.
	[ adj ] = readBPCHSingle(strcat(geosdir,'OptData/gctm.gdt.01'),'C_IJ_GDE','CO2bal', ...
	strcat(geosdir,'tracerinfo.dat'),strcat(geosdir,'diaginfo.dat'),true,true,false);

	adj=reshape(adj,[],1);

	% Only keep land cells
	amodel = adj(landall);
	
