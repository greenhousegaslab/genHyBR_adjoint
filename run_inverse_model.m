%----------------------------------------------------------------------------------------
% run_inverse_model.m
%
% This script will launch an inverse model using genHyBR. 
%	The script is set up for inverse modeling using GEOS-Chem.
%
% Note from Scot: For now, this script is a skeleton script that we can customize with the
%	individual inverse modeling setup that we plan to use for GEOS-Chem.
%
% initial version T.Cho, Jun. 10, 2021
% updated by J. Chung Aug 12, 2021
% updated by S. Miller, Dec 29, 2022
% updated for GenHyBR by S. Miller, Sept. 2023
%----------------------------------------------------------------------------------------


%-----------------%
% Other notes:    %
%-----------------%

	% In papers from Anna Michalak's group, they often set "m" as the length of s and "n" as the length of Z.
	% Julianne Chung and Arvind Saibaba use the opposite notation -- "n" is the length of s and "m" the length of Z.
	% In this script, I've used Julianne and Arvind's convention.


%--------------------------------%
% Set the run name and path      %
%--------------------------------%

	% Set the name of the case study
	% Here, I've created a variable called "casestudy" where you can set a unique
	% name for the case study you're running. Below, you can have GEOS-Chem run different
	% options depending on the name of the case study.
	casestudy = 'case1';

	% Path to the inverse modeling inputs (e.g., prior emissions, etc.)
	switch casestudy
	    case 'case1'
	        inpath = '';
	    case 'case2'
	        inpath = '';
	end;


%----------------------------------------------------------------------------------------%
% Check to see if the workspace has already been created. If so, read in that workspace. %
%----------------------------------------------------------------------------------------%

	workspacename = strcat(inpath,'inversion_inputs_',casestudy,'.mat');

	if exist(workspacename) == 2;
	
	load(workspacename);
	
	
%------------------------------------------------------------%
% If not, create all input variables for the inverse model.  %
%------------------------------------------------------------%

	else;


%-------------------------------%
% Set additional run options    %
%-------------------------------%

	% maximum iteration of hybrid method
	iter = 50;

	% Choose regularization option : 'optimal', 'dp'
	% 	The 'optimal' method will produce the best results but only works with 
	% 	synthetic problems. I think the 'dp' method was the preferred approach by Jiahua
	%	for real-world inverse problems.
	RegOptions = 'dp';
	nRO = length(RegOptions);

	% Set the number of time periods in the inverse model
	% (i.e., the number of time periods for which emissions get estimated)
	ntimes = 31;	


%--------------------------------------%
% Set covariance matrix parameters     %
%--------------------------------------%

	% theta(1): sigma_R. This parameter defines R.
	%	Units on theta(1) are standard deviation or ppb.
		
	% The parameters below apply to s1, not to s2.
	% theta(2): sigma_Q. This parameter defines the diagonal elements of Q.
	%	The units on theta(2) are standard deviation (same units as the fluxes).
	%   Note: we usually set this value to 1 and then solve for the value in the inverse
	%   model via the regularization parameter.
	% theta(3): Decorrelation length (in km).
	% theta(4): Decorrelation time (in days).

	switch casestudy
	    case 'case1'
	        theta = [ 2.000 1 555.420 9.854 ];
	    case 'case2'
	        theta = [ 2.000 1 585.680 12.366 ];
	end;

	disp('Covariance matrix parameters');
	disp(theta);
	
	
%----------------------------------------%
% Define the X matrix                    %
%----------------------------------------%

	X = [ set the X matrix here however you'd like ];
	
	n = size(X,1);
	p = size(X,2);


%------------------------------%
% Read in the observations     %
%------------------------------%

	% The order of the observations here must match the order that they're read in using
	% the observation operator in GEOS-Chem.
	
	Z = [ read in the observations here ];
	m = length(Z);
		
	% Create the variable "b"
	% For the unknown mean case:
	% 	Set b = Z, unless you're setting a non-zero prior on beta (see Taewon's paper)
	% For the known mean case:
	% 	To create b, run the prior emissions estimate through GEOS-Chem, and subtract that 
	% 	estimate from Z. These GEOS-Chem runs can either be done offline and read in here, or they can
	% 	be done online here. If the prior estimate is zero, you can skip this step and set b = Z;
		
	b = Z;


%-------------------------------------%
% Create the R covariance matrix      %
%-------------------------------------%

	% Note: if using GEOS-Chem, the observation operator in the adjoint will likely already
	%	handle calculations with R. If that is the case, then value of R set here doesn't matter
	%	and won't actually be used in any of the inverse modeling calculations.

	disp('Create R');
	R = theta(1)^2*ones(n,1);
	R = spdiags(R,0,n,n);


%-----------------------%
% Load 'distmat.mat'    %
%-----------------------%

	% This matrix is a square, symmetric matrix that gives the distance between
	% each pair of grid boxes in the model domain. The units on this distance matrix
	% can be anything you'd like, but the units must match the decorrelation length used
	% to define Q.

	deltamat = load(strcat(inpath,'distmat.mat'));
	m1 = size(deltamat,1);


%----------------------------------------------------------------%
% Create E -- spatial covariances for the Q covariance matrix    %
%----------------------------------------------------------------%

	% Create E
	% Spherical covariance model
	E  = 1 - 1.5 .* (deltamat ./theta(3))  + 0.5 .* (deltamat.^3 ./ theta(3).^3);
	E(deltamat > theta(3)) = 0;


%---------------------------------------------------------------%
% Create D -- temporal covariances for the Q covariance matrix  %
%---------------------------------------------------------------%

	% Create D
	% Create time distance matrix
	days = 1:ntimes;
	days = days’;
	days = days * ones(1,length(days));
	days = abs(days - days’);

	% Spherical covariance model
	D = 1 - 1.5 .* (days   ./theta(4))  + 0.5 .* (days.^3   ./ theta(4).^3);
	D(days   > theta(4)) = 0;
	
	sigmaQ = theta(2);
	D = (sigmaQ*sigmaQ') .* D;
	

%------------------------%
% Create augmented Q     %
%------------------------%

	% "kronMat" is an object class that Julianne created for matrices that are Kronecker
	%	products. This object class will override the standard Matlab matrix math
	%	operations and will instead implement more efficient operations that exploit the
	%	structure of the Kronecker product.
	% 'QtilMat' is also an object class that Julianne created for the Q tilda structure 
	%	used in the unknown mean case.
	
	disp('Q-kron');

	Q = kronMat(D,E);

	% Note: you can set the value of gamma her to whatever value you'd like.
	% GenHybr includes a prior pdf on beta. We almost always assume that the prior mean of 
	% beta is zero. Here, gamma defines the diagonal covariance matrix for beta.
	% In this case, gamma defines the prior covariance on beta.
	gamma = 10;
	Qbeta = (1/gamma^2) * eye(p);

	% Additional, required math for the unknown mean case:
	Q = QtilMat(Qbeta,Q,X);


%---------------------------------------------------------------%
% Save workspace in case run crashes or needs to be restarted   %
%---------------------------------------------------------------%

	save(workspacename);
	
	end; % End workspacename if statement


%--------------------------------------%
%  Run genHyBR to estimate emissions   %
%--------------------------------------%

	disp('Run_genHyBR');
	thr = 1;

	% Set genHyBR options
	input = HyBR_lsmrset('InSolv', 'tikhonov', 'RegPar',RegOptions, 'Iter', iter,'thr', thr, 'Unknownmean','on');
	
	% Launch genHyBR
	[x_out, output] = genHyBR_adjoint(b, Q, R, input);
	

%---------------------------------------------------%
% Extract s (and beta) from the output of genHyBR   %
%---------------------------------------------------%

	% The 'output' variable contains a lot of summary information on the run.

	% For the known mean case, s = x_out

	s    = x_out(1:n,:);
	beta = x_out((n+1):length(x_out),:);

	
%---------------------------------------------%
%  Save the inverse modeling outputs to file  %
%---------------------------------------------%

	save(strcat(inpath,'emissions_estimate_',casestudy,'.mat'),'s','beta');


%----------------------------------------------------------------------------------------
% END OF SCRIPT
%----------------------------------------------------------------------------------------
