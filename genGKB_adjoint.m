function [U, B, V, QV] = genGKB_adjoint(Q, R, U, B, V, QV, options)
%
%     [U, B, V, QV] = genGKB_adjoint(Q, R, U, B, V, QV, options)
%
%  Performs one step of generalized Golub-Kahan bidiagonalization.
%
% Input:
%       Q, R - covariance matrices
%       U, V - accumulation of vectors
%          B - bidiagonal matrix
%         QV - accumulation of Qv vectors
%    options - structure from HyBR (see HyBRset)
%
% Output:
%       U, V - updated "orthogonal" matrix
%          B - updated bidiagonal matrix
%         QV - matrix of vectors Q*v_j (used for sampling)
%
%  Refs:
%   Arioli. "Generalized Golub-Kahan bidiagonalization and stopping
%       criteria. SIMAX, 2013.
%   Chung and Saibaba. "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems", submitted 2016
%
%   J.Chung and A.Saibaba, 2016
%   J. Chung, modified in 2022 for GEOS-Chem problem

% Determine if we need to do reorthogonalization or not.
reorth = strcmp(HyBR_lsmrget(options,'Reorth'), {'on'});

% Determine whether we're solving the unknown mean case
unmean = strcmp(HyBR_lsmrget(options,'Unknownmean'), {'on'});

% Set dimensions
k = size(B,2)+1;
if unmean;  % Find length of s vector
    n2 = sizen(Q); 
    p  = sizep(Q); 
else;
    n2 = size(Q,1);
end;   

if k == 1
    % do an extra forward/adjoint to get things started
%    v = A'*(R\U(:));
    v = afunadj(U(:),R); % SMM: The GEOS-Chem adjoint already divides by R, so no need to inlcude here.
    if unmean; v = [v; zeros(p,1)]; end;

    alpha = normM(v,Q);
    V(:,1) = v / alpha;
    B = alpha;
    return
end
if reorth % Need reorthogonalization
    v = V(:,k-1);
    Qv = Q*v;
    [h,h2] = afun(Qv(1:n2,:),R); % Only use the first n elements of Qv (remaining elements not needed, even for unknown mean case)

    % For unknown mean case, we need to add zeros to the end of h2 to get dimensions to work
    if unmean; h2 = [h2; zeros(p,1)]; end;

    u = h - B(k-1,k-1)*U(:,k-1); % B(k-1,k-1) is alpha from the previous iteration

    % Reorthogonalize U
    for j = 1:k-1
        temp = (U(:,j)'*(R\u));
        u = u - temp*U(:,j);
    end
    beta = normM(u, @(x)R\x);
    u = u / beta;
    U = [U, u];

    %     v = A'*(R\u) - beta*v;
    scalar = (B(k-1,k-1)^2 + beta^2)/ beta; % Scaler = (beta_k-1^2 + beta_k)/beta_k
    if k == 2
        v = h2/beta - scalar*v;
    else
        scalar2 = (B(k-1,k-1) * B(k-1,k-2))/ beta; 
        v = h2/beta - scalar*v - scalar2*V(:,k-2);
    end

    % Reorthogonalize V 
    for j = 1:k-1
        temp = (V(:,j)'*(Q*v));
        v = v - temp*V(:,j);
    end
    alpha = normM(v,Q);
    v = v / alpha;

else % Do not need reorthogonalization, save on storage
    v = V(:,k-1);

    Qv = Q*v;
    [h,h2] = afun(Qv(1:n2,:),R);

    % For unknown mean case, we need to add zeros to the end of h2 to get dimensions to work
    if unmean; h2 = [h2; zeros(p,1)]; end;

    u = h - B(k-1,k-1)*U;

    beta = normM(u, @(x)R\x);
    U = u / beta; % override U

    scalar = (B(k-1,k-1)^2 + beta^2)/ beta; % (alpha^2 + beta^2)/beta
    if k == 2
        v = h2/beta - scalar*v;
    else
        scalar2 = (B(k-1,k-1) * B(k-1,k-2))/ beta;
        v = h2/beta - scalar*v - scalar2*V(:,k-2);
    end
    alpha = normM(v,Q);
    v = v / alpha;
end

V = [V, v];
B = [B, zeros(k-1,1); [zeros(1,k-2), beta, alpha]];

if nargout > 3
    QV = [QV, Qv];
end
end


function nrm = normM(v, M)
if isa(M, 'function_handle')
    Mv = M(v);
else
    Mv = M*v;
end
nrm = sqrt(v'*Mv);
end
