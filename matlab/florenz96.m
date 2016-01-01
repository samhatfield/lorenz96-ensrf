
% ===================================================
% Disretized transition operator of the Lorenz model
% ===================================================

function [fX] = florenz96(Nx,X,F)

% Inputs
% ------
% Nx = System dimension.
% X  = model state.
% F  = Lorenz model parameter.

% Outputs
% -------
% fX = f*X.


% Program
% -------

%% Set fX(1).
fX(1) = (X(2) - X(Nx-1))*X(Nx) - X(1) + F;

%% Set fX(2).
fX(2) = (X(3) - X(Nx))*X(1) - X(2) + F;

%% Set fX(3:Nx-1).
for i = 3:Nx-1
  fX(i) = (X(i+1) - X(i-2))*X(i-1) - X(i) + F;
end

%% Set fX(Nx).
fX(Nx) = (X(1) - X(Nx-2))*X(Nx-1) - X(Nx) + F;
