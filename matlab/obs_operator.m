
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine represents the observation operator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HX] = obs_operator(X,caso)

% Inputs
% ------
% X = array of state vectors.
% caso = observation vetor case.

% Outputs
% -------
% HX = H*X.


% Program
% -------

% All variables are observed.
if strcmp(caso,'all') == 1
  HX = X;
end

% Every other variable is observed.
if strcmp(caso,'half') == 1
  HX = X(1:2:end,:);
end

% Every other 3 variables is observed.
if strcmp(caso,'quarter') == 1
  HX = X(1:3:end,:);
end