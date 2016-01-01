
% ===================================================================
% Discretized Lorenz96 model with the fourth-order Runge-Kutta method
% The model is integrated for one time step dt starting from X0
% ===================================================================

function [X1] = lorenz96(Nx,X0,dt,F)


% Inputs
% ------
% Nx = System dimension.
% X0 = initial state.
% dt = time step.
% F  = Lorenz96 model parameter.

% Outputs
% -------
% X1 = state vector after one step integration with the Lorenz-96 model.


% Program
% -------

K1 = florenz96(Nx,X0,F);
K2 = florenz96(Nx,X0+0.5*dt*K1,F);
K3 = florenz96(Nx,X0+0.5*dt*K2,F);
K4 = florenz96(Nx,X0+dt*K3,F);

X1 = X0 + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);

