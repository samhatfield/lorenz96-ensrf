
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %
% This routine applies the Ensemble Kalman Filter to        %
% the Lorenz96 model -                                      %
% Aneesh C. S. (Code modified from Ibrahim Hoteit, 2007)    %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%% Reset the Gaussian generator to its 1000-th state.
randn('state',1000);


% =============================
% Set Lorenz96 Model Parameters
% =============================

Nx = 40;            % Nbre of variables in the Lorenz model.
F  = 8;             % Lorenz model parameter.
dt = 0.05;          % Time step ( = 6 hours).
t1 = dt;            % Initial integration time.
t2 = 20;            % Final integration time  (= 30 days = 4*0.05*30).
% t2 = 73;            % Final integration time  (= 1 year = 4*0.05*365).

% Set total number of model steps
nmsteps = length(t1:dt:t2);

% Initialize solution arrays
ref_states     = zeros(nmsteps,Nx);
free_states    = zeros(nmsteps,Nx);
filter_states = zeros(nmsteps,Nx);
forecast_states = zeros(nmsteps,Nx);



% =========================
% Simulate Reference States
% =========================

%%% Set initial state.
% The initial state is obtaind by integrating the model N0 steps starting from 
N0 = 5000;
Xn = F*ones(1,Nx);  
Xn(20) = (1+0.001)*F;

% Integrate the model.
for nt = 1:N0
  Xn = lorenz96(Nx,Xn,dt,F);
end
X0 = Xn;

%%% Set model error --> N(0,varx).
varx   = 1*ones(1,Nx)   ;
sigmax = sqrt(varx)     ;
alfa   = 0.0            ;       % Set to 0 for no model error.

%%% Sinmulate reference states.
% Integrate the model between t1 and t2.
Xn = X0;
for nt = 1:nmsteps
  % Save trajectory.
  ref_states(nt,:) = Xn;
  % Integrate and perturb states.
  Xn = lorenz96(Nx,Xn,dt,F) + alfa*sqrt(dt)*sigmax.*randn(1,Nx);
end
clear Xn;



% ====================
% Extract Observations
% ====================

%%% Which Observation system?
% 'all' stands for all variables are observed.
% 'half' stands for every other variable is observed.
%cas = input('Which observation system? all, half,','s');
cas = 'all';
% cas = 'half';
%cas = 'quarter';

%%% Extract observations.
obs = obs_operator(ref_states',cas);

%%% Set observation system parameters.
% Set dimension of the observation vector.
No = size(obs,1);                     
% Set frequency of the observations.
dto = dt;       % Observation every model time step (= 6hours).
% dto = 2*dt;    % Observation every 12 hours.
% dto = 4*dt;    % Observation every day.
%dto = 8*dt;    % Observation every 2 days.
%dto = 16*dt;    % Observation every 4 days.

% Set observation error covariance matrix: N(0,varo).
varo   = 0.01*ones(1,No);       
sigmao = sqrt(varo);
Ro     = diag(varo);

%%% Perturb observations.
for it = 1:size(obs,2)
  obs_states(it,:) = obs(:,it)';% + sigmao.*randn(1,No);
end



% ==========
% Run Filter
% ==========

% Set filter parameters
% =====================

%%% Set the size of the ensemble.
Ne = 500;

%%% Filter step (in model steps).
fstep = dto/dt;

%%% Set total number of filter steps.
nfsteps = length(t1:dto:t2)-1;

%%% Set inflation factor
rho = 1.2;


% Initilization Step
% ==================

%%% Set the filter initial state --> N(Xa0,var0).
%Xa0    = zeros(1,Nx);     % Start from the origin of the state space
%Xa0    = X0 - 0.5*X0;     % Start from model initial state + perturbation.
Xa0    = mean(ref_states,1); % Start from the mean of the reference states.
var0   = 3*ones(1,Nx);
sigma0 = sqrt(var0);

%%% Generate the initial ensemble members.
for i = 1:Ne
  Ensa(:,i) = Xa0' + sigma0'.*randn(Nx,1);
end

%%% Save initial filter estimate.
forecast_states(1,:)  = mean(Ensa,2)';
filter_states(1,:)  = mean(Ensa,2)';


% -------------
% Free-run
% -------------

%%% Integrate the model between t1 and t2 starting from the filter initial state.
Xn = Xa0;
for nt = 1:nmsteps
  % Save trajectory.
  free_states(nt,:) = Xn;
  % Integrate and perturb states.
  Xn = lorenz96(Nx,Xn,dt,F) + alfa*sqrt(dt)*sigmax.*randn(1,Nx);
end
clear Xn;


% -------------
% Filter's Loop
% -------------

for nf = 1:nfsteps
  
  disp(['Step-' num2str(nf)]);
  
  % Forecast Step
  % =============
  
  %%% Run model for each analysis member to generate forecast ensemble.
  Ensf = Ensa;
  for l = 1:fstep
    for i = 1:Ne
      Xn(i,:) = lorenz96(Nx,Ensf(:,i)',dt,F) + alfa*sqrt(dt)*sigmax.*randn(1,Nx);
      Ensf(:,i) = Xn(i,:)';
    end
    % Save filter estimate.
%     filter_states((nf-1)*fstep+l+1,:) = mean(Ensf,2)';
    forecast_states(nf*fstep+1,:) = mean(Ensf,2)';
  end
  
  %%% Compute forecast state.
  Ensfm = mean(Ensf,2);
  
  disp(['  End Forecast  ']);

  
  % Analysis Step
  % =============
  
  %%% Get observation vector at this time and perturb it.
  for i = 1:Ne
    Obstab(:,i) = obs_states((nf-1)*fstep+1,:)' + sigmao'.*randn(No,1);
  end
  
  %%% Compute innovation = Y - H*Ensf.
  Innovtab = Obstab - obs_operator(Ensf,cas);
  
  %%% Compute Ensfp = Ensf - mean(Ensf).
  Ensfp = Ensf - Ensfm(:,ones(1,Ne));
  Varf_states(nf,:) = diag((1/(Ne-1))*(Ensfp*Ensfp'));
            
  %%% Compute H*Ensfp.
  HEnsfp = obs_operator(Ensfp,cas);
  
  %%% Compute covf*H(t).
  covfHt = (rho/(Ne-1))*Ensfp*HEnsfp';
  
  %%% Compute H*covf*H(t) + R.
  HcovfHt = (rho/(Ne-1))*HEnsfp*HEnsfp' + Ro;
  
  %%% Compute Gain.
  Gain = covfHt*inv(HcovfHt);
  
  gainmax(nf*fstep+1) = max(abs(Gain(:)));
  
  %%% Compute analysis ensemble.
  Ensa = Ensf + Gain*Innovtab;
  
  %%% Save filter estimate.
  filter_states(nf*fstep+1,:) = mean(Ensa,2)';
  
  disp(['  End Analysis  ']);
  
end



% =======
% Results
% =======

%%% Set times of filter steps;
intf = [t1:5*dto:5*t2];
  

%%% Compute and Plot RMS.
% Compute RMS
for nf = 1:nfsteps+1
  rmse(nf) = (1/sqrt(Nx))*norm(ref_states((nf-1)*fstep+1,:)-free_states((nf-1)*fstep+1,:),2);
  rmsf(nf) = (1/sqrt(Nx))*norm(ref_states((nf-1)*fstep+1,:)-forecast_states((nf-1)*fstep+1,:),2);
  rmsa(nf) = (1/sqrt(Nx))*norm(ref_states((nf-1)*fstep+1,:)-filter_states((nf-1)*fstep+1,:),2);
end


% Plot RMS
figure;
hold on; 
plot(intf,rmse,'k');
plot(intf,rmsf,'r');
plot(intf,rmsa,'g');
legend('RMS - Free-Run','RMS - Forecast','RMS - Analysis');
title(['ENKF: ' num2str(Ne) ' Members'],'fontsize',14);


%%% Display RMS for the whole state vector and for each variable.
disp('')
disp('================');
disp(['Mean Free-Run RMS = ' num2str(mean(rmse))]);
disp(['Mean Forecast RMS = ' num2str(mean(rmsf))]);
disp(['Mean Analysis RMS = ' num2str(mean(rmsa))]);
disp('================');
for i = 1:Nx
  disp(['RMS VAR-' num2str(i) ': Free-Run = ' num2str((1/sqrt(nfsteps+1))*norm(ref_states(1:fstep:end,i)-free_states(1:fstep:end,i))) ' || Forecast = ' num2str((1/sqrt(nfsteps+1))*norm(ref_states(1:fstep:end,i)-forecast_states(1:fstep:end,i))) ' || Analysis = ' num2str((1/sqrt(nfsteps+1))*norm(ref_states(1:fstep:end,i)-filter_states(1:fstep:end,i)))]);
end


%%% Plot some variables.
% VAR-1
figure
subplot(2,1,1); hold on; 
plot(intf,ref_states(1:fstep:end,1),'k'); 
plot(intf,forecast_states(1:fstep:end,1),'r');
plot(intf,filter_states(1:fstep:end,1),'g');
title('VAR-1');
legend('Reference','Forecast','Analysis');
subplot(2,1,2); hold on; 
plot(intf,abs(ref_states(1:fstep:end,1)-free_states(1:fstep:end,1)),'k');
plot(intf,abs(ref_states(1:fstep:end,1)-forecast_states(1:fstep:end,1)),'r');
plot(intf,abs(ref_states(1:fstep:end,1)-filter_states(1:fstep:end,1)),'g');
title('Diff Var-1');
%suptitle(['ENKF: ' num2str(Ne) ' Members'],12);
legend('Free-Run','Forecast','Analysis');
% VAR-10
figure
subplot(2,1,1); hold on; 
plot(intf,ref_states(1:fstep:end,10),'k'); 
plot(intf,forecast_states(1:fstep:end,10),'r');
plot(intf,filter_states(1:fstep:end,10),'g');
title('VAR-10');
legend('Reference','Forecast','Analysis');
subplot(2,1,2); hold on;
plot(intf,abs(ref_states(1:fstep:end,10)-free_states(1:fstep:end,10)),'k');
plot(intf,abs(ref_states(1:fstep:end,10)-forecast_states(1:fstep:end,10)),'r');
plot(intf,abs(ref_states(1:fstep:end,10)-filter_states(1:fstep:end,10)),'g');
title('Diff Var-10');
%suptitle(['ENKF: ' num2str(Ne) ' Members'],12);
legend('Free-Run','Forecast','Analysis');
% VAR-22
figure
subplot(2,1,1); hold on; 
plot(intf,ref_states(1:fstep:end,22),'k'); 
plot(intf,forecast_states(1:fstep:end,22),'r');
plot(intf,filter_states(1:fstep:end,22),'g');
title('VAR-22');
legend('Reference','Forecast','Analysis');
subplot(2,1,2); hold on; 
plot(intf,abs(ref_states(1:fstep:end,22)-free_states(1:fstep:end,22)),'k');
plot(intf,abs(ref_states(1:fstep:end,22)-forecast_states(1:fstep:end,22)),'r');
plot(intf,abs(ref_states(1:fstep:end,22)-filter_states(1:fstep:end,22)),'g');
title('Diff Var-22'); 
%suptitle(['ENKF: ' num2str(Ne) ' Members'],12);
legend('Free-Run','Forecast','Analysis');
% VAR-37
figure
subplot(2,1,1); hold on; 
plot(intf,ref_states(1:fstep:end,37),'k'); 
plot(intf,forecast_states(1:fstep:end,37),'r');
plot(intf,filter_states(1:fstep:end,37),'g');
title('VAR-37');
legend('Reference','Forecast','Analysis');
subplot(2,1,2); hold on; 
plot(intf,abs(ref_states(1:fstep:end,37)-free_states(1:fstep:end,37)),'k');
plot(intf,abs(ref_states(1:fstep:end,37)-forecast_states(1:fstep:end,37)),'r');
plot(intf,abs(ref_states(1:fstep:end,37)-filter_states(1:fstep:end,37)),'g');
title('Diff Var-37'); 
%suptitle(['ENKF: ' num2str(Ne) ' Members'],12);
legend('Free-Run','Forecast','Analysis');


%%% Save solution.
eval(['save Results/enkf=obs_' cas num2str(dto) '_N_' num2str(Ne) ...
      '_rho_' num2str(rho) '.mat intf Nx fstep nfsteps Ne rmse rmsf' ...
      ' rmsa ref_states free_states forecast_states filter_states;']);  

