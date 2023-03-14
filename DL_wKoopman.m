%% 
% This code will generate all the desired matrices as well as controllers 
% for running the dual-loop simulation.
% First, run this code, and then the simulink model to simulate dynamics.
%
% The Koopman part of the current code is inspired from the available code
% by Milan Korda at: https://github.com/MilanKorda/KoopmanMPC/raw/master/KoopmanMPC.zip

clear all; clc
% addpath('./Resources')
% rng(2141444)
%% ---------- Writing System Dynamics ------------------

f = @DCMotor_Shifted;

% Discretizating the dynamics using Runge-Kutta 4
deltaT = 0.01;
k1 = @(t,x,u) (f(t,x,u));
k2 = @(t,x,u) (f(t,x + k1(t,x,u)*deltaT/2,u));
k3 = @(t,x,u) (f(t,x + k2(t,x,u)*deltaT/2,u));
k4 = @(t,x,u) (f(t,x + k1(t,x,u)*deltaT,u));
f_d = @(t,x,u) (x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)));

%% ---------- Defining variables and matrices --------------
n = 2; % number of states
m = 1; % number of control inputs
nd = 1; % noise dimension

noise_intensity = 1; 

% System matrices - Including noise and performance matrices
B1 = 0.01*ones(n,1);
C1 = [0 0; 0 1; 0 0];
D12 = [0;0;1];
nz = size(C1,1);
C2 = [0 1];
ny = size(C2,1);
D21 = 1;

%% ------ Setting up the basis functions -------------
basisFunction = 'rbf';
% RBF centers
Nrbf = 100;
cent = rand(n,Nrbf)*2 - 1;
rbf_type = 'thinplate'; 
% Lifting mapping - RBFs + the state itself
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n;

%% ----- Data genertion for System Identification ---------
Nsim = 100; % Total no. of simulations for generating data set
Ntraj = 200; % No. of trajectories for each simulation

% Control Input - Random forcing
Ubig = 4*(2*rand([Nsim Ntraj]) - 1); % Input varies between - 4 and 4

% Noise matrix
Wnbig = noise_intensity*(2*rand([Nsim Ntraj]) - 1);

% Random initial conditions
Xcurrent = (rand(n,Ntraj)*2 - 1);
Xcurrent_n = Xcurrent;

X = []; Y = []; U = []; Wn = []; Y_out = []; Z = [];
for i = 1:Nsim
    Xnext = f_d(0,Xcurrent,Ubig(i,:)) + B1*Wnbig(i,:);
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Wn = [Wn Wnbig(i,:)];
    Z = [Z C1*Xcurrent + D12*Ubig(i,:)];
    Y_out = [Y_out C2*Xcurrent + D21*Wnbig(i,:)]; 
    Xcurrent = Xnext;
end

%% ------- Lifting the dynamics ---------

% Lifting the dynamics
Xlift = liftFun(X);
Ylift = liftFun(Y);

%% --------- System Matrix calculation ----------

% Matrix structure 
% [A B1 B2; C1 D11 D12; C2 D21 D22] where D11 = D22 = 0

W = [Ylift; Z; Y_out];
V = [Xlift; Wn; U];
VVt = V*V';
WVt = W*V';
M = WVt * pinv(VVt); 

A_lift = M(1:Nlift,1:Nlift);
B1_lift = M(1:Nlift,Nlift+1:Nlift+nd);
B2_lift = M(1:Nlift,Nlift+nd+1:end);

C1_lift = M(Nlift+1:Nlift+nz,1:Nlift);
D12_lift = M(Nlift+1:Nlift+nz,Nlift+nd+1:end);

C2_lift = M(Nlift+nz+1:end,1:Nlift);
D21_lift = M(Nlift+nz+1:end,Nlift+1:Nlift+nd);


% Checking the stability of the identified dynamics
if abs(max(eig(A_lift))) > 1
    disp('Identified System is Unstable')
else
    disp('Identified System is Stable')
end

%% *********************** Predictor comparison ***************************

Tmax = 5;
Nsim = Tmax/deltaT;
u_dt = @(i) (0);

% Initial condition
x0 = [5;30];
x_true = x0;

% Lifted initial condition
xlift = liftFun(x0);

% Simulate
for i = 0:Nsim-1
    % Simulating Koopman dynamics
    xlift = [xlift, A_lift*xlift(:,end) + B2_lift*u_dt(i)]; % Lifted dynamics
  
    % Simulating True dynamics
    x_true = [x_true, f_d(0,x_true(:,end),u_dt(i)) ];
    
end

figure
plot([0:Nsim]*deltaT,x_true(1,:),'-b','linewidth', 2); hold on
plot([0:Nsim]*deltaT,xlift(1,:), '--r','linewidth',2)
xlabel('Time')
ylabel('Output')
legend('True','Koopman')
set(gca,'FontSize',25)

figure
plot([0:Nsim]*deltaT,x_true(2,:),'-b','linewidth', 2); hold on
plot([0:Nsim]*deltaT,xlift(2,:), '--r','linewidth',2)
xlabel('Time')
ylabel('Output')
legend('True','Koopman')
set(gca,'FontSize',25)

%% ------- Designing Feedback control ----------------

% Checking the Controllability
PErank = rank(ctrb(A_lift,B2_lift));
T = eig(A_lift) - A_lift;

if PErank == size(A_lift,1)
    disp('System is Controllable')
elseif PErank < size(A_lift,1)
    disp('System is not Controllable...')
    if rank([T B2_lift]) == size(A_lift,1)
        disp('... BUT Stabilizable')
    end
else
    error('System is not Stabilization. Cannot design the control')
end

%------ Designing feedback gain using LQR -----------
A_new = A_lift;
B1_new = B1_lift;
B2_new = B2_lift;
C1_new = C1_lift;
C2_new = C2_lift;
D12_new = D12_lift;
D21_new = D21_lift;

Q = (C2_lift'*C2_lift);
R = 1;
K = dlqr(A_new,B2_new,Q,R);

%--------- Solving Kalman Filter equation to estimate L ---------
Qn = eye(Nlift);
Rn = eye(1);

Kf = (dlqr(A_new',C2_new',Qn,Rn))';   % alternatively, possible to design using "LQR" code
sysKF = ss(A_new-Kf*C2_new,[B2_new Kf],eye(Nlift),0*[B2_new Kf]);  % Kalman filter estimator


% Estimation gains based on Koopman matrices
K_lqg = -K;
L_lqg = -Kf;

nx = Nlift;

% Augmented matrices given in eqn. 17 of the paper
A_aug = [A_new + B2_new*K_lqg -B2_new*K_lqg; zeros(nx,nx) A_new+L_lqg*C2_new];
B1_aug = [B1_new;B1_new + L_lqg*D21_new];
B2_aug = [B2_new;zeros(nx,size(B2_new,2))];
C1_aug = [C1_new + D12_new*K_lqg -D12_new*K_lqg];
C2_aug = [zeros(size(C2_new,1),nx) -C2_new];
D12_aug = D12_new;
D21_aug = -D21_new;
r11_aug = size(C1_new,1);
c11_aug = size(B1_new,2);
D11_aug = zeros(r11_aug,c11_aug);
r22_aug = size(C2_new,1);
c22_aug = size(B2_new,2);
D22_aug = zeros(r22_aug,c22_aug);
D_aug = [D11_aug D12_aug; D21_aug D22_aug];

% Creating state-space for Hinf control design
P = ss(A_aug,[B1_aug B2_aug],[C1_aug;C2_aug],D_aug,Nsim);

%------------ Designing Hinf control ---------------
nmeas = 1;
ncont = 1;
[K_hinf_aug,CL_hinf_aug,GAM_hinf_aug,INFO_hinf_aug] = hinfsyn(P,nmeas,ncont); % Hinf controller


% Defining initial conditions to be used for simulating dynamics in Simulink
IC_original = [x0' zeros(1,Nrbf)];
IC_augmented = [x0' zeros(1,Nlift+Nrbf)];

%% All the matrices and control obtained will now be used for running the Simulink model 