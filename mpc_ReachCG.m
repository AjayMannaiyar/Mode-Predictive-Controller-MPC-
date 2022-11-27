clear all; close all;

%% Control of time delay system using inverse optimal constrained control (MPC)

% PRIOR TO RUNNING

% 1) Fill in missing data in simulink_ReachCG.slx

% 2) Add YALMIP-master and all subfolders to Matlab path : (set Matlab
% working directory to folder directly above YALIMIP-master; then right
% click on YALMIP-master folder in Matlab folder view and select 'add to
% path -> Selected folder and sub-folders'

% 3) Fill in missing code lines hilghilghted by '%???' below


%% Simulink_ReachCG.slx  Simulink model parameters

Ts = 5;        % Sampling interval [hours]
Tdelay = 25;   % time delay [hours]
simtime = 130; % simulation time [samples]


baseflow = 0.6;       % baseflow at Gowangardie weir
delta_baseflow = 0.3; % limit to change of baseflow

Tstep = 105;    % time at which water starts to be withdrawn at Gowangardie [hours]               
QGdraw = 0.4;   % demanded withdrawal flow rate at Gowangardie for t >= Tstep [m^3/s]
QGdrawinit = 0; % demanded flow rate at Gowangardie for t < Tstep [m^3/s] 

Ki = 0.05;     % Favourite controller: I-controller, K(z) = Ki/(z-1)


%% State Space Plant model
A =[0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1;
    0 0 0 0 0];
B = [0; 0; 0; 0; 1]; 
C = [1 0 0 0 0];
D = 0;
dplant = ss(A,B,C,D);

% Augment with disturbance model
Ad = [A B;zeros(1,length(A)) 1]; 
Bd = [B;0];
Cd = [C 0]; 
Dd = 0; 
dplantd = ss(Ad,Bd,Cd,Dd);

%% Controller Model
% Favourite controller: I-controller, K(z) = Ki/(z-1)
Ki = 0.05; 
Ak = 1; 
Bk = -Ki;   %Negative feedback is considered
Ck = 1; 
Dk = 0; 

% Closed loop A matrix between Am and Ak
Acld = [Ad Bd*Ck; Bk*Cd Ak]; 

nd = length(Ad);
nkd = length(Ak);


%% Calculating matrix T
[Uud, Tud] =eig(Acld);                         % Computing matrix U (Eqn. 20 of Foo & Weyer (2011))
Usubud = [Uud(:,1) Uud(:,3:7)];               % Constructing left partitioned matrix U by removing eigenvector associated with the fastest real eigenvalue.
Ttud = Usubud(nd+1: nd+nkd, :)/Usubud(1:nd,:);  % Calculating matrix T 
Ttud = real(Ttud);                            % If T is a complex number, consider only the real part

%% Calculating the MPC weight matrices
Qtud = (-Ck*Ttud)'*(-Ck*Ttud); 
Rtud = 1; 
Stud = (-Ck*Ttud)';

Kcud =(-Ck*Ttud);        % Computing state feedback gain K
Kfud = pinv(Ttud)*Bk;    % Computing observer gain L

% Obtaining state space representation of the full order observer (i.e. Eqns 8-9 of Foo & Weyer (2011))
Akud = Ad -Bd*Kcud - Kfud*Cd; 
Bkud = Kfud; 
Ckud = Kcud;  
Dkud = 0; 

% Compute the transfer function from the full order observer state space model
[numud,denud] = ss2tf(Akud,Bkud,Ckud,Dkud,1);
EsGdud =  minreal(tf(numud,denud,Ts));

%% Setting Up MPC Simulation
% Prediction Horizon
Hp = 5; 

% % Modelling disturbance matrix
Bdist = -1;
ndoff = size(Bdist,2);

% Water Offtakes Withdrawal Demand
Doff = 1*[QGdrawinit*ones(1,round(Tstep/Ts -1)) QGdraw*ones(1,simtime)];
Doff1 = 0*Doff;
Cstviol = Doff1;

%% Setting up YALMIP parameters
nxmd = size(Ad,1);
numd = size(Bd,2);
nymd = size(Cd,1);

% Defining Water Withdrawal Matrices
Boffd = -1;
ndoffd = size(Boffd,2);

% Defining YALMIP variables for initial condition.
u0d = sdpvar(numd,1); % ==>u(k-1)
dud = sdpvar(numd,Hp);
x0d = sdpvar(nxmd,1); % ==>x(0)
doff = sdpvar(ndoffd,Hp);

% Defining YALMIP variables for soft constraint
scud = sdpvar(nymd,Hp);
scld = sdpvar(nymd,Hp);

objd = 0; % Initialising cost function
cstd = []; % Initialising constraint matrix

xd = x0d; % Initialising state variable
xhd = x0d; % Initialising observer state variable
ud = u0d; % Initialising input variable

% Initial value for PLANT
X0d = [0.0 * ones(nxmd,1)];
Xd = zeros(nxmd,simtime);
Xd(:,1) = X0d;
Ud = 0.0 * ones(numd,simtime);

% Initial value for OBSERVER
Xh0d = [0.0*ones(nxmd,1)];
Xhd = zeros(nxmd,simtime);
Xhd(:,1) = Xh0d;

% Collecting variable for prediction horizon******************************
for k = 1:Hp
    % Full Order Observer with Doff 
    ud = ud + dud(:,k);            
    xd = Ad*xd + Bd*ud;
    yd = Cd*xd + (-1)*doff(:,k);  
    
    % MPC Cost Function
    objd = objd + xd'*Qtud*xd + ud'*Rtud*ud + xd'*Stud*ud + ud'*Stud'*xd... % MPC cost function
        + 10*scud(:,k) + 10*scld(:,k);                                      % Soft constraint
    
    % Specifying constraint
%     cstd = [cstd, (ud + dud(:,k)) <= XX];                 % Specifying input constraint
   cstd = [cstd, -delta_baseflow-scld(:,k) <= yd <= delta_baseflow+scud(:,k)]; % Specifying output constraint
%    UNCOMMENT THE ABOVE LINE TO IMPOSE OUTPUT FLOW CONSTRAINT
    cstd = [cstd, scud(:,k) >= 0];
    cstd = [cstd, scld(:,k) >= 0];
end

% Solving optimisation using YALMIP
for t = (1):(simtime-Hp)
    % Full Order Observer with augmented disturbance model
    optisolved = solvesdp([cstd, x0d == Xhd(:,t),...
        u0d == Ud(:,t),...
        doff == (Doff(:,t:t+Hp-1))],...
        objd,sdpsettings('solver','quadprog','verbose',0)); % The built-in QP solver in MATLAB is quadprog
    if (optisolved.problem == 1)
        yalmiperror(optisolved.problem)
        return;
    end
    
    % Obtaining solution from the optimisation*****************
    duoptd = double(dud(:,1));
    Ud(:,t+1) = Ud(:,t) + duoptd;
    Ud(:,t+1) = Ud(:,t+1);
    
    % Plant
    Xd(:,t+1) = Ad*Xd(:,t) + Bd*Ud(:,t+1);
    Yd(:,t) = Cd*Xd(:,t) + (-1)*Doff(:,t);
    
    % Impose hard limit on water withdrawal(MPC solution)
    if (Yd(:,t)<-delta_baseflow)
      Cstviol(:,t) = Yd(:,t)+delta_baseflow;
      Doff(:,t) = Doff(:,t)+Cstviol(:,t);
      Yd(:,t) = Cd*Xd(:,t) + (-1)*Doff(:,t);
    end
    
    % Implementing Full Order Observer Method
    Xhd(:,t+1) = (Ad - Kfud*Cd - Bd*Kcud)*Xhd(:,t) + Kfud*(Yd(:,t));
    
    display(t)
end

%% Plotting results
% Generated from simulink file 'simulink_ReachCG.slx'

sim('simulink_ReachCG.slx')

load reachCG_U.mat
load reachCG_Y.mat

tsim = (0:simtime-Hp-1)*Ts;

figure(5);
subplot 211;
plot(rU(1,1:simtime-Hp),rU(2,1:simtime-Hp) + baseflow,'k-','LineWidth',1);
hold on;
plot(tsim,Ud(1,1:end-Hp) + baseflow,'r--','LineWidth',1);
grid on;
title('Control action, u(k), Flow at Casey''s Weir');
legend('K_{fav} (I-Controller)','Full Order Observer','Location','Best');
xlabel('Time (hours)');
ylabel('m^3/s');
axis([1000/60 3e4/60 0.55 1.1]);

subplot 212;
plot(rY(1,2:simtime-Hp),rY(2,2:simtime-Hp) + baseflow,'k-','LineWidth',1);
hold on;
plot(tsim(2:end),(Yd(1,1:end-1)) + baseflow,'r--','LineWidth',1);
grid on;
title('Output, y(k), Flow at Gowangardie Weir');
legend('K_{fav} (I-Controller)','Full Order Observer','Location','Best');
xlabel('Time (hours)');
ylabel('m^3/s');
axis([1000/60 3e4/60 0.1 0.7]);