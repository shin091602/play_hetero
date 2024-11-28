function [DF,Ud] = DF_qpoms(Z,Ud,p,FR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the error Jacobian in calculating QPO by using GMOS method
%% By:Soi Yamaguchi
%%% input
% Z :current state variables [X;T;rho;w]
% Ud :torus solution dictionary
% p :parameter dictionary
% FR :fourier matrices

%%% output
% F :error vector
% Ud :updated torus solution package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY OPEN
d = p("d");
M = p("M");
N = p("N");

% current state variables
T = Z(d*N*M+1);
rho = Z(d*N*M+2);

% current frequencies
w0 = Z(end-1);
w1 = Z(end);

% fourier matrices
Fr = FR{1,2};
IFr = FR{2,2};
DFr = FR{3,2};
Q = FR{4,2};
Q = rotation_matrix(Q,rho,p);
R = IFr*Q*Fr;

% initialize Jacobian
DF = zeros(d*N*M+p("pcon")+2,d*N*M+4);

%% INVARIANCE CONDITION PARTIALS
% wrt U0
DF(1:d*N,1:d*N) = -eye(d*N,d*N);
DF(1:d*N,end-(d*N)-3:end-4) = R*Ud{6,2}(:,:,end);

% wrt T
DF(1:d*N,d*N*M+1) = R*Ud{7,2};

% wrt rho
DF(1:d*N,d*N*M+2) = -DFr*R*Ud{5,2}(:,end);

%% MULTIPLE SHOOTING PARTIALS
for i=2:M
    k=d*N;
    % row idx
    ridx = (i*k-(k-1)):i*k;
    % column idx
    cidx = (i*k-(k-1))-k:i*k-k;
    % partials
    DF(ridx,cidx) = Ud{6,2}(:,:,i-1);
    DF(ridx,ridx) = -eye(k);
end

%% PHASE CONDITION PARTIALS
DP = [(1/N).*Ud{1,2}';(1/N).*Ud{2,2}'];
DF(d*N*M+1:d*N*M+2,1:d*N) = DP(1:2,:);

%% CONSISTENCY PARTIALS
% c0
DF(end-1,d*N*M+1) = w0;
DF(end-1,d*N*M+3) = T;

% c1
DF(end,d*N*M+1) = w1;
DF(end,d*N*M+2) = -1;
DF(end,d*N*M+4) = T;

%% JACOBI INTEGRAL PARTIALS
Dsi = [zeros(1,d*M*N),1,0,0,0];

DF = [DF;Dsi];

end