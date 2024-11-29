%% calculate initial conditions for stable and unstable manifolds of a periodic orbit which has two pair of eigenvalue which is not one.
function [XS_left, XS_right, XU_left, XU_right, Y] = fun_manifold_cr3bp(mu, x0, t0, N, xpert, options)
% mu      : mass ratio of the primaries
% x0      : initial state of a periodic orbit
% t0      : the period of a periodic orbit
% N       : the number of points on a periodic orbit
% xpert   : weight for a position of the eigenvector
% options : options for ode

% calculate eigenvalues of the monodromy matrix
% Propogate orbit for one period and gather data at N points
tspan = linspace(0, t0, N);
X0 = [x0; reshape(eye(6), 36, 1)];
[~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),tspan,X0,options);

% Monodromy matrix
M = reshape(Y(end, 7:42), 6, 6);

% Eigenvector and eigenvalue analysis
[vec, val] = eig(M);
val = diag(val);
val_t = val;

[~, index_s_1] = min(abs(val_t)); % stable
[~, index_u_1] = max(abs(val_t)); % unstable

val_t(index_s_1) = NaN;
val_t(index_u_1) = NaN;

[~, index_s_2] = min(abs(val_t)); % stable
[~, index_u_2] = max(abs(val_t)); % unstable

if (imag(val(index_s)) ~= 0) || (imag(val(index_u)) ~= 0)
    error('Imag eigenvalues are dominant');
end

% Stable and unstable eigenvectors
vector_stable_1 = vec(:,index_s_1);
if vector_stable_1(1) < 0
    vector_stable_1 = - vector_stable_1;
end
vector_unstable_1 = vec(:,index_u_1);
if vector_unstable_1(1) < 0
    vector_unstable_1 = - vector_unstable_1;
end
vector_stable_2 = vec(:,index_s_2);
if vector_stable_2(1) < 0
    vector_stable_2 = - vector_stable_2;
end
vector_unstable_2 = vec(:,index_u_2);
if vector_unstable_2(1) < 0
    vector_unstable_2 = - vector_unstable_2;
end
% calculate manifolds
XS_left = zeros(6, N);
XS_right = zeros(6, N);
XU_left = zeros(6, N);
XU_right = zeros(6, N);

% Apply perturbations to each of N points of a periodic orbit
for iteration = 1:N
    % Grab state transition matrix at the fixed point
    Phi = reshape(Y(iteration, 7:42), 6, 6);

    % Grab state at the fixed point
    x_star = Y(iteration, 1:6)';

    % Map stable and unstable vectors forward
    S_1 = Phi*vector_stable_1;
    S_1 = S_1/norm(S_1);
    U_1 = Phi*vector_unstable_1;
    U_1 = U_1/norm(U_1);
    S_2 = Phi*vector_stable_2;
    S_2 = S_2/norm(S_2);
    U_2 = Phi*vector_unstable_2;
    U_2 = U_2/norm(U_2);

    % Create perturbation vector
    % vpert = xpert*norm(x_star(4:6))/norm(x_star(1:3));
    % pert = [ones(3,1).*xpert; ones(3,1).*vpert];
    pert = ones(6,1).*xpert;

    % Perturb conditions
    XS_left(:,iteration)  = x_star - S_1.*pert;
    XS_right(:,iteration) = x_star + S.*pert;

    XU_left(:,iteration)  = x_star - U.*pert;
    XU_right(:,iteration) = x_star + U.*pert;
end





end
