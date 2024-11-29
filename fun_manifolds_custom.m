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
vector_stable = vec(:,index_s);
if vector_stable(1) < 0
    vector_stable = - vector_stable;
end
vector_unstable = vec(:,index_u);
if vector_unstable(1) < 0
    vector_unstable = - vector_unstable;
end





end
