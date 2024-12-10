%% Start of script
close all;  %close all figures
clear;      %clear all variables
clc;        %clear the command terminal
format long
%warning off

% line width
set(0, 'DefaultLineLineWidth', 0.8) % default 0.5pt
set(0, 'DefaultAxesLineWidth', 0.8)
set(0, 'DefaultTextLineWidth', 0.8)

% font size
set(0, 'DefaultTextFontSize', 13)
set(0, 'DefaultAxesFontSize', 13)

% font name
set(0, 'DefaultTextFontName', 'Times New Roman')
set(0, 'DefaultAxesFontName', 'Times New Roman')
set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')

% figure color
set(0, 'DefaultFigureWindowStyle', 'docked');
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardCopy', 'off');

close

current_pass = pwd;
addpath("Functions_base/");

%% initial setting
[mu, ~, ~, ~]      = parameter(2); % Earth-Moon
[L1, L2, L3, L4, L5] = librationPoints(mu);

count_max     = 8000;
plot_interval = 100;

Jacobi_min = 3.0;
Jacobi_max = 3.1723;

color      = jet;
Jacobi_lim = linspace(Jacobi_min, Jacobi_max, size(color,1));
Interp_c   = griddedInterpolant(Jacobi_lim, color);

Jacobi       = zeros(count_max/plot_interval, 1);
x0_corrected = zeros(count_max/plot_interval, 6);
t0_corrected = zeros(count_max/plot_interval, 1);

%% オプション設定
options_ODE_1   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14, 'Events', @(t,x) odestop_hetero_1(t,x,mu));
options_ODE_2   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14, 'Events', @(t,x) odestop_hetero_2(t,x,mu));
options_ODE_3   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14, 'Events', @(t,x) odestop_hetero_3(t,x,mu));
options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);

%% ヤコビ定数設定
C_xn = 3.12;
% C_xn = Jacobi_const([L2;zeros(3,1)],mu);
C_error = 1e-13;
disp(C_xn);

%continuation-----------------------------------------------------------------------------------
iteration_max = 100;
threshold     = 1e-10;
delta         = 4e-5;
count         = 0;
count_max     = 8000;
detect_period_orbit_1 = zeros(8,1);

x0_1 = [0.836900082907655, 0, 0, 0, 1.770874936727959e-06, 0];
t0_1 = 1.3458;

%Loop for first orbit------------------------------------------------------------------------------------------------
while 1
    count = count + 1;

    % Differential correction-----------------------------------------------
    for iteration = 1:iteration_max
        [x_n_1, t_n_1, C] = fun_differential_correction_cr3bp(x0_1, t0_1, mu, options_ODE);

        tspan = [0 2*t_n_1];
        [t_corrected_1, x_corrected_1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan, x_n_1, options_ODE);

        x_error_1 = norm(x_corrected_1(end, :) - x_corrected_1(1, :));

        if x_error_1 < threshold
            break;
        end

        if x_error_1 > 1e+3
            disp('calculation diverged');
            return;
        end

        if iteration == iteration_max
            count = count - 1;
            x0_1(1) = x0_1(1) - delta;
            delta = delta/10;
        end

        x0_1 = x_n_1;
        t0_1 = t_n_1;
    end

    if C < C_xn
        x0_1(1) = x0_1(1) - delta;
        delta = delta/5;
        disp(strcat('delta changed : count = ', num2str(count)));

    end

    if(count >= 100)&&(mod(count,100)==0)
        disp(strcat('count = ', num2str(count)));
    end

    if abs(C-C_xn)<C_error
        additional = [C; t_n_1; x_n_1(:,1)];
        detect_period_orbit_1 = additional;
        break;
    end

    x0_1(1) = x0_1(1) + delta;

    if count == count_max - 1
        break
    end
end

disp('L1-sucucessfuly finished');

%% 周期軌道計算---------------------------------
x0_1=detect_period_orbit_1(3:8);
t0_1=detect_period_orbit_1(2);

tspan_1 = [0 2*t0_1];
[t_1, x_1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan_1, x0_1, options_ODE);

%% 面外方向に多様体を伝播
xpert = 1e-4;
N = 100;

[XS_left, XS_right, XU_left, XU_right, Y_2] = fun_manifolds_custom(mu, x0_1, 2*t0_1, N, xpert, options_ODE);
tf = 15;
tspan_s = [tf_2 0];
tspan_u = [0 tf_2];
% プロット

figure();
hold on

for i = 1:N
    [~,ys_left] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_s, XS_left(:,i), options_ODE);
    ys_left = flipud(ys_left);
    f1_p1 = plot3(ys_left(:,1),ys_left(:,2),ys_left(:,3),'g');
end
for i = 1:N
    [~,ys_right] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_s, XS_right(:,i), options_ODE);
    ys_right = flipud(ys_right);
    f1_p2 = plot3(ys_right(:,1),ys_right(:,2),ys_right(:,3),'Color','#77AC30');
end
% for i = 1:N
%     [~,yu_left] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_u, XU_left(:,i), options_ODE);
%     f1_p3 = plot3(yu_left(:,1),yu_left(:,2),yu_left(:,3),'m');
% end
% for i = 1:N
%     [~,yu_right] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_u, XU_right(:,i), options_ODE);
%     f1_p4 = plot3(yu_right(:,1),yu_right(:,2),yu_right(:,3),'r');
% end
plot3(x_1(:,1), x_1(:,2),x_1(:,3),'k');
plot3(L1(1),L1(2),L1(3),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);

view_angle = [40 30];

view(gca,view_angle(1),view_angle(2));
% xlim([0.96 1.02]);
% ylim([-0.02 0.02]);
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
grid on
hold off

%% 面外，面内両方の成分を持つ多様体をプロット
[XS_left, XS_right, XU_left, XU_right, Y_2] = fun_manifolds_double(mu, x0_1, 2*t0_1, N, xpert, options_ODE);
tf = 15;
tspan_s = [tf_2 0];
tspan_u = [0 tf_2];
% プロット

figure();
hold on

for i = 1:N
    [~,ys_left] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_s, XS_left(:,i), options_ODE);
    ys_left = flipud(ys_left);
    f1_p1 = plot3(ys_left(:,1),ys_left(:,2),ys_left(:,3),'g');
end
for i = 1:N
    [~,ys_right] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_s, XS_right(:,i), options_ODE);
    ys_right = flipud(ys_right);
    f1_p2 = plot3(ys_right(:,1),ys_right(:,2),ys_right(:,3),'Color','#77AC30');
end
% for i = 1:N
%     [~,yu_left] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_u, XU_left(:,i), options_ODE);
%     f1_p3 = plot3(yu_left(:,1),yu_left(:,2),yu_left(:,3),'m');
% end
% for i = 1:N
%     [~,yu_right] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_u, XU_right(:,i), options_ODE);
%     f1_p4 = plot3(yu_right(:,1),yu_right(:,2),yu_right(:,3),'r');
% end
plot3(x_1(:,1), x_1(:,2),x_1(:,3),'k');
plot3(L1(1),L1(2),L1(3),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);

view_angle = [40 30];

view(gca,view_angle(1),view_angle(2));
xlim([0.7, 1.3]);
ylim([-0.2 0.2]);
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
grid on
hold off