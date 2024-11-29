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

options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
iteration_max = 100;
threshold     = 1e-10;
delta         = 2e-5;
count         = 0;
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

%% initial values
x0 = [1.155693450711950, 0, 0, 0, 1.666612527554842e-06, 0];
t0 = 1.68665;

%% calc Lyapunov family
f1 = figure();

while 1
    count = count + 1;

    % differential correction
    for iteration = 1:iteration_max
        [x_n, t_n, C] = fun_differential_correction_cr3bp(x0, t0, mu, options_ODE);

        tspan = [0 2*t_n];
        [t_corrected, x_corrected] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan, x_n, options_ODE);

        x_error = norm(x_corrected(end, :) - x_corrected(1, :));

        if x_error < threshold
            break;
        end

        if x_error > 1e+3
            disp('calculation diverged');
            return;
        end

        if iteration == iteration_max
            disp('do not finish');
            return;
        end

        x0 = x_n;
        t0 = t_n;
    end

    % plot orbit
    if (count >= plot_interval) && (mod(count, plot_interval) == 0)
        Jacobi(count/plot_interval, 1) = C;
        x0_corrected(count/plot_interval, :) = x_corrected(1, :);
        t0_corrected(count/plot_interval, 1) = t_corrected(end);
        rgb = Interp_c(C);
        plot(x_corrected(:, 1), x_corrected(:, 2), 'Color', rgb);
        hold on
        disp(strcat('count = ', num2str(count)));
    end

    x0(1) = x0(1) - delta;

    if count == count_max
        break
    end
end

f1_p1 = plot(L2(1), L2(2), '*', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
f1_p2 = plot((1-mu), 0, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
axis equal
colormap jet;
c = colorbar;
ylabel(c, 'Jacobi constant [-]', 'FontSize', 15);
caxis([Jacobi_min Jacobi_max]);
c.Ticks = linspace(Jacobi_min, Jacobi_max, 6);
xlabel('$x$[km]');
ylabel('$y$[km]');
grid on
legend([f1_p1, f1_p2], {'$L_2$', 'Earth'}, 'Location', 'northwest');
hold off

%% check the eigenvalue
video = VideoWriter('video/lyapunov_eigenvalue_1.mp4', 'MPEG-4');
video.FrameRate = 5;
open(video);
for i = 1:count_max/plot_interval
    X0 = [x0_corrected(i,:),reshape(eye(6),1,[])];
    tspan = [0,t0_corrected(i)];
    [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),tspan,X0,options_ODE);

    monodromy = reshape(Y(end,7:end),6,6);
    [V,D] = eig(monodromy);

    theta = linspace(0, 2*pi, 1001)';
    x_circle = cos(theta);
    y_circle = sin(theta);

    figure;
    hold on
    plot(x_circle,y_circle,'k');
    plot(diag(D),'o','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',10);
    axis equal
    xlabel('real part');
    ylabel('imaginary part');
    grid on
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    text(0.1, 0.9, ['jacobi constant: ' num2str(Jacobi(i))], 'Units', 'normalized','FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
    hold off

    frame = getframe(gcf);
    writeVideo(video, frame);

    close(gcf);

end

close(video)