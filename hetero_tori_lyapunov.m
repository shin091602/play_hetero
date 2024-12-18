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
addpath("funvtion_tori/");

%% 初期設定
[mu, ~, ~, ~]      = parameter(2); % Earth-Moon
[L1, L2, L3, L4, L5] = librationPoints(mu);

% 定数の定義
p = dictionary();
p("Ms") = 5.9724e+24;%mass of Sun　→ earth
p("M1") = 7.3458e+22;%mass of Earth　→ moon
p("mu") = 1.21536E-02;%mu--EM
p("G") = 6.67300000000000e-11;%gravitational constant
p("chara_length_CR3BP") = 3.844e+8;%[m]
p("chara_mass_CR3BP") = 5.9724e+24;%[kg]
p("chara_time_CR3BP") = 3.7748e+5;%[s]
p("N_CR3BP_SE") = 1.6645e-05;

% トーラス計算のための定数の定義
p("d") = 6; % dimension of variables
p("N") = 25; % number of points in circle　NとMは一緒で奇数, defalut: 15
p("M") = 25; %number of circle --the number of N has to be equal to the number of M.
p("K") = 2e-2; %radius of circle これ何？
p("null_acr") = p("K")*1e-1; %null accuracy
p("Iteration") = 50; %iteration of GMOS
p("Threshold") = 1e-8; %converge threshold
p("phsflg") = true;%flag of phase condition --ON
p("ctswitch") = true;%center tangenet switch
p("pcon") = 2;%phase condition index
p("congflg") = false;%converged flag
p("HLflg") = true;%Halo:true, Lyapunov:false
p("xpert") = 1e-4;%deorbit velocity
%** for manifold **%
p("snap_ini_time") = 0.8*pi;%initial time of snapshots %default:pi/2
p("snap_fin_time") = 1.2*pi;%final time of snapshots %default:0.8pi
p("snap_span") = 2;%number of snapshots defalut:6
p("num_iter") = 200;%number of grid points for interpolation
snap_time = linspace(p("snap_ini_time"),p("snap_fin_time"),p("snap_span"));

% トーラス計算の際のPACのための定数の定義
pacqp = dictionary();
pacqp("n") = 5;%iteration of PAC --if one need to calculate the family, change its value.
pacqp("tol") = 1e-7;%error tolerance default = 1e-7
pacqp("itmax") = 100;%max iteration
pacqp("optit") = 5;%maximum step length
pacqp("smax") = 1e-5;%maximum step length
pacqp("jcr") = p("d")*p("N")*p("M")+4;%Jacobian rows to consider for initial family tangent
pacqp("issparse") = true;%sparse Jacobian flag
pacqp("N") = p("N")*p("M");%number of continuation function points
pacqp("fpidx") = p("d")*p("N")*p("M");%continuation function point indices
pacqp("plotflg") = false;%plot flag
pacqp("coflg") = false;%correction only flag
pacqp("fdcheck") = false;%finite difference check

% リアプノフ軌道計算のための定数の定義
lypp = dictionary();
lypp("C_error") = 1e-13; %ヤコビ定数の精度
lypp("iteration_max") = 100; %一回の計算の限界値
lypp("threshold") = 1e-10; %閾値
lypp("count_max") = 8000; %計算の限界値
lypp("xpert") = 1e-4; %deorbit verocity
lypp("N") = 100; % 軌道の分割数
lypp("tf") = 15; % 多様体の伝播時間

%% オプション設定
options_ODE_1   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14, 'Events', @(t,x) odestop_hetero_1(t,x,mu));
options_ODE_2   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14, 'Events', @(t,x) odestop_hetero_2(t,x,mu));
options_ODE   = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);

%% トーラスの計算

% L2-ハロー軌道の初期値
x0_t = [1.116913489266650; 0; 0.021776489613409; 0; 0.186114812512697; 0];
T_t = 3.407900091010860 * 2;

% トーラスを作成する周期軌道とfix pointをプロット
z0 = [x0_t;T_t];
[zpo,Xd] = po_mulshoot_initialization(z0,p);

% フーリエ変換に関する処理
[R,IR,DR] = fourier_matrix(p);
FR = cell(5,2);
FR{1,1} = "Fr";
FR{1,2} = R;
FR{2,1} = "IFr";
FR{2,2} = IR;
FR{3,1} = "DFr";
FR{3,2} = DR;
FR{4,1} = "Q";
FR{4,2} = zeros(p("d")*p("N"),p("d")*p("N"));
FR{5,1} = "R";
FR{5,2} = zeros(p("d")*p("N"),p("d")*p("N"));

% 初期点,一周期後の点，ρ回転させた点
[z0,phi0,Ud0,rho,Ec] = fun_center_manifold_cr3bp([zpo(1:p("d"));zpo(end)],p,FR);%(Ref:(4)-(5))

% PSEUDO-ARCLENGTH CONTINUATION によるトーラスの計算
sol_qpos = PAC_qpoms(z0,zpo,1e-5,Ud0,phi0,p,pacqp,FR);%(Ref:(19)-(23))

% トーラスの多様体を計算
% トーラスを指定
nm = 5;%1~pacqp("n")
solc_qpos_m = sol_qpos{nm,1};
fin_qpos_m = solc_qpos_m{1,2};

% calculate the directions of the designated manifolds
del_w_us_inter_full = direction_manifold(fin_qpos_m,p,FR);

% vector arrows
U0 = zeros(p("d")*p("M")*p('N'),1);

%% 多様体の計算のために初期値を不安定方向にずらす
check_jacobi = zeros(p("N"),p("M"));
for i=1:p("M")
    X0 = fin_qpos_m;
    del_w_us = del_w_us_inter_full(:,:,i);
    for j=1:p("N")%designate torus point
        % Create perturbation vector
        vpert = p("xpert")*norm(X0(p("d")*j-2:p("d")*j))/norm(X0(p("d")*j-5:p("d")*j-3));
        pert = [ones(3,1).*p("xpert"); ones(3,1).*vpert];

        % Initial points of torus unstable manifolds　odeに投げる初期値
        U0(p("d")*((i-1)*p("M")+j)-5:p("d")*((i-1)*p("M")+j)) = X0(p("d")*((i-1)*p("M")+j)-5:p("d")*((i-1)*p("M")+j))-pert.*del_w_us(:,j);
        check_jacobi(i,j) = Jacobi_const(U0(p("d")*((i-1)*p("M")+j)-5:p("d")*((i-1)*p("M")+j)), mu);
    end
end

%% リアプノフ軌道の計算

% ターゲットとなるトーラスのヤコビ定数を取得
C_xn = mean(check_jacobi(:));
delta = 4e-5;
count = 0;
detect_lyp = zeros(8,1);

% L1-リアプノフ軌道の初期値
x0_1 = [0.836900082907655, 0, 0, 0, 1.770874936727959e-06, 0];
t0_1 = 1.3458;

%Loop for first orbit------------------------------------------------------------------------------------------------
while 1
    count = count + 1;

    % Differential correction-----------------------------------------------
    for iteration = 1:lypp("iteration_max")
        [x_n_1, t_n_1, C] = fun_differential_correction_cr3bp(x0_1, t0_1, mu, options_ODE);

        tspan = [0 2*t_n_1];
        [t_corrected_1, x_corrected_1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan, x_n_1, options_ODE);

        x_error_1 = norm(x_corrected_1(end, :) - x_corrected_1(1, :));

        if x_error_1 < lypp("threshold")
            break;
        end

        if x_error_1 > 1e+3
            disp('calculation diverged');
            return;
        end

        if iteration == lypp("iteration_max")
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

    if (count >= 100) && (mod(count,100)==0)
        disp(strcat('count = ', num2str(count)));
    end

    if abs(C-C_xn) < lypp("C_error")
        additional = [C; t_n_1; x_n_1(:,1)];
        detect_lyp = additional;
        break;
    end

    x0_1(1) = x0_1(1) + delta;

    if count == lypp("count_max") - 1
        break
    end
end
disp('L1-sucucessfuly finished');

% 周期軌道計算
x0_1=detect_lyp(3:8);
t0_1=detect_lyp(2);

tspan_1 = [0 2*t0_1];
[t_1, x_1] = ode113(@(t, x) fun_cr3bp(t, x, mu), tspan_1, x0_1, options_ODE);

[XS_left, XS_right, XU_left, XU_right, Y_2] = fun_manifolds_double(mu, x0_1, 2*t0_1, lypp("N"), lypp("xpert"), options_ODE);
tspan_s = [lypp("tf") 0];
tspan_u = [0 lypp("tf")];

%% プロット
hinterm = figure();
hold on
grid on
box on
xlim([0.7, 1.3]);
ylim([-0.2 0.2]);
zlim([-0.06 0.06])
view([7 25]);
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');

% 元のトーラスをプロット
solc_qpos = sol_qpos{nm,1};
fin_qpos = solc_qpos{1,2};
%meshgrid
tht0 = linspace(-2*pi,4*pi,p("N")*3);
tht1 = linspace(-2*pi,4*pi,p("M")*3)';
%interpolated grid
tht0_n = linspace(-2*pi,4*pi,p("num_iter"));
tht1_n = linspace(-2*pi,4*pi,p("num_iter"));
% interpolate
[ri,~] = interpolate_torus(fin_qpos,tht0,tht1,tht0_n,tht1_n,p);
%color setting
CO = zeros(p("num_iter"),p("num_iter"),3);
CO(:,:,1) = 0.3010.*ones(p("num_iter")); % red
CO(:,:,2) = 0.5450.*ones(p("num_iter")); % green
CO(:,:,3) = 0.7.*ones(p("num_iter")); % blue
hs = surf(ri{1,1},ri{2,1},ri{3,1},CO);
% hold on

% 多様体をプロット
xe_tori = zeros(6,p("M"),p("N"));
for i=1:p("N")
    for j=1:p("M")
        % [~,X_sample] = ode113(@(t,x) fun_cr3bp(t,x,p("mu")),[0 8],U0(p("d")*((i-1)*p("M")+j)-5:p("d")*((i-1)*p("M")+j)));
        [~,X_sample,te,xe,~] = ode113(@(t,x) fun_cr3bp(t,x,p("mu")),[0 8],U0(p("d")*((i-1)*p("M")+j)-5:p("d")*((i-1)*p("M")+j)),options_ODE_2);
        if xe
            xe_tori(:,j,i) = xe;
        end
        plot3(X_sample(:,1),X_sample(:,2),X_sample(:,3),"r","LineWidth",0.5);
        hold on
    end
end

% 面を補完
[rim,~,~,~] = interpolate_manifold(fin_qpos_m,del_w_us_inter_full,tht0,tht1,tht0_n,tht1_n,p);

%color
CO_un = zeros(p("num_iter"),p("num_iter"),3);
CO_un(:,:,1) = 0.7010.*ones(p("num_iter")); % red

% for visualization
shading interp
lightangle(27,36)
lightangle(27,36)

% リアプノフ軌道のプロット
xes_right = zeros(6,100,lypp("N"));
for i = 1:lypp("N")
    for j = 1:100
        [~,ys_right,~,xe,~] = ode113(@(t,x) fun_cr3bp(t,x,mu), tspan_s, XS_right(:,j,i), options_ODE_1);
        ys_right = flipud(ys_right);
        if xe
            xes_right(:,j,i) = xe;
        end
        f1_p2 = plot3(ys_right(:,1),ys_right(:,2),ys_right(:,3),'Color','#77AC30','LineWidth',0.2);
    end
end

plot3(x_1(:,1), x_1(:,2),x_1(:,3),'k');
plot3(L1(1),L1(2),L1(3),'*','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10);

hold off

%%
figure();
hold on
for i = 1:15
    for j = 1:15
        plot3(xe_tori(2,j,i),xe_tori(3,j,i),xe_tori(4,j,i),'ro')
    end
end
for i = 1:lypp("N")
    for j = 1:100
        plot3(xes_right(2,j,i),xes_right(3,j,i),xes_right(4,j,i),'bo')
    end
end
xlabel('y');
ylabel('z');
zlabel('dx');
grid on
hold off

%%
figure();
hold on
for i = 1:15
    for j = 1:15
        plot3(xe_tori(2,j,i),xe_tori(3,j,i),xe_tori(6,j,i),'ro')
    end
end
for i = 1:lypp("N")
    for j = 1:100
        plot3(xes_right(2,j,i),xes_right(3,j,i),xes_right(6,j,i),'bo')
    end
end
xlabel('y');
ylabel('z');
zlabel('dx');
xlim([0 0.1]);
ylim([-0.025, 0.01]);
zlim([-0.0001,0]);
grid on
hold off