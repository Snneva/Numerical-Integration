%% Setup
clc; clear; close all;

% Rentang integrasi
a = 0; 
b = 5; 
N = 100;               % jumlah partisi (harus genap untuk Simpson)
dx = (b - a)/N;       

% Fungsi
y = @(t) -(cos(2*pi*t))/(2*pi) - (sin(3*pi*t))/(3*pi);

% Integral eksak
I_exact = integral(y, a, b);

% Titik partisi
x = linspace(a,b,N+1);
y_vals = y(x);

t = linspace(a,b,1000);

%% 1. Riemann Left, Mid, Right
% --- Left
x_left = x(1:end-1);
R_left = sum(y(x_left))*dx;

% --- Right
x_right = x(2:end);
R_right = sum(y(x_right))*dx;

% --- Mid
x_mid = (x(1:end-1)+x(2:end))/2;
R_mid = sum(y(x_mid))*dx;

% --- Plot Left, Mid, Right
figure;
subplot(3,1,1);
plot(t,y(t),'b','LineWidth',2); hold on; grid on;
for i=1:N
    fill([x(i) x(i+1) x(i+1) x(i)], [0 0 y(x(i)) y(x(i))], ...
        'r','FaceAlpha',0.3,'EdgeColor','r');
end
title(['Riemann Left = ', num2str(R_left)]);

subplot(3,1,2);
plot(t,y(t),'b','LineWidth',2); hold on; grid on;
for i=1:N
    xm = (x(i)+x(i+1))/2;
    fill([x(i) x(i+1) x(i+1) x(i)], [0 0 y(xm) y(xm)], ...
        'g','FaceAlpha',0.3,'EdgeColor','g');
end
title(['Riemann Mid = ', num2str(R_mid)]);

subplot(3,1,3);
plot(t,y(t),'b','LineWidth',2); hold on; grid on;
for i=1:N
    fill([x(i) x(i+1) x(i+1) x(i)], [0 0 y(x(i+1)) y(x(i+1))], ...
        'm','FaceAlpha',0.3,'EdgeColor','m');
end
title(['Riemann Right = ', num2str(R_right)]);

%% 2. Trapezoidal Rule
I_trap = (dx/2)*(y_vals(1) + 2*sum(y_vals(2:end-1)) + y_vals(end));

figure;
plot(t,y(t),'b','LineWidth',2); hold on; grid on;
for i=1:N
    fill([x(i) x(i+1) x(i+1) x(i)], [0 0 y_vals(i+1) y_vals(i)], ...
        'c','FaceAlpha',0.3,'EdgeColor','k');
    plot([x(i) x(i+1)], [y_vals(i) y_vals(i+1)], 'k', 'LineWidth',1);
end
title(['Trapezoidal Rule = ', num2str(I_trap)]);

%% 3. Simpson's Rule (N harus genap)
I_simp = (dx/3)*(y_vals(1) + ...
    4*sum(y_vals(2:2:end-1)) + ...
    2*sum(y_vals(3:2:end-2)) + ...
    y_vals(end));

figure;
plot(t,y(t),'b','LineWidth',2); hold on; grid on;
for i=1:2:N
    xx = linspace(x(i),x(i+2),50);
    % Interpolasi kuadratik (Lagrange)
    L0 = ((xx-x(i+1)).*(xx-x(i+2)))/((x(i)-x(i+1))*(x(i)-x(i+2)));
    L1 = ((xx-x(i)).*(xx-x(i+2)))/((x(i+1)-x(i))*(x(i+1)-x(i+2)));
    L2 = ((xx-x(i)).*(xx-x(i+1)))/((x(i+2)-x(i))*(x(i+2)-x(i+1)));
    yy = y_vals(i)*L0 + y_vals(i+1)*L1 + y_vals(i+2)*L2;
    fill([xx fliplr(xx)], [yy zeros(size(xx))], 'y','FaceAlpha',0.3,'EdgeColor','none');
end
title(['Simpson Rule = ', num2str(I_simp)]);

%% 4. Monte Carlo
N_mc = 3300;
x_rand = a + (b-a)*rand(1,N_mc);
y_rand = y(x_rand);
I_mc = (b-a)*mean(y_rand);

figure;
plot(t,y(t),'b','LineWidth',2); hold on; grid on;
scatter(x_rand(1:3300), y_rand(1:3300), 20, 'r','filled');
title(['Monte Carlo = ', num2str(I_mc)]);

%% 5. Tabel Error
methods = {'Riemann Left'; 'Riemann Mid'; 'Riemann Right'; ...
           'Trapezoidal'; 'Simpson'; 'Monte Carlo'};
values  = [R_left; R_mid; R_right; I_trap; I_simp; I_mc];
errors  = 100*abs(values - I_exact)/abs(I_exact);

T = table(methods, values, errors, ...
    'VariableNames', {'Metode','Nilai_Integral','Error_Persen'});

disp('================ Komparasi ================');
disp(['Nilai Integral = ', num2str(I_exact)]);
disp(T);
