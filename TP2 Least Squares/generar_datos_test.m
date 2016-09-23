% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales I
% Trabajo Práctico 2: 
%   - Estimación de parámetros utilizando LS -
% 1° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Script para generar datos de test - MATLAB
% -------------------------------------------------------------------------- %

% close all;
% clear all;


% -------------------------------------------------------------------------- %
% Posición inicial
initial_position = [0, 0];

% Velocidad inicial
initial_velocity = [1, 3];

% Varianza del ruido de los acelerómetros
acel_variance = [0.25, 0.64];

% Constante universal de aceleración de la gravedad
g = 9.8; % [m/s^2]
% -------------------------------------------------------------------------- %
Ekx = -0.0301;
Eky = 0.01;
Esx = 0.0767;
Esy = -0.0175;


% -------------------------------------------------------------------------- %
% --- Ensayo de acelerómetros, para estimar los errores de escala y sesgo ---%
% -------------------------------------------------------------------------- %

N = 20001;
tita = linspace(0, 2*pi, N)';

A_real_x = -g*sin(tita);
A_real_y = -g*cos(tita);
Vx = normrnd(0, acel_variance(1), [N, 1]);
Vy = normrnd(0, acel_variance(2), [N, 1]);
datos(:,1) = A_real_x * (1 + Ekx) + Esx + Vx;
datos(:,2) = A_real_y * (1 + Eky) + Esy + Vy;

save('ensayo_test.mat', 'tita', 'datos');


% -------------------------------------------------------------------------- %
% ---- Ensayo de trayectoria ----
% -------------------------------------------------------------------------- %

N = 6000;
t = linspace(0, 60, N)';
T = t(2) - t(1);

% Px = linspace(0,10, N)';
Px = t;
% Py = -Px.*(Px-10);
Py = Px.*(Px-10).*(Px-30).*(Px-60);
% Py = -t.*(t-60);

Vex = diff(Px)/T;
Vey = diff(Py)/T;

Ax = diff(Vex)/T;
Ay = diff(Vey)/T;

N = length(Ax);
Vx = normrnd(0, acel_variance(1), [N, 1]);
Vy = normrnd(0, acel_variance(2), [N, 1]);

Aerr(:,1) = Ax * (1 + Ekx) + Esx + Vx;
Aerr(:,2) = Ay * (1 + Eky) + Esy + Vy;

P0 = [Px(1), Py(1)];

V0 = [Vex(1), Vey(1)];

save('acel_test.mat', 't', 'Aerr');
save('pos_real_test.mat', 'Px', 'Py');
save('initial_test.mat', 'P0', 'V0');
save('acel_real_test.mat', 'Ax', 'Ay');
