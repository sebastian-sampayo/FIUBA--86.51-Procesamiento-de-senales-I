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
% Script principal - MATLAB
% -------------------------------------------------------------------------- %

close all;
clear all;

% -------------------------------------------------------------------------- %
% ---- Parámetros iniciales ----
% -------------------------------------------------------------------------- %
% ## TEST ## (Para usar datos generados con "generar_datos_test.m")
test = true; % false: Sin test, true: Con test
if (test)
    generar_datos_test
end
% ##
% Nombre del archivo de datos del ensayo realizado con los acelerómetros
% Contiene: 
%   tita:   los valores del ángulo rotado θ ($\theta$). Nx1
%   datos:  los valores de aceleraciones medidas en x e y, respectivamente, 
%           correspondientes a cada ángulo θ ($\theta$). Nx2
test_file_name = 'ensayo.mat';
% ## TEST ##
if (test)
    test_file_name = 'ensayo_test.mat';
end
% ##
load(test_file_name);

% Nombre del archivo de los puntos A, B, C y D
% Contiene:
%   A, B, C y D. 
points_file_name = 'puntos.mat';
load(points_file_name);

% Nombre del archivo de la aceleración medida durante la trayectoria del vehículo
% Contiene:
%   t:      tiempo en segundos
%   Aerr:   la aceleración en los ejes x e y, respectivamente.
acel_file_name = 'acel.mat';
% ## TEST ##
if (test)
    acel_file_name = 'acel_test.mat';
end
% ##
load(acel_file_name);

% Posición inicial
initial_position = [0, 0];
% ## TEST ##
if (test)
    load('initial_test.mat');
    initial_position = P0;
end
% ##

% Velocidad inicial
initial_velocity = [1, 3];
% ## TEST ##
if (test)
    initial_velocity = V0;
end
% ##

% Varianza del ruido de los acelerómetros
acel_variance = [0.25, 0.64];

% Constante universal de aceleración de la gravedad
g = 9.8; % [m/s^2]
% -------------------------------------------------------------------------- %


% -------------------------------------------------------------------------- %
% ---- Ejercicio 1) ----
% -------------------------------------------------------------------------- %
% Indique cómo haría para, a partir de los datos del ensayo, armar un modelo 
% (y = Ax + η) que le permita estimar los sesgos y factores de escala para 
% cada uno de los acelerómetros.

% y = Ax + η
% b = Hc + v
% b := Aceleración medida, Nx1
% H := Matriz de 2 columnas. 
%        Columna 1: Aceleración Real. Columna 2: todos 1's.
% c := Vector columna de 2 componentes. 
%        Componente 1: (1+Error de Escala). Componente 2: Error de sesgo
% v := Ruido, Nx1.



% -------------------------------------------------------------------------- %
% ---- Ejercicio 2) ----
% -------------------------------------------------------------------------- %
% Estime los valores de los sesgos y factores de escala, a partir de los datos 
% del ensayo (archivo: ensayo.mat) que se le suministraron. Calcule la varianza 
% del estimador.

% [x_hat, cov_x_hat] = mvue(y, H, Rv) 

% Cantidad de muestras:
N = length(tita);
% Se asume Ruido Blanco Gaussiano, por lo tanto, la matriz de correlación es:
% Teóricamente la matriz identidad multiplicada por la varianza del ruido.
% Sin embargo, como esto utiliza demasiada memoria innecesariamente, lo reduzco
% a un escalar igual a la varianza. Todo da igual, ver las fórmulas.
Rv_x = acel_variance(1); % * diag(ones(N,1)); 
Rv_y = acel_variance(2); % * diag(ones(N,1));

H = [-g*sin(tita), ones(N,1)];
[c_hat_x, cov_c_hat_x] = mvue(datos(:,1), H, Rv_x);
clear H;

H = [-g*cos(tita), ones(N,1)];
[c_hat_y, cov_c_hat_y] = mvue(datos(:,2), H, Rv_y);
clear H;

% Pasaje de variables de errores obtenidos:
scale_error_x = c_hat_x(1) - 1
scale_error_y = c_hat_y(1) - 1
bias_error_x = c_hat_x(2)
bias_error_y = c_hat_y(2)

scale_error_x_variance = cov_c_hat_x(1,1)
scale_error_y_variance = cov_c_hat_y(1,1)
bias_error_x_variance = cov_c_hat_x(2,2)
bias_error_y_variance = cov_c_hat_y(2,2)
% -------------------------------------------------------------------------- %


% -------------------------------------------------------------------------- %
% ---- Ejercicio 3) ----
% -------------------------------------------------------------------------- %
% Calcule la trayectoria del vehículo. De las posiciones de los cuatro puntos 
% suministrados A,B,C o D (archivo: puntos.mat) ¿a cuál de ellos llega el
% vehículo?.

% Estimador de Aceleración real:
A_real_x = (Aerr(:,1) - bias_error_x)/(1 + scale_error_x);
A_real_y = (Aerr(:,2) - bias_error_y)/(1 + scale_error_y);

% Calculo la velocidad y la posición a partir de la aceleración medida 
% aproximando las integrales con el método del trapecio.
T_step = t(2)-t(1);
integrate = @(x) (T_step*cumtrapz(x));
velocity_x_measured = initial_velocity(1) + integrate(Aerr(:,1));
velocity_y_measured = initial_velocity(2) + integrate(Aerr(:,2));
position_x_measured = initial_position(1) + integrate(velocity_x_measured);
position_y_measured = initial_position(2) + integrate(velocity_y_measured);

% Calculo la velocidad y la posición a partir del estimador de aceleración real
velocity_x_real = initial_velocity(1) + integrate(A_real_x);
velocity_y_real = initial_velocity(2) + integrate(A_real_y);
position_x_real = initial_position(1) + integrate(velocity_x_real);
position_y_real = initial_position(2) + integrate(velocity_y_real);

figure
hold all
plot(position_x_measured, position_y_measured, 'k--', 'Linewidth', 3)
plot(position_x_real, position_y_real, 'Linewidth', 3)
if (test == false)
    plot(A(1), A(2), 'o', 'MarkerSize', 10)
    plot(B(1), B(2), 'o', 'MarkerSize', 10)
    plot(C(1), C(2), 'o', 'MarkerSize', 10)
    plot(D(1), D(2), 'o', 'MarkerSize', 10)
end
legend('Posición medida', 'Posición real estimada', 'A', 'B', 'C', 'D',...
    'Location', 'NorthWest');
title('Trayectoria del vehículo');
grid on;
xlabel('x');
ylabel('y');

% print('-dpng', 'trayectoria.png');

% ## TEST ##
if (test)
    load('pos_real_test.mat');
    plot(Px, Py, 'r--')
%     axis([0 max(Px) min(Py) max(Py)]);
%     print('-dpng', 'trayectoria_test.png');
%     load('acel_real_test.mat');
%     figure
%     hold all
%     plot(Ax)
%     plot(A_real_x)
%     plot(Aerr(:,1))
%     legend('Original', 'Estimada');
%     
%     figure
%     hold all
%     plot(Ay)
%     plot(A_real_y)
%     plot(Aerr(:,2))
%     legend('Original', 'Estimada');
end
% ##
