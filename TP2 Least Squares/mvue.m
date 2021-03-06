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
% Función MVUE (Minimum Variance Unbiased Estimator) - MATLAB
% -------------------------------------------------------------------------- %
%
%% Función MVUE (Minimum Variance Unbiased Estimator)
%
%   Calcula el estimador insesgado de mínima varianza para un modelo lineal:
%
%        $$  y = Hx + v  $$
%
%   donde 'x' es determinístico pero desconocido mientras que 'v' es ruido
%   aleatorio, y H es conocida.
%
%   Uso:
%     [x_hat, cov_x_hat] = mvue(y, H, Rv)
%
%   donde:
%     y:          es un vector de Nx1
%     H:          es una matriz de NxP
%     Rv:         es la matriz de correlación del ruido 'v', de NxN
%     x_hat:      es el estimador de 'x', de Px1
%     cov_x_hat:  es la matriz de covarianza del estimador 'x_hat', de PxP
%
%   Fórmulas utilizadas:
%     
%     $$  x
%     $$  cov_x
%     $$  Rv

function [x_hat, cov_x_hat] = mvue(y, H, Rv)
  x_hat = inv(H'*inv(Rv)*H) * H'*inv(Rv)*y;
  cov_x_hat = inv(H'*inv(Rv)*H);
end
