% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales I
% Trabajo Práctico 1: Diseño de filtro con fase lineal
%   - Transformador de Hilbert para codificación de voz -
% 1° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Función WLS - MATLAB
% -------------------------------------------------------------------------- %
%
%% Función WLS - Weighted Least Squares - %
%
%   Calcula los coeficientes de un filtro FIR de fase lineal generalizada
% de tipo y orden genérico (I, II, III o IV), por el método de cuadrados
% mínimos ponderados.
% 
%   Uso:
%     h = WLS(M, f, A, W, flg_type)
%  
%   donde:
%     M es el orden del filtro a diseñar
%     f es un vector de pares de puntos de frecuencia en el rango [0,1],
%       donde 1 corresponde a la frecuencia de Nyquist.
%     W es la función peso evaluada en cada banda correspondiente a cada par
%       de valores en el vector _f_.
%     A es la función de amplitud de la respuesta en frecuencia deseada evaluada
%       en los puntos especificados en _f_.
%
%      $$ H \left( e^{j\omega} \right) = A(\omega) e^{-j(\frac{\omega M}{2} - \beta)} $$
%
%      $$ A(\omega) = F(\omega) \sum_{k=0}^{r} g[k]cos(\omega k) $$
%
%     flg_type es el Tipo de filtro a diseñar (I, II, III o IV = {1,2,3,4})
%       Tipo I:
%         $$ r = \frac{M}{2} $$
%         $$ \beta = 0 $$
%         $$ F(\omega) = 1 $$
%       Tipo II:
%         $$ r = \frac{M-1}{2} $$
%         $$ \beta = 0 $$
%         $$ F(\omega) = cos(\omega/2) $$
%       Tipo III:
%         $$ r = \frac{M-2}{2} $$
%         $$ \beta = -frac{\pi}{2} $$
%         $$ F(\omega) = sin(\omega) $$
%       Tipo IV:
%         $$ r = \frac{M-1}{2} $$
%         $$ \beta = -frac{\pi}{2} $$
%         $$ F(\omega) = sin(\omega/2) $$
%
%     h son los M+1 coeficientes del filtro diseñado. En caso de argumentos
%       inválidos, devuelve NaN.

function h = WLS(M, f_bands, A_bands, W_bands, flg_type)
  % ---- Validaciones ----
  % Si la dimensión de f_bands es distinta a la de A_bands, error
  if (length(f_bands)~=length(A_bands))
    h = NaN;
    return;
  % Si dicha dimensión no es par, error
  elseif (mod(length(A_bands),2))
    h = NaN;
    return;
  % Si la dimensión de W_bands no es la mitad que la de f_bands y A_bands, error
  elseif (length(W_bands)~=(length(A_bands)/2))
    h = NaN;
    return;
  end %if
  % -----------------------

  alfa = 20;
  w_points = alfa * M;
  % Construyo un vector w de "w_points" puntos entre 0 y pi
  w = linspace(0, pi, w_points)';

  % Construyo los vectores Amplitud(A) y Peso (W) evaluados en cada w entre 0 y pi.
  A = zeros(w_points,1);
  W = zeros(w_points,1);
  for i = 1:2:length(A_bands)
    % Índice izquierdo de la banda i-ésima:
    np_1 = floor(f_bands(i)*(w_points-1)) + 1; 
    % Índice derecho de la banda i-ésima:
    np_2 = floor(f_bands(i+1)*(w_points-1)) + 1; 
    % Interpolo linealmente las amplitudes de los extremos de la banda:
    A(np_1:np_2) = linspace(A_bands(i), A_bands(i+1), np_2-np_1+1)'; 
    % Construyo la función peso:
    W(np_1:np_2) = W_bands(i/2 + .5)*ones(np_2-np_1+1,1); 
  end

  % Configuro los valores de "r" y "F" según el tipo de filtro pedido
  % Asimismo valido que los argumentos sean coherentes.
  if (flg_type == 1)
    if (mod(M,2)) 
      h = NaN;
      return; % M debe ser PAR
    end
    r = M/2;
    F = ones(w_points,1);

  elseif (flg_type == 2)
    if (~mod(M,2)) 
      h = NaN;
      return; % M debe ser IMPAR
    end
    r = (M-1)/2;
    F = cos(w/2);

  elseif (flg_type == 3)
    if (mod(M,2)) 
      h = NaN;
      return; % M debe ser PAR
    end
    r = (M-2)/2;
    F = sin(w);

  elseif (flg_type == 4)
    if (~mod(M,2)) 
      h = NaN;
      return; % M debe ser IMPAR
    end
    r = (M-1)/2;
    F = sin(w/2);

  else
    h = NaN;
    return; % El tipo de filtro solo puede ser {1,2,3,4}
  end %if

  W = diag(W);
  F = diag(F);
  C = F * cos(w * (0:r));
  % Queremos proyectar el vector A sobre Col(C),
  % utilizando el producto interno definido por W.
  % Las coordenadas de dicha proyección en Col(C) son "g".
  % donde Col(C) es el subespacio de cosenos multiplicados por F(\omega).
  % W*A = (W*C)*g
  % (W*C)' * W*A = (W*C)'*(W*C)*g
  %        => g^ = (C*W)^+ * (W*A)
  % donde (·)^+ es la operación pseudo-inversa.
  
  g = ((W*C) \ (W*A));

  % Existe una relación lineal entre "g(n)" y "h(n)", obtenida analíticamente.
  % Mapeo de coeficientes "g" a "h"
  h = zeros(M+1,1);

  if (flg_type == 1)
    n = (1:r)';
    h(n) = .5 * g(r+2 - n);
    h(r+1) = g(1);
    n = [n; r+1];
    h(M+2-n) = h(n); % El +2 es debido a que Matlab no permite índice '0'.

  elseif (flg_type == 2)
    n = (2:r)';
    h(1) = .5 * g(r+1);
    h(n) = .25 * ( g(r+2 - n) + g(r+3 - n) );
    h(r+1) = .5 * ( g(1) + .5*g(2) );
    n = [1; n; r+1];
    h(M+2-n) = h(n); % El +2 es debido a que Matlab no permite índice '0'.

  elseif (flg_type == 3)
    n = (3:r)';
    h(1) = .25 * g(r+1);
    h(2) = .25 * g(r);
    h(n) = .25 * ( g(r+2 - n) - g(r+4 - n) );
    h(r+1) = .25 * ( 2*g(1) - g(3) );
    n = [1;2; n; r+1];
    h(M+2-n) = -h(n); % El +2 es debido a que Matlab no permite índice '0'.
    % El mapeo anterior es para beta = pi/2
    % Me conviene hacerlo para beta = -pi/2:
    h = -h; % Esto es equivalente a negar las expresiones anteriores directamente

  elseif (flg_type == 4)
    n = (2:r)';
    h(1) = .25 * g(r+1);
    h(n) = .25 * ( g(r+2 - n) - g(r+3 - n) );
    h(r+1) = .25 * ( 2*g(1) - g(2) );
    n = [1; n; r+1];
    h(M+2-n) = -h(n); % El +2 es debido a que Matlab no permite índice '0'.
    % El mapeo anterior es para beta = pi/2
    % Me conviene hacerlo para beta = -pi/2:
    h = -h; % Esto es equivalente a negar las expresiones anteriores directamente

  end %if

end %function WLS
