% FIUBA
% Materia: Procesamiento de Señales I
% TP
% Sampayo, Sebastián Lucas 
% 
%% Función WLS - Weighted Least Squares -
% 
%   Uso:
%     h = WLS(W, H, M, flg_type)
%  
%   donde:
%     W es la función peso evaluada entre 0 y $\pi$
%     A es la función de amplitud de la respuesta en frecuencia deseada evaluada entre 0 y $\pi$
%        $$ H \left( e^{j\omega} \rigth) = A(\omega) e^{-j(frac{\omega M}{2} - \beta)} $$
%     M es el orden del filtro a diseñar
%     flg_type es el Tipo de filtro a diseñar (I, II, III o IV = {1,2,3,4})
%     h son los coeficientes del filtro diseñado

function h = WLS(M, A, W, flg_type)
  w_points = length(W);
  % Si la dimensión de W es distinta a la de A, error
  if (w_points~=length(A))
    h = NaN;
    return;
  end

  % Construyo un vector w de "w_points" puntos entre 0 y pi
  w = linspace(0, pi, w_points)';
  % Para esta implementación, agregar "_bands" a los nombres de los argumentos
  % y agregar f_bands como argumento, que es las frecuencias de cada banda normalizadas en pi
  % Además checkear que A_bands y f_bands tengan la misma dimensión, y que W_bands tenga la mitad de dicha dimensión
  % alfa = 20;
  % w_points = alfa * M;
  % Construyo los vectores Amplitud(A) y Peso (W) evaluados en cada w entre 0 y pi.
  % A = zeros(w_points,1);
  % W = zeros(w_points,1);
  % for i = 1:2:length(A_bands)
  %   np_1 = floor(f_bands(i)*w_points); % Índice izquierdo de la banda i-ésima
  %   np_2 = floor(f_bands(i+1)*w_points); % Índice derecho de la banda i-ésima
  %   A(np_1:np_2) = linspace(A_bands(i), A_bands(i+1), np_2-np_1+1)'; % Interpolo linealmente las amplitudes de los extremos de la banda
  %   W(np_1:np_2) = W_bands(i/2)*ones(np_2-np_1+1,1); % Construyo la función peso
  % end

%  r = M/2 * (flg_type==1) + (M-1)/2 * (flg_type==2 || flg_type==4) + (M-2)/2 * (flg_type==3);
%  F = ones(w_points,1) * (flg_type==1) + cos(w/2) * (flg_type==2) + sin(w) * (flg_type==3) + sin(w/2) * (flg_type==4);

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
  end

  W = diag(W);
  F = diag(F);
  C = F * cos(w * (0:r));
  % Queremos proyectar el vector A sobre Col(C),
  % utilizando el producto interno definido por W.
  % Las coordenadas de dicha proyección en Col(C) son "g".
  % W*A = (W*C)*g
  % (W*C)' * W*A = (W*C)'*(W*C)*g
  % C'*W^2*A = C'*W^2*C * g
  % => g^ = ( C'*W^2*C )^-1 * C' *W*A
  %    g^ = (C*W)'' * (W*A)
  
  g = ((W*C) \ (W*A));

  % Existe una relación lineal entre "g(n)" y "h(n)"
  % Mapeo de coeficientes "g" a "h"
  h = zeros(M+1,1);

  if (flg_type == 1)
    n = (1:r)';
    h(n) = .5 * g(r+2 - n);
    h(r+1) = g(1);
    h(M+2-n) = h(n); % El +2 es debido a que Matlab no permite índice '0'.

  elseif (flg_type == 2)
    n = (2:r)';
    h(1) = .5 * g(r+1);
    h(n) = .25 * ( g(r+2 - n) + g(r+3 - n) );
    h(r+1) = .5 * ( g(1) + .5*g(2) );
    h(M+2-n) = h(n); % El +2 es debido a que Matlab no permite índice '0'.

  elseif (flg_type == 3)
    n = (3:r)';
    h(1) = .25 * g(r+1);
    h(2) = .25 * g(r);
    h(n) = .25 * ( g(r+2 - n) - g(r+4 - n) );
    h(r+1) = .25 * ( 3*g(1) - g(3) );
    h(M+2-n) = -h(n); % El +2 es debido a que Matlab no permite índice '0'.

  elseif (flg_type == 4)
    n = (2:r)';
    h(1) = .25 * g(r+1);
    h(n) = .25 * ( g(r+2 - n) - g(r+3 - n) );
    h(r+1) = .25 * ( 2*g(1) - g(2) );
    h(M+2-n) = -h(n); % El +2 es debido a que Matlab no permite índice '0'.

  end

end
