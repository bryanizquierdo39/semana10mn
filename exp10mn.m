clear; clc;

fprintf('LABORATORIO 10 METODOS NUMERICOS');

%% 1. INGRESO DE DATOS

fprintf('Ingresa f(x) como string.\n');
fprintf('Ejemplos validos:\n');
fprintf('  exp(-x.^2/2)/sqrt(2*pi)\n');
fprintf('  cos(x.^2)\n');
fprintf('  x.^5 - 2*x.^3 + 4\n\n');

func_str = input('f(x) = ','s');

a_str = input('Limite inferior a = ','s');
b_str = input('Limite superior b = ','s');

a = str2double(a_str);
b = str2double(b_str);

if isnan(a) || isnan(b)
    error('Los limites deben ser numericos.');
end

%% 2. CREAR FUNCION DE FORMA SEGURA

try
    f = str2func(['@(x) ' func_str]);

    % prueba de evaluacion
    f(0);

catch
    error(['Error en la funcion ingresada.' newline ...
           'Ejemplo valido: x.^5 - 2*x.^3 + 4']);
end

%% 3. VALOR DE REFERENCIA

valor_exacto = integral(f,a,b,'AbsTol',1e-12,'RelTol',1e-12);

%% 4. PARAMETROS

N = 12;       % 12 subintervalos -> 13 puntos
n_gauss = 5;  % valor por defecto para funciones no polinomicas

%% 5. DETECCION DE POLINOMIO

es_polinomio = isempty(regexp(func_str,...
    'exp|cos|sin|tan|log|sqrt|asin|acos|atan|sinh|cosh|tanh',...
    'once'));

if es_polinomio

    grado = grado_polinomio(func_str);

    n_optimo = ceil((grado + 1)/2);

    fprintf('\nDetectado polinomio grado %d.\n',grado);
    fprintf('Gauss-Legendre usara n=%d para error teorico cero.\n',n_optimo);

else

    n_optimo = n_gauss;

    fprintf('\nFuncion no polinomica.\n');
    fprintf('Gauss-Legendre usara n=%d.\n',n_optimo);

end

%% 6. CALCULO DE LOS METODOS

[trap_aprox, trap_p]     = trapecio_compuesto(f,a,b,N);
[simp13_aprox, simp13_p] = simpson_13_compuesto(f,a,b,N);
[simp38_aprox, simp38_p] = simpson_38_compuesto(f,a,b,N);
[gauss_aprox, gauss_p]   = gauss_legendre(f,a,b,n_optimo);

%% 7. TABLA DE RESULTADOS

fprintf('\n');
fprintf('Metodo                | Puntos F(x) | Aproximacion    | Error Absoluto\n');
fprintf('-----------------------------------------------------------------------\n');

imprimir_fila('Trapecio',trap_p,trap_aprox,valor_exacto);
imprimir_fila('Simpson 1/3',simp13_p,simp13_aprox,valor_exacto);
imprimir_fila('Simpson 3/8',simp38_p,simp38_aprox,valor_exacto);
imprimir_fila('Gauss-Legendre',gauss_p,gauss_aprox,valor_exacto);

fprintf('-----------------------------------------------------------------------\n');
fprintf('Valor exacto (integral) = %.12f\n',valor_exacto);

%% 8. ANALISIS

err_trap  = abs(trap_aprox  - valor_exacto);
err_simp  = abs(simp13_aprox - valor_exacto);
err_gauss = abs(gauss_aprox - valor_exacto);

fprintf('\n================ ANALISIS =================\n\n');

fprintf('1) EFICIENCIA COMPUTACIONAL\n');
fprintf('Trapecio y Simpson usan %d evaluaciones.\n',trap_p);
fprintf('Gauss-Legendre usa %d evaluaciones.\n',gauss_p);
fprintf('Menos evaluaciones significa menor costo computacional.\n\n');

fprintf('2) COMPARACION DE ERRORES\n');
fprintf('Error Trapecio      : %.3e\n',err_trap);
fprintf('Error Simpson 1/3  : %.3e\n',err_simp);
fprintf('Error Gauss        : %.3e\n\n',err_gauss);

if es_polinomio
    fprintf('Para un polinomio grado %d,\n',grado);
    fprintf('Gauss con n=%d integra exactamente hasta grado %d.\n',...
        n_optimo,2*n_optimo-1);
    fprintf('Por ello el error es practicamente cero.\n\n');
end

fprintf('3) OBSERVACION\n');
fprintf('Los nodos de Gauss no son equidistantes.\n');
fprintf('Su ubicacion optimiza la precision de la integral.\n');

%% FUNCIONES LOCALES

function [I,n_eval] = trapecio_compuesto(f,a,b,N)

    h = (b-a)/N;

    x = a:h:b;

    fx = f(x);

    I = h*(0.5*fx(1) + sum(fx(2:end-1)) + 0.5*fx(end));

    n_eval = length(x);

end

function [I,n_eval] = simpson_13_compuesto(f,a,b,N)

    if mod(N,2) ~= 0
        N = N + 1;
    end

    h = (b-a)/N;

    x = a:h:b;

    fx = f(x);

    I = h/3 * ( ...
        fx(1) + ...
        fx(end) + ...
        4*sum(fx(2:2:end-1)) + ...
        2*sum(fx(3:2:end-2)));

    n_eval = length(x);

end


function [I,n_eval] = simpson_38_compuesto(f,a,b,N)

    if mod(N,3) ~= 0
        N = N + (3 - mod(N,3));
    end

    h = (b-a)/N;

    x = a:h:b;

    fx = f(x);

    suma = fx(1) + fx(end);

    for i = 2:length(x)-1

        if mod(i-1,3) == 0
            suma = suma + 2*fx(i);
        else
            suma = suma + 3*fx(i);
        end

    end

    I = 3*h/8 * suma;

    n_eval = length(x);

end

function [I,n_eval] = gauss_legendre(f,a,b,n)

    [t,w] = gauss_nodes_weights(n);

    t = t(:);
    w = w(:);

    x = (b-a)/2 .* t + (b+a)/2;

    I = (b-a)/2 * sum(w .* f(x));

    n_eval = n;

end

function [x,w] = gauss_nodes_weights(n)

    beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^(-2));

    T = diag(beta,1) + diag(beta,-1);

    [V,D] = eig(T);

    x = diag(D);

    [x,idx] = sort(x);

    w = 2*(V(1,idx).^2);

    x = x(:);
    w = w(:);

end

function imprimir_fila(metodo,puntos,aprox,exacto)

    err = abs(aprox - exacto);

    fprintf('%-21s | %-11d | %-14.8f | %.3e\n',...
        metodo,puntos,aprox,err);

end

function g = grado_polinomio(str)

    tokens = regexp(str,'x(?:\.\^|\^)(\d+)','tokens');

    if isempty(tokens)

        if contains(str,'x')
            g = 1;
        else
            g = 0;
        end

    else

        grados = zeros(1,length(tokens));

        for k = 1:length(tokens)
            grados(k) = str2double(tokens{k}{1});
        end

        g = max(grados);

    end

end