function distribution()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROGRAMA DISTRIBUCIÓN DE PARTÍCULAS EN COMPARTIMENTOS
    %
    % Versión en Octave del código Fortran original.
    % Simula la distribución aleatoria de N partículas en P compartimentos
    % y compara con la solución analítica (para P=2).
    %
    % Uso:
    %   1. Ejecutar en Octave: `distribution()`
    %   2. Ingresar parámetros cuando se solicite.
    %   3. Resultados se guardan en 'output.dat' y 'analytical.dat'.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% --- CONSTANTES Y VARIABLES GLOBALES ---
    Maxp = 10;       % Máximo número de compartimentos
    Maxn = 100000;   % Máximo número de partículas

    % Matriz de distribución (P x N+1) y contadores de partículas por compartimento
    global Dist Nop;
    Dist = zeros(Maxp, Maxn + 1);  % +1 para incluir el índice 0
    Nop = zeros(1, Maxp);

    %% --- INICIALIZACIÓN DEL GENERADOR ALEATORIO ---
    % (Equivalente a Sstmm + Genrand en Fortran)
    rand_seed = 0.001 * mod((10 + 10 * floor(now * 86400)), 1000);  % Usa el tiempo actual como semilla
    rand_seed = max(0.001, min(0.999, rand_seed));  % Asegura rango [0.001, 0.999]
    rand('state', rand_seed);  % Fija semilla (Octave)

    %% --- LECTURA DE PARÁMETROS ---
    N = input('Número de partículas (N)? ');
    P = input('Número de compartimentos (P)? ');
    Ncycle = input('Número de ciclos de simulación? ');

    % Validación de parámetros
    if (P < 2 || P > Maxp || N < 2 || N > Maxn)
        error('Error: Parámetros fuera de rango. P debe estar entre 2 y %d, N entre 2 y %d.', Maxp, Maxn);
    end

    %% --- SIMULACIÓN MONTE CARLO ---
    for I = 1:Ncycle
        for Kkk = 1:1000
            % 1. Distribuir N partículas aleatoriamente en P compartimentos
            Nop = zeros(1, P);  % Reiniciar contadores
            for part = 1:N
                compartment = 1 + floor(P * rand());  % Entero aleatorio entre 1 y P
                Nop(compartment) += 1;
            end

            % 2. Registrar en el histograma
            for J = 1:P
                if (Nop(J) <= Maxn)
                    Dist(J, Nop(J) + 1) += 1;  % +1 porque Octave indexa desde 1
                end
            end
        end
    end

    %% --- NORMALIZACIÓN DE RESULTADOS ---
    Dist = Dist / (Ncycle * 1000);  % Convierte conteos a probabilidades

    %% --- GUARDAR RESULTADOS NUMÉRICOS ---
    % Distribución simulada
    fh = fopen('output.dat', 'w');
    for Jj = 0:N
        fprintf(fh, '%d ', Jj);  % Número de partículas
        for J = 1:P
            fprintf(fh, '%f ', Dist(J, Jj + 1));  % Probabilidad para cada compartimento
        end
        fprintf(fh, '\n');
    end
    fclose(fh);

    %% --- SOLUCIÓN ANALÍTICA (PARA P=2) ---
    if (P == 2)
        fh = fopen('analytical.dat', 'w');
        for Jj = 0:N
            % Fórmula: P(n) = exp[log(N!) - log(n!) - log((N-n)!) - N*log(2)]
            log_prob = faculty(N) - faculty(Jj) - faculty(N - Jj) - N * log(2);
            prob = exp(log_prob);
            fprintf(fh, '%d %f\n', Jj, prob);
        end
        fclose(fh);
    else
        % Crear archivo vacío si P != 2
        fh = fopen('analytical.dat', 'w');
        fclose(fh);
    end

    disp('Simulación completada. Resultados en output.dat y analytical.dat.');
end

%% --- FUNCIÓN AUXILIAR: LOG-FACTORIAL ---
function result = faculty(n)
    % Calcula log(n!) = log(1) + log(2) + ... + log(n)
    % Evita overflow para n grande.
    if (n < 0)
        error('n debe ser no negativo.');
    elseif (n == 0 || n == 1)
        result = 0;
    else
        result = sum(log(1:n));
    end
end
