close all;
clear all;
clc;

%% Parte experimental
% Lendo o sinal
[sinal, fs] = audioread('antarctica.wav'); %wavread

%% 1
% Cortando sinal e obtendo os parâmetros do modelo LPC com 10 coeficientes
sinal_cortado = sinal(200:439);

[ak, sig10]= lpc(sinal_cortado.*hamming(240), 10);

%% 2
% Calculando o ganho G
qsi = sig10;
G = sqrt(qsi);

%GH = G*H;
figure(5);
freqz(G,ak, 512);
% title('Resposta em frequência do filtro G.H(z)');

%% 3
% Calculando o periodograma
hold on;
periodogram(sinal_cortado,[],512);

legend('Filtro GH(z)', 'Transformada do trecho de sinal');

%% 4
% Aumentando a ordem do modelo LPC
figure(6);

% n = 20
[ak20, sig20]= lpc(sinal_cortado.*hamming(240), 20);
freqz(sqrt(sig20), ak20, 512);

%n = 40
hold on;
[ak40, sig40]= lpc(sinal_cortado.*hamming(240), 40);
freqz(sqrt(sig40), ak40, 512);

%n = 80
hold on;
[ak80, sig80]= lpc(sinal_cortado.*hamming(240), 80);
freqz(sqrt(sig80), ak80, 512);

hold on;
periodogram(sinal_cortado,[],512);

legend('Ordem 20', 'Ordem 40', 'Ordem 80', 'Periodograma');

%% 5
% plotando o Yaapt

% Figuras 1, 2, 3 e 4
pitch = yaapt(sinal, fs, 1, [], 1, 1);
pitch = [0 pitch];

%% 6 
% Testando o codificador

% Calculando os parâmetros do modelo LPC
aks = zeros(11,46);
sigs = zeros(46,1);
Gs = zeros(138,1);
gerado = zeros(138*80,1);

% Gera coefiientes LPC para cada 30ms 
for i = 1:46
    [aks(:,i), sigs(i)] = lpc(sinal(((i-1)*240 + 1):(i*240)).*hamming(240), 10);
end

% Calcula o valor de G para cada 10ms
for i = 1:138
    if (pitch(i) ~= 0) % Sinal sonoro
        Gs(i) = sqrt((8000/pitch(i))*sigs(ceil(i/3)));
    else % Sinal surdo
        Gs(i) = sqrt(sigs(ceil(i/3)));
    end
end

for i = 1:138
    if (pitch(i) ~= 0) % Sinal sonoro
        gerado((80*(i-1))+1:round(8000/pitch(i)):80*i) = 1;
        gerado((80*(i-1))+1:80*i) = filter(Gs(i), aks(:,ceil(i/3)), gerado((80*(i-1))+1:80*i));
    else % Sinal surdo
        aleatorio = randn(80,1);
        gerado(80*(i-1)+1:80*i) = Gs(i)*aleatorio*sqrt(var(aleatorio));
    end
end

exc_sinal_cortado = randn(80,1);
saida_sinal_cortado = G*exc_sinal_cortado*sqrt(var(exc_sinal_cortado));

gerado = [saida_sinal_cortado; gerado];

%Gerando gráfico

figure(7);
plot(sinal);
hold on;
plot(gerado);
legend(["Sinal original", "Sinal de saída"]);
title("Comparação entre os sianis original e de saída");

%% G

% Quantizando os coeficientes

Ba = 7;
Bg = 5;

Q_aks = zeros(11,46);

for i = 1:46
    Q_aks(:,i) = quantize3(aks(:,i), Ba);
end

Q_Gs = quantize3(Gs, Bg);
Q_pitch = quantize3(pitch, Bg);

Q_gerado = zeros(138*80,1);

% Repetino a codificação

for i = 1:138
    if (pitch(i) ~= 0) % Sinal sonoro
        Q_gerado((80*(i-1))+1:round(8000/Q_pitch(i)):80*i) = 1;
        Q_gerado((80*(i-1))+1:80*i) = filter(Q_Gs(i), Q_aks(:,ceil(i/3)), Q_gerado((80*(i-1))+1:80*i));
    else % Sinal surdo
        aleatorio = randn(80,1);
        Q_gerado(80*(i-1)+1:80*i) = Q_Gs(i)*aleatorio*sqrt(var(aleatorio));
    end
end

exc_sinal_cortado = randn(80,1);
Q_saida_sinal_cortado = quantize3(G, Bg)*exc_sinal_cortado*sqrt(var(exc_sinal_cortado));

Q_gerado = [Q_saida_sinal_cortado; Q_gerado];

figure(8);
plot(sinal);
hold on;
plot(Q_gerado);
legend(["Sinal original", "Sinal de saída"]);
title("Comparação entre os sianis original e de saída com coeficientes quantizados");


%% RESULTADO
%Para ouvir use:
% sound(gerado, fs)
% sounda(Q_gerado, fs)

%gerando arquivos .wav
audiowrite("antarctice_Qgerado.wav", Q_gerado, fs);
audiowrite("antarctica_gerado.wav", gerado, fs);






