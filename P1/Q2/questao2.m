%% Questão 2

%% 1
close all;

[sinal, fs] = audioread('antarctica.wav');

sinal_2x = sinal(1:2:end);


%%O áudio fica demasiado rápido e impossível de ser entendido.

%% 2



%clear;


% Definições
L = 240;
N = L;
Q = 512;
K = 2;
gama = 1;

vel = .8;
%[sinal, fs] = audioread('antarctica.wav'); 

%% 1
% Cortando sinal em trechos de 240 amostras (30 ms)
% trechos = zeros(L,46);
% 
% for i = 1:46
%         %trechos(:,46) = sinal(240*(i-1)+1:240*i);
% end

%% 2
% Base com 512 funções aleatórias
fnc_base = randn(N, Q);

%% 3
% Condição inicial do filtro de trato vocal
zs = zeros(1, 10); 

%% 4

gerado = [];

for i = 1:40
    % A
    % Determinando os coeficientes LPC do quadro atual
    trecho = sinal(L*(ceil(i)-1)+1:L*ceil(i));
    [aq, sig]= lpc(trecho.*hamming(L), 10);
    
    subquadro = trecho(1:L);

    % C
    % Filtrando as 512 sequ"Encias da base pelo filtro
    fnc_base_filt = zeros(N, Q);
    fnc_base_filt = filter(1, aq, fnc_base);

    % D
    [y0, zs] = filter(1, aq, zeros(N,1));

    % E
    e0 = subquadro - y0;

    % F
    [ganhos, indices] = find_Nbest_components(e0, fnc_base_filt, K);
    
    % G
    d = fnc_base(:,indices)*ganhos;
    
    % H

    if( vel == 0.8)
        [gerado_trecho, zs] = filter(1, aq, [d(1:2:96);d(1:end)], zs);   
    else
        [gerado_trecho, zs] = filter(1, aq, d(120:end), zs); 
    end

    
    gerado = [gerado; gerado_trecho];
    

end

figure(1);
plot(sinal);
hold on
plot(gerado);
legend("Sinal original", "CELP");
title("Comparação do sinal original com o reconstruído (codificador CELP)");

audiowrite("antarctida_CELP.wav", gerado, 8000);


