clear;
close all;

load desejado2.mat
load entrada2.mat
%%

siz = length(u2);
N = siz;

% R = xcorr(u2, 1, 'biased'); % acima de 2 coeficientes já fica muito abaixo de 1.
% R = [R(2), R(1); R(1), R(2)];
% 
% 
% p = xcorr(d2, u2, 'biased'); % verificaão matriz p
% p = [p(siz); p(siz+1)];
% % 
%  wO = R\p;
% 
% J_min = var(d2) - transpose(wO)*p;
% % 
%  lambda = eig(R);
%  
%  mi_max = 2/max(lambda);

[W2,erro]= lms(u2, d2, 200, N, 5e-3);

pot_erro = mag2db(var(erro)); %% n sei se está certo isso aqui

y = d2 - transpose(erro);


ERLE = erle(d2, erro, 1024, 8000);


figure(2); 
plot(W2(:,1)); 
hold on 
plot(W2(:,2)) 
legend('w1', 'w2')

%W2 = [W2(end,1) ; W2(end,2)]
%ERLE = 10*log10(ˆ2)


freqz(W2(end, :), 1)