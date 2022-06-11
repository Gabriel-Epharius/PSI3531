clear;
close all;

load desejado1.mat
load entrada1.mat
%%

N = 1e4;

R = xcorr(u1, 1, 'biased'); % acima de 2 coeficientes já fica muito abaixo de 1.
R = [R(2), R(1); R(1), R(2)];


p = xcorr(d1, u1, 'biased'); % verificaão matriz p
p = [p(1e4); p(1e4+1)];
% 
 wO = R\p;

J_min = var(d1) - transpose(wO)*p;
% 
 lambda = eig(R);
 
 mi_max = 2/max(lambda);
 

[W,erro]= lms(u1, d1, 200, N, 0.0005);
 
figure(1); 
plot(W(:,1)); 
hold on 
plot(W(:,2)) 
legend('w1', 'w2')

pot_erro = mag2db(var(erro));


W1 = [W(end,1) ; W(end,2)]

freqz(W(end, :), 1)

% 
%não está convergindo.
