function [W,erro]=lms(x, d, M, N, mu)

% [W,erro] = lms(x, d, M, N, mu)
% Algoritmo LMS 
%
% Par?metros de sa?da
% W -> matriz dos coeficientes do filtro transversal
% erro -> sinal de erro (d-y)
%
% Par?metros de entrada
% x -> sequ?ncia de observa??o
% d -> sequ?ncia desejada
% M -> n?mero de coeficientes do filtro adaptativo
% N -> n?mero de itera??es
% mu -> passo de adapta??o

% MTMS 14/10/1999

X=zeros(M,1);
W=zeros(N+1,M);
erro = zeros(1,N);
y = zeros(1,N);
for n=1:N
   X=[x(n);X(1:M-1)];
   y=W(n,:)*X;
   erro(n)=d(n)-y(n); 
   W(n+1,:)=W(n,:)+mu*erro(n)*X';
end;