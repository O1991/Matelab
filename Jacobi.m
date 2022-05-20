function [x, iter]=Jacobi(A,b,tol,kmax)
%      Método iterativo de Jacobi
%      [x, iter]=Jacobi(A,b,tol,kmax) determina uma aproximação para
%      a solução do sistema Ax=b usando o método iterativo de
%      Jacobi tomando como aproximação inicial o vector nulo
%
%INPUTS:   A - matriz do sistema
%          b - vector coluna com n componentes
%          tol - tolerância para o erro
%          kmax - número máximo de iterações
%
%OUTPUTS:  x - solução (vector coluna)
%          iter - número de iterações


% Com A=D+L+U, D diagonal e L e U estritamente triangulares, inferior e
% superior respectivamente, temos a matriz de iteração de Jacobi
% BJ=-inv(D)*(L+U); o método converge SSE o raio espectral de BJ < 1.
% Para mais detalhes consulte por exemplo o livro do Prof. Heitor Pina,
% pag.343.


n=length(b);
x=zeros(n,1);   % aproximação inicial (vector nulo)
erro=tol+1;     % aqui qualquer valor maior do que tol serve
iter=0;         % inicializa um contador de iterações

while erro>tol & iter < kmax
     y=x;        % armazena os valores da iteraçao anterior
     for i=1:n
         h(i)=(b(i)-dot(A(i,:),y))/A(i,i); % implementação da fórmula
         x(i)=x(i)+h(i);                  % 8.2.4 no livro referido
     end
     erro=norm(h,inf);
     iter=iter+1;
     iter, x, pause
end

if erro<=tol
     disp(['Atingiu-se a precisão desejada em ',num2str(iter-1),'iterações']);
else
     disp(['O método diverge ou converge mas kmax é insuficiente para permitir alcançar a precisão desejada']);
end