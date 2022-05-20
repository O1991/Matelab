function [x, k]=jacobi_geral(A,b,tol,kmax)
 % Encontra a solu��o  relacionada com a distribui��o  de temperatura que se 
 % em estado de equilibrio em  uma placa quadrada, 
 % para condi��es de contorno espec�ficas
  % utiliza  iter��es de jacobi_geral
  % tol � a toler�ncia para o crit�rio de paragem
  %  A matriz tridiagonal(unidimensional) ,ou uma matriz compostas por cinco diagonai(bidimensional)
  % b vetor dos termos independente
  % kmax n�mero maximo de itera��es 
  % x solu��o aproximada a solu��o axata do problema
  % K n�mero de itera��es para obter a solu��o aproximada
   % Na grade , usou-se i para linha (ascendemte),  j para coluna (� esquerda-diereita) 
  % Ver Algoritmo geral para m�todos iterativos na disserta��o p�gina 14
  
  n=length(b);     % Comprimento do vetor independente
                         % 'n' � o n�mero de pontos de grade em cada dire��o 
                         % (incluindo limite pontos) ou de uma recta ou eixo.
 
  x=zeros(n,1);   % aproxima��o inicial (vector nulo)
  erro=tol+1;      % aqui qualquer valor maior do que tol serve
  k=0;                % inicializa um contador de itera��es
  h=zeros(n,1);   % pre-alocalocar para maior rapidez
  
while erro>tol && k < kmax  % crit�rio de paragem
     y=x;                                           % armazena os valores da itera��o anterior
     for i=1:n
         h(i)=(b(i)-dot(A(i,:),y))/A(i,i);  % implementa��o da f�rmula (3.7)
         x(i)=x(i)+h(i);                         %  ver na p�gina 16, da disserta��o
     end
     erro=norm(h,inf)/norm(x,inf);         % erro relativo 
     k=k+1;                                        % contador para o n�mero de itera��es 
     
     %iter, u, pause
end
   if erro<=tol
      disp(['Atingiu-se a precis�o desejada em ',num2str(k),' itera��es.']);
   else
      disp('O m�todo diverge ou converge mas kmax � insuficiente para permitir alcan�ar a precis�o desejada.');
   end