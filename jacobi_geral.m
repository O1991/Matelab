function [x, k]=jacobi_geral(A,b,tol,kmax)
 % Encontra a solução  relacionada com a distribuição  de temperatura que se 
 % em estado de equilibrio em  uma placa quadrada, 
 % para condições de contorno específicas
  % utiliza  iterções de jacobi_geral
  % tol é a tolerância para o critério de paragem
  %  A matriz tridiagonal(unidimensional) ,ou uma matriz compostas por cinco diagonai(bidimensional)
  % b vetor dos termos independente
  % kmax número maximo de iterações 
  % x solução aproximada a solução axata do problema
  % K número de iterações para obter a solução aproximada
   % Na grade , usou-se i para linha (ascendemte),  j para coluna (á esquerda-diereita) 
  % Ver Algoritmo geral para métodos iterativos na dissertação página 14
  
  n=length(b);     % Comprimento do vetor independente
                         % 'n' é o número de pontos de grade em cada direção 
                         % (incluindo limite pontos) ou de uma recta ou eixo.
 
  x=zeros(n,1);   % aproximação inicial (vector nulo)
  erro=tol+1;      % aqui qualquer valor maior do que tol serve
  k=0;                % inicializa um contador de iterações
  h=zeros(n,1);   % pre-alocalocar para maior rapidez
  
while erro>tol && k < kmax  % critério de paragem
     y=x;                                           % armazena os valores da iteração anterior
     for i=1:n
         h(i)=(b(i)-dot(A(i,:),y))/A(i,i);  % implementação da fórmula (3.7)
         x(i)=x(i)+h(i);                         %  ver na página 16, da dissertação
     end
     erro=norm(h,inf)/norm(x,inf);         % erro relativo 
     k=k+1;                                        % contador para o número de iterações 
     
     %iter, u, pause
end
   if erro<=tol
      disp(['Atingiu-se a precisão desejada em ',num2str(k),' iterações.']);
   else
      disp('O método diverge ou converge mas kmax é insuficiente para permitir alcançar a precisão desejada.');
   end