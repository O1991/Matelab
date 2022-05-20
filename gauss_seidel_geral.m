function [x, k]=gauss_seidel_geral(A, b, tol,kmax)
  
  % Encontra a solução  relacionada com a distribuição  de temperatura que se 
  % em estado de equilibrio em  uma placa quadrada, 
  % para condições de contorno específicas
  % utiliza  iterções de gauss_seidel_geral
  % tol é a tolerância para o critério de paragem
  % kmax número maximo de iterações
  %  A matriz tridiagonal(unidimensional) ,ou uma matriz compostas por cinco diagonai(bidimensional)
  % b vetor dos termos independente 
  % x solução aproximada a solução axata do problema
  % K número de iterações para obter a solução aproximada
  % Na grade , usou-se i para linha (ascendemte),  j para coluna (á esquerda-diereita) 
  % Ver Algoritmo geral para métodos iterativos na dissertação página 14

   n=length(b);             %  Comprimento do vetor independente
   x=zeros(n,1);            % aproximação inicial (vector nulo)
   erro=tol+1;               % qualquer valor maior do que tol serve
   k=0;                        % inicializa um contador de iteracoes
   h=zeros(n,1);            % pre-alocalocar para maior rapidez

   while erro>tol && k<kmax   % critério de paragem 
     for i=1:n
       h(i)=(b(i)-dot(A(i,1:n),x))/A(i,i);  % implementação da fórmula (3.15)
       x(i)=x(i)+h(i);                            % ver na página 17 da dissertação
     end
     erro=norm(h,inf)/norm(x,inf);  % erro relativo 
     k=k+1;                                % contador para o número de iterações 
   end

   if erro<=tol
      disp(['Atingiu-se a precisão desejada em ',num2str(k),' iterações.']);
   else
      disp('O método diverge ou converge mas kmax é insuficiente para permitir alcançar a precisão desejada.');
   end