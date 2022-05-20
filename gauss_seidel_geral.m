function [x, k]=gauss_seidel_geral(A, b, tol,kmax)
  
  % Encontra a solu��o  relacionada com a distribui��o  de temperatura que se 
  % em estado de equilibrio em  uma placa quadrada, 
  % para condi��es de contorno espec�ficas
  % utiliza  iter��es de gauss_seidel_geral
  % tol � a toler�ncia para o crit�rio de paragem
  % kmax n�mero maximo de itera��es
  %  A matriz tridiagonal(unidimensional) ,ou uma matriz compostas por cinco diagonai(bidimensional)
  % b vetor dos termos independente 
  % x solu��o aproximada a solu��o axata do problema
  % K n�mero de itera��es para obter a solu��o aproximada
  % Na grade , usou-se i para linha (ascendemte),  j para coluna (� esquerda-diereita) 
  % Ver Algoritmo geral para m�todos iterativos na disserta��o p�gina 14

   n=length(b);             %  Comprimento do vetor independente
   x=zeros(n,1);            % aproxima��o inicial (vector nulo)
   erro=tol+1;               % qualquer valor maior do que tol serve
   k=0;                        % inicializa um contador de iteracoes
   h=zeros(n,1);            % pre-alocalocar para maior rapidez

   while erro>tol && k<kmax   % crit�rio de paragem 
     for i=1:n
       h(i)=(b(i)-dot(A(i,1:n),x))/A(i,i);  % implementa��o da f�rmula (3.15)
       x(i)=x(i)+h(i);                            % ver na p�gina 17 da disserta��o
     end
     erro=norm(h,inf)/norm(x,inf);  % erro relativo 
     k=k+1;                                % contador para o n�mero de itera��es 
   end

   if erro<=tol
      disp(['Atingiu-se a precis�o desejada em ',num2str(k),' itera��es.']);
   else
      disp('O m�todo diverge ou converge mas kmax � insuficiente para permitir alcan�ar a precis�o desejada.');
   end