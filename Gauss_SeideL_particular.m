function [u,ITER]=Gauss_SeideL_particular(n, tol, kmax)
% Encontra a solu��o  relacionada com a distribui��o  de temperatura que se 
% em estado de Equil�brio em  uma placa quadrada, 
% para condi��es de contorno espec�ficas
% utiliza  iter��es de Gauss_SeideL_particula
% n � o n�mero de pontos de grade em cada dire��o ( sem incluir pontos limite )
% tol � a toler�ncia para o crit�rio de paragem
% kmax n�mero maximo de itera��es
%Na grade , usou-se i para linha (ascendemte),  j para coluna (� esquerda-diereita) 
% conjunto de valores limite, De temperatura � zero 
% no bordo superior e 100 nos restantes bordos

N=n+2;   % N � o n�mero de pontos de grade em cada dire��o (incluindo pontos limite)

u(1,1:N)=100;  % limite inferior
u(1:N,1)=100;  % limite esquerdo
u(1:N,N)=100;  % limite direito
u(N,1:N)=0;     % limite acima

w(1,1:N)=100;  % limite inferior
w(1:N,1)=100;  % limite esquerdo
w(1:N,N)=100;  % limite direito 
w(N,1:N)=0;     % limite acima

% aproxima��o inicial para pontos interiores
u(2:N-1,2:N-1)=50;

% Calcular a solu��o do estado est�vel(equil�brio)
DIFF=tol+1;      % qualquer valor maior do que tol serve
ITER=0;            % inicializa um contador de iteracoes

while DIFF>tol && ITER<kmax        % crit�rio de paragem 
    for i=2:N-1
        for j=2:N-1
            w(i,j)=(w(i-1,j)+w(i,j-1)+u(i,j+1)+u(i+1,j))/4;   % implementa��o da f�rmula (4.22)
        end                                                                   % ver na p�gina 29 da disserta��o
    end
    DIFF=max(max(abs(w-u)))/max(max(abs(w)));   % erro relativo 
    u=w;                  % armazena os valores da itera��o recente
    ITER=ITER+1;      % contador para o n�mero de itera��es 
end
   if DIFF<=tol
      disp(['Atingiu-se a precis�o desejada em ',num2str(ITER),' itera��es.']);
   else
      disp('O m�todo diverge ou converge mas kmax � insuficiente para permitir alcan�ar a precis�o desejada.');
   end
% EXIBE UM MAPA COLORIDO DE DISTRIBUI��ES
 u=flipud(u);% image(u)






