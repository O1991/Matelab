function [u,ITER]=Gauss_SeideL_particular(n, tol, kmax)
% Encontra a solução  relacionada com a distribuição  de temperatura que se 
% em estado de Equilíbrio em  uma placa quadrada, 
% para condições de contorno específicas
% utiliza  iterções de Gauss_SeideL_particula
% n é o número de pontos de grade em cada direção ( sem incluir pontos limite )
% tol é a tolerância para o critério de paragem
% kmax número maximo de iterações
%Na grade , usou-se i para linha (ascendemte),  j para coluna (á esquerda-diereita) 
% conjunto de valores limite, De temperatura é zero 
% no bordo superior e 100 nos restantes bordos

N=n+2;   % N é o número de pontos de grade em cada direção (incluindo pontos limite)

u(1,1:N)=100;  % limite inferior
u(1:N,1)=100;  % limite esquerdo
u(1:N,N)=100;  % limite direito
u(N,1:N)=0;     % limite acima

w(1,1:N)=100;  % limite inferior
w(1:N,1)=100;  % limite esquerdo
w(1:N,N)=100;  % limite direito 
w(N,1:N)=0;     % limite acima

% aproximação inicial para pontos interiores
u(2:N-1,2:N-1)=50;

% Calcular a solução do estado estável(equilíbrio)
DIFF=tol+1;      % qualquer valor maior do que tol serve
ITER=0;            % inicializa um contador de iteracoes

while DIFF>tol && ITER<kmax        % critério de paragem 
    for i=2:N-1
        for j=2:N-1
            w(i,j)=(w(i-1,j)+w(i,j-1)+u(i,j+1)+u(i+1,j))/4;   % implementação da fórmula (4.22)
        end                                                                   % ver na página 29 da dissertação
    end
    DIFF=max(max(abs(w-u)))/max(max(abs(w)));   % erro relativo 
    u=w;                  % armazena os valores da iteração recente
    ITER=ITER+1;      % contador para o número de iterações 
end
   if DIFF<=tol
      disp(['Atingiu-se a precisão desejada em ',num2str(ITER),' iterações.']);
   else
      disp('O método diverge ou converge mas kmax é insuficiente para permitir alcançar a precisão desejada.');
   end
% EXIBE UM MAPA COLORIDO DE DISTRIBUIÇÕES
 u=flipud(u);% image(u)






