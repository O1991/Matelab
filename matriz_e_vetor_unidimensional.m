function [A ,b]=matriz_e_vetor_unidimensional(n)
  % 'n' é o número de pontos da reta ou eixo  (incluindo pontos no limite)
  % 'b' vetor dos termos independente
  % 'A' Matriz tridiagonal
  
  alfa=100;  % valor definido no primeiro  ponto limite (valor da temperatura)
  beta=alfa; % valor definido no ultimo  ponto limite
  d=-2*ones(1,n);  % vetor com elementos todos 2
  c=ones(1,n-1);    % vetor com elementos todos 2
  A=diag(d)+diag(c,1)+diag(c,-1);  % matriz tridiagonal, com diagonal principal 2 
                                                  % e duas  diagonais secundarias com 1
                                                  
  b(1)=-alfa;    % primeira componete do vetor dos termos independente
  b(2:n-1)=0;   % Segunda até a penúltima componete do vetor dos termos independente
  b(n)=-beta;   % ultima componete do vetor dos termos independente                
  end