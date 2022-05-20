function [T ,b1]=matriz_e_vetor_bidimensional(n)
  % 'n' é o número de pontos de grade em cada direção (sem incluir os pontos limite)
  % 'b1' vetor dos termos independente
  % 'T' Matriz composta por cinco diagonais
   
 
  alfa=100;              % conjunto de valores limite De temperatura é zero 
  delta=alfa;            % no bordo superior e 100 nos restantes bordos
  beta=alfa;
  super_gama=0;
  alfa(1:n)=alfa;               % limite no bordo inferior
  beta(1:n)=beta;             % limite no bordo  esquerdo
  gama(1:n)=super_gama;  % limite no bordo direito
  delta(1:n)=delta;            % limite no bordo superior 
  
  v1=[alfa(1)+delta(1) delta(2:n-1) beta(1)+delta(n)];   %vetor em que as suas componete são os valores 
                                                                              % do bordo inferior 
  m=numel(2:n-1);
  v=[alfa(2), zeros(1,n-2), beta(2)]; % vetor em que a primeira e a ultima componete é 100 e as restantes são zeros
  p=ones(1,m);               % vetor unitário 
  v2=kron(p,v);               % vetor em que as suas componete são as n's linhas internas da malha
 
 
  v3=[gama(1)+delta(n), gama(2:n-1), gama(n)+beta(n)]; % vetor em que as suas componentes são os valores
                                                                                 % do bordo superior  
                                                        
  fr=[v1,v2, v3];                % vetor dos valores das temperatura nas arestas
  b1=-fr;                          %  vetor dos termos independentes do problema
  
  id=eye(n);                       % matriz identidade
  d=-2*ones(1,n);               % vetor com elementos todos igual a 2
  c=ones(1,n-1);                 % vetor com elementos todos igual a 1
  A=diag(d)+diag(c,1)+diag(c,-1);   % matriz tridiagonal, com diagonal principal 2 
                                                  % e duas  diagonais secundarias com 1
                                                  
  T=kron(id,A)+kron(A,id);     % Matriz composta por cinco diagonais onde a diagonal principal e igual a -4, 
                                         %duais diagonais secundarias igual a 1  e outras duas composta por zero e 1. 
  end