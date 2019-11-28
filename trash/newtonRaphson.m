%estrutura das equacoes: x'Ax+b'x+c = 0

%erro admissivel
ea = 1e-10;

%maximo numero de iteracoes
mi = 1000;

%numero de equacoes
e = 10;
%numero de variaveis
v = 10;

%gerando aleatoriamente os coeficientes
A = cell(0);
b = cell(0);
c = cell(0);
for i=1:e
    A{end+1} = 5*(rand(v)-0.5);
    b{end+1} = 5*(rand(v,1)-0.5);
    c{end+1} = 5*(rand-0.5);
end

%solucao inicial
x = ones(v,1);
err0 = inf;%erro de referencia
it = 0;%numero de iteracoes

tic
while true

    %testando uma solucao
    f = zeros(e,1);
    for i=1:e
        f(i) = x.'*A{i}*x+b{i}.'*x+c{i};
    end

    err = sum(abs(f))/e;
    if err>=err0 || it>mi
        x = NaN*x;%nao convergiu
        break;
    end

    if err<=ea %condicao de parada: erro muito baixo
        break;%convergiu
    end

    %calculando a proxima solucao--------

    %Jacobiano
    J = zeros(v,e);
    for i=1:e
        %gradiente
        g = (A{i}.'+A{i})*x+b{i};
        %compondo o jacobiano
        J(i,:) = g.';
    end

    %proxima solucao
    x = x - J\f;

    it=it+1;
end
toc
disp(err);
disp(x);
