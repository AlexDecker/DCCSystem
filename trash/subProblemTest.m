%este script testa as formulas do subproblema. No caso, roda com diferentes instancias
%e para so quando encontra um erro
clear all;

rng('shuffle');

err = 1e-9;%tolerancia

timeList = [];
convergiu = 0;
nconvergiu = 0;

while true
    nr = round(4*rand)+1;
    nt = round(4*rand)+1;

    R = 3*diag(rand(nt+nr,1));%resitencia
    
    MT = rand(nt);MR = rand(nr);%acoplamento TX-TX e RX-RX
    
    ZT = R(1:nt,1:nt)-1i*(MT+MT'-2*diag(diag(MT)))/2;%impedancia TX
    ZR = R(nt+1:end,nt+1:end)-1i*(MR+MR'-2*diag(diag(MR)))/2;%impedancia RX
    M  = -1i*rand(nr,nt);%acoplamento TX-RX

    Z = [ZT, M.';
         M, ZR];

    H = -ZR\M;%canal TX/RX
    ZTT = ZT+M.'*H;

    %testando a formula para abs(iR(k))^2
    iT = 2*(rand(nt,1)-0.5)+(2i)*(rand(nt,1)-0.5);%um vetor de conrrente de transmissao qualquer
    iR = H*iT;%calculando a corrente de recebimento correspondente
    x = [real(iT);imag(iT)];
    for k=1:nr
        Hk = H(k,:)'*H(k,:);
        iRk = x.'*[real(Hk), -imag(Hk); imag(Hk), real(Hk)]*x;
        if abs(iRk-abs(iR(k))^2)/max(1,abs(iR(k))^2)>err
            error('corrente alvo');
        end
    end
    
    %testando a geracao das correntes de transmissao validas (para obter tensoes reais)
    L = [imag(ZTT),real(ZTT)];
    x = 10*(eye(2*nt)-pinv(L)*L)*(rand(2*nt,1)-0.5);
    iT = x(1:nt)+(1i)*x(nt+1:end);
    vT = ZTT*iT;
    if abs(imag(vT))/abs(vT)>err
        error('tensao real');
    end

    %testando a formula da potencia ativa (utilizando o iT do item anterior)
    p = x.'*[real(ZTT), -imag(ZTT);zeros(nt), zeros(nt)]*x;
    ref = real(iT'*ZTT*iT);
    if abs(p-ref)/max(1,abs(ref))>err
        error('potencia ativa');
    end

    %gerando uma instancia que com certeza tem pelo menos solucao x
    miT2 = real(iT).^2+imag(iT).^2;%maxima amplitude para a corrente de transmissao
    %p eh a potencia ativa maxima
    miR2 = abs(H*iT).^2;%corrente alvo
    
    disp('Testes preliminares executados com sucesso');
    
    %Usando a Newton-Raphson com apenas igualdades---------------------------------
    
    
    ttl1 = 100;%solucoes iniciais distintas
    tic
    while(true)
        x1 = rand(2*nt,1);%solucao inicial
        ttl2 = 100;%iteracoes restantes
        err0 = inf;
        while true
            %testando uma solucao
            f = zeros(nr+2*nt+1,1);
            %restricoes de corrente alvo
            for k=1:nr
                Hk = H(k,:)'*H(k,:);
                f(k) = x1.'*[real(Hk), -imag(Hk); imag(Hk), real(Hk)]*x1 - miR2(k);
            end
            %restricoes de tensao real
            f(nr+1:nr+nt) = [imag(ZTT),real(ZTT)]*x1;
            %restricoes de limite de corrente de transmissao
            f(nr+nt+1:nr+2*nt) = x1(1:nt).^2+x1(nt+1:end).^2-abs(iT).^2;
            %restricao de potencia
            f(end) = x1.'*[real(ZTT), -imag(ZTT); zeros(nt), zeros(nt)]*x1-p;

            err1 = mean(abs(f));
            if err1>=err0 || ttl2==1
                x = NaN*x;%solucao nao encontrada
                break;
            end

            if err1<=err
                break;%solucao encontrada
            end

            %calculando a proxima solucao------------

            %Jacobiano
            J = zeros(nr+2*nt+1,2*nt);
            %restricoes de corrente alvo
            for k=1:nr
                Hk = H(k,:)'*H(k,:);
                A = [real(Hk), -imag(Hk); imag(Hk), real(Hk)];
                J(k,:) = x1.'*(A+A.');
            end
            %restricoes de tensao real
            J(nr+1:nr+nt,:) = [imag(ZTT),real(ZTT)];
            %restricoes de limite de corrente de transmissao
            J(nr+nt+1:nr+2*nt,:) = 2*[diag(x1(1:nt)),diag(x1(nt+1:end))];
            %restricao de potencia
            J(end,:) = x1.'*([2*real(ZTT), -imag(ZTT); -imag(ZTT).', zeros(nt)]).';

            %proxima solucao
            x1 = x1 - J\f;

            ttl2 = ttl2-1;
        end
        if err1<=err
            convergiu = convergiu+1;
            break;
        end
        if ttl1==0 %acabaram as chances
            nconvergiu = nconvergiu+1;
            break;
        end
        ttl1 = ttl1-1;
    end
    timeList = [timeList,toc];
    toc
    disp('Instancia do sistema de equacoes executada com sucesso');

    tic
    toc
    disp('Instancia do sistema de inequacoes executada com sucesso');
end
