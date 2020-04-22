nt = randi(4)+1;
nr = nt + randi(5);

ZTr = diag(rand(nt,1));
ZTi = rand(nt)-0.5; ZTi = ZTi+ZTi.';
ZT = ZTr + (1i)*ZTi;

ZRr = diag(rand(nr,1));
ZRi = rand(nr)-0.5; ZRi = ZRi+ZRi.';
ZR = ZRr + (1i)*ZRi;

%%%%%%%%%%%%%%%%%%%%%%%%%%

M_i = rand(nr,nt)-0.5;
M = (1i)*M_i;

Z = [ZT, M.'; M, ZR];

VTrefr = rand(nt,1)-0.5;
VTrefi = rand(nt,1)-0.5;

VTref = VTrefr + (1i)*VTrefi;

Vref = [VTref; zeros(nr,1)];

I = Z\Vref;

IT = I(1:nt);
IR = I(nt+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%
if (nt<=nr)
	disp('Teste do sistema superdeterminado com nr maior que nt');
	M*IT+ZR*IR
	ZT*IT+M.'*IR-VTref

	M(1:nt,1:nt)*IT+ZR(1:nt,1:nr)*IR

	IT+M(1:nt,1:nt)\(ZR(1:nt,1:nr)*IR)
end
%%%%%%%%%%%%%%%%%%%%%%%%%
if (nt==nr)
	disp('Teste da formula de tensao imaginaria');
	(-ZT*(M\ZR)+M.')*IR-VTref
	real(-ZT*(M\ZR)+M.')*imag(IR)+imag(-ZT*(M\ZR)+M.')*real(IR)-imag(VTref)
	
	real(-(ZTr+(1i)*ZTi)*(eye(nt)/M)*(ZRr+(1i)*ZRi))-real(-ZT*(M\ZR)+M.')
	
	-ZTr*(eye(nt)/M)*((1i)*ZRi)...
	-((1i)*ZTi)*(eye(nt)/M)*ZRr...
	-real(-ZT*(M\ZR)+M.')
	
	imag(-(ZTr+(1i)*ZTi)*(eye(nt)/M)*(ZRr+(1i)*ZRi)+M.')-imag(-ZT*(M\ZR)+M.')
	
	imag(-ZTr*(eye(nt)/M)*ZRr)...
	+imag(ZTi*(eye(nt)/M)*ZRi)...
	+M_i.'...
	-imag(-ZT*(M\ZR)+M.')
	
	(1i)*ZTr*(eye(nt)/M)*ZRr...
	-(1i)*ZTi*(eye(nt)/M)*ZRi...
	-(1i)*M.'...
	-imag(-ZT*(M\ZR)+M.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nt==nr)
	disp('Teste da formula de geração da matrix de impedância');
	x = rand(nr,1)-0.5;
	y = rand(nr,1)-0.5;
	
	izt = (-M.'*y*pinv(x)-ZTr*(eye(nt)/M)*ZRi+ZTr*(eye(nt)/M)*ZRr*y*pinv(x))/((eye(nt)/M)*(ZRr+ZRi*y*pinv(x)));
	
	imag((-(ZTr+(1i)*izt)*(M\ZR)+M.')*((1i)*x+y))
	
	W = 10*(1i)*(rand(nt)-0.5);
	L = -(1i)*ZTr*(eye(nt)/M)*ZRi-(1i)*izt*(eye(nt)/M)*ZRr;
	S = (1i)*ZTr*(eye(nt)/M)*ZRr-(1i)*izt*(eye(nt)/M)*ZRi-(1i)*M.';
	
	(L+W*(eye(nt)-x*pinv(x)))*x+S*y
	
	izt = (-M.'*y*pinv(x) + ZTr*(eye(nt)/M)*ZRr*y*pinv(x)+W*(eye(nt)-x*pinv(x))...
			-ZTr*(eye(nt)/M)*ZRi...
			)/((eye(nt)/M)*(ZRr+ZRi*y*pinv(x)));
			
	L = -(1i)*ZTr*(eye(nt)/M)*ZRi-(1i)*izt*(eye(nt)/M)*ZRr;
	S = (1i)*ZTr*(eye(nt)/M)*ZRr-(1i)*izt*(eye(nt)/M)*ZRi-(1i)*M.';
	
	L*x+S*y
	(L+W*(eye(nt)-x*pinv(x)))*x+S*y
	
	L-real(-(ZTr+(1i)*izt)*(M\ZR)+M.')
	
	S-imag(-(ZTr+(1i)*izt)*(M\ZR)+M.')
	
	imag((-(ZTr+(1i)*izt)*(M\ZR)+M.')*((1i)*x+y))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nt<nr)
	disp('Teste da formula de geração da matrix de impedância com nr maior que nt');
	
	iM = (M.'*M)\(M.');
	
	M*IT+ZR*IR
	ZT*IT+M.'*IR-VTref
	
	IT+iM*ZR*IR
	(M.'-ZT*iM*ZR)*IR-VTref
	
	(1i)*ZTr*iM*ZRr...
	-(1i)*ZTi*iM*ZRi...
	-(1i)*M.'...
	- imag(M.'-ZT*iM*ZR)
	
	-ZTr*iM*((1i)*ZRi)...
	-((1i)*ZTi)*iM*ZRr...
	-real(M.'-ZT*iM*ZR)
	
	%gerando um M válido para x e y
	
	x = rand(nr,1)-0.5;
	y = rand(nr,1)-0.5;
	
	%M*it = -ZR*ir
	%M(1:nt,:)*it = -ZR(1:nt,:)*((1i)*x+y)
	%it = -(eye(nt)/M(1:nt,:))*ZR(1:nt,:)*((1i)*x+y)
	m1 = (1i)*(rand(nt,nt)-0.5);
	ir = (1i)*x+y;
	it = -(eye(nt)/m1)*ZR(1:nt,:)*ir;
	
	%M = [m1;m2]
	%m2*it =-ZR(nt+1:end,:)*((1i)*x+y)
	m2 = zeros(nr-nt,nt);
	for i = 1:nr-nt
		A = [real(it).';imag(it).']
		b = [-(1i)*imag(ZR(nt+i,:)*ir);(1i)*real(ZR(nt+i,:)*ir)]
		m = pinv(A)*b;
		m.'*it+ZR(nt+i,:)*ir
		m2(i,:) = m.';
	end
	
	M = [m1;m2];
	
	M*it + ZR*ir
	
	M.'*M*it + M.'*ZR*ir
	
	it + (M.'*M)\(M.')*ZR*ir
	
	iM = (M.'*M)\(M.');
	
	it+iM*ZR*ir
	
	izt = rand(nt)-0.5;
	v = (ZTr+(1i)*izt)*it+(M.')*ir;
	
	(ZTr+(1i)*izt)*it+(M.')*ir - v
	
	izt*real(it)+(ZTr*imag(it)+imag((M.')*ir - v))
	
	izt = -(ZTr*imag(it)+imag((M.')*ir - v))*pinv(real(it));
	
	izt*real(it)+(ZTr*imag(it)+imag((M.')*ir - v))
	%%%%%%%%%%%
	
	M*it + ZR*ir
	
	%imag(ZT) = -(1i)*MT
	%real(ZT) = diag(RT)
	RT = -ones(nt,1);
	MT = (1i)*(rand(nt)-0.5);
	%enquanto tiver resistencia negativa
	while true
		
		v = (1i)*MT*real(it)-imag(M.'*ir);
		
		RT = v./imag(it);
		
		if max(RT<=0)==0
			break;
		end
		
		[~,k] = min(RT);
		
		dm = (1i)*v(k)/real(it(k)) - (1i)*sign(imag(it(k)))*sign(real(it(k)))*gamrnd(1,2);
		
		MT(k,k) = MT(k,k) + dm;
	end
	ZT = diag(RT) + MT;
	imag(ZT*it+M.'*ir)
end

