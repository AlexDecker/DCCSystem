%generates a set of files for the simpleTracking mobility model
calcInd = true;
plotAnim = false;
nFiles = 10;
X = 2;
Y = 3;
R = 3;
D = 0.4;
mVel = 0.1;
nFrames = 15;
dir = 'mobData2';
M0 = -ones(X*Y+R);
for i=1:nFiles
    M = simpleTracking(X,Y,R,D,plotAnim,calcInd,[dir,'/',num2str(i),'.mat'],M0,mVel,nFrames)
	if calcInd
    	M0(1:X*Y,1:X*Y) = M(1:X*Y,1:X*Y);
    	M0 = M0-diag(diag(M0))+diag(diag(M));
	end
end
