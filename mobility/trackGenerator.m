%generates a set of files for the simpleTracking mobility model
nFiles = 10;
X = 2;
Y = 3;
R = 3;
D = 1;
mVel = 0.2;
nFrames = 15;
dir = 'mobData1';
M0 = -ones(X*Y+R);
for i=1:nFiles
    M = simpleTracking(X,Y,R,D,false,true,[dir,'/',num2str(i),'.mat'],M0,mVel,nFrames)
    M0(1:X*Y,1:X*Y) = M(1:X*Y,1:X*Y);
    M0 = M0-diag(diag(M0))+diag(diag(M))
end
