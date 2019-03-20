%creates a set of coils disposed along a XxY grid for the transmitting setup
%and R receiving independent coils that perform random walk.
%-plotAnimation = boolean, show the animation of the coils in 3D?
%-evalMutualCoupling = boolean, calculate the interactions between the coils and
%save it in a file?
%-D = acceptable box side length (considering the center of the TX setup)
%-file = output filename
%-M0 = preprocessed mutual coupling matrix with no magnetic permeability constant.
%use -1 for certain position to be calculated when evalMutualCoupling==true.i
%-mVel = maximum velocity (meters per frame for each direction)
%-nFrames = number of frames for this simulation

function M = simpleTracking(X, Y, R, D, plotAnimation, evalMutualCoupling,...
    file,M0,mVel,nFrames)
    
    M = [];

    w = 1e+5;%(dummie)
    mi = pi*4e-7;%(dummie)

    %Dimensions (transmitting coils)
    R2_tx = 0.1262;%external radius, in order to the total area to be 0.05m2
    N_tx = 17;%number of turns
    wire_radius_tx = 0.0015875;%wire thickness(m) diam = 1/8''
    R1_tx = R2_tx-4*N_tx*wire_radius_tx;%internal radius
    
    %Dimensions (receiving coils)
    R1_rx = 0.001368;%inner radius, tunned
    N_rx = 21.9649;%number of turns, tunned
    wire_radius_rx = 0.00079375;%wire radius (m) diam = 1/16''
    R2_rx = R1_rx+2*N_rx*wire_radius_rx;%external radius
    A_rx=0.011272;B_rx=0.00068937;%inner rectangle dimensions, tunned

    pts_tx = 750;%resolution of each coils (number of points)
    pts_rx = 750;

    stx = 0.04;%espacing between transmitters (acording with the illustration
    %of the paper. For generating an area of 0.3822m2, it must be 0.0
    
    coilPrototypeRX = QiRXCoil(R1_rx,R2_rx,N_rx,A_rx,B_rx,wire_radius_rx,pts_rx);
    coilPrototypeTX = SpiralPlanarCoil(R2_tx,R1_tx,N_tx,wire_radius_tx,pts_tx);
    
    %centroid of the tx setup (used for positioning the box)
    centerX = 0;
    centerY = 0;
    txGroups = [];
    for x=0:X-1
        for y=0:Y-1
            group = struct('coils', struct('obj',translateCoil(...
                            coilPrototypeTX,x*(2*R2_tx+stx),y*(2*R2_tx+stx),0)),...
                    'R', -1, 'C', -1);
            centerX = centerX + x*(R2_tx+stx);
            centerY = centerY + y*(R2_tx+stx);
            txGroups = [txGroups, group];
        end
    end
    centerX = centerX/(X*Y);
    centerY = centerY/(X*Y);
    %endpoints of the box
    X0 = centerX-D/2; X1 = centerX+D/2;
    Y0 = centerY-D/2; Y1 = centerY+D/2;
    
    %the initial pose for the receivers
    rxGroups = [];
    for r = 1:R
        %getting valid coordinates
        while(true)
            x = rand*(X1-X0)+X0;
            y = rand*(Y1-Y0)+Y0;
            z = rand*D;
            %testing for conflicts
            conflict = false;
            for i=1:length(txGroups)
                zVar = abs(txGroups(i).coils.obj.Z-z);
                xyDist = sqrt((txGroups(i).coils.obj.X-x)^2+(txGroups(i).coils.obj.Y-y)^2);
                if((zVar<=wire_radius_tx+wire_radius_rx)&&(xyDist<=R2_tx+R2_rx+max(A_rx,B_rx)/2))
                    conflict = true;
                    break;
                end
            end
            if ~conflict
                for i=1:length(rxGroups)
                    zVar = abs(rxGroups(i).coils.obj.Z-z);
                    xyDist = sqrt((rxGroups(i).coils.obj.X-x)^2+(rxGroups(i).coils.obj.Y-y)^2);
                    if((zVar<=2*wire_radius_rx)&&(xyDist<=2*R2_rx+max(A_rx,B_rx)))
                        conflict = true;
                        break;
                    end
                end
                if ~conflict
                    break;
                end
            end
        end 
        group = struct('coils', struct('obj',translateCoil(coilPrototypeRX,x,y,z)),...
                    'R', -1, 'C', -1);
        rxGroups = [rxGroups, group];
    end
    groupList = [txGroups, rxGroups];
    environment = Environment(groupList,w,mi);
    envList = [environment];
    %random walk
    for n = 2:nFrames
        newRxGroups = [];%the rx coils after being iterated
        for j=1:length(rxGroups)
            %getting valid coordinates
            while(true)
                %new possible displacement values
                dx = rand*mVel-mVel/2;
                dy = rand*mVel-mVel/2;
                dz = rand*mVel-mVel/2;
                %testing for conflicts
                conflict = false;
                for i=1:length(txGroups)
                    zVar = abs(txGroups(i).coils.obj.Z-rxGroups(j).coils.obj.Z+dz);
                    xyDist = sqrt((txGroups(i).coils.obj.X-rxGroups(j).coils.obj.X+dx)^2+...
                             (txGroups(i).coils.obj.Y-rxGroups(j).coils.obj.Y+dy)^2);
                    if((zVar<=wire_radius_tx+wire_radius_rx)&&...
                        (xyDist<=R2_tx+R2_rx+max(A_rx,B_rx)/2))
                        conflict = true;
                        break;
                    end
                end
                if ~conflict
                    for i=1:j-1 %verify only the already updated coils
                        zVar = abs(rxGroups(i).coils.obj.Z-rxGroups(j).coils.obj.Z+dz);
                        xyDist = sqrt((rxGroups(i).coils.obj.X-rxGroups(j).coils.obj.X+dx)^2+...
                                (rxGroups(i).coils.obj.Y-rxGroups(j).coils.obj.Y+dy)^2);
                        if((zVar<=2*wire_radius_rx)&&...
                            xyDist<=2*R2_rx+max(A_rx,B_rx))
                            conflict = true;
                            break;
                        end
                    end
                    if ~conflict
                        break;
                    end
                end
            end
            group = struct('coils', struct(...
                    'obj',translateCoil(rxGroups(j).coils.obj,-dx,-dy,-dz)),...
                    'R', -1, 'C', -1);
            newRxGroups = [newRxGroups, group];
        end
        rxGroups = newRxGroups;
        groupList = [txGroups, rxGroups];
        environment = Environment(groupList,w,mi);
        envList = [envList,environment];
    end

    ok = true;
    for i=1:length(envList)
        ok = ok && check(envList(i));
    end

    if(ok)
        if evalMutualCoupling
            disp('Starting coupling evaluation...');
            envList(1) = evalM(envList(1),M0);%computing first frame
            M = envList(1).M;
            %the self inductances and the transmitting-part sub-matrix must not be
            %recalculated
            M0(1:X*Y,1:X*Y) = M(1:X*Y,1:X*Y);
            M0 = M0 - diag(diag(M0)) + diag(diag(M));
            %calculating the rest of the frames
            for i=1:length(envList)
                envList(i) = evalM(envList(i),M0);
            end
            save(file,'envList');
        end

        if plotAnimation
            figure;
            plotCoil(coilPrototypeRX);
            figure;
            plotCoil(coilPrototypeTX);
            figure;
            animation(envList,0.4,0);
        end
        disp('Calculations finished.');
    else
        error('Something is wrong with the environments.')
    end
end
