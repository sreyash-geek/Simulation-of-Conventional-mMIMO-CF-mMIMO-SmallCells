clc;
close all;
clear all;


%% Simulation Parameters

squareLength = 1000; %Length of the square coverage area in meters
nbrAPs = 64; %Total number of APs in the area
numberOfAntennas = 1; %Range of number of antennas
nbrAPsPerDim = sqrt(nbrAPs); %APs per dimension
SNRedge = 1/8; %Set the cell edge SNR to -9 dB, which corresponds to 0 dB with an 9 dBi antenna
Kdims = 20; %Number of UE locations in horizontal/vertical dimension
K = Kdims*Kdims; %Number of UEs per cell
alpha = 4; %Pathloss exponent
maxSE = 8; %Maximum spectral efficiency
bandwidth = 10; %Bandwidth in MHz

%Distance between APs in horizontal/vertical dimension
interSiteDistance = squareLength/nbrAPsPerDim;
interSiteDistanceCellular = squareLength/2;
%Put out the APs on a square grid
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[nbrAPs 1]);

%Prepare to put out the UEs in the cells
UEpositions = zeros(K,nbrAPs);
perBS = zeros(nbrAPs,1);
%Prepare to compute SEs
SEs_cellular = zeros(K,nbrAPs,length(numberOfAntennas));
SEs_cellfree = zeros(K,nbrAPs,length(numberOfAntennas));

%% Main Simulation Loop

%Go through all the cells
for l = 1:nbrAPs
    %Put out K users in each cell
    while min(perBS(l))<K
        posTmp = repmat(linspace(0,1,Kdims),[Kdims 1]);
        posTmp2 = posTmp';
        posX = interSiteDistance*posTmp(:) + real(BSpositions(l)) - interSiteDistance/2;
        posY = interSiteDistance*posTmp2(:) + imag(BSpositions(l)) - interSiteDistance/2;
        UEpositions(:,l) = posX + 1i*posY;
        perBS(l) = K;
    end

    %Compute the distance from the UEs to AP l
    distancesSquaredBSl = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(l,:),[K 1]) ),[],2);
    %Compute the signal power for each UE
    signalPower = ((interSiteDistanceCellular/2)./distancesSquaredBSl).^(alpha);
    %Compute the interference power for each UE
    interferencePower = zeros(size(distancesSquaredBSl));
    for j = 1:nbrAPs
        %Compute the distance from the users to UE j
        distancesSquaredBSj = min(abs( repmat(UEpositions(:,l),[1 size(APpositionsWrapped,2)]) - repmat(APpositionsWrapped(j,:),[K 1]) ),[],2);
        %Add interference from non-serving APs
        if j~=l
            interferencePower = interferencePower + ((interSiteDistanceCellular/2)./distancesSquaredBSj).^(alpha);
        end
    end
    %Compute SE for different number of antennas, assuming i.i.d. Rayleigh fading and perfect CSI
    for m = 1:length(numberOfAntennas)
        SEs_cellfree(:,l,m) = log2(1+SNRedge*numberOfAntennas(m)*(signalPower + interferencePower));
        SEs_cellfree(SEs_cellfree(:,l,m)>maxSE,l,m) = maxSE;
        SEs_cellular(:,l,m) = log2(1+numberOfAntennas(m)*signalPower./(interferencePower+1/SNRedge));
        SEs_cellular(SEs_cellular(:,l,m)>maxSE,l,m) = maxSE;
    end
end

%% Plotting Results

%Plot the simulation results
figure;
hold on; box on; grid on;
for l = 1:nbrAPs
    surf(reshape(real(UEpositions(:,l)),[Kdims Kdims]),reshape(imag(UEpositions(:,l)),[Kdims Kdims]),bandwidth*reshape(SEs_cellfree(:,l,m),[Kdims Kdims]));
end
xlim([0 1000]);
ylim([0 1000]);
zlim([0 bandwidth*maxSE]);
xlabel('Position [m]');
ylabel('Position [m]');
zlabel('Data rate [Mbit/s/user]');
clim([0,bandwidth*maxSE]);
view(-52,28);
set(gca,'fontsize',16);

%Plot the simulation results
figure;
hold on; box on; grid on;
for l = 1:nbrAPs
    surf(reshape(real(UEpositions(:,l)),[Kdims Kdims]),reshape(imag(UEpositions(:,l)),[Kdims Kdims]),bandwidth*reshape(SEs_cellular(:,l,m),[Kdims Kdims]));
end
xlim([0 1000]);
ylim([0 1000]);
zlim([0 bandwidth*maxSE]);
xlabel('Position [m]');
ylabel('Position [m]');
zlabel('Data rate [Mbit/s/user]');
clim([0,bandwidth*maxSE]);
view(-52,28);
set(gca,'fontsize',16);