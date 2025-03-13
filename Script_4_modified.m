clc;
close all;
clear all;

%% Simulation Parameters

squareLength = 1000;      % Length of the square coverage area in meters
nbrAPs = 64;              % Total number of APs in the area
numberOfAntennas = 1;     % Range of number of antennas
nbrAPsPerDim = sqrt(nbrAPs);  % APs per dimension
SNRedge = 1/8;            % Cell edge SNR (linear scale)
Kdims = 50;               % Number of UE locations in horizontal/vertical dimension
K = Kdims*Kdims;          % Number of UEs per cell
alpha = 4;                % Pathloss exponent
maxSE = 10;                % Maximum spectral efficiency
bandwidth = 10;           % Bandwidth in MHz

% Distance between APs in horizontal/vertical dimension
interSiteDistance = squareLength / nbrAPsPerDim;
interSiteDistanceCellular = squareLength / 2;

% Put out the APs on a square grid (centers of cells)
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2, [nbrAPsPerDim 1]);
locationsGridVertical   = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i * locationsGridVertical(:);

% Compute alternative AP locations by using wrap-around
wrapHorizontal = repmat([-squareLength 0 squareLength], [3, 1]);
wrapVertical   = wrapHorizontal';
wrapLocations  = wrapHorizontal(:)' + 1i * wrapVertical(:)';
APpositionsWrapped = repmat(BSpositions, [1, length(wrapLocations)]) + repmat(wrapLocations, [nbrAPs, 1]);

%% Power Consumption Models (same for all realizations)
% Cellular mMIMO (centralized)
P0_mMIMO = 15;       % [W] Fixed base power
P_ant_mMIMO = 0.1;   % [W] per antenna
P_total_mMIMO = P0_mMIMO + nbrAPs * P_ant_mMIMO;

% Cell-Free mMIMO (distributed)
P0_cellfree = 0.3;   % [W] Base power per AP
P_ant_cellfree = 0.05; % [W] Additional power per AP
P_total_cellfree = nbrAPs * (P0_cellfree + P_ant_cellfree);

%% Main Simulation Loop 
nbrSim = 100;   % Number of simulation realizations
EE_cellular_all = zeros(nbrSim,1);
EE_cellfree_all   = zeros(nbrSim,1);

% Aggregate SE values from all realizations
aggregatedSEs_cellular = [];
aggregatedSEs_cellfree   = [];

for sim = 1:nbrSim
    % For each realization, generate UE positions for all cells.
    % Use the same grid-based approach but add a random offset per cell.
    UEpositions = zeros(K, nbrAPs);
    for l = 1:nbrAPs
        posTmp = repmat(linspace(0,1,Kdims), [Kdims, 1]);
        posTmp2 = posTmp';
        shiftX = rand() * interSiteDistance;
        shiftY = rand() * interSiteDistance;
        posX = interSiteDistance * posTmp(:) + real(BSpositions(l)) - interSiteDistance/2 + shiftX;
        posY = interSiteDistance * posTmp2(:) + imag(BSpositions(l)) - interSiteDistance/2 + shiftY;
        UEpositions(:,l) = posX + 1i * posY;
    end
    
    % Prepare arrays for SE in this realization (per UE)
    SEs_cellular = zeros(K, nbrAPs, length(numberOfAntennas));
    SEs_cellfree   = zeros(K, nbrAPs, length(numberOfAntennas));
    
    % Compute SE for each cell
    for l = 1:nbrAPs
        distancesBSl = min( abs( repmat(UEpositions(:,l), [1, size(APpositionsWrapped,2)]) - ...
                               repmat(APpositionsWrapped(l,:), [K, 1]) ), [], 2 );
        signalPower = ((interSiteDistanceCellular/2) ./ distancesBSl).^alpha;
        interferencePower = zeros(size(distancesBSl));
        for j = 1:nbrAPs
            distancesBSj = min( abs( repmat(UEpositions(:,l), [1, size(APpositionsWrapped,2)]) - ...
                                 repmat(APpositionsWrapped(j,:), [K, 1]) ), [], 2 );
            if j ~= l
                interferencePower = interferencePower + ((interSiteDistanceCellular/2) ./ distancesBSj).^alpha;
            end
        end
        for m = 1:length(numberOfAntennas)
            SEs_cellfree(:,l,m) = log2(1 + SNRedge * numberOfAntennas(m) * (signalPower + interferencePower));
            SEs_cellfree(SEs_cellfree(:,l,m) > maxSE, l, m) = maxSE;
            SEs_cellular(:,l,m) = log2(1 + numberOfAntennas(m) * signalPower ./ (interferencePower + 1/SNRedge));
            SEs_cellular(SEs_cellular(:,l,m) > maxSE, l, m) = maxSE;
        end
    end
    
    % Aggregate the SE values across all UEs for this realization
    aggregatedSEs_cellular = [aggregatedSEs_cellular; SEs_cellular(:)];
    aggregatedSEs_cellfree   = [aggregatedSEs_cellfree; SEs_cellfree(:)];
    
    % Compute overall average SE (averaging over all cells and UEs) for this realization
    avgSE_cellular = mean(SEs_cellular(:));
    avgSE_cellfree   = mean(SEs_cellfree(:));
    
    % Compute overall Energy Efficiency for the realization
    EE_cellular = avgSE_cellular / P_total_mMIMO;
    EE_cellfree   = avgSE_cellfree / P_total_cellfree;
    
    EE_cellular_all(sim) = EE_cellular;
    EE_cellfree_all(sim)   = EE_cellfree;
end

%% Compute Overall Mean EE (across realizations)
meanEE_cellular = mean(EE_cellular_all);
meanEE_cellfree   = mean(EE_cellfree_all);

%% Plotting Results

% Energy Efficiency Analysis
figure;
bar([meanEE_cellular, meanEE_cellfree], 'FaceColor','flat');
set(gca, 'XTickLabel', {'Cellular mMIMO','Cell-Free mMIMO'});
ylabel('Energy Efficiency (bit/s/Hz per Watt)');
title('Energy Efficiency (Cell-Free vs Cellular)');
grid on;
set(gca, 'FontSize', 16);
xtips = get(gca, 'XTick');
ytips = [meanEE_cellular, meanEE_cellfree];
for i = 1:length(ytips)
    text(xtips(i), ytips(i), sprintf('%.2f', ytips(i)), 'HorizontalAlignment','center',...
         'VerticalAlignment','bottom', 'FontSize', 14);
end

% CDF Plot: Spectral Efficiency Distribution
figure;
cdfplot(aggregatedSEs_cellular);
hold on;
cdfplot(aggregatedSEs_cellfree);
xlabel('Spectral Efficiency (bit/s/Hz)');
ylabel('CDF');
title('CDF of Spectral Efficiency (Cell-Free vs Cellular)');
legend('Cellular mMIMO','Cell-Free mMIMO','Location','best');
grid on;
set(gca, 'FontSize', 16);

