close all; clear all; clc

% Stanley Kogin
% BIOEN 217 Final Project

% File containing elements, their crystal structures, and lattice
% parameters.
small_f = readcell('CullityInfo');

% File containing coefficients to calculate the intensity factors for each
% element.
big_f = readcell('StructureFactorCoefficients');

% Wavelength 'lamda' is assumed to be that of K_alpha peak off a target of
% selected material.
target = input('Select target element; choices are Mo, Cu, Co, Fe, & Cr: ', 's');

if target == 'Mo'
    lamda = 0.711;
elseif target == 'Cu'
    lamda = 1.542;
elseif target == 'Co'
    lamda = 1.790;
elseif target == 'Fe'
    lamda = 1.937;
elseif target == 'Cr'
    lamda = 2.291;
else 
    disp('Target is unavailable, please retry with another');    
end

element = input('Select element/ion pattern: ', 's');
ind_1 = 0;
ind_2 = 0;

for i = 2:length(small_f(:,1))
    p = cast(small_f(i,1), 'char');
    if length(element) == length(p) && sum(element == p) == length(element)
            ind_1 = i;
            break;
    end
end

for i = 2:length(big_f(:,1))
    l = cast(big_f(i,1), 'char');
    if length(element) == length(l) && sum(element == l) == length(element)
            ind_2 = i;
            break;
    end
end

if ind_1 == 0
    disp('Element/ion is unavailable, please retry with another element');
end

structure = cast(small_f(ind_1, 2), 'char');

% Set of simple cubic plane Miller coordinates in increasing h^2 + k^2 + l^2 order.
uvw = [1 0 0; 1 1 0; 1 1 1; 2 0 0; 2 1 0; 2 1 1; 2 2 0; 3 0 0; 3 1 0; 
    3 1 1; 2 2 2; 3 2 0; 3 2 1; 4 0 0; 4 1 0; 4 1 1; 3 3 1; 4 2 0;
    4 2 1; 3 3 2; 4 2 2; 5 0 0; 5 1 0; 5 1 1; 5 2 0; 5 2 1; 4 4 0];

% Using default simple cubic planar coordinates to find which ones are
% acceptable for each crystal structure.
check = 1;

% Specific Miller coordinates for the selected crystal structure.
hkl = zeros(1,3);

for i = 1:length(uvw)
    if structure == 'SCC'
        hkl = uvw;
        check = check + 1;
    elseif mod(sum(uvw(i,:)), 2) == 0 & structure == 'BCC'
        hkl(check,:) = uvw(i,:);
        check = check + 1;
    elseif (sum(mod(uvw(i,:), 2)) == 3 || sum(mod(uvw(i,:), 2)) == 0) & structure == 'FCC'
        hkl(check,:) = uvw(i,:);
        check = check + 1;
    elseif mod((uvw(i,1) + 2*uvw(i,2)), 3) ~= 0 & mod(uvw(i,3),2) == 0 & structure == 'HCP'
        hkl(check,:) = uvw(i,:);
        check = check + 1;       
    end
end

% Very important matrix that will serve as our XRD table, values from left
% to right are; (1) planar spacing, (2) theta, (3) q-factor, (4) intensity factor,
% (5) absorption factors, (6) multiplicity, (7) Lorentz-Polarization factor, (8) intensity,
% and finally (9) intensity percentage. Each value is calculated off the last
% value in the table.
xrd = zeros(length(hkl(:,1)),9);

xrd(:,1) = cell2mat(small_f(ind_1,3)) ./ sqrt(hkl(:,1).^2 + hkl(:,2).^2 + hkl(:,3).^2);
xrd(:,2) = asind(lamda ./ (2*xrd(:,1)));
xrd(:,3) = sind(xrd(:,2)) / lamda;

% Coefficients need to converted from 'cell' to 'char' then to doubles.
coeffs = str2num(cast(big_f(ind_2,[2:10]), 'char'));

% Calculating the structure factors.
for i = 1:2:7
    xrd(:,4) = xrd(:,4) + coeffs(i) .* exp(-coeffs(i + 1) .* xrd(:,3).^2) + coeffs(9);
end

if structure == 'BCC'
    xrd(:,5) = 4*xrd(:,4).^2;
elseif structure == 'FCC'
    xrd(:,5) = 16*xrd(:,4).^2;
else 
    xrd(:,5) = xrd(:,4).^2;
end

% Checking multiplicity conditions.
for i = 1:length(hkl(:,1))
    if hkl(i,1) ~= hkl(i,2) && hkl(i,1) ~= hkl(i,3) && hkl(i,2) ~= hkl(i,3)
        if hkl(i,1) ~= 0 && hkl(i,2) ~= 0 && hkl(i,3) ~= 0
            xrd(i,6) = 48;
        else
            xrd(i,6) = 24;
        end
        
    elseif hkl(i,1) == hkl(i,2) && hkl(i,1) == hkl(i,3)
        xrd(i,6) = 8;
    elseif sum(hkl(i,:)) == hkl(i,1) || sum(hkl(i,:)) == hkl(i,2) || sum(hkl(i,:)) == hkl(i,3)        
        xrd(i,6) = 6;     
    elseif sum(hkl(i,[1 2])) == sum(hkl(i,:)) || sum(hkl(i,[2 3])) == sum(hkl(i,:)) || sum(hkl(i,[1:2:3])) == sum(hkl(i,:))
        xrd(i,6) = 12;
    else
        xrd(i,6) = 24;
    end
end

% Calculating Lorentz-Polarization Factors
xrd(:,7) = (1 + cosd(2.*xrd(:,2)).^2) ./  (cosd(xrd(:,2)).*sind(xrd(:,2)).^2);

% Calculating intensities and intensity percentages.
xrd(:,8) = xrd(:,5) .* xrd(:,6) .* xrd(:,7);
xrd(:,9) = xrd(:,8) ./ max(xrd(:,8));

% Plotting each XRD peak.
ind_3 = 1;
figure(1) 
for i = 1:length(xrd(:,1))
    if xrd(i,2) < 90
        plot(2*xrd(i,2)*ones(2), [0 xrd(i,9)], 'b');
        xlabel('2-theta'); ylabel('Intensity %');
        xlim([0 180]);
        text(2*xrd(i,2), xrd(i,9), mat2str(hkl(i,:)));
        hold on;
        pause(0.2);
    end
end


% Citations [IEEE]
% B. D. Cullity and S. R. Stock, Elements of X-ray diffraction. Miejsce nieznane: Pearson India Education Services, 2015.
% P. Hadley, ?Atomic form factors,? Atomic form factors. [Online]. Available: http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php. [Accessed: 03-Dec-2019].