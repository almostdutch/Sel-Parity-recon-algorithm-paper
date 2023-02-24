function [kspaceOP_spiritRecon, imgOP_spiritRecon, kspaceEP_spiritRecon, imgEP_spiritRecon, ...
    kspaceOP_spRecon, imgOP_spRecon, kspaceEP_spRecon, imgEP_spRecon, cpmgPhaseMap, residuals] = selpar_recon(ParsAndDataForRecon)
% function [kspaceOP_spiritRecon, imgOP_spiritRecon, kspaceEP_spiritRecon, imgEP_spiritRecon, ...
%     kspaceOP_spRecon, imgOP_spRecon, kspaceEP_spRecon, imgEP_spRecon, cpmgPhaseMap, residuals] = selpar_recon(ParsAndDataForRecon)
%
% Script implementing selpar recon algorithm #1
%

% kspace dims
[Nfe, Npe, Ncha] = size(ParsAndDataForRecon.kspaceOP);

if size(ParsAndDataForRecon.phaseMap, 3) ~= Ncha
    errorMsg = ['Dimensions mismatch, both (CPMG phase map and kspaceOP/kspaceEP) should have the same number of channels.'...
        'Hint: perhaps, something to do with the coil compression, aborting...\n'];
    error(errorMsg);
end

% Weighting scheme to control kspace signal decay and image PSF
wts_func_handle = @(a, b, t) 1 - a.*exp(-b.*t);
a_coef = ParsAndDataForRecon.a_coef;
b_coef = ParsAndDataForRecon.b_coef;
t = 1:Npe;
wts2 = wts_func_handle(a_coef, b_coef, t);
wts2_OP = wts2(ParsAndDataForRecon.echoParityArray == 1);
wts2_EP = wts2(ParsAndDataForRecon.echoParityArray == 0);
wts = ones(1, Npe);
wts(ParsAndDataForRecon.indxLinesOP) = wts2_OP;
wts(ParsAndDataForRecon.indxLinesEP) = wts2_EP;
wts = repmat(wts, [Nfe, 1]);

% Apply weighting scheme to kspace
ParsAndDataForRecon.kspaceOP = bsxfun(@times, ParsAndDataForRecon.kspaceOP, wts);
ParsAndDataForRecon.kspaceEP = bsxfun(@times, ParsAndDataForRecon.kspaceEP, wts);

% Simplified selective parity recon for b0 data
if strcmp(ParsAndDataForRecon.selparReconFlag, 'no')
    % Effectively just SPIRIT recon

    % No phase difference, hence safe to combine OP and EP kspace data
    ParsAndDataForRecon.kspaceOP = ParsAndDataForRecon.kspaceOP + ParsAndDataForRecon.kspaceEP;
    ParsAndDataForRecon.kspaceEP = ParsAndDataForRecon.kspaceOP;
    
    % Estimate the missing kspace data
    [kspaceOP_spiritRecon, kspaceEP_spiritRecon] = estimate_missing_kspace_lines(ParsAndDataForRecon);
    imgOP_spiritRecon = ifftshift(ifft2(fftshift(kspaceOP_spiritRecon)));
    imgEP_spiritRecon = ifftshift(ifft2(fftshift(kspaceEP_spiritRecon)));
    
    % Selective parity recon
    kspaceOP_spRecon = kspaceOP_spiritRecon;
    kspaceEP_spRecon = kspaceEP_spiritRecon;
    imgOP_spRecon = imgOP_spiritRecon;
    imgEP_spRecon = imgEP_spiritRecon;
    
    % CPMG phase map
    cpmgPhaseMap = angle(imgOP_spiritRecon);
    
    % Residuals
    residuals =  0;
    
    return;
end

% Two stage selective recon for bxxxx data
correctionFactor = 2;
% step 1 
% Estimate the missing lines
[kspaceOP_spiritRecon, kspaceEP_spiritRecon] = estimate_missing_kspace_lines(ParsAndDataForRecon);
imgOP_spiritRecon = ifftshift(ifft2(fftshift(kspaceOP_spiritRecon)));
imgEP_spiritRecon = ifftshift(ifft2(fftshift(kspaceEP_spiritRecon)));

% Stage 1
% Match the power of EP to that of OP by using the external CPMG phase map
% (from b0 data) to reduce the degrees of freedom
kspaceOP_spRecon = kspaceOP_spiritRecon;
kspaceEP_spRecon = kspaceEP_spiritRecon;
imgOP_spRecon = zeros(Nfe, Npe, Ncha);
imgEP_spRecon = zeros(Nfe, Npe, Ncha);
residuals = zeros(1, ParsAndDataForRecon.spNiter);
for iterNo = 1:ParsAndDataForRecon.spNiter
    for chaNo = 1:Ncha
        % step 2: transform to image domain
        imgOP_spReconTemp = ifftshift(ifft2(fftshift(kspaceOP_spiritRecon(:, :, chaNo)))); % fixed
        imgEP_spReconTemp = ifftshift(ifft2(fftshift(kspaceEP_spRecon(:, :, chaNo)))); % updated on each iter to increase the power
        
        % step 3: CPMG phase map is fixed (from b0 data)
        cpmgPhaseMap = ParsAndDataForRecon.phaseMap(:,:,chaNo);
        
        % step 4: generate pseudo images for OP and EP
        imgOP_pseudo = conj(imgEP_spReconTemp.*exp(-1i.*cpmgPhaseMap)).*exp(1i.*cpmgPhaseMap);
        imgEP_pseudo = conj(imgOP_spReconTemp.*exp(-1i.*cpmgPhaseMap)).*exp(1i.*cpmgPhaseMap);
        
        % step 5: force information sharing between either parity and its pseudo equivalent
        imgOP_new = imgOP_spReconTemp + (imgOP_pseudo - imgOP_spReconTemp) / correctionFactor;
        imgEP_new = imgEP_spReconTemp + (imgEP_pseudo - imgEP_spReconTemp) / correctionFactor;
        
        % step 6: transform to kspace
        kspaceOP_new = fftshift(fft2(ifftshift(imgOP_new)));
        kspaceEP_new = fftshift(fft2(ifftshift(imgEP_new)));
        
        % step 7: force data consistency
        kspaceOP_new(:,ParsAndDataForRecon.indxLinesOP) = ParsAndDataForRecon.kspaceOP(:, ParsAndDataForRecon.indxLinesOP, chaNo);
        kspaceEP_new(:,ParsAndDataForRecon.indxLinesEP) = ParsAndDataForRecon.kspaceEP(:, ParsAndDataForRecon.indxLinesEP, chaNo);
        
        % step 8: update selpar kspace data for subsequent iterations
        kspaceOP_spRecon(:, :, chaNo) = kspaceOP_new;
        kspaceEP_spRecon(:, :, chaNo) = kspaceEP_new;
        
        % Keep track of residuals
        imgOP_spRecon(:, :, chaNo) = ifftshift(ifft2(fftshift(kspaceOP_new)));
        imgEP_spRecon(:, :, chaNo) = ifftshift(ifft2(fftshift(kspaceEP_new)));
    end
    % residuals based on SoS images
    residuals(iterNo) = norm(sos(imgOP_spRecon, 3) - sos(imgEP_spRecon, 3), 'fro') / norm(sos(imgOP_spRecon, 3), 'fro');
end

% Stage 2
% Do the actual selective parity recon
kspaceOP_spRecon = kspaceOP_spiritRecon;
kspaceEP_spRecon = kspaceEP_spRecon; % from % Stage 1
imgOP_spRecon = zeros(Nfe, Npe, Ncha);
imgEP_spRecon = zeros(Nfe, Npe, Ncha);
residuals = zeros(1, ParsAndDataForRecon.spNiter);
for iterNo = 1:ParsAndDataForRecon.spNiter
    for chaNo = 1:Ncha
        % step 2: transform to image domain
        imgOP_spReconTemp = ifftshift(ifft2(fftshift(kspaceOP_spRecon(:, :, chaNo))));
        imgEP_spReconTemp = ifftshift(ifft2(fftshift(kspaceEP_spRecon(:, :, chaNo))));
        
        % step 3: CPMG phase map is fixed (from b0 data)
        cpmgPhaseMap = ParsAndDataForRecon.phaseMap(:,:,chaNo);
        
        % step 4: generate pseudo images for OP and EP
        imgOP_pseudo = conj(imgEP_spReconTemp.*exp(-1i.*cpmgPhaseMap)).*exp(1i.*cpmgPhaseMap);
        imgEP_pseudo = conj(imgOP_spReconTemp.*exp(-1i.*cpmgPhaseMap)).*exp(1i.*cpmgPhaseMap);
        
        % step 5: force information sharing between either parity and its pseudo equivalent
        imgOP_new = imgOP_spReconTemp + (imgOP_pseudo - imgOP_spReconTemp) / correctionFactor;
        imgEP_new = imgEP_spReconTemp + (imgEP_pseudo - imgEP_spReconTemp) / correctionFactor;
        
        % step 6: transform to kspace
        kspaceOP_new = fftshift(fft2(ifftshift(imgOP_new)));
        kspaceEP_new = fftshift(fft2(ifftshift(imgEP_new)));
        
        % step 7: force data consistency
        kspaceOP_new(:,ParsAndDataForRecon.indxLinesOP) = ParsAndDataForRecon.kspaceOP(:, ParsAndDataForRecon.indxLinesOP, chaNo);
        kspaceEP_new(:,ParsAndDataForRecon.indxLinesEP) = ParsAndDataForRecon.kspaceEP(:, ParsAndDataForRecon.indxLinesEP, chaNo);
        
        % step 8: update selpar kspace data for subsequent iterations
        kspaceOP_spRecon(:, :, chaNo) = kspaceOP_new;
        kspaceEP_spRecon(:, :, chaNo) = kspaceEP_new;
        
        % Keep track of residuals
        imgOP_spRecon(:, :, chaNo) = ifftshift(ifft2(fftshift(kspaceOP_new)));
        imgEP_spRecon(:, :, chaNo) = ifftshift(ifft2(fftshift(kspaceEP_new)));
    end
    % residuals based on SoS images
    residuals(iterNo) = norm(sos(imgOP_spRecon, 3) - sos(imgEP_spRecon, 3), 'fro') / norm(sos(imgOP_spRecon, 3), 'fro');
end
cpmgPhaseMap = ParsAndDataForRecon.phaseMap; % all channels
end

