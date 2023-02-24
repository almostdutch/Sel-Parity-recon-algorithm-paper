function [kspaceOP_spiritRecon, imgOP_spiritRecon, kspaceEP_spiritRecon, imgEP_spiritRecon, ...
    kspaceOP_spRecon, imgOP_spRecon, kspaceEP_spRecon, imgEP_spRecon, cpmgPhaseMap, residuals] = selpar_recon(ParsAndDataForRecon)
% function [kspaceOP_spiritRecon, imgOP_spiritRecon, kspaceEP_spiritRecon, imgEP_spiritRecon, ...
%     kspaceOP_spRecon, imgOP_spRecon, kspaceEP_spRecon, imgEP_spRecon, cpmgPhaseMap, residuals] = selpar_recon(ParsAndDataForRecon)
%
% Script implementing selpar recon algorithm #2
%

% kspace dims
[Nfe, Npe, Ncha] = size(ParsAndDataForRecon.kspaceOP);

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

% Apply weighting coefficients to kspace
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

% Selective parity recon
% Create kspace masks
maskOP = zeros(Nfe, Npe);
maskOP(:, Npe / 2 + 1 : end) = 1;
maskOP_pseudo = double(~maskOP); % complementary mask
if ismember(Npe / 2 + 1, ParsAndDataForRecon.indxLinesEP)
    maskOP_pseudo(:, Npe / 2 + 1, :) = 1; % if EP does have k0 line, set it to 1
end
maskOP_avg = maskOP + maskOP_pseudo; % taking care of kspace data overlap
maskOP_avg(maskOP_avg == 0) = 1; % no gaps!

% Flip kspaceEP left to right, OP and EP are originate from the same half of k-space
ParsAndDataForRecon.kspaceEP = circshift(fliplr(ParsAndDataForRecon.kspaceEP), [0 1]); 

% step 1 
% Estimate the missing lines
[kspaceOP_spiritRecon, kspaceEP_spiritRecon] = estimate_missing_kspace_lines(ParsAndDataForRecon);
imgOP_spiritRecon = ifftshift(ifft2(fftshift(kspaceOP_spiritRecon)));
imgEP_spiritRecon = ifftshift(ifft2(fftshift(kspaceEP_spiritRecon)));

% Memory preallocation
imgOP_spRecon = zeros(Nfe, Npe, Ncha);
imgEP_spRecon = zeros(Nfe, Npe, Ncha);
ParsAndDataForRecon.spNiter = 1; % non-iterative recon!
for chaNo = 1:Ncha    
    % step 2  
    % Correct for CPMG angle
    
    % CPMG phase map is fixed (from b0 data)
    cpmgPhaseMapCorrTerm = exp(-1i.*ParsAndDataForRecon.phaseMap(:,:,chaNo));
    
    % OP
    imgOP_corr = ifftshift(ifft2(fftshift(reshape(kspaceOP_spiritRecon(:, :, chaNo), [Nfe, Npe]))));
    imgOP_corr = imgOP_corr.*cpmgPhaseMapCorrTerm;
    
    % EP
    imgEP_corr = ifftshift(ifft2(fftshift(reshape(kspaceEP_spiritRecon(:, :, chaNo), [Nfe, Npe]))));
    imgEP_corr = imgEP_corr.*cpmgPhaseMapCorrTerm;
    
    % step 3
    % Generate pseudo-OP
    imgOP_corr_pseudo = conj(imgEP_corr);
    
    % step 4
    % Combine OP and pseudo-OP data   
    kspaceOP = fftshift(fft2(ifftshift(imgOP_corr)));
    kspaceOP_pseudo = fftshift(fft2(ifftshift(imgOP_corr_pseudo)));
    
    kspaceOP = maskOP.*kspaceOP + maskOP_pseudo.*kspaceOP_pseudo;
    kspaceOP = kspaceOP./maskOP_avg; 
    imageOP = ifftshift(ifft2(fftshift(kspaceOP)));
    
    imageEP = imageOP; % to keep consistent interface
    imgOP_spRecon(:, :, chaNo) = imageOP;
    imgEP_spRecon(:, :, chaNo) = imageEP;
end
residuals = norm(sos(imgOP_spRecon, 3) - sos(imgEP_spRecon, 3), 'fro') / norm(sos(imgOP_spRecon, 3), 'fro');
kspaceOP_spRecon = fftshift(fft2(ifftshift(imgOP_spRecon)));
kspaceEP_spRecon = fftshift(fft2(ifftshift(imgEP_spRecon)));
cpmgPhaseMap = ParsAndDataForRecon.phaseMap; % all channels
end
