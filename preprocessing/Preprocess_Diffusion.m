function Preprocess_Diffusion(vol_B0uw, M, qmat, bvals, gwinfo, dcminfo, outdir, isFOCUS, isARTPro, pe_dim, params)

if exist('params', 'var')
  EddyCorrFlag = params.EddyCorrFlag;
  GradWarpFlag = params.GradWarpFlag;
  MotionCorrFlag = params.MotionCorrFlag;
else
  EddyCorrFlag = 1;
  GradWarpFlag = 1;
  MotionCorrFlag = 1;
end


% Eddy current correction
if EddyCorrFlag && ~isempty(qmat)
    fprintf('%s -- %s.m:    Running Eddy Current Correction...\n',datestr(now),mfilename);
    vol_ecc = QD_Eddy(vol_B0uw, qmat, bvals, pe_dim);
    fname_ecc = sprintf('%s/DWI_vol_preprocessed.mgz', outdir);
    QD_save_mgh(single(vol_ecc), fname_ecc, M);
else
    if isempty(qmat)
       fprintf('%s -- %s.m:    WARNING: No qmat, skipping eddy current correction...\n',datestr(now),mfilename);
    end
    vol_ecc = vol_B0uw; 
    fname_ecc = sprintf('%s/DWI_vol_preprocessed.mgz', outdir);
    QD_save_mgh(single(vol_ecc), fname_ecc, M);
end


% Gradient Unwarping
if GradWarpFlag && ~isempty(gwinfo) && isfield(gwinfo, 'gwtype')
  fprintf('%s -- %s.m:    Gradient Unwarping...\n',datestr(now),mfilename);
  vol_uw = zeros(size(vol_ecc));

  manufacturer = dcminfo.Manufacturer;
  match_GE = ~isempty(regexpi(manufacturer, 'ge'));
  match_Siemens = ~isempty(regexpi(manufacturer, 'siemens'));

  if match_GE
    if isFOCUS || isARTPro
      gradwarp_flag = 1;
      disp('Applying gradwarp correction through-plane only');
    else 
      % Otherwise, assume research sequence with on-scanner gradwarp correction disabled 
      gradwarp_flag = 0;
      disp('Applying gradwarp correction in 3D');
    end
  elseif match_Siemens
    gradwarp_flag = gwinfo.unwarpflag;
    disp('Applying gradwarp correction according to mmil_get_gradwarpinfo');
  else
    gradwarp_flag = gwinfo.unwarpflag;
    disp('Applying gradwarp correction according to mmil_get_gradwarpinfo');
  end

  for i = 1:size(vol_uw,4)
    voltmp = mgh2ctx(vol_ecc(:,:,:,i),M);
    voltmp_uw = ctx_unwarp_grad(voltmp, gwinfo.gwtype, gradwarp_flag, gwinfo.isoctrflag);
    vol_uw(:,:,:,i) = voltmp_uw.imgs;
  end

else
  fprintf('%s -- %s.m:    WARNING: No gradient unwarping performed...\n',datestr(now),mfilename);
  vol_uw = vol_ecc;
end
fname_uw = sprintf('%s/DWI_vol_preprocessed.mgz', outdir);
QD_save_mgh(single(vol_uw), fname_uw, M);


% Motion Correction 
if MotionCorrFlag && ~isempty(qmat)
    fprintf('%s -- %s.m:    Running Motion Correction...\n',datestr(now),mfilename);
    [vol_ecc_mc, qmat_mc, Mreg_mat] = QD_MotionCorr(vol_uw, M, qmat, bvals);
    fname_mc = sprintf('%s/DWI_vol_preprocessed.mgz', outdir);
    QD_save_mgh(single(vol_ecc_mc), fname_mc, M);
    save(sprintf('%s/qmat_corrected.mat', outdir), 'qmat_mc');
else
    if isempty(qmat)
      fprintf('%s -- %s.m:    WARNING: No qmat, skipping motion correction...\n',datestr(now),mfilename);
    end
    fname_mc = fname_uw;
    qmat_mc = qmat;
    save(sprintf('%s/qmat_corrected.mat', outdir), 'qmat_mc');
end


fprintf('%s -- %s.m:    Finished Successfully.\n',datestr(now),mfilename);

end
