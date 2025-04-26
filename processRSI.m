function status = processRSI(RSI_path, RSI_path_rev, T2_path, DCE_path, output_path, params)
% processRSI - Process RSI data set
%
% Input Arguments
%   RSI_path: Path to directory containing DICOM files for a single RSI series
%   RSI_path_rev (optional): Path to directory containing DICOM files for a single RSI series with phase-encoding direction opposite to that of RSI_path
%   T2_path (optional): Path to directory containing DICOM files for a reference anatomical T2w volume
%   DCE_path (optional): Path to directory (or directories) containing DICOM files for the raw (unsubtracted) images from a single DCE acquisition
%                        Type: Cell array containing one or more character vectors; {'', '', ''} 
%   output_path: Path to directory where results will be saved
%   params: Structure of input parameters


%% ============== Parse Inputs ============== %%

status = 1;

% Make output directory
if ~exist(output_path,'dir')
    mkdir(output_path);
end

% Check if T2 needs to be ignored
if ~isempty(T2_path)
  T2_flag = 1;
else
  T2_flag = 0;
end

% Swap fwd and rev scans, if applicable, s.t. the fwd scan stretches the anatomy (hopefully)
if ~isempty(RSI_path_rev)
  if length(dir(RSI_path)) == length(dir(RSI_path_rev))
    file_fwd = dir(RSI_path);
    file_fwd = fullfile(file_fwd(3).folder, file_fwd(3).name);
    info_fwd = dicominfo(file_fwd);
    Manufacturer = info_fwd.Manufacturer;

    if strcmpi(info_fwd.Manufacturer, 'ge medical systems')
      seq_name_fwd = info_fwd.Private_0019_109c;
      tmp_path1 = RSI_path;
      tmp_path2 = RSI_path_rev;
      if any(strcmp(seq_name_fwd, {'epi2alt', 'epi2_ART', 'epi2_pepolarFOCUSFLEX', 'epi2_pepolarFLEX'}))
	RSI_path = tmp_path2;
	RSI_path_rev = tmp_path1;
      end
    end

  end
end


%% ============== load RSI Data ============== %%

% check dicom data
impax_flag = 0;

fprintf('%s -- %s.m:    Loading RSI data...\n',datestr(now),mfilename);
if iscell(RSI_path)
    flist = RSI_path;
else
    flist = recursive_dir(RSI_path);
end
for i = 1:length(flist)
    try
        dcminfo = dicominfo(flist{i}); 
	dcminfo = fix_impax_dcm_tags(dcminfo);
	break
    catch
        continue
    end
end

if ~exist('dcminfo','var')
    error('No RSI DICOMs found.'); 
end

Manufacturer = dcminfo.Manufacturer;

% Read vendor-specific tag labels
if strcmpi(Manufacturer, 'ge medical systems')
   disp('Reading GE vendor tags');
   dicomdict('set', 'gems-dicom-dict.txt');
   dcminfo_GE = dicominfo(flist{i});
   if T2_flag
     flist_T2 = recursive_dir(T2_path);
     dcminfo_T2_GE = dicominfo(flist_T2{i});
   end
   if ~isempty(RSI_path_rev)
     flist_rev = recursive_dir(RSI_path_rev);
     dcminfo_rev_GE = dicominfo(flist_rev{i});
   end
   dicomdict('factory');
end

if isfield(dcminfo, 'Private_0043_1039') == 1 && strcmp(class(dcminfo.Private_0043_1039),'char') == 1
   impax_flag = 1;
   params.B0DISCO = 0;
   params.EddyCorrFlag = 0;
   params.MotionCorrFlag = 0;
end

isARTPro = 0;
switch lower(Manufacturer)
       
  case 'ge medical systems' % GE --------------------------------------------------------------------

    % Check if data is FOCUS or not
    isFOCUS = 0;
    if isfield(dcminfo, 'Private_0019_109c') == 1
      if ischar(dcminfo.Private_0019_109c) == 1
	if isempty(strfind(lower(dcminfo.Private_0019_109c), 'focus')) == 0
	  isFOCUS = 1;
	end
      end
    end

    if isempty(strfind(lower(dcminfo.SeriesDescription), 'focus')) == 0
      isFOCUS = 1;
    end

    if isfield(params, 'ForceFOCUS') == 1
      if params.ForceFOCUS == 1
	isFOCUS = 1;
      elseif params.ForceFOCUS == 0
	isFOCUS = 0;
      end
    end

    % Check for ART-Pro protocol
    if isfield(dcminfo, 'Private_0019_109c') 
      if ~isempty(strfind(lower(dcminfo.Private_0019_109c), 'art')) || ~isempty(strfind(lower(dcminfo.Private_0019_109c), 'alt'))
	isARTPro = 1;
      end
    end

    % grab tensor.dat file
    if isfield(params, 'TensorFile')
      tensorfile = params.TensorFile;
    elseif isfield(dcminfo, 'Private_0019_10b6') == 1
      tensor_ID = dcminfo.Private_0019_10b6;
      tensorfile = sprintf('tensor%d.dat',tensor_ID);
      disp(['Tensorfile: ' tensorfile]);
    else
      tensorfile = 'tensor0.dat';
    end

    % load dicom data
    if isempty(which(tensorfile)) || strcmp(tensorfile, 'tensor0.dat')

      disp('Assuming GE product DWI or acquisition on DV29+ platform')
      [rsidat, M, qmat, bvals, gwinfo_rsi, dcminfo] = ReadDicomDiffusionData(RSI_path);
      if numel(unique(bvals)) == 2 % DV29.0 with valid qmat, but only b=0 and bmax listed in DICOM header
	 bmax = max(bvals);
	 bvals = bmax*sum(qmat.^2,2);
	 bvals(bvals<400) = round(bvals(bvals<400));
	 bvals(bvals>400&bvals<1000) = round(bvals(bvals>400&bvals<1000), -1);
	 bvals(bvals>=1000) = round(bvals(bvals>=1000), -2);
      end

      if ~isempty(RSI_path_rev)
	[rsidat_rev, M_rev, qmat_rev, bvals_rev, gwinfo_rsi_rev, dcminfo_rev] = ReadDicomDiffusionData(RSI_path_rev);
	if numel(unique(bvals_rev)) == 2 % DV29.0 with valid qmat, but only b=0 and bmax listed in DICOM header
	   bmax_rev = max(bvals_rev);
	   bvals_rev = bmax_rev*sum(qmat_rev.^2,2);
	   bvals_rev(bvals_rev<400) = round(bvals_rev(bvals_rev<400));
	   bvals_rev(bvals_rev>400&bvals_rev<1000) = round(bvals_rev(bvals_rev>400&bvals_rev<1000), -1);
	   bvals_rev(bvals_rev>=1000) = round(bvals_rev(bvals_rev>=1000), -2);
	end
      end

    else

      disp('Using custom tensor file')
      [rsidat, M, qmat, bvals, gwinfo_rsi, dcminfo] = ReadDicomDiffusionData(RSI_path);
      bmax = max(bvals);
      numt2s = dcminfo.Private_0019_10df;
      numdirs = dcminfo.Private_0019_10e0;
      qmat = ReadTensorFile(numdirs,numt2s,tensorfile);
      bvals = bmax*sum(qmat.^2,2);
      
      bvals(bvals<400) = round(bvals(bvals<400));
      bvals(bvals>400&bvals<1000) = round(bvals(bvals>400&bvals<1000), -1);
      bvals(bvals>=1000) = round(bvals(bvals>=1000), -2);

      if ~isempty(RSI_path_rev)
        [rsidat_rev, M_rev, qmat_rev, bvals_rev, gwinfo_rsi_rev, dcminfo_rev] = ReadDicomDiffusionData(RSI_path_rev);
	bmax_rev = max(bvals_rev);
        numt2s = dcminfo_rev.Private_0019_10df;
	numdirs = dcminfo_rev.Private_0019_10e0;
        qmat_rev = ReadTensorFile(numdirs,numt2s,tensorfile);
        bvals_rev = bmax_rev*sum(qmat_rev.^2,2);
        
	bvals_rev(bvals_rev<400) = round(bvals_rev(bvals_rev<400));
	bvals_rev(bvals_rev>400&bvals_rev<1000) = round(bvals_rev(bvals_rev>400&bvals_rev<1000), -1);
	bvals_rev(bvals_rev>=1000) = round(bvals_rev(bvals_rev>=1000), -2);
      end

    end

    % Check if sequence is integrated or not
    integrated_rev_flag = 1;
    if isfield(dcminfo, 'Private_0019_109c')
      pulse_seq_name = lower(dcminfo.Private_0019_109c);
      if ~isempty( strfind(pulse_seq_name, 'pepolar') )
	integrated_rev_flag = 1;
      elseif ~isempty(strfind(pulse_seq_name, 'art')) || ~isempty(strfind(pulse_seq_name, 'alt'))
	integrated_rev_flag = 0;
      elseif strcmp(pulse_seq_name, 'epi2')
	integrated_rev_flag = 0;
      end
    end

    if ~isempty(RSI_path_rev)
      integrated_rev_flag_rev = 1;
      if isfield(dcminfo_rev, 'Private_0019_109c')
        pulse_seq_name_rev = lower(dcminfo_rev.Private_0019_109c);
        if ~isempty( strfind(pulse_seq_name_rev, 'pepolar') )
          integrated_rev_flag_rev = 1;
        elseif ~isempty(strfind(pulse_seq_name_rev, 'art')) || ~isempty(strfind(pulse_seq_name_rev, 'alt'))
          integrated_rev_flag_rev = 0;
        elseif strcmp(pulse_seq_name_rev, 'epi2')
	  integrated_rev_flag_rev = 0;
	end
      end
    end
    

  case {'siemens', 'siemens healthineers'}  % Siemens --------------------------------------------------------------------

    isFOCUS = 0;
    if ~isempty(regexpi(dcminfo.SeriesDescription, 'ZOOMit'))
      isFOCUS = 1;
    end

    integrated_rev_flag = 0;
    integrated_rev_flag_rev = 0;

    [rsidat, M, qmat, bvals, dcminfo] = QD_Read_DICOM_Diffusion_Directory_Siemens(RSI_path);
    deratefac = sqrt(bvals/max(bvals));
    qmat = qmat.*repmat(colvec(deratefac), [1 3]);
    [gwinfo_rsi, ~] = mmil_get_gradwarpinfo(dcminfo);

    if ~isempty(RSI_path_rev)
      [rsidat_rev, M_rev, qmat_rev, bvals_rev, dcminfo_rev] = QD_Read_DICOM_Diffusion_Directory_Siemens(RSI_path_rev);
      if max(bvals_rev) > 0
	deratefac = sqrt(bvals_rev/max(bvals_rev));
	qmat_rev = qmat_rev.*repmat(colvec(deratefac), [1 3]);
      end
      [gwinfo_rsi_rev, ~] = mmil_get_gradwarpinfo(dcminfo_rev);
    end


  case {'philips', 'philips healthcare'}  % Philips --------------------------------------------------------------------

    isFOCUS = 0;
    if ~isempty(regexpi(dcminfo.SeriesDescription, 'ZOOM'))
      isFOCUS = 1;
    end

    integrated_rev_flag = 0;
    integrated_rev_flag_rev = 0;

    [rsidat, M, qmat, bvals, gwinfo_rsi, dcminfo] = ReadDicomDiffusionDataPhilips(RSI_path);
    if bvals(end)>0 && all(qmat(end,:)==[0 0 0]) % Remove synthesized image if it exists
      fprintf('Chopping off synthesized volume\n');
      qmat = qmat(1:end-1,:);
      bvals = bvals(1:end-1);
      rsidat = rsidat(:,:,:,1:end-1);
    end

    if ~isempty(RSI_path_rev)
      [rsidat_rev, M_rev, qmat_rev, bvals_rev, gwinfo_rsi_rev, dcminfo_rev] = ReadDicomDiffusionDataPhilips(RSI_path_rev);
      if bvals_rev(end)>0 && all(qmat_rev(end,:)==[0 0 0]) % Remove synthesized image if it exists
	fprintf('Chopping off synthesized volume\n');
	qmat_rev = qmat_rev(1:end-1,:);
	bvals_rev = bvals_rev(1:end-1);
	rsidat_rev = rsidat_rev(:,:,:,1:end-1);
      end
    end


  otherwise
    error('Unsupported scanner manufacturer')

end


% Save unprocessed DWI data
fprintf('%s -- %s.m:    Saving unprocessed DWI data...\n', datestr(now), mfilename);
if strcmpi(Manufacturer, 'ge medical systems')
   disp('Saving GE vendor tags');
   dcminfo_tmp = dcminfo;
   dcminfo = dcminfo_GE;
   save(fullfile(output_path, 'dcminfo.mat'), 'dcminfo');
   dcminfo = dcminfo_tmp;
else 
   save(fullfile(output_path, 'dcminfo.mat'), 'dcminfo');
end

if integrated_rev_flag == 1
  QD_save_mgh( rsidat, fullfile(output_path,'DWI_vol_unprocessed_concat_fwd_rev.mgz'), M );
else
  QD_save_mgh( rsidat, fullfile(output_path,'DWI_vol_unprocessed_fwd.mgz'), M );
end

if exist('rsidat_rev', 'var')
  if integrated_rev_flag_rev == 1
    QD_save_mgh( rsidat_rev, fullfile(output_path,'REV_DWI_vol_unprocessed_concat_fwd_rev.mgz'), M_rev );
  elseif integrated_rev_flag_rev == 0
    QD_save_mgh( rsidat_rev, fullfile(output_path,'REV_DWI_vol_unprocessed_fwd.mgz'), M_rev );
  end
end


% Separate forward and reverse PE polarity volumes
pedim = 2; % column direction - A/P
if strcmp(dcminfo.InPlanePhaseEncodingDirection,'ROW')
  fprintf('%s -- %s.m:    WARNING: Left/Right Phase encode direction!\n', datestr(now), mfilename);
  pedim = 1;
  params.B0DISCO = 0;
end

if strcmpi(Manufacturer,'ge medical systems')

  % Flip "reverse" volumes (which are now the "forward" volumes)
  if exist('rsidat_rev', 'var') && (integrated_rev_flag==1)
    rsidat = flip(rsidat, pedim);
  end
  
  b0inds = find(bvals == 0);
  
  if integrated_rev_flag == 0
    fprintf('%s -- %s.m:    This is NOT an integrated sequence...\n', datestr(now), mfilename);
    b0_vol_fwd = mgh2ctx( mean(rsidat(:,:,:,b0inds),4), M );
    if exist('rsidat_rev', 'var')
      b0inds_rev = find(bvals_rev == 0);
      b0_vol_rev = mgh2ctx( mean(rsidat_rev(:,:,:,b0inds_rev),4), M_rev );
    else
      fprintf('%s -- %s.m:    Only one b=0 volume found, so no B0 distortion correction...\n', datestr(now), mfilename);
      b0_vol_rev = [];
      params.B0DISCO = 0;
    end
    
  elseif integrated_rev_flag == 1
    fprintf('%s -- %s.m:    This IS an integrated sequence...\n', datestr(now), mfilename);
    b0_vol_rev = mgh2ctx( flip(rsidat(:,:,:,b0inds(1)),pedim), M );
    b0_vol_fwd = mgh2ctx( mean(rsidat(:,:,:,b0inds(2:end)),4), M );

  end

else

  b0inds = find(bvals == 0);
  b0_vol_fwd = mgh2ctx( mean(rsidat(:,:,:,b0inds),4), M );

  if exist('rsidat_rev', 'var')
    b0inds_rev = find(bvals_rev == 0);
    b0_vol_rev = mgh2ctx( mean(rsidat_rev(:,:,:,b0inds_rev),4), M_rev );
  else
    fprintf('%s -- %s.m:    Only one b=0 volume found, so no B0 distortion correction...\n', datestr(now), mfilename);
    b0_vol_rev = [];
    params.B0DISCO = 0;
  end

end

fname_b0_fwd = fullfile(output_path, 'b0_vol_fwd.mgz');
QD_ctx_save_mgh( b0_vol_fwd, fname_b0_fwd );
if ~isempty(b0_vol_rev)
   fname_b0_rev = fullfile(output_path, 'b0_vol_rev.mgz');
   QD_ctx_save_mgh( b0_vol_rev, fname_b0_rev );
end


% Remove reverse PE polarity volume from data, if it exists
if integrated_rev_flag == 1 
   rsidat = rsidat(:,:,:,2:end);
   bvals = bvals(2:end);
   qmat = qmat(2:end,:);
end
QD_save_mgh( rsidat, fullfile(output_path,'DWI_vol_unprocessed_fwd.mgz'), M );
save(fullfile(output_path, 'bvals.mat'), 'bvals');
save(fullfile(output_path, 'qmat.mat'), 'qmat');

if exist('rsidat_rev', 'var')
  if integrated_rev_flag_rev == 1
    rsidat_rev = rsidat_rev(:,:,:,2:end);
    bvals_rev = bvals_rev(2:end);
    qmat_rev = qmat_rev(2:end,:);
  end
  QD_save_mgh( rsidat_rev, fullfile(output_path,'REV_DWI_vol_unprocessed_fwd.mgz'), M_rev );
  save(fullfile(output_path, 'REV_bvals.mat'), 'bvals_rev');
  save(fullfile(output_path, 'REV_qmat.mat'), 'qmat_rev');
end

uqvecs = unique(qmat, 'rows');
if size(uqvecs, 1) == 1 || any(isnan(uqvecs(:)))
   disp('WARNING: Bad qmat, disabling eddy current correction and motion correction');
   params.EddyCorrFlag = 0;
   params.MotionCorrFlag = 0;
end


% Compute and save averages of unprocessed DWI data 
ubvals = unique(bvals);
rsidat_averaged = zeros( size(rsidat,1), size(rsidat,2), size(rsidat,3), length(ubvals) );
for b = 1:length(ubvals)
    b_inds = find(bvals == ubvals(b));
    rsidat_averaged(:,:,:,b) = mean( rsidat(:,:,:,b_inds), 4 );
end
fprintf('%s -- %s.m:    Saving averages of unprocessed DWI data...\n', datestr(now), mfilename);
QD_save_mgh( rsidat_averaged, fullfile(output_path,'DWI_vol_unprocessed_averaged.mgz'), M );

if exist('rsidat_rev', 'var')
   ubvals_rev = unique(bvals_rev);
   rsidat_averaged_rev = zeros( size(rsidat_rev,1), size(rsidat_rev,2), size(rsidat_rev,3), length(ubvals_rev) );
   for b = 1:length(ubvals_rev)
     b_inds = find(bvals_rev == ubvals_rev(b));
     rsidat_averaged_rev(:,:,:,b) = mean( rsidat_rev(:,:,:,b_inds), 4 );
   end
   fprintf('%s -- %s.m:    Saving averages of unprocessed DWI data (reverse PE polarity)...\n', datestr(now), mfilename);
   QD_save_mgh( rsidat_averaged_rev, fullfile(output_path,'REV_DWI_vol_unprocessed_averaged.mgz'), M_rev );
end


% Collect some RSI protocol parameters
nFrames = length(bvals);
TE = dcminfo.EchoTime;
TR = dcminfo.RepetitionTime;
ModelName = dcminfo.ManufacturerModelName;

fprintf('\n\n');
fprintf('Protocol Parameters:\n');
fprintf('\t Scanner Model: %s - %s\n',Manufacturer,ModelName);
fprintf('\t EchoTime: %.1f ms\n',TE);
fprintf('\t RepetitionTime: %.1f ms\n',TR);
fprintf('\t Number of frames: %d \n',nFrames);
fprintf('\t b-values: %.1f mms-2\n',ubvals);
fprintf('\n\n');

if exist('bvals_rev', 'var')
   fprintf('Protocol Parameters (reverse PE polarity):\n');
   fprintf('\t Scanner Model: %s - %s\n',Manufacturer,ModelName);
   fprintf('\t EchoTime: %.1f ms\n',TE);
   fprintf('\t RepetitionTime: %.1f ms\n',TR);
   fprintf('\t Number of frames: %d \n', length(bvals_rev));
   fprintf('\t b-values: %.1f mms-2\n', ubvals_rev);
   fprintf('\n\n');
end

[rows, cols, slices, ~] = size(rsidat);
if isfield(dcminfo, 'AcquisitionMatrix')
  matrix = double(dcminfo.AcquisitionMatrix);
  matrix = matrix(matrix~=0);
else
  matrix = [rows cols]';
end
pixel_spacing = dcminfo.PixelSpacing;
FOV_row = rows*pixel_spacing(1);
real_vxl_size_row = FOV_row./matrix(1);

if exist('rsidat_rev', 'var') && (integrated_rev_flag == 0)
  [rows_rev, cols_rev, slices_rev, ~] = size(rsidat_rev);
  if isfield(dcminfo_rev, 'AcquisitionMatrix')
    matrix_rev = double(dcminfo_rev.AcquisitionMatrix);
    matrix_rev = matrix_rev(matrix_rev~=0);
  else
    matrix_rev = [rows_rev cols_rev]';
  end
  pixel_spacing_rev = dcminfo_rev.PixelSpacing;
  FOV_row_rev = rows_rev*pixel_spacing_rev(1);
  real_vxl_size_row_rev = FOV_row_rev./matrix_rev(1);
end


% Check if forward and reverse volumes are compatible
flag_compute_avgs = 1;
if exist('rsidat_rev', 'var') && (integrated_rev_flag == 0)

   if ~isequal(size(rsidat_averaged(:,:,:,1)), size(rsidat_averaged_rev(:,:,:,1)))
      fprintf('WARNING: Forward and reverse volumes are different sizes!\n');
      params.B0DISCO = 0;
      flag_compute_avgs = 0;
   end

   diff = abs(M(:,4) - M_rev(:,4));
   if any(diff > real_vxl_size_row)
      fprintf('WARNING: Scan coverage of forward and reverse volumes are different!\n');
      params.B0DISCO = 0;
      flag_compute_avgs = 0;
   end

end


%% ============================= T2 Data ============================== %%
if T2_flag ~= 0
    
  fprintf('%s -- %s.m:    Loading T2 data...\n',datestr(now),mfilename);
  fprintf('T2 series: %s\n', T2_path);
  [volT2, gwarpInfoT2, dcminfo_T2] = ReadDicomT2Data(T2_path); volT2 = volT2{1};
  dcminfo_T2 = dcminfo_T2(1);
  [rows_T2, cols_T2, slices_T2, ~] = size(volT2.imgs);
  matrix_T2 = double(dcminfo_T2.AcquisitionMatrix);
  matrix_T2 = matrix_T2(matrix_T2~=0);
  pixel_spacing_T2 = dcminfo_T2.PixelSpacing;
  FOV_row_T2 = rows_T2*pixel_spacing_T2(1);

  % Unwarp T2
  fprintf('%s -- %s.m:    Unwarping T2...\n',datestr(now),mfilename);
  if params.GradWarpFlag && isfield(gwarpInfoT2, 'gwtype') && isfield(gwarpInfoT2, 'unwarpflag') && isfield(gwarpInfoT2, 'isoctrflag')
    volT2 = ctx_unwarp_grad(volT2, gwarpInfoT2.gwtype, gwarpInfoT2.unwarpflag, gwarpInfoT2.isoctrflag);
  else
    fprintf('WARNING: T2 gradwarp correction disabled\n');
  end

  fprintf('%s -- %s.m:    Saving T2 data...\n',datestr(now),mfilename);
  QD_ctx_save_mgh( volT2,  fullfile(output_path, 'T2_corrected_GUW.mgz') );

else

  fprintf('%s -- %s.m:    Not using any T2 data...\n',datestr(now),mfilename);

end


%% ========================== Start Processing ========================== %%
% B0 distortion correction
if params.B0DISCO

   if ~exist('ubvals_rev', 'var')
     ubvals_rev = NaN;
   end

   if any(isnan(ubvals_rev)) || (numel(ubvals_rev)<2)
     fprintf('%s -- %s.m:    Performing B0 distortion correction (using RPG)...\n', datestr(now), mfilename);   
     disp_field = RPG_estimate_displacements(fname_b0_fwd, fname_b0_rev, params.RPG);
     ctx_rsidat = mgh2ctx(rsidat, M);
     ctx_rsidat_B0uw = apply_displacement_field(ctx_rsidat, disp_field);
     QD_ctx_save_mgh( ctx_rsidat_B0uw, fullfile(output_path, 'DWI_vol_B0disCo_only.mgz') );
     rsidat = ctx_rsidat_B0uw.imgs;
     
     if exist('rsidat_rev', 'var')
       disp_field_rev = disp_field;
       disp_field_rev.imgs = -disp_field.imgs;
       ctx_rsidat_rev = mgh2ctx(rsidat_rev, M_rev);
       ctx_rsidat_rev_B0uw = apply_displacement_field(ctx_rsidat_rev, disp_field_rev);
       QD_ctx_save_mgh( ctx_rsidat_rev_B0uw, fullfile(output_path, 'REV_DWI_vol_B0disCo_only.mgz') );
       rsidat_rev = ctx_rsidat_rev_B0uw.imgs;
     end
   end

   if exist('rsidat_averaged_rev', 'var') && (numel(ubvals_rev)>=2)
      fprintf('%s -- %s.m:    Performing B0 distortion correction (by registering normalized forward and reverse scans)...\n', datestr(now), mfilename);

      % For now, enforce identical b-values between "forward" and "reverse" scans
      bvals_common = intersect(ubvals, ubvals_rev);
      inds_fwd = ismember(ubvals, bvals_common);
      inds_rev = ismember(ubvals_rev, bvals_common);
      tmp_vol_fwd = rsidat_averaged(:,:,:,inds_fwd);
      tmp_vol_rev = rsidat_averaged_rev(:,:,:,inds_rev);

      ctx_rsidat_averaged = mgh2ctx(tmp_vol_fwd, M);
      ctx_rsidat_averaged_rev = mgh2ctx(tmp_vol_rev, M_rev);
      disp_field = estimate_B0shift_imreg(ctx_rsidat_averaged, ctx_rsidat_averaged_rev);
      disp_field_rev = disp_field;
      disp_field_rev.imgs = -disp_field.imgs;

      ctx_rsidat = mgh2ctx(rsidat, M);
      ctx_rsidat_rev = mgh2ctx(rsidat_rev, M_rev);
      ctx_rsidat_B0uw = apply_displacement_field(ctx_rsidat, disp_field);
      ctx_rsidat_rev_B0uw = apply_displacement_field(ctx_rsidat_rev, disp_field_rev);
      QD_ctx_save_mgh( ctx_rsidat_B0uw, fullfile(output_path, 'DWI_vol_B0disCo_only.mgz') );
      QD_ctx_save_mgh( ctx_rsidat_rev_B0uw, fullfile(output_path, 'REV_DWI_vol_B0disCo_only.mgz') );
      rsidat = ctx_rsidat_B0uw.imgs;
      rsidat_rev = ctx_rsidat_rev_B0uw.imgs;
   end

   QD_ctx_save_mgh( disp_field, fullfile(output_path, 'B0_corr_shift_map.mgz') );
      
end


% Other corrections: Eddy current, gradient warping, motion ---------------------------------------------
if exist('ubvals_rev', 'var') && (numel(ubvals_rev)>=2)
   fprintf('%s -- %s.m:    Performing other corrections on reverse PE polarity volumes...\n', datestr(now), mfilename);
   Preprocess_Diffusion(rsidat_rev, M_rev, qmat_rev, bvals_rev, gwinfo_rsi_rev, dcminfo_rev, output_path, isFOCUS, isARTPro, pedim, params);
   fname_corr_rev = sprintf('%s/DWI_vol_preprocessed.mgz', output_path);
   [rsidat_rev, M_rev] = QD_load_mgh(fname_corr_rev);
   qmat_rev = load(sprintf('%s/qmat_corrected.mat',output_path)); qmat_rev = qmat_rev.qmat_mc;
   b0inds_rev = find(prod(qmat_rev,2)==0);
   volb0_rev = mgh2ctx(mean(rsidat_rev(:,:,:,b0inds_rev),4), M_rev);
   volb0_rev.imgs(find(isnan(volb0_rev.imgs))) = 0;
   volb0_rev.imgs(find(isinf(volb0_rev.imgs))) = 0;
   movefile(fullfile(output_path,'DWI_vol_preprocessed.mgz'), fullfile(output_path,'REV_DWI_vol_preprocessed.mgz'));
   movefile(fullfile(output_path,'qmat_corrected.mat'), fullfile(output_path,'REV_qmat_corrected.mat'));
end

fprintf('%s -- %s.m:    Performing other corrections...\n', datestr(now), mfilename);
Preprocess_Diffusion(rsidat, M, qmat, bvals, gwinfo_rsi, dcminfo, output_path, isFOCUS, isARTPro, pedim, params);
fname_corr = sprintf('%s/DWI_vol_preprocessed.mgz', output_path);
[rsidat, M] = QD_load_mgh(fname_corr);
qmat = load(sprintf('%s/qmat_corrected.mat',output_path)); qmat = qmat.qmat_mc;
b0inds = find(prod(qmat,2)==0);
volb0 = mgh2ctx(mean(rsidat(:,:,:,b0inds),4), M);
volb0.imgs(find(isnan(volb0.imgs))) = 0;
volb0.imgs(find(isinf(volb0.imgs))) = 0;


% Segment prostate from T2 volume -----------------------------------------
if (params.ProstateSeg == 1) && (T2_flag ~= 0)

  % Check if prostate has already been contoured for another RSI series
  contour_exists = 0;
  path_date = fullfile(output_path, '../');
  cmd = sprintf('find %s -name "prostate_contour_T2_space.mgz"', path_date);
  [~, cmdout] = system(cmd);
  if ~isempty(cmdout)
    fprintf('%s -- %s.m:    Copying prostate contour from another RSI series...\n',datestr(now),mfilename);
    cmdout = split(cmdout);
    match_contour = find(~cellfun(@isempty, cmdout));
    path_contour = cmdout{match_contour};
    copyfile(path_contour, output_path);
    contour = QD_ctx_load_mgh(fullfile(output_path, 'prostate_contour_T2_space.mgz'));
    target = QD_ctx_load_mgh(fname_corr);
    contour_dwi_space = vol_resample(contour, target, eye(4));
    QD_ctx_save_mgh( contour_dwi_space, fullfile(output_path, 'prostate_contour_DWI_space.mgz') );
    contour_exists = 1;
  end

  if contour_exists == 0

    if strcmpi(params.ProstateSegVendor, 'cmig')
      fprintf('%s -- %s.m:    Segmenting prostate from T2 volume using CMIG software...\n',datestr(now),mfilename);
      contour = contour_prostate_cmig( fullfile(output_path, 'T2_corrected_GUW.mgz'), params.ProstateSegContainer );

    elseif strcmpi(params.ProstateSegVendor, 'cortechs')
      fprintf('%s -- %s.m:    Segmenting prostate from T2 volume using CorTechs'' software...\n',datestr(now),mfilename);
      contour = contour_prostate_cortechs( fullfile(output_path, 'T2_corrected_GUW.mgz'), params.ProstateSegContainer );

    else
      error('Prostate segmentation was enabled, but ProstateSegVendor parameter was not recognized');
    end

    if isempty(contour)
      params.ProstateSeg = 0;
    else
      cmd = ['mv ' fullfile(output_path,'T2_corrected_GUW_seg.mgz') ' ' fullfile(output_path,'prostate_contour_T2_space.mgz')];
      system(cmd);
      fprintf('%s -- %s.m:    Resampling prostate contour into DWI space...\n',datestr(now),mfilename);
      target = QD_ctx_load_mgh(fname_corr);
      contour_dwi_space = vol_resample(contour, target, eye(4));
      QD_ctx_save_mgh( contour_dwi_space, fullfile(output_path, 'prostate_contour_DWI_space.mgz') );
    end

  end

end


% Compute average DWI volume and save it -----------------------------------
DWI_vol_averaged = zeros( size(rsidat,1), size(rsidat,2), size(rsidat,3), length(ubvals) );
for b = 1:length(ubvals)
    b_inds = find(bvals == ubvals(b));
    DWI_vol_averaged(:,:,:,b) = mean( rsidat(:,:,:,b_inds), 4 );
end
fprintf('%s -- %s.m:    Saving hardware artifact-corrected DWI data...\n', datestr(now), mfilename);
QD_save_mgh( DWI_vol_averaged, fullfile(output_path,'DWI_vol_preprocessed_averaged.mgz'), M );

if exist('rsidat_rev', 'var')
  DWI_vol_averaged_rev = zeros( size(rsidat_rev,1), size(rsidat_rev,2), size(rsidat_rev,3), length(ubvals_rev) );
  for b = 1:length(ubvals_rev)
    b_inds = find(bvals_rev == ubvals_rev(b));
    DWI_vol_averaged_rev(:,:,:,b) = mean( rsidat_rev(:,:,:,b_inds), 4 );
  end
  fprintf('%s -- %s.m:    Saving hardware artifact-corrected DWI data...\n', datestr(now), mfilename);
  QD_save_mgh( DWI_vol_averaged_rev, fullfile(output_path,'REV_DWI_vol_preprocessed_averaged.mgz'), M_rev );
end


% Compute conventional ADC from average DWI volume w/out noise correction --
bvals_for_ADC = ubvals <= 1000;
conv_ADC_vol = compute_ADCs( DWI_vol_averaged(:,:,:,bvals_for_ADC), ubvals(bvals_for_ADC) );
fprintf('%s -- %s.m:    Saving conventional ADC maps...\n', datestr(now), mfilename);
QD_save_mgh( conv_ADC_vol, fullfile(output_path,'conventional_ADC_map.mgz'), M );

if exist('DWI_vol_averaged_rev', 'var') && (numel(ubvals_rev)>=2)
   bvals_for_ADC_rev = ubvals_rev <= 1000;
   conv_ADC_vol_rev = compute_ADCs( DWI_vol_averaged_rev(:,:,:,bvals_for_ADC_rev), ubvals_rev(bvals_for_ADC_rev) );
   fprintf('%s -- %s.m:    Saving conventional ADC maps...\n', datestr(now), mfilename);
   QD_save_mgh( conv_ADC_vol_rev, fullfile(output_path,'REV_conventional_ADC_map.mgz'), M_rev );
end


% Fit RSI model to data without noise correction ------------------------------
fprintf('%s -- %s.m:    Fitting RSI model to data without noise correction...\n', datestr(now), mfilename);
Sobs = zeros( length(ubvals), numel(DWI_vol_averaged(:,:,:,1)) );
for b = 1:length(ubvals)
    vol_b = DWI_vol_averaged(:,:,:,b);
    Sobs(b,:) = vol_b(:)';
end
temp_data = struct('Sobs', Sobs, 'bvals', ubvals', 'tol', 0);

nonneg_flag = 0;
[~, ~, C_mat_legacy] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);
nonneg_flag = 1;
[~, ~, C_mat] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);

C_vol_legacy = zeros(size(DWI_vol_averaged,1), size(DWI_vol_averaged,2), size(DWI_vol_averaged,3), size(C_mat,1));
C_vol = zeros(size(DWI_vol_averaged,1), size(DWI_vol_averaged,2), size(DWI_vol_averaged,3), size(C_mat,1));
for c = 1:size(C_mat, 1)
    C_vol_legacy(:,:,:,c) = reshape( C_mat_legacy(c,:), size(C_vol_legacy,1), size(C_vol_legacy,2), size(C_vol_legacy,3) );
    C_vol(:,:,:,c) = reshape( C_mat(c,:), size(C_vol,1), size(C_vol,2), size(C_vol,3) );
end
fprintf('%s -- %s.m:    Saving RSI signal-contribution (C) maps...\n', datestr(now), mfilename);
QD_save_mgh( C_vol, fullfile(output_path,'RSI_C_vol_noNC.mgz'), M );
QD_save_mgh( C_vol_legacy, fullfile(output_path,'RSI_C_vol_legacy_noNC.mgz'), M );

if exist('DWI_vol_averaged_rev', 'var') && (numel(ubvals_rev)>=2)
  Sobs = zeros( length(ubvals_rev), numel(DWI_vol_averaged_rev(:,:,:,1)) );
  for b = 1:length(ubvals_rev)
    vol_b = DWI_vol_averaged_rev(:,:,:,b);
    Sobs(b,:) = vol_b(:)';
  end
  temp_data = struct('Sobs', Sobs, 'bvals', ubvals_rev', 'tol', 0);

  nonneg_flag = 0;
  [~, ~, C_mat_legacy] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);
  nonneg_flag = 1;
  [~, ~, C_mat] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);

  C_vol_legacy_rev = zeros(size(DWI_vol_averaged_rev,1), size(DWI_vol_averaged_rev,2), size(DWI_vol_averaged_rev,3), size(C_mat,1));
  C_vol_rev = zeros(size(DWI_vol_averaged_rev,1), size(DWI_vol_averaged_rev,2), size(DWI_vol_averaged_rev,3), size(C_mat,1));
  for c = 1:size(C_mat, 1)
    C_vol_legacy_rev(:,:,:,c) = reshape( C_mat_legacy(c,:), size(C_vol_rev,1), size(C_vol_rev,2), size(C_vol_rev,3) );
    C_vol_rev(:,:,:,c) = reshape( C_mat(c,:), size(C_vol_rev,1), size(C_vol_rev,2), size(C_vol_rev,3) );
  end
  fprintf('%s -- %s.m:    Saving RSI signal-contribution (C) maps...\n', datestr(now), mfilename);
  QD_save_mgh( C_vol_rev, fullfile(output_path,'REV_RSI_C_vol_noNC.mgz'), M_rev );
  QD_save_mgh( C_vol_legacy_rev, fullfile(output_path,'REV_RSI_C_vol_legacy_noNC.mgz'), M_rev );
end


% Apply Anders' noise correction ------------------------------------------
if params.CorrectNoise == 1

  smf = 100;
  DWI_vol_nc = correct_noise_amd(rsidat, bvals, smf);
  DWI_avg_nc = zeros( size(DWI_vol_nc,1), size(DWI_vol_nc,2), size(DWI_vol_nc,3), length(ubvals) );
  for b = 1:length(ubvals)
    b_inds = find(bvals == ubvals(b));
    DWI_avg_nc(:,:,:,b) = mean( DWI_vol_nc(:,:,:,b_inds), 4 );
  end

  % Save noise-corrected volumes
  fprintf('%s -- %s.m:    Saving noise-corrected DWI data...\n', datestr(now), mfilename);
  QD_save_mgh( DWI_vol_nc, fullfile(output_path,'DWI_vol_NC.mgz'), M );
  QD_save_mgh( DWI_avg_nc, fullfile(output_path,'DWI_vol_NC_averaged.mgz'), M );

else

  DWI_vol_nc = rsidat;
  DWI_avg_nc = DWI_vol_averaged;

end

if exist('rsidat_rev', 'var') && (numel(ubvals_rev)>=2)
  if params.CorrectNoise == 1

    DWI_vol_nc_rev = correct_noise_amd(rsidat_rev, bvals_rev, smf);
    DWI_avg_nc_rev = zeros( size(DWI_vol_nc_rev,1), size(DWI_vol_nc_rev,2), size(DWI_vol_nc_rev,3), length(ubvals_rev) );
    for b = 1:length(ubvals_rev)
      b_inds = find(bvals_rev == ubvals_rev(b));
      DWI_avg_nc_rev(:,:,:,b) = mean( DWI_vol_nc_rev(:,:,:,b_inds), 4 );
    end

    % Save noise-corrected volumes
    fprintf('%s -- %s.m:    Saving noise-corrected DWI data...\n', datestr(now), mfilename);
    QD_save_mgh( DWI_vol_nc_rev, fullfile(output_path,'REV_DWI_vol_NC.mgz'), M_rev );
    QD_save_mgh( DWI_avg_nc_rev, fullfile(output_path,'REV_DWI_vol_NC_averaged.mgz'), M_rev );

  else

    DWI_vol_nc_rev = rsidat_rev;
    DWI_avg_nc_rev = DWI_vol_averaged_rev;

  end
end


% Fit RSI model to noise-corrected data ------------------------------------------------
if params.CorrectNoise == 1

  fprintf('%s -- %s.m:    Fitting RSI model to noise-corrected data...\n', datestr(now), mfilename);
  Sobs = zeros( length(ubvals), numel(DWI_avg_nc(:,:,:,1)) );
  for b = 1:length(ubvals)
    vol_b = DWI_avg_nc(:,:,:,b);
    Sobs(b,:) = vol_b(:)';
  end
  temp_data = struct('Sobs', Sobs, 'bvals', ubvals', 'tol', 0);
  obs_fwd = Sobs;

  nonneg_flag = 0;
  [~, ~, C_mat_legacy] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);
  nonneg_flag = 1;
  [pred_fwd, ~, C_mat] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);

  C_vol_legacy = zeros(size(DWI_vol_averaged,1), size(DWI_vol_averaged,2), size(DWI_vol_averaged,3), size(C_mat,1));
  C_vol = zeros(size(DWI_vol_averaged,1), size(DWI_vol_averaged,2), size(DWI_vol_averaged,3), size(C_mat,1));
  for c = 1:size(C_mat, 1)
    C_vol_legacy(:,:,:,c) = reshape( C_mat_legacy(c,:), size(C_vol_legacy,1), size(C_vol_legacy,2), size(C_vol_legacy,3) );
    C_vol(:,:,:,c) = reshape( C_mat(c,:), size(C_vol,1), size(C_vol,2), size(C_vol,3) );
  end
  fprintf('%s -- %s.m:    Saving RSI signal-contribution (C) maps...\n', datestr(now), mfilename);
  QD_save_mgh( C_vol, fullfile(output_path,'RSI_C_vol_NC.mgz'), M );
  QD_save_mgh( C_vol_legacy, fullfile(output_path,'RSI_C_vol_legacy_NC.mgz'), M );

  if exist('DWI_avg_nc_rev', 'var') && (numel(ubvals_rev)>=2)
    fprintf('%s -- %s.m:    Fitting RSI model to noise-corrected data...\n', datestr(now), mfilename);
    Sobs = zeros( length(ubvals_rev), numel(DWI_avg_nc_rev(:,:,:,1)) );
    for b = 1:length(ubvals_rev)
      vol_b = DWI_avg_nc_rev(:,:,:,b);
      Sobs(b,:) = vol_b(:)';
    end
    temp_data = struct('Sobs', Sobs, 'bvals', ubvals_rev', 'tol', 0);
    obs_rev = Sobs;

    nonneg_flag = 0;
    [~, ~, C_mat_legacy] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);
    nonneg_flag = 1;
    [pred_rev, ~, C_mat] = RK_breastRSI_ga_Spred_amd(params.ModelADCs, temp_data, nonneg_flag);

    C_vol_legacy_rev = zeros(size(DWI_vol_averaged_rev,1), size(DWI_vol_averaged_rev,2), size(DWI_vol_averaged_rev,3), size(C_mat,1));
    C_vol_rev = zeros(size(DWI_vol_averaged_rev,1), size(DWI_vol_averaged_rev,2), size(DWI_vol_averaged_rev,3), size(C_mat,1));
    for c = 1:size(C_mat, 1)
      C_vol_legacy_rev(:,:,:,c) = reshape( C_mat_legacy(c,:), size(C_vol_rev,1), size(C_vol_rev,2), size(C_vol_rev,3) );
      C_vol_rev(:,:,:,c) = reshape( C_mat(c,:), size(C_vol_rev,1), size(C_vol_rev,2), size(C_vol_rev,3) );
    end
    fprintf('%s -- %s.m:    Saving RSI signal-contribution (C) maps...\n', datestr(now), mfilename);
    QD_save_mgh( C_vol_rev, fullfile(output_path,'REV_RSI_C_vol_NC.mgz'), M_rev );
    QD_save_mgh( C_vol_legacy_rev, fullfile(output_path,'REV_RSI_C_vol_legacy_NC.mgz'), M_rev );
  end

end


% Get median signal within the bladder at b=0 ----------------------
vol_b0 = DWI_avg_nc(:,:,:,1);
C3_vol = C_vol(:,:,:,3);
urine_mask = get_bladder_mask(C3_vol, 0.75);
k = median( vol_b0(urine_mask) );
QD_save_mgh( urine_mask, fullfile(output_path,'urine_mask.mgz'), M );
save( fullfile(output_path, 'urine_norm_scalar.mat'), 'k' );

if exist('DWI_avg_nc_rev', 'var')
  vol_b0_rev = DWI_avg_nc_rev(:,:,:,1);
  C3_vol_rev = C_vol_rev(:,:,:,3);
  urine_mask_rev = get_bladder_mask(C3_vol_rev, 0.75);
  k = median( vol_b0_rev(urine_mask_rev) );
  QD_save_mgh( urine_mask_rev, fullfile(output_path,'REV_urine_mask.mgz'), M_rev );
  save( fullfile(output_path, 'REV_urine_norm_scalar.mat'), 'k' );
end


% Compute mb0 value ------------------------------------------------
if exist('contour_dwi_space', 'var')

  k = median( vol_b0(logical(contour_dwi_space.imgs)) );
  save( fullfile(output_path, 'mb0_scalar.mat'), 'k' );
  mb0 = k;

  if exist('vol_b0_rev', 'var')
     k = median( vol_b0_rev(logical(contour_dwi_space.imgs)) );
     save( fullfile(output_path, 'REV_mb0_scalar.mat'), 'k' );
     mb0_rev = k;
  end

end


% Compute RSIrs -----------------------------------------------------
if exist('mb0', 'var')

  RSIrs_fwd = 1000 .* (C_vol(:,:,:,1)./mb0);
  QD_save_mgh( RSIrs_fwd, fullfile(output_path,'RSIrs.mgz'), M );

  if exist('C_vol_rev', 'var') && exist('mb0_rev', 'var')
    RSIrs_rev = 1000 .* (C_vol_rev(:,:,:,1)./mb0_rev);
    QD_save_mgh( RSIrs_rev, fullfile(output_path,'REV_RSIrs.mgz'), M_rev );
  end

end


% Average forward and reverse volumes -------------------------------
% DWI
if exist('DWI_vol_averaged_rev', 'var') && (numel(ubvals)==numel(ubvals_rev)) && flag_compute_avgs
  mean_DWI_vol_averaged = (DWI_vol_averaged + DWI_vol_averaged_rev) ./ 2;
else
  mean_DWI_vol_averaged = DWI_vol_averaged;
end
QD_save_mgh( mean_DWI_vol_averaged, fullfile(output_path,'meanFwdRev_DWI_vol_preprocessed_averaged.mgz'), M );
mean_DWI_ctx_avg = mgh2ctx(mean_DWI_vol_averaged, M);

% DWI NC (only used for Anders' QC, for now)
if exist('DWI_avg_nc_rev', 'var') && (numel(ubvals)==numel(ubvals_rev)) && flag_compute_avgs
  mean_DWI_nc = (DWI_avg_nc + DWI_avg_nc_rev) ./ 2;
else
  mean_DWI_nc = DWI_avg_nc;
end

% ADC
if exist('conv_ADC_vol_rev', 'var') && flag_compute_avgs
  conv_ADC_vol_avg = (conv_ADC_vol + conv_ADC_vol_rev) ./ 2;
else
  conv_ADC_vol_avg = conv_ADC_vol;
end
QD_save_mgh( conv_ADC_vol_avg, fullfile(output_path,'meanFwdRev_conventional_ADC_map.mgz'), M );

% C maps
if exist('C_vol_rev', 'var') && (numel(ubvals)==numel(ubvals_rev)) && flag_compute_avgs
  C_vol_avg = (C_vol + C_vol_rev) ./ 2;
else
  C_vol_avg = C_vol;
end
QD_save_mgh( C_vol_avg, fullfile(output_path,'meanFwdRev_RSI_C_vol_NC.mgz'), M );

% RSIrs
if exist('RSIrs_fwd', 'var')
  if exist('RSIrs_rev', 'var') && (numel(ubvals)==numel(ubvals_rev)) && flag_compute_avgs
    RSIrs_avg = (RSIrs_fwd + RSIrs_rev) ./ 2;
  else
    RSIrs_avg = RSIrs_fwd;
  end
  QD_save_mgh( RSIrs_avg, fullfile(output_path,'meanFwdRev_RSIrs.mgz'), M );
  ctx_RSIrs_avg = mgh2ctx(RSIrs_avg, M);
end


% Identify prostate lesions ------------------------------------------------
if exist('contour_dwi_space', 'var')
  vol_lesions = lesion_id(contour_dwi_space, ctx_RSIrs_avg);
end


%% ========================== Create color overlay ============================= %%
if exist('RSIrs_avg', 'var') && T2_flag
   ctx_RSIrs_avg = vol_resample(ctx_RSIrs_avg, volT2, eye(4));
   color_range = [80 180];
   vol_overlay = create_color_overlay(ctx_RSIrs_avg, volT2, color_range);
end


%% ========================== Check protocol compliance ========================== %%
if params.CheckProtocolCompliance && T2_flag && exist('DWI_vol_averaged_rev', 'var') 

   scan_info_dwi_fwd.dcminfo = dcminfo;
   scan_info_dwi_fwd.dcminfo.FOV = FOV_row;
   scan_info_dwi_fwd.slices = volb0.dimd;
   scan_info_dwi_fwd.M = volb0.Mvxl2lph;
   scan_info_dwi_fwd.bvals = ubvals;

   scan_info_dwi_rev.dcminfo = dcminfo_rev;
   scan_info_dwi_rev.dcminfo.FOV = FOV_row_rev;
   scan_info_dwi_rev.slices = volb0_rev.dimd;
   scan_info_dwi_rev.M = volb0_rev.Mvxl2lph;
   scan_info_dwi_rev.bvals = ubvals_rev;

   scan_info_T2.dcminfo = dcminfo_T2;
   scan_info_T2.dcminfo.FOV = FOV_row_T2;
   scan_info_T2.slices = volT2.dimd;
   scan_info_T2.M = volT2.Mvxl2lph;
   
   load(which(params.ProtocolReference));
   problem_flag = 0;
   try
     problem_flag = check_RSI_protocol(output_path, [RSI_path '/..'], mandatory_series_list, ref_info_dwi_fwd, scan_info_dwi_fwd, ref_info_T2, scan_info_T2, ref_info_dwi_rev, scan_info_dwi_rev);
   catch ME
     fprintf('ERROR: %s\n', ME.message);
   end

   if problem_flag == 0
      delete(fullfile(output_path, 'protocol_compliance_report.txt'));
   end

end


%% ========================== Anders' QC metrics ========================== %%
% Check predicted vs. observed framewise dMRI forward and reverse => model fit and uncorrected motion / eddy currents
if exist('contour_dwi_space', 'var')
  vol_dist = bwdist(contour_dwi_space.imgs>0.5);
  vol_mask = vol_dist<5;
  ivec_mask = find(vol_mask); % Need version that considers anisotropic voxel size
  for i = 1:size(obs_fwd,1)
    corrvec_fit_fwd(i) = corr(colvec(pred_fwd(i,ivec_mask)),colvec(obs_fwd(i,ivec_mask)));
    if exist('DWI_avg_nc_rev', 'var') && (numel(ubvals)==numel(ubvals_rev))
      corrvec_fit_rev(i) = corr(colvec(pred_rev(i,ivec_mask)),colvec(obs_rev(i,ivec_mask)));
    else
      corrvec_fit_rev = [];
    end
  end
else
  corrvec_fit_fwd = [];
  corrvec_fit_rev = [];
end

% Check dMRI forward and reverse => distortion correction accuracy (should also compute the predicted Jacobian, and adjust intensities?)
if params.B0DISCO && exist('DWI_avg_nc_rev', 'var') && flag_compute_avgs && exist('contour_dwi_space', 'var')
  for fi = 1:size(tmp_vol_fwd,4)
    vol_fwd = tmp_vol_fwd(:,:,:,fi);
    vol_rev = tmp_vol_rev(:,:,:,fi);
    corrvec_fwd_rev(fi) = corr(colvec(vol_fwd(ivec_mask)),colvec(vol_rev(ivec_mask)));
  end
else
  corrvec_fwd_rev = [];
end

% Check T2 vs dMRI b=0, compile joint probability map (in vicinity of prostate) => motion between scans
if exist('contour_dwi_space', 'var')
  volT2_dwi_ctx = vol_resample(volT2,contour_dwi_space,eye(4));
  k = median(volT2_dwi_ctx.imgs(logical(contour_dwi_space.imgs)));
  volT2_dwi = volT2_dwi_ctx.imgs/k;
  volb0_dwi = mean_DWI_nc(:,:,:,1)/mb0;
  hv_t2 = linspace(0,10,100);
  hv_b0 = linspace(0,10,110);
  hc_t2_b0 = hist3([volT2_dwi(ivec_mask) volb0_dwi(ivec_mask)],{hv_t2 hv_b0});
else
  hc_t2_b0 = [];
end

% Should save corrvec_fit_fwd corrvec_fit_rev corrvec_fwd_rev hcc_t2_b0 to file
save(fullfile(output_path,'QC_metrics.mat'),'corrvec_fit_fwd','corrvec_fit_rev','corrvec_fwd_rev','hc_t2_b0');


%% ========================== Load DCE Data ========================== %%
if ~isempty(DCE_path{1})

  % Check if DCE data has already been loaded for another RSI series
  dce_exists = 0;
  path_date = fullfile(output_path, '../');
  cmd = sprintf('find %s -name "DCE.mgz"', path_date);
  [~, cmdout] = system(cmd);
  if ~isempty(cmdout)
    fprintf('%s -- %s.m:    Copying DCE data from another RSI series...\n',datestr(now),mfilename);
    cmdout = split(cmdout);
    match_dce = find(~cellfun(@isempty, cmdout));
    path_dce = cmdout{match_dce};
    copyfile(path_dce, output_path);
    [fpath_dce, ~, ~] = fileparts(path_dce);
    copyfile(fullfile(fpath_dce, 't.mat'), output_path);
    ctx_dce = QD_ctx_load_mgh(fullfile(output_path, 'DCE.mgz'));
    t = load(fullfile(output_path, 't.mat'));
    t = t.t;
    dce_exists = 1;
  end

  if dce_exists == 0
    fprintf('%s -- %s.m:    Loading DCE data...\n', datestr(now), mfilename);
    try
      path_parts = regexp(DCE_path{1}, filesep, 'split');
      acq_name = path_parts{end};
      acq_name = nixify(acq_name);
      fname_out_opt = fullfile(output_path, 'DCE.mgz');
      fname_out_t = fullfile(output_path, 't.mat');

      if length(DCE_path) > 1
	path_tmp_in = fullfile(output_path, 'tmp_dce');
	mkdir(path_tmp_in);
	fprintf('%s -- %s.m:    Consolidating DCE files...\n', datestr(now), mfilename);
	for i = 1:length(DCE_path)
	  cmd = sprintf('cp --backup=t ''%s''/* ''%s''/', DCE_path{i}, path_tmp_in);
	  system(cmd);
	end
	DCE_path_consolidated = path_tmp_in;
      else
	DCE_path_consolidated = DCE_path{1};
      end

      [ctx_dce, t] = load_dce_data(DCE_path_consolidated);
      QD_ctx_save_mgh(ctx_dce, fname_out_opt);
      save(fname_out_t, 't');

      if exist('path_tmp_in', 'var')
	rmdir(path_tmp_in, 's');
      end

    catch ME
      fprintf('%s\n', ME.message);
    end

  end

end


%% ========================== Write Out DICOMs ================================= %%
if params.WriteDICOMS == 1

dcm_write_dir = fullfile(output_path, 'DICOM_Output');
mkdir(dcm_write_dir);


% T2
if T2_flag ~= 0 && params.SelectDICOMS.T2
  fprintf('%s -- %s.m:    Writing T2 DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
  T2_dcm_label = 'RSI_anatomic_T2W';
  dcm_hdr_struct = dcminfo_T2;
  dcm_hdr_struct.SeriesDescription = T2_dcm_label;
  dcm_hdr_struct.SeriesNumber = 6969;
  outpath_T2 = fullfile(dcm_write_dir, T2_dcm_label);
  dicomwrite_cmig(volT2, outpath_T2, dcm_hdr_struct)
end


% DCE
if exist('ctx_dce', 'var') && params.SelectDICOMS.DCE
   fprintf('%s -- %s.m:    Writing DCE DICOMs ---------------------------------------------------\n', datestr(now), mfilename);

   dce_dcm_label = 'DCE_subtraction';
   outpath_dcm = fullfile(dcm_write_dir, dce_dcm_label);
   dce_dcm_seriesnum = 333;

   vol_enhance = mean(ctx_dce.imgs(:,:,:,2:10), 4) - ctx_dce.imgs(:,:,:,1);
   vol_enhance(vol_enhance<0) = 0;
   ctx_dce.imgs = vol_enhance;

   if exist('volT2', 'var')
     vol_dce_t2 = zeros(rows_T2, cols_T2, slices_T2, size(ctx_dce.imgs,4));
     tmp = ctx_dce;
     for i = 1:size(ctx_dce.imgs,4)
       tmp.imgs = ctx_dce.imgs(:,:,:,i);
       ctx_resamp = vol_resample(tmp, volT2, eye(4));
       vol_dce_t2(:,:,:,i) = ctx_resamp.imgs;
     end
     ctx2write = volT2;
     ctx2write.imgs = vol_dce_t2;
     dcm_hdr_struct = dcminfo_T2;
     dcm_hdr_struct.SeriesDescription = dce_dcm_label;
     dcm_hdr_struct.SeriesNumber = dce_dcm_seriesnum;
     dicomwrite_cmig(ctx2write, outpath_dcm, dcm_hdr_struct);
   else
     fname_dce = dir(DCE_path{1});
     fname_dce = fullfile(fname_dce(3).folder, fname_dce(3).name);
     dcm_hdr_struct = dicominfo(fname_dce);
     dcm_hdr_struct.SeriesDescription = dce_dcm_label;
     dcm_hdr_struct.SeriesNumber = dce_dcm_seriesnum;
     dicomwrite_cmig(ctx_dce, outpath_dcm, dcm_hdr_struct);
   end

end


% DWI
if params.SelectDICOMS.DWI
  fprintf('%s -- %s.m:    Writing DWI DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
  DWI_dcm_label = 'RSI_DWI_averages';
  DWI_dcm_seriesnum = 6970;
     
  for b = 1:size(mean_DWI_ctx_avg.imgs, 4)
    ctx2write = mean_DWI_ctx_avg;
    ctx2write.imgs = mean_DWI_ctx_avg.imgs(:,:,:,b);

    dwi_dcm_dir = fullfile(dcm_write_dir, DWI_dcm_label, ['b-value_' num2str(ubvals(b))]);
    dwi_label = [DWI_dcm_label '_b-value_' num2str(ubvals(b))];

    if exist('volT2', 'var')
      ctx2write = vol_resample(ctx2write, volT2, eye(4));
      dcm_hdr_struct = dcminfo_T2;
      dcm_hdr_struct.SeriesDescription = dwi_label;
      dcm_hdr_struct.SeriesNumber = DWI_dcm_seriesnum;
      dicomwrite_cmig(ctx2write, dwi_dcm_dir, dcm_hdr_struct);
    else 
      dcm_hdr_struct = dcminfo;
      dcm_hdr_struct.SeriesDescription = dwi_label;
      dcm_hdr_struct.SeriesNumber = DWI_dcm_seriesnum;
      dicomwrite_cmig(ctx2write, dwi_dcm_dir, dcm_hdr_struct);
    end

    DWI_dcm_seriesnum = DWI_dcm_seriesnum + 1;

  end
end


% ADC
if params.SelectDICOMS.ADC
  fprintf('%s -- %s.m:    Writing ADC DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
  ADC_dcm_label = 'ADC';
  ADC_dcm_seriesnum = 8008;
     
  ctx2write = mean_DWI_ctx_avg;
  ctx2write.imgs = 1e6 .*  conv_ADC_vol_avg;

  RescaleSlope = 1/1e6;
  RescaleIntercept = 0;
  WindowWidth = '';
  WindowCenter = '';

  if exist('volT2', 'var')
    ctx2write = vol_resample(ctx2write, volT2, eye(4));
    dcm_hdr_struct = dcminfo_T2;
    dcm_hdr_struct.SeriesDescription = ADC_dcm_label;
    dcm_hdr_struct.SeriesNumber = ADC_dcm_seriesnum;
    dcm_hdr_struct.RescaleSlope = RescaleSlope;
    dcm_hdr_struct.RescaleIntercept = RescaleIntercept;
    dcm_hdr_struct.WindowWidth = WindowWidth;
    dcm_hdr_struct.WindowCenter = WindowCenter;
    dicomwrite_cmig(ctx2write, fullfile(dcm_write_dir, ADC_dcm_label), dcm_hdr_struct);
  else 
    dcm_hdr_struct = dcminfo;
    dcm_hdr_struct.SeriesDescription = ADC_dcm_label;
    dcm_hdr_struct.SeriesNumber = ADC_dcm_seriesnum;
    dcm_hdr_struct.RescaleSlope = RescaleSlope;
    dcm_hdr_struct.RescaleIntercept = RescaleIntercept;
    dcm_hdr_struct.WindowWidth = WindowWidth;
    dcm_hdr_struct.WindowCenter = WindowCenter;
    dicomwrite_cmig(ctx2write, fullfile(dcm_write_dir, ADC_dcm_label), dcm_hdr_struct);
  end

end


% RSI C maps
if params.SelectDICOMS.RSICmaps
  fprintf('%s -- %s.m:    Writing RSI C map DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
  RSI_C_label = 'RSI_C';
  RSI_C_seriesnum = 4200;
  
  ctx = mean_DWI_ctx_avg;
  for c = 1:size(C_vol_avg, 4)
    ctx2write = ctx;
    ctx2write.imgs = C_vol_avg(:,:,:,c);

    rsi_dcm_dir = fullfile(dcm_write_dir, [RSI_C_label '_maps'], [RSI_C_label num2str(c)]);
    rsi_label = [RSI_C_label num2str(c)];

    if exist('volT2', 'var')
      ctx2write = vol_resample(ctx2write, volT2, eye(4));
      dcm_hdr_struct = dcminfo_T2;
      dcm_hdr_struct.SeriesDescription = rsi_label;
      dcm_hdr_struct.SeriesNumber = RSI_C_seriesnum;
      dicomwrite_cmig(ctx2write, rsi_dcm_dir, dcm_hdr_struct);
    else
      dcm_hdr_struct = dcminfo;
      dcm_hdr_struct.SeriesDescription = rsi_label;
      dcm_hdr_struct.SeriesNumber = RSI_C_seriesnum;
      dicomwrite_cmig(ctx2write, rsi_dcm_dir, dcm_hdr_struct);
    end

    RSI_C_seriesnum = RSI_C_seriesnum + 1;

  end
end


% RSI biomarker
if exist('mb0', 'var') && params.SelectDICOMS.RSIrs
   fprintf('%s -- %s.m:    Writing RSI biomarker DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
   label = 'RSIrs_Experimental';
   seriesnum = 666;

   ctx2write = mean_DWI_ctx_avg;
   ctx2write.imgs = RSIrs_avg;
   ctx2write = vol_resample(ctx2write, volT2, eye(4));

   dcm_hdr_struct = dcminfo_T2;
   dcm_hdr_struct.SeriesDescription = label;
   dcm_hdr_struct.SeriesNumber = seriesnum;
   dcm_hdr_struct.WindowWidth = 100;
   dcm_hdr_struct.WindowCenter = 130;
   dicomwrite_cmig(ctx2write, fullfile(dcm_write_dir, label), dcm_hdr_struct);
end


% RSIrs color overlay
if exist('vol_overlay', 'var') && params.SelectDICOMS.RSIrsOverlay
   fprintf('%s -- %s.m:    Writing RSI biomarker color overlay ---------------------------------------------------\n', datestr(now), mfilename);
   label = 'RSIrs_Overlay_Experimental';
   seriesnum = 667;

   ctx2write = vol_overlay;
   ctx2write.imgs = 255 .* ctx2write.imgs;
   dcm_hdr_struct = dcminfo_T2;
   dcm_hdr_struct.ColorType = 'truecolor';
   dcm_hdr_struct.SamplesPerPixel = 3;
   dcm_hdr_struct.PhotometricInterpretation = 'RGB';
   dcm_hdr_struct.SmallestImagePixelValue = 0;
   dcm_hdr_struct.LargestImagePixelValue = 255;
   dcm_hdr_struct.BitDepth = 8;
   dicomwrite_cmig(ctx2write, fullfile(dcm_write_dir, label), dcm_hdr_struct);
end


% Prostate mask RT
if exist('contour_dwi_space', 'var') && params.SelectDICOMS.ProstateAutoSeg_RT
   fprintf('%s -- %s.m:    Writing prostate mask RTstruct DICOM ---------------------------------------------------\n', datestr(now), mfilename);
   segSTRUCT.number = 1; 
   segSTRUCT.name = 'Prostate_Autoseg'; 

   segSTRUCT.seg = contour.imgs; 

   if exist(fullfile(dcm_write_dir, 'RSIrs_Experimental'), 'dir')
     path_ref = fullfile(dcm_write_dir, 'RSIrs_Experimental');
   elseif exist(fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1'), 'dir')
     path_ref = fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1');
   elseif exist(fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0'), 'dir')
     path_ref = fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0');
   end
   RK_write_segSTRUCT(segSTRUCT, path_ref, fullfile(dcm_write_dir, 'Prostate_Mask_RTst_DWI'), 'Prostate_Mask_AutoSeg', 999); 

end


% Prostate mask SEG
if exist('contour_dwi_space', 'var') && params.SelectDICOMS.ProstateAutoSeg_SEG
   fprintf('%s -- %s.m:    Writing prostate mask SEG DICOM ---------------------------------------------------\n', datestr(now), mfilename);

   if exist(fullfile(dcm_write_dir, 'RSIrs_Experimental'), 'dir')
     path_ref = fullfile(dcm_write_dir, 'RSIrs_Experimental');
   elseif exist(fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1'), 'dir')
     path_ref = fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1');
   elseif exist(fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0'), 'dir')
     path_ref = fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0');
   end

   seg_fields.data = contour;
   seg_fields.source_images = path_ref;

   seg_fields.algorithm.name = 'Cortechs AI Segmentation';
   seg_fields.algorithm.version = 'v69';
   seg_fields.algorithm.family = 'ArtificialIntelligence';

   seg_fields.description.label = 'Prostate';
   seg_fields.description.type = 'Organ';
   seg_fields.description.tracking_id = 'Cortechs Prostate Autoseg';

   seg_fields.metadata.series_description = 'SEG_prostate';
   seg_fields.metadata.series_number = 808;
   seg_fields.metadata.instance_number = 1;
   seg_fields.metadata.manufacturer = 'Cortechs.AI';
   seg_fields.metadata.manufacturer_model_name = 'Cortechs Autoseg';
   seg_fields.metadata.software_versions = 'v69';
   seg_fields.metadata.device_serial_number = '69';

   write_dicom_seg(seg_fields, fullfile(dcm_write_dir, 'Prostate_Mask_SEG_DWI'), params.PythonVEnv);

end


% Lesion masks RT
if exist('vol_lesions', 'var')
  if ~isempty(vol_lesions) && params.SelectDICOMS.LesionAutoSeg_RT
    fprintf('%s -- %s.m:    Writing lesion mask RTstruct DICOM ---------------------------------------------------\n', datestr(now), mfilename);

    for i = 1:length(vol_lesions)
      segSTRUCT(i).number = i; 
      segSTRUCT(i).name = ['Lesion_Autoseg_' num2str(i)]; 
      ctx_lesion = mgh2ctx(vol_lesions{i}, M);
      ctx_lesion_t2 = vol_resample(ctx_lesion, volT2, eye(4));
      segSTRUCT(i).seg = ctx_lesion_t2.imgs; 
    end

    if exist(fullfile(dcm_write_dir, 'RSIrs_Experimental'), 'dir')
      path_ref = fullfile(dcm_write_dir, 'RSIrs_Experimental');
    elseif exist(fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1'), 'dir')
      path_ref = fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1');
    elseif exist(fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0'), 'dir')
      path_ref = fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0');
    end
    RK_write_segSTRUCT(segSTRUCT, path_ref, fullfile(dcm_write_dir, 'Lesion_Mask_RTst_DWI'), 'Lesion_Mask_AutoSeg', 1066); 
    
  end
end


% Lesion masks SEG
if exist('vol_lesions', 'var')
  if ~isempty(vol_lesions) && params.SelectDICOMS.LesionAutoSeg_SEG
    fprintf('%s -- %s.m:    Writing lesion mask SEG DICOM ---------------------------------------------------\n', datestr(now), mfilename);

    if exist(fullfile(dcm_write_dir, 'RSIrs_Experimental'), 'dir')
      path_ref = fullfile(dcm_write_dir, 'RSIrs_Experimental');
    elseif exist(fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1'), 'dir')
      path_ref = fullfile(dcm_write_dir, 'RSI_C_maps', 'RSI_C1');
    elseif exist(fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0'), 'dir')
      path_ref = fullfile(dcm_write_dir, 'RSI_DWI_averages', 'b-value_0');
    end

    lesion_series_num = 809;
    seg_fields.metadata.series_description = 'SEG_lesions';
    seg_fields.metadata.series_number = lesion_series_num;
    seg_fields.metadata.instance_number = 1;
    seg_fields.metadata.manufacturer = 'CMIG';
    seg_fields.metadata.manufacturer_model_name = 'CMIG RSI AI';
    seg_fields.metadata.software_versions = 'v69';
    seg_fields.metadata.device_serial_number = '69';

    seg_fields.source_images = path_ref;
    
    seg_fields.algorithm.name = 'GMIG RSI AI';
    seg_fields.algorithm.version = 'v69';
    seg_fields.algorithm.family = 'ArtificialIntelligence';

    seg_fields.data = volT2;
    seg_fields.data.imgs = [];
    for i = 1:length(vol_lesions)
      ctx_lesion = mgh2ctx(vol_lesions{i}, M);
      ctx_lesion_t2 = vol_resample(ctx_lesion, volT2, eye(4));
      seg_fields.data.imgs = cat(4, seg_fields.data.imgs, ctx_lesion_t2.imgs);
      seg_fields.description(i).type = 'Abnormal';
      seg_fields.description(i).tracking_id = 'RSI-bright lesion';
      seg_fields.description(i).label = ['Lesion_' num2str(i)];
    end

    write_dicom_seg(seg_fields, fullfile(dcm_write_dir, 'Lesion_Mask_SEG_DWI'), params.PythonVEnv);
    
  end
end

end


%% ========================== Finish ===================== %%

fprintf('%s -- %s.m:    Cleaning up...\n',datestr(now),mfilename);

if params.DebugFlag == 0
    
   if exist(fullfile(output_path, 'avgEIP.mgz'), 'file') ~= 0
     delete( fullfile(output_path, 'avgEIP.mgz') );
     delete( fullfile(output_path, 'difEIP.mgz') );
     delete( fullfile(output_path, 'fJ.mgz') );
     delete( fullfile(output_path, 'rJ.mgz') );
   end

end

status = 0;

fprintf('%s -- %s.m:    Finished!\n\n\n',datestr(now),mfilename);

end
