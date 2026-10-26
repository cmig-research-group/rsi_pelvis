function [ctx_T2, dcminfo_T2, implant_detector_output, prostate_detector_output, contour_prostate, contour_pz, contour_urethra] = process_T2(path_dcm_dir, path_output, params)

% NOTE: Hip implant and prostate detectors are currently disabled due to unreliable performance
%       implant_detector_output is hard-coded to 0 (assumes no hip implant)
%       prostate_detector_output is hard-coded to 1 (assumes prostate is present)


% Set default values for parameters
defaults.ProstateSeg = 1;
defaults.UrethraSeg = 1;
defaults.ProstateSegContainer = 'singularity';
defaults.PathSingularitySif = '/space/bil-nas01/1/cmig_bil/containers/autoseg_prostate/autoseg_prostate.sif';
defaults.WriteDICOMS = 1;

% Load user-defined parameters
switch class(params)
  case 'struct' % Already in MATLAB 
  case {'char', 'string'}
    [~, ~, ext] = fileparts(params);
    if strcmp(ext, '.mat')
      load(params);
    else
      params = eval_file(params);
      fprintf('\n\n');
    end
  otherwise
    fprintf('WARNING: Unable to parse parameter file, using default options\n');
    params = struct;
end

% Insert default value for missing parameters
fields_default = fieldnames(defaults);
for i = 1:length(fields_default)
  if ~isfield(params, fields_default{i})
    params.(fields_default{i}) = defaults.(fields_default{i});
  end
end

 
if ~exist(path_output, 'dir')
  mkdir(path_output);
end


% Unset LD_PRELOAD environment variable for system calls to singularity
% Otherwise need to parse glibc warnings
orig_LD_PRELOAD = getenv('LD_PRELOAD');
setenv('LD_PRELOAD', '');


% Load T2 volume -------------------------------------------------------
fprintf('%s -- %s.m:    Loading T2 data...\n', datestr(now), mfilename);
[ctx_T2, dcminfo_T2] = QD_read_dicomdir(path_dcm_dir);
dcminfo_T2 = dcminfo_T2(1).fields;
[gwarpInfoT2, ~] = mmil_get_gradwarpinfo(dcminfo_T2);


% Unwarp T2 ------------------------------------------------------------
fprintf('%s -- %s.m:    Unwarping T2...\n',datestr(now),mfilename);
if isfield(gwarpInfoT2, 'gwtype') && isfield(gwarpInfoT2, 'unwarpflag') && isfield(gwarpInfoT2, 'isoctrflag')
  ctx_T2 = ctx_unwarp_grad(ctx_T2, gwarpInfoT2.gwtype, gwarpInfoT2.unwarpflag, gwarpInfoT2.isoctrflag);
else
  fprintf('WARNING: T2 gradwarp correction disabled\n');
end

fprintf('%s -- %s.m:    Saving T2 data...\n',datestr(now),mfilename);
QD_ctx_save_mgh( ctx_T2,  fullfile(path_output, 'T2_corrected_GUW.mgz') );


% Check for hip implant -------------------------------------------------
% Disable for now, until we have a more reliable solution
% implant_detector_output = detect_hip_implant( fullfile(path_output, 'T2_corrected_GUW.mgz'), params.ProstateSegContainer );
implant_detector_output = 0;
if implant_detector_output
  disp('WARNING: Patient may have hip implant');
end


% Segment prostate from T2 volume -----------------------------------------
if params.ProstateSeg
  fprintf('%s -- %s.m:    Segmenting prostate from T2 volume using CMIG software...\n',datestr(now),mfilename);
  [contour_prostate, prostate_detector_output] = contour_prostate_cmig( fullfile(path_output, 'T2_corrected_GUW.mgz'), params );

  if ~isempty(contour_prostate)
    cmd = ['mv ' fullfile(path_output,'T2_corrected_GUW_seg.mgz') ' ' fullfile(path_output,'prostate_contour_T2_space.mgz')];
    system(cmd);
  end

else
  contour_prostate = [];
  prostate_detector_output = [];
end


% Segment PZ from T2 volume -----------------------------------------
if params.ProstateSeg
  fprintf('%s -- %s.m:    Segmenting PZ from T2 volume using CMIG software...\n',datestr(now),mfilename);
  contour_pz = contour_pz_cmig( fullfile(path_output, 'T2_corrected_GUW.mgz'), params );

  if ~isempty(contour_pz)
    cmd = ['mv ' fullfile(path_output,'T2_corrected_GUW_seg.mgz') ' ' fullfile(path_output,'pz_contour_T2_space.mgz')];
    system(cmd);
  end

else
  contour_pz = [];
end


% Segment urethra from T2 volume -----------------------------------------
if params.UrethraSeg
  fprintf('%s -- %s.m:    Segmenting urethra from T2 volume using CMIG software...\n',datestr(now),mfilename);
  contour_urethra = contour_urethra_cmig( fullfile(path_output, 'T2_corrected_GUW.mgz'), params );

  if ~isempty(contour_urethra)
    cmd = ['mv ' fullfile(path_output,'T2_corrected_GUW_seg.mgz') ' ' fullfile(path_output,'urethra_contour_T2_space.mgz')];
    system(cmd);
  end

else 
  contour_urethra = [];
end


% Return LD_PRELOAD to its original value
setenv('LD_PRELOAD', orig_LD_PRELOAD);


% Write DICOM files ------------------------------------------------------
if ~params.WriteDICOMS
  return
end

dcm_write_dir = fullfile(path_output, 'DICOM_Output');
mkdir(dcm_write_dir);

% T2
fprintf('%s -- %s.m:    Writing T2 DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
dcm_hdr_struct = dcminfo_T2;
dcm_hdr_struct.SeriesDescription = sprintf('PPro_%s', dcm_hdr_struct.SeriesDescription);
sernum = num2str(dcm_hdr_struct.SeriesNumber);
sernum_new = str2double([sernum '00']);
dcm_hdr_struct.SeriesNumber = sernum_new;
outpath_T2 = fullfile(dcm_write_dir, 'PPro_T2W');
dicomwrite_cmig(ctx_T2, outpath_T2, dcm_hdr_struct);

% RT DICOM contours
struct_ind = 1;
if ~isempty(contour_prostate) 
  segSTRUCT(struct_ind).number = struct_ind; 
  segSTRUCT(struct_ind).name = 'Prostate_Autoseg'; 
  segSTRUCT(struct_ind).seg = contour_prostate.imgs; 
  struct_ind = struct_ind + 1;
end
if ~isempty(contour_pz) 
  segSTRUCT(struct_ind).number = struct_ind; 
  segSTRUCT(struct_ind).name = 'PZ_Autoseg'; 
  segSTRUCT(struct_ind).seg = contour_pz.imgs; 
  struct_ind = struct_ind + 1;
end
if ~isempty(contour_urethra)
  segSTRUCT(struct_ind).number = struct_ind;
  segSTRUCT(struct_ind).name = 'Urethra_Autoseg';
  segSTRUCT(struct_ind).seg = contour_urethra.imgs;
  struct_ind = struct_ind + 1;
end

path_ref = outpath_T2;

if exist('segSTRUCT', 'var')
  fprintf('%s -- %s.m:    Writing RTstruct DICOMs ---------------------------------------------------\n', datestr(now), mfilename);
  RK_write_segSTRUCT(segSTRUCT, path_ref, fullfile(dcm_write_dir, 'RT_DICOM'), 'PPro_contours', sernum_new+1); 
end


end
