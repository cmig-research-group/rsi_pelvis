function [prostate_mask, prostate_detector_output] = contour_prostate_cmig(path_to_axT2_mgz, container)

% Unset LD_PRELOAD environment variable for system calls to singularity
% Otherwise need to parse glibc warnings
orig_LD_PRELOAD = getenv('LD_PRELOAD');
setenv('LD_PRELOAD', '');

[filepath, name, ext] = fileparts(path_to_axT2_mgz);
path_input = fullfile(filepath, 'seg');
path_output = filepath;
mkdir(path_input);


% First, check to see if there is a prostate in the volume
container_path_in = sprintf('%s/%s%s', path_output, name, ext);

if strcmpi(container, 'docker')
     cmd = sprintf('sudo docker run --ipc="host" --mount type=bind,src=%s,dst=%s --entrypoint=/app/miniconda3/bin/conda localhost/autoseg_prostate run -n nnUNet python3 -Wignore /app/3D_inference_prostate_detector.py %s', path_output, path_output, container_path_in);

elseif strcmpi(container, 'singularity')
  path_sif = '/space/bil-syn01/1/cmig_bil/containers/autoseg_prostate/autoseg_prostate.sif';
  path_tmp = fullfile(path_output, 'tmp');
  mkdir(path_tmp);
  cmd = sprintf('singularity exec -B %s:%s -B %s:/app/tmp %s python3 -Wignore /app/3D_inference_prostate_detector.py %s', path_output, path_output, path_tmp, path_sif, container_path_in);
end

disp(['Command: ' cmd]);
[status, cmdout] = system(cmd);
cmdout = split(cmdout);
match_notempty = find(~cellfun(@isempty, cmdout));
has_prostate = str2double(cmdout{match_notempty});

if ~has_prostate
  disp('WARNING: Patient may not have prostate');
  prostate_detector_output = 0;
else
  prostate_detector_output = 1;
end


% Segment prostate from volume
fname_nifti = fullfile(path_input, 'prostate_900_0000.nii.gz');

ctx_t2 = QD_ctx_load_mgh(path_to_axT2_mgz);
ctx_save_nifti(ctx_t2, fname_nifti);

if strcmpi(container, 'docker')
  cmd = sprintf('sudo docker run --ipc="host" --mount type=bind,src=%s,dst=/data_in --mount type=bind,src=%s,dst=/data_out localhost/autoseg_prostate', path_input, path_output);

elseif strcmpi(container, 'singularity')
  path_sif = '/space/bil-syn01/1/cmig_bil/containers/autoseg_prostate/autoseg_prostate.sif';
  path_tmp = fullfile(filepath, 'tmp');
  mkdir(path_tmp);
  cmd = sprintf('singularity run -B %s:/data_in -B %s:/data_out_predict -B %s:/data_out %s', path_input, path_tmp, path_output, path_sif); 
end
disp(['Command: ' cmd]);
system(cmd);


rmdir(path_input, 's');
if exist('path_tmp', 'var')
  rmdir(path_tmp, 's');
end

vol_seg_nifti = niftiread(fullfile(filepath, 'prostate_900.nii.gz'));
nifti_info = niftiinfo(fullfile(filepath, 'prostate_900.nii.gz'));

% Check to see if prostate segmentation failed (created mask of all zeros)
if ~any(vol_seg_nifti(:)) % Prostate mask is all zeros

  disp('WARNING: Prostate segmentation returned blank mask');
  delete( fullfile(filepath, 'prostate_900.nii.gz') );
  prostate_mask = [];

else

  % Make sure segmentation has same slice ordering as T2 volume
  M_RAS = nifti_info.Transform.T';
  M_LPH = M_RAS_TO_LPH * M_RAS;

  slice_ordering_T2 = round(ctx_t2.Mvxl2lph(3,3));
  slice_ordering_seg = round(M_LPH(3,3));
  if slice_ordering_seg ~= slice_ordering_T2
    contour_flipped = ctx_t2;
    contour_flipped.imgs = flip(vol_seg_nifti, 3);
    contour = contour_flipped;
  else
    contour = ctx_t2;
    contour.imgs = vol_seg_nifti;
  end

  contour.imgs = double(contour.imgs);
  QD_ctx_save_mgh( contour, fullfile(filepath, [name '_seg.mgz']) );
  prostate_mask = contour;

end

delete(fullfile(filepath, 'prostate_900.nii.gz'));
delete(fullfile(filepath, '*.json'));

% Return LD_PRELOAD to its original value
setenv('LD_PRELOAD', orig_LD_PRELOAD);

end
