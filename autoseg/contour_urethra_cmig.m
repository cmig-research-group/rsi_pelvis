function urethra_mask = contour_urethra_cmig(path_to_axT2_mgz, container, params)

[filepath, name, ~] = fileparts(path_to_axT2_mgz);
path_input = fullfile(filepath, 'seg');
path_output = filepath;
mkdir(path_input);
fname_nifti = fullfile(path_input, 'prostate_900_0000.nii.gz');

ctx_t2 = QD_ctx_load_mgh(path_to_axT2_mgz);
ctx_save_nifti(ctx_t2, fname_nifti)

if strcmpi(container, 'docker')
     cmd = sprintf('sudo docker run --ipc="host" --mount type=bind,src=%s,dst=/data_in --mount type=bind,src=%s,dst=/data_out --entrypoint=/app/miniconda3/bin/conda localhost/autoseg_prostate run -n nnUNet /bin/bash -c /app/run_seg_urethra.sh', path_input, path_output);

elseif strcmpi(container, 'singularity')
  path_sif = '/space/bil-syn01/1/cmig_bil/containers/autoseg_prostate/autoseg_prostate.sif';
  cmd = sprintf('singularity exec -B %s:/data_in -B %s:/data_out %s /app/miniconda3/bin/conda run -n nnUNet /bin/bash -c /app/run_seg_urethra.sh', path_input, path_output, path_sif);
end
disp(['Command: ' cmd]);
system(cmd);

rmdir(path_input, 's');
if exist('path_tmp', 'var')
  rmdir(path_tmp, 's');
end

vol_seg_nifti = niftiread(fullfile(filepath, 'prostate_900.nii.gz'));
nifti_info = niftiinfo(fullfile(filepath, 'prostate_900.nii.gz'));

% Check to see if urethra segmentation failed (created mask of all zeros)
if ~any(vol_seg_nifti(:)) % Urethra mask is all zeros

  disp('WARNING: Urethra segmentation returned blank mask');
  delete( fullfile(filepath, 'prostate_900.nii.gz') );
  urethra_mask = [];

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
  urethra_mask = contour;

end

delete(fullfile(filepath, 'prostate_900.nii.gz'));
delete(fullfile(filepath, '*.json'));

end
