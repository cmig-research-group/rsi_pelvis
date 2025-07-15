function [prostate_mask, prostate_detector_output] = contour_prostate_cortechs(path_to_axT2_mgz, container)

[filepath, name, ~] = fileparts(path_to_axT2_mgz);
name_seg = [name '_seg.mgz'];

mkdir(fullfile(filepath, 'seg'));
copyfile( path_to_axT2_mgz, fullfile(filepath, 'seg') );

if ~exist('container', 'var')
  container = 'singularity';
end
disp(['Container: ' container]);

container = lower(container);

if ~isdeployed
  if strcmp(container, 'docker')
    cmd = ['sudo docker run -i --rm -u 0:0 -v ' fullfile(filepath, 'seg') ':/input:ro -v ' fullfile(filepath, 'seg') ':/output 696021503801.dkr.ecr.us-east-1.amazonaws.com/prostate-segmentation:dev'];
  else
    mkdir(fullfile(filepath, 'seg_out'));
    cmd = ['singularity run -B ' fullfile(filepath, 'seg') ':/app/input -B ' fullfile(filepath, 'seg_out') ':/app/output /space/bil-syn01/1/cmig_bil/Cortechs/prostate_seg.sif'];
  end
  disp(['Command: ' cmd]);
  system(cmd);

else
  cmd = ['sudo ' which('call_docker.sh') ' '  '''' fullfile(filepath, 'seg') ''''];
  disp(['Command: ' cmd]);
  system(cmd);

end

if strcmp(container, 'docker')
  copyfile( fullfile(filepath, 'seg', name_seg), filepath );
  rmdir(fullfile(filepath, 'seg'), 's');
else
  copyfile( fullfile(filepath, 'seg_out', name_seg), filepath );
  rmdir(fullfile(filepath, 'seg'), 's');
  rmdir(fullfile(filepath, 'seg_out'), 's');
end

contour = QD_ctx_load_mgh(fullfile(filepath, name_seg));
volT2 = QD_ctx_load_mgh(path_to_axT2_mgz);

% Check to see if prostate segmentation failed (created mask of all zeros)
if ~any(contour.imgs(:)) % Prostate mask is all zeros

  disp('WARNING: Prostate segmentation returned blank mask');
  delete( fullfile(filepath, name_seg) );
  prostate_mask = [];
  prostate_detector_output = 0;

else

  prostate_detector_output = 1;

  % Make sure segmentation has same slice ordering as T2 volume
  slice_ordering_T2 = round(volT2.Mvxl2lph(3,3));
  slice_ordering_seg = round(contour.Mvxl2lph(3,3));
  if slice_ordering_seg ~= slice_ordering_T2
    contour_flipped = volT2;
    contour_flipped.imgs = flip(contour.imgs, 3);
    contour = contour_flipped;
    QD_ctx_save_mgh( contour, fullfile(filepath, name_seg) );
  end
  prostate_mask = contour;

end

end
