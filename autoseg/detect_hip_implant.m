function implant_detector_output = detect_hip_implant(path_to_axT2_mgz, container)

% Unset LD_PRELOAD environment variable for system calls to singularity
% Otherwise need to parse glibc warnings
orig_LD_PRELOAD = getenv('LD_PRELOAD');
setenv('LD_PRELOAD', '');

[filepath, name, ext] = fileparts(path_to_axT2_mgz);

container_path_in = sprintf('%s/%s%s', filepath, name, ext);

if strcmpi(container, 'docker')
  cmd = sprintf('sudo docker run --ipc="host" --mount type=bind,src=%s,dst=%s --entrypoint=/app/miniconda3/bin/conda localhost/autoseg_prostate run -n nnUNet python3 -Wignore /app/3D_inference_hip_implant_detector.py %s', filepath, filepath, container_path_in);

elseif strcmpi(container, 'singularity')
  path_sif = '/space/bil-syn01/1/cmig_bil/containers/autoseg_prostate/autoseg_prostate.sif';
  path_tmp = fullfile(filepath, 'tmp');
  mkdir(path_tmp);
  cmd = sprintf('singularity exec -B %s:%s -B %s:/app/tmp %s python3 -Wignore /app/3D_inference_hip_implant_detector.py %s', filepath, filepath, path_tmp, path_sif, container_path_in);
end

disp(['Command: ' cmd]);
[status, cmdout] = system(cmd);
if status ~= 0
  fprintf('Error encountered during hip implant detection\n')
  implant_detector_output = [];
else
  cmdout = split(cmdout);
  match_notempty = find(~cellfun(@isempty, cmdout));
  implant_detector_output = str2double(cmdout{match_notempty});
end

if exist('path_tmp', 'var')
  rmdir(path_tmp, 's');
end

% Return LD_PRELOAD to its original value
setenv('LD_PRELOAD', orig_LD_PRELOAD);

end
