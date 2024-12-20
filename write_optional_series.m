function write_optional_series(path_struct, exam_output_dir)

% Axial T2
if ~isempty(path_struct.T2_ax{1})
  for i = 1:length(path_struct.T2_ax)
    try
      path_parts = regexp(path_struct.T2_ax{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
	mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'T2_ax.mgz'); 
      vol = QD_read_dicomdir(path_struct.T2_ax{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Sagittal T2 cube
if ~isempty(path_struct.T2_sag_cube{1})
  for i = 1:length(path_struct.T2_sag_cube)
    try
      path_parts = regexp(path_struct.T2_sag_cube{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'T2_sag_cube.mgz');
      vol = QD_read_dicomdir(path_struct.T2_sag_cube{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Conventional DWI
if ~isempty(path_struct.DWI_conventional{1})
  for i = 1:length(path_struct.DWI_conventional)
    try
      path_parts = regexp(path_struct.DWI_conventional{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'DWI_conventional.mgz');
      fname_out_qmat = fullfile(path_out_opt, 'qmat.mat');
      fname_out_bvals = fullfile(path_out_opt, 'bvals.mat');
      fname_out_gwinfo = fullfile(path_out_opt, 'gwinfo.mat');
      fname_out_dcminfo = fullfile(path_out_opt, 'dcminfo.mat');
      [vol, M, qmat, bvals, gwinfo, dcminfo] = ReadDicomDiffusionData(path_struct.DWI_conventional{i});
      QD_save_mgh(vol, fname_out_opt, M);
      save(fname_out_qmat, 'qmat');
      save(fname_out_bvals, 'bvals');
      save(fname_out_gwinfo, 'gwinfo');
      save(fname_out_dcminfo, 'dcminfo');
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Conventional ADC
if ~isempty(path_struct.ADC_conventional{1})
  for i = 1:length(path_struct.ADC_conventional)
    try
      path_parts = regexp(path_struct.ADC_conventional{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'ADC_conventional.mgz');
      vol = QD_read_dicomdir(path_struct.ADC_conventional{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Dixon in-phase
if ~isempty(path_struct.Dixon_in{1})
  for i = 1:length(path_struct.Dixon_in)
    try
      path_parts = regexp(path_struct.Dixon_in{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'Dixon_in.mgz');
      vol = QD_read_dicomdir(path_struct.Dixon_in{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Dixon out-phase
if ~isempty(path_struct.Dixon_out{1})
  for i = 1:length(path_struct.Dixon_out)
    try
      path_parts = regexp(path_struct.Dixon_out{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'Dixon_out.mgz');
      vol = QD_read_dicomdir(path_struct.Dixon_out{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Dixon fat
if ~isempty(path_struct.Dixon_fat{1})
  for i = 1:length(path_struct.Dixon_fat)
    try
      path_parts = regexp(path_struct.Dixon_fat{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'Dixon_fat.mgz');
      vol = QD_read_dicomdir(path_struct.Dixon_fat{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

% Dixon water
if ~isempty(path_struct.Dixon_water{1})
  for i = 1:length(path_struct.Dixon_water)
    try
      path_parts = regexp(path_struct.Dixon_water{i}, filesep, 'split');
      acq_name = path_parts{end};
      fprintf('Loading %s\n', acq_name);
      acq_name = nixify(acq_name);
      path_out_opt = fullfile(exam_output_dir, acq_name);
      if ~exist(path_out_opt, 'dir')
        mkdir(path_out_opt);
      end
      fname_out_opt = fullfile(path_out_opt, 'Dixon_water.mgz');
      vol = QD_read_dicomdir(path_struct.Dixon_water{i});
      QD_ctx_save_mgh(vol, fname_out_opt);
    catch ME
      fprintf('%s\n', ME.message);
    end
  end
end

end
