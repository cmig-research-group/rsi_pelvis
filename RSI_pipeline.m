function series_output_list = RSI_pipeline(input_dir, output_dir, params)

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
    error('Unable to parse parameter file');
end

if ~exist(input_dir, 'dir')
   fprintf('Input directory does not exist: %s\n', input_dir);
   return
end

if strcmp(input_dir(end), '/')
   input_dir = input_dir(1:end-1);
end

indx = regexp(input_dir, '\/');
bottom_dir = input_dir(indx(end)+1:end);
contents_bottom_dir = dir(input_dir);
contents_bottom_dir = {contents_bottom_dir.name};

% Check if contents of input_dir are series folders, which means input_dir is a date subfolder, not a top-level patient folder 
match_series = any(~cellfun(@isempty, regexp(contents_bottom_dir, '^Series_')));

if match_series
  patient_name = input_dir(indx(end-1)+1:indx(end)-1);
  pat_dir.name = bottom_dir;
  pat_dir.folder = input_dir(indx(1):indx(end)-1);
else
  patient_name = input_dir(indx(end)+1:end);
  pat_dir = dir(input_dir);
end

exclude = {'.', '..', '.DS_Store'};

series_output_list = {};
for v = 1:length(pat_dir)
  if any(strcmp(pat_dir(v).name, exclude))
    continue
  end

  exam_dir = fullfile(pat_dir(v).folder, pat_dir(v).name);
  disp('Fetching paths to RSI data...');
  paths = fetch_RSI_paths(exam_dir);
  disp(paths);

  if (length(paths.RSI)==1) && isempty(paths.RSI{1})
    disp(['No RSI data found in ' exam_dir]);
    continue
  end

  exam_output_dir = fullfile(output_dir, patient_name, pat_dir(v).name);
  exam_output_dir = nixify(exam_output_dir);

  % Load optional series ----------------------------------------------------------------------------
  if params.WriteOptionalSeries
    try
      write_optional_series(paths, exam_output_dir);
    catch ME
      fprintf('WARNING: %s\n', ME.message);
    end
  end

  % RSI processing ----------------------------------------------------------------------------------
  for i = 1:length(paths.RSI)
    if ~isempty(paths.RSI{i})

      RSI_path_parts = regexp(paths.RSI{i}, filesep, 'split');
      acq_name = RSI_path_parts{end};
      acq_name = nixify(acq_name);
      
      series_output_dir = fullfile(exam_output_dir,acq_name);
      series_output_list{end+1} = series_output_dir;

      disp(['Processing ' paths.RSI{i}]);
      try

	path_selected_T2 = select_T2w_series(paths.T2_ax); 
	processRSI( paths.RSI{i}, paths.RSI_rev{i}, path_selected_T2, paths.DCE, series_output_dir, params );

      catch ME
	fprintf('ERROR: %s\n', ME.message);    
        fprintf('ERROR: Failed to process %s, moving on...\n\n', paths.RSI{i});
      end

    end
  end

end

end
