function RSI_pipeline(input_dir, output_dir, params)

if ~exist(input_dir, 'dir')
   fprintf('Input directory does not exist: %s\n', input_dir);
   return
end

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

if strcmp(input_dir(end), '/')
   input_dir = input_dir(1:end-1);
end

indx = regexp(input_dir, '\/');
bottom_dir = input_dir(indx(end)+1:end);
match_date = ~isempty(regexp(bottom_dir, '\d{8}$'));
match_date2 = ~isempty(regexp(bottom_dir, '\d{8}_Study*'));

if match_date || match_date2
  patient_name = input_dir(indx(end-1)+1:indx(end)-1);
  pat_dir.name = bottom_dir;
  pat_dir.folder = input_dir(indx(1):indx(end)-1);
else
  patient_name = input_dir(indx(end)+1:end);
  pat_dir = dir(input_dir);
end

exclude = {'.', '..', '.DS_Store'};

for v = 1:length(pat_dir)
  if ~any(strcmp(pat_dir(v).name, exclude))

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

	disp(['Processing ' paths.RSI{i}]);
	try

	  processRSI( paths.RSI{i}, paths.RSI_rev{i}, paths.T2_ax{end}, paths.DCE, fullfile(exam_output_dir,acq_name), params );

	catch ME
	  fprintf('ERROR: %s\n', ME.message);    
          fprintf('ERROR: Failed to process %s, moving on...\n\n', paths.RSI{i});
	end

      end
    end

  end
end

end
