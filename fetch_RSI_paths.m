function paths = fetch_RSI_paths(exam_dir)

% Enumerate data to be collected -----------------------------
paths = struct;

% Mandatory
paths.RSI = {''};

% Optional
paths.T2_ax = {''};
paths.T2_sag_cube = {''};
paths.DWI_conventional = {''};
paths.ADC_conventional = {''};
paths.DCE = {''};
paths.Dixon_in = {''};
paths.Dixon_out = {''};
paths.Dixon_fat = {''};
paths.Dixon_water = {''};


% Patterns used to ID data -----------------------------------
% Ignore
global_exclude = {'.', '..', '.DS_Store'};

% RSI 
seqID_RSI_GE = {'epi2_pepolarFOCUSFLEX', 'epi2_pepolarFLEX', 'epi2_ART', 'epi2_revART', 'epi2alt', 'epi2altoff'};
patterns_seqID_RSI_Siemens = {'ep_b', 'ez_b', 'WI_b'};
patterns_exclude_RSI = {'trace', 'adc', 'color'};
patterns_exclude_RSI_siemens = {'trace', 'adc', 'color', 'fa', 'tensor'};

% Axial T2 for RSI overlay
patterns_T2_ax_plane = {'ax', 'tra'};
patterns_T2_ax_contrast = {'T2', 'FSE'};
patterns_exclude_T2_ax = {'water', 'fat', 'flex', 'sag', 'cor', 'reformat', 'cube', 'T1', 'bifurcation', 'prop', '3d', 'space'};

% Sag Cube
seqID_cube_GE = 'Cube';
patterns_T2_sag_cube = 'Sag';

% Conventional DWI
seqID_DWI_GE = 'epi2';
patterns_DWI = 'DWI';

% Conventional ADC
patterns_ADC = 'ADC';

% DCE
patterns_DCE = {'DISCO', 'TRICKS'};
patterns_exclude_DCE = 'SUB';

% Dixon
patterns_Dixon_in = 'InPhase:[\w\s-]*LAVA[\w\s-]*Flex';
patterns_Dixon_out = 'OutPhase:[\w\s-]*LAVA[\w\s-]*Flex';
patterns_Dixon_fat = 'FAT:[\w\s-]*LAVA[\w\s-]*Flex';
patterns_Dixon_water = 'WATER:[\w\s-]*LAVA[\w\s-]*Flex';


% -----------------------------------------------------------
RSI_path_num = 1;
T2_ax_path_num = 1;
T2_sag_cube_path_num = 1;
DWI_conventional_path_num = 1;
ADC_conventional_path_num = 1;
DCE_path_num = 1;
in_path_num = 1;
out_path_num = 1;
fat_path_num = 1;
water_path_num = 1;

acqs = dir(exam_dir);
for i = 1:length(acqs)
  if ~any(strcmp(acqs(i).name, global_exclude))
    acq_path = fullfile(exam_dir, acqs(i).name);

    ims = dir(acq_path);
    try
      im = ims(3).name;
      info = dicominfo(fullfile(acq_path, im));
      info = fix_impax_dcm_tags(info);
    catch ME
      fprintf('%s\n', ME.message)
      continue
    end

    if ~isfield(info, 'ImageType')
       continue
    end
    im_type = info.ImageType;
    
    SeriesDescription = info.SeriesDescription;
    manufacturer = info.Manufacturer;
    if strcmpi(manufacturer, 'ge medical systems')
      manufacturer = 'ge';
    elseif any(strcmpi(manufacturer, {'siemens', 'siemens healthineers'}))
      manufacturer = 'siemens';
    end

    if ~isfield(info,'Private_0019_109c')
      info.Private_0019_109c = '';
    end
    
    if ~isfield(info,'SequenceName')
      info.SequenceName = '';
    end

    if strcmp(manufacturer, 'ge')
      seq_name = info.Private_0019_109c;
    elseif strcmp(manufacturer, 'siemens')
      seq_name = info.SequenceName;
    end

    % Check for RSI data
    if strcmp(manufacturer, 'ge')
       match_RSI_seq = any(~cellfun(@isempty, regexpi(seq_name, seqID_RSI_GE)));
       match_RSI_name = ~isempty(regexpi(SeriesDescription, 'RSI'));
       match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_RSI)));
       if (match_RSI_seq || match_RSI_name) && ~match_exclude
	  paths.RSI{RSI_path_num} = acq_path;
	  series_description_list{RSI_path_num} = SeriesDescription; % Save for later filtering of reverse acquisitions
	  RSI_path_num = RSI_path_num + 1;
       end
    end
    if strcmp(manufacturer, 'siemens')
       match_RSI_seq = any(~cellfun(@isempty, regexpi(seq_name, patterns_seqID_RSI_Siemens)));
       match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_RSI_siemens)));
       if match_RSI_seq && ~match_exclude
	  paths.RSI{RSI_path_num} = acq_path;
	  series_description_list{RSI_path_num} = SeriesDescription; % Save for later filtering of reverse acquisitions
	  RSI_path_num = RSI_path_num + 1;
       end
    end

    % Check for anatomical T2 DICOMs
    match_T2_plane = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_T2_ax_plane)));
    match_T2_contrast = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_T2_ax_contrast)));
    match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_T2_ax)));
    if match_T2_plane && match_T2_contrast && ~match_exclude
      paths.T2_ax{T2_ax_path_num} = acq_path;
      T2_ax_path_num = T2_ax_path_num + 1;
    end

    % Check for sag T2 cube
    if strcmp(seq_name, seqID_cube_GE) && ~isempty(regexpi(SeriesDescription, patterns_T2_sag_cube))
       paths.T2_sag_cube{T2_sag_cube_path_num} = acq_path;
       T2_sag_cube_path_num = T2_sag_cube_path_num + 1;
    end

    % Check for conventional DWI
    if strcmp(seq_name, seqID_DWI_GE) && ~isempty(regexpi(SeriesDescription, patterns_DWI))
       paths.DWI_conventional{DWI_conventional_path_num} = acq_path;
       DWI_conventional_path_num = DWI_conventional_path_num + 1;
    end

    % Check for conventional ADC
    if strcmp(seq_name, seqID_DWI_GE) && ~isempty(regexpi(SeriesDescription, patterns_ADC))
       paths.ADC_conventional{ADC_conventional_path_num} = acq_path;
       ADC_conventional_path_num = ADC_conventional_path_num + 1;
    end

    % Check for DCE
    match_name_dce = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_DCE)));
    match_exclude_dce = isempty(regexpi(SeriesDescription, patterns_exclude_DCE));
    match_orig_dce = ~isempty(regexpi(im_type, 'ORIGINAL'));
    if match_name_dce && match_exclude_dce && match_orig_dce
       paths.DCE{DCE_path_num} = acq_path;
       DCE_path_num = DCE_path_num + 1;
    end

    % Check for Dixon
    if ~isempty(regexpi(SeriesDescription, patterns_Dixon_in)) % In-phase
       paths.Dixon_in{in_path_num} = acq_path;
       in_path_num = in_path_num + 1;
    end
    if ~isempty(regexpi(SeriesDescription, patterns_Dixon_out)) % Out-phase
       paths.Dixon_out{out_path_num} = acq_path;
       out_path_num = out_path_num + 1;
    end
    if ~isempty(regexpi(SeriesDescription, patterns_Dixon_fat)) % Fat
       paths.Dixon_fat{fat_path_num} = acq_path;
       fat_path_num = fat_path_num + 1;
    end
    if ~isempty(regexpi(SeriesDescription, patterns_Dixon_water)) % Water
       paths.Dixon_water{water_path_num} = acq_path;
       water_path_num = water_path_num + 1;
    end

  end
end


if isempty(paths.RSI{1})
   return
end


% Separate reverse PE volumes, if any
paths.RSI_rev = cell(size(paths.RSI));
for i = 1:length(paths.RSI)
  series1 = series_description_list{i};
  for j = 1:length(paths.RSI)
    series2 = series_description_list{j};
    if strcmp(series2, [series1 '_rev'])
      paths.RSI_rev{i} = paths.RSI{j};
      paths.RSI{j} = [];
    end
  end
end


end
