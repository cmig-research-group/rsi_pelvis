function subject_info = unpack_data(data_dir, sorted_dir, sort_variable)

if ~exist('sort_variable', 'var')
  sort_variable = 'PatientName';
end

t = datetime;
t = posixtime(t);
t = num2str(t);
t = strrep(t, '.', '');

fname_sh = fullfile(sorted_dir, sprintf('mv_cmds_%s.sh',t) );
subject_info = struct('subject_ids', [], 'scan_dates', [], 'study_uids', [], 'fingerprints', []);
subject_info = generate_cmds(data_dir, sorted_dir, fname_sh, subject_info, sort_variable);
subject_info = rmfield(subject_info, 'fingerprints');

disp('Moving all those files...');
cmd1 = ['chmod 777 ' fname_sh];
system(cmd1);
cmd2 = fname_sh;
status = system(cmd2);

if status == 0
  delete(fname_sh);
else
  error('Sorting failed');
end

end


function subject_info = generate_cmds(data_dir, sorted_dir, fname_sh, subject_info, sort_variable)

exclude_exact = {'.', '..', '.DS_Store'};
exclude_regexp = 'XX_'; % Ignore Philips metadata files for now  
contents = dir(data_dir);

if ~exist(sorted_dir, 'dir')
   mkdir(sorted_dir);
end

fID = fopen(fname_sh, 'a');

disp(['Sorting ' data_dir]);
for i = 1:length(contents)

  match_exclude_exact = any(strcmp(contents(i).name, exclude_exact));
  match_exclude_regexp = ~isempty(regexp(contents(i).name, exclude_regexp));
  if match_exclude_exact | match_exclude_regexp
    continue
  end

  item = fullfile(contents(i).folder, contents(i).name);

  if exist(item, 'dir')
    subject_info = generate_cmds(item, sorted_dir, fname_sh, subject_info, sort_variable);
  end

  try

    if exist(item, 'dir') == 7
      continue
    end

    info = dicominfo(item);

    study_uid = info.StudyInstanceUID;
    subject_id = info.PatientID;
    subject_name = info.PatientName.FamilyName;
    subject_name = nixify(subject_name);
    study_date = info.StudyDate; 
    series_num = info.SeriesNumber;
    if ~isfield(info, 'SeriesDescription')
      info.SeriesDescription = ['series' num2str(series_num)];
    end
    series_name = info.SeriesDescription;
    series_name = nixify(series_name);
    series_name = strrep(series_name, '/', '_');

    if strcmp(subject_id, '')
      subject_id = 'UnknownID';
    end
    if strcmp(subject_name, '')
      subject_name = 'UnknownName';
    end
    if strcmp(study_date, '')
      study_date = 'UnknownDate';
    end

    % Top-level sorting variable; can be patient's MRN or name
    if strcmp(sort_variable, 'PatientID')
      id = subject_id;
    elseif strcmp(sort_variable, 'PatientName')
      id = subject_name;
    end
    
    subject_dir = fullfile(sorted_dir, id);
    if ~exist(subject_dir, 'dir')
      mkdir(subject_dir);
    end

    fingerprint = sprintf('%s_%s_%s', id, study_date, study_uid);
    fingerprint_list = subject_info.fingerprints;
    match_fingerprint = find(strcmp(fingerprint, fingerprint_list));
    if isempty(match_fingerprint)
       subject_info.fingerprints{end+1} = fingerprint;
       match_fingerprint = length(subject_info.fingerprints);
       subject_info.subject_ids{end+1} = id;
       subject_info.scan_dates{end+1} = study_date;
       subject_info.study_uids{end+1} = study_uid;
    end

    % Consider multiple scans on a single date
    inds_id = strcmp(subject_info.subject_ids, id);
    inds_date = strcmp(subject_info.scan_dates, study_date);
    matches_id_date = find(inds_id & inds_date);
    first_match = min(matches_id_date);
    study_num = (match_fingerprint - first_match) + 1;

    if study_num == 1
      visit_dir = fullfile(subject_dir, study_date);
    elseif study_num > 1
      study_date_and_num = sprintf('%s_Study%s', study_date, num2str(study_num));
      visit_dir = fullfile(subject_dir, study_date_and_num);
    end
    if ~exist(visit_dir, 'dir')
      mkdir(visit_dir);
    end

    series_dir = fullfile(visit_dir, ['Series_' num2str(series_num) '__' series_name]);
    if ~exist(series_dir, 'dir')
      mkdir(series_dir);
    end

    fprintf(fID, 'mv ''%s'' ''%s%c%s''\n', item, series_dir, filesep, contents(i).name);
    
  catch ME
    disp(ME.message);
  end

end

fclose(fID);

end
