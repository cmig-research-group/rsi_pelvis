function unpack_data(data_dir, sorted_dir)

t = datetime;
t = posixtime(t);
t = num2str(t);
t = strrep(t, '.', '');

fname_sh = fullfile(sorted_dir, sprintf('mv_cmds_%s.sh',t) );
study_uids = struct;
study_uids = generate_cmds(data_dir, sorted_dir, fname_sh, study_uids);

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


function study_uids = generate_cmds(data_dir, sorted_dir, fname_sh, study_uids)

exclude = {'.', '..', '.DS_Store'};
contents = dir(data_dir);

if ~exist(sorted_dir, 'dir')
   mkdir(sorted_dir);
end

fID = fopen(fname_sh, 'a');

disp(['Sorting ' data_dir]);
for i = 1:length(contents)
  if ~any(strcmp(contents(i).name, exclude))

    item = fullfile(contents(i).folder, contents(i).name);

    if exist(item, 'dir')
      study_uids = generate_cmds(item, sorted_dir, fname_sh, study_uids);
    end

    try
      info = dicominfo(item);

      study_uid = info.StudyInstanceUID;
      subject_name = info.PatientName.FamilyName;
      subject_name = nixify(subject_name);
      pad_flag = 0;
      num_inds = regexp(subject_name, '\d');
      if ~isempty(num_inds) && (num_inds(1) == 1)
	 pad_flag = 1;
	 subject_name = ['PAD' subject_name];
      end

      study_date = info.StudyDate;      
      name_and_date = [subject_name '_' study_date];
      if ~isfield(study_uids, name_and_date)
	study_uids.(name_and_date) = {};
      end
      uids_list = study_uids.(name_and_date);
      study_num = find(strcmp(study_uid, uids_list));
      if isempty(study_num)
         uids_list{end+1} = study_uid;
         study_num = length(uids_list);
	 study_uids.(name_and_date) = uids_list;
      end

      subject_ID = info.PatientID;
      
      series_num = info.SeriesNumber;
      if ~isfield(info, 'SeriesDescription')
	 info.SeriesDescription = ['series' num2str(series_num)];
      end
      series_name = info.SeriesDescription;
      series_name = nixify(series_name);
      series_name = strrep(series_name, '/', '_');

      if strcmp(study_date, '')
	 study_date = 'UnknownDate';
      end

      if strcmp(subject_ID, '')
	subject_ID = 'UnknownID';
      end

      if strcmp(subject_name, '')
	subject_name = 'UnknownName';
      end

      if pad_flag == 0
	subject_dir = fullfile(sorted_dir, subject_name);
      else
	subject_dir = fullfile(sorted_dir, subject_name(4:end));
      end
      if ~exist(subject_dir, 'dir')
	mkdir(subject_dir);
      end

      % Consider multiple scans on a single date
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
      
    catch
      continue
    end

  end
end

fclose(fID);

end
