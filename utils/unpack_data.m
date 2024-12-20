function unpack_data(data_dir, sorted_dir)

t = datetime;
t = posixtime(t);
t = num2str(t);
t = strrep(t, '.', '');

fname_sh = fullfile(sorted_dir, sprintf('mv_cmds_%s.sh',t) );
generate_cmds(data_dir, sorted_dir, fname_sh);

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


function generate_cmds(data_dir, sorted_dir, fname_sh)

contents = dir(data_dir);

if ~exist(sorted_dir, 'dir')
   mkdir(sorted_dir);
end

fID = fopen(fname_sh, 'a');

disp(['Sorting ' data_dir]);
for i = 1:length(contents)
  if strcmp(contents(i).name,'.') == 0 && strcmp(contents(i).name,'..') == 0

    item = fullfile(contents(i).folder, contents(i).name);

    if exist(item, 'dir')
      generate_cmds(item, sorted_dir, fname_sh);
    end

    try
      info = dicominfo(item);
      study_date = info.StudyDate;
      subject_name = info.PatientName.FamilyName;
      subject_name = nixify(subject_name);
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

      subject_dir = fullfile(sorted_dir, subject_name);
      if ~exist(subject_dir, 'dir')
	mkdir(subject_dir);
      end

      visit_dir = fullfile(subject_dir, study_date);
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
