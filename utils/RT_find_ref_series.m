function path_series = RT_find_ref_series(ref_UID, paths_candidates)

path_series = '';

for i = 1:length(paths_candidates)
  path_candidate = paths_candidates{i};
  files = dir(path_candidate);
  fname_ref = fullfile(files(3).folder, files(3).name);
  dcminfo = dicominfo(fname_ref);
  StudyInstanceUID = dcminfo.StudyInstanceUID;
  if strcmp(StudyInstanceUID, ref_UID)
     path_series = path_candidate;
     return
  end
end

end
