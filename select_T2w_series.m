function selected_T2 = select_T2w_series(candidate_T2_series)

if length(candidate_T2_series) == 1
  selected_T2 = candidate_T2_series{1};
  return
end

series_nums = zeros(size(candidate_T2_series));
for i = 1:length(candidate_T2_series)
    ims = dir(candidate_T2_series{i});
    im = fullfile(ims(3).folder, ims(3).name);
    info_T2 = dicominfo(im);
    series_nums(i) = info_T2.SeriesNumber;
end

ind_highest_series = find(series_nums==max(series_nums));

selected_T2 = candidate_T2_series{ind_highest_series};

end
