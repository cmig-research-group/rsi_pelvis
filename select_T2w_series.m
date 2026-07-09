function selected_T2 = select_T2w_series(candidate_T2_series, RSI_series)

if length(candidate_T2_series) == 1
  selected_T2 = candidate_T2_series{1};
  return
end

% Get series numbers
series_nums = zeros(size(candidate_T2_series));
for i = 1:length(candidate_T2_series)
  ims = dir(candidate_T2_series{i});
  im = fullfile(ims(3).folder, ims(3).name);
  info_T2 = dicominfo(im);
  series_nums(i) = info_T2.SeriesNumber;
end

% Get scan coverage
manufacturer = find_manufacturer(info_T2);
if strcmp(manufacturer, 'ge')
  [rsidat, M] = ReadDicomDiffusionData(RSI_series);
elseif strcmp(manufacturer, 'philips')
  [rsidat, M] = ReadDicomDiffusionDataPhilips(RSI_series);
elseif strcmp(manufacturer, 'siemens')
  [rsidat, M] = ReadDicomDiffusionDataSiemens(RSI_series);
end
ctx_rsi = mgh2ctx(rsidat, M);
slices_rsi = size(ctx_rsi.imgs, 3);
center_rsi = ctx_rsi.lphcent';
corner1_rsi = ctx_rsi.Mvxl2lph(1:3,4)'; 

series_slices = zeros(size(candidate_T2_series));
series_centers = zeros(numel(candidate_T2_series), 3);
series_corner1s = zeros(numel(candidate_T2_series), 3);
for i = 1:length(candidate_T2_series)
  ctx_t2 = QD_read_dicomdir(candidate_T2_series{i});
  series_slices(i) = size(ctx_t2.imgs, 3);
  series_centers(i,:) = ctx_t2.lphcent';
  series_corner1s(i,:) = ctx_t2.Mvxl2lph(1:3,4)';
end

% Compare scan coverage
diff_slices = abs(series_slices - slices_rsi);
diff_centers = sqrt( sum( (series_centers - center_rsi).^2, 2) );
diff_corner1s = sqrt( sum( (series_corner1s - corner1_rsi).^2, 2) );

best_match_slices = diff_slices == min(diff_slices);
best_match_center = (diff_centers == min(diff_centers))';
best_match_corner1 = (diff_corner1s == min(diff_corner1s))';
match_score = best_match_slices + best_match_center + best_match_corner1;
best_total_match = find(match_score == max(match_score));

if length(best_total_match) < 2
   selected_T2 = candidate_T2_series{best_total_match}
else
  best_series_nums = series_nums(best_total_match);
  ind_highest_series = find(best_series_nums == max(best_series_nums));
  selected_T2 = candidate_T2_series{best_total_match(ind_highest_series)};
end
  
end
