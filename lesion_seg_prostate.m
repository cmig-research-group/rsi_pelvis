function [vol_lesions, vec_rsirs_max] = lesion_seg_prostate(ctx_rsirs, ctx_prostmask)

vol_rsirs = ctx_rsirs.imgs;
mask = logical(ctx_prostmask.imgs);
size_prostate = sum(mask(:));

% Identify connected components of potential lesions
threshvol = vol_rsirs > 80;
CC = bwconncomp(threshvol);

relevant_lesions = struct;
lesion_cntr = 1;
for i = 1:length(CC.PixelIdxList)

  vol_blob = false(size(vol_rsirs));
  vol_blob(CC.PixelIdxList{i}) = true;
  size_blob = numel(CC.PixelIdxList{i});

  % Lesion-selection criteria
  vec_rsirs_blob = vol_rsirs(vol_blob);
  at_least_100 = any(vec_rsirs_blob>=100);
  bigger_than_prostate = size_blob > size_prostate;
  too_high = any(vec_rsirs_blob>=1000);

  quarter_blob = round(size_blob/4);
  intersect_prostate = mask & vol_blob;
  size_intersection = sum(intersect_prostate(:));
  quarter_within = size_intersection >= quarter_blob;
  vec_rsirs_intersect = vol_rsirs(intersect_prostate);
  rsirs_max_total = max(vec_rsirs_blob);
  rsirs_max_intersect = max(vec_rsirs_intersect);
  bright_in_prostate = rsirs_max_intersect >= (0.9*rsirs_max_total);

  if at_least_100 && quarter_within && bright_in_prostate && ~bigger_than_prostate && ~too_high
    relevant_lesions(lesion_cntr,1).rsirs_max = rsirs_max_total;
    relevant_lesions(lesion_cntr,1).indxs = CC.PixelIdxList{i};
    lesion_cntr = lesion_cntr + 1;
  end

end

if ~isfield(relevant_lesions, 'indxs')
  vol_lesions = [];
  vec_rsirs_max = [];
  return
end

% Assemble binary volume of potential lesion masks
vol_lesions = false([size(vol_rsirs) length(relevant_lesions)]);
for i = 1:length(relevant_lesions)
  vol_tmp = vol_lesions(:,:,:,i);
  vol_tmp(relevant_lesions(i).indxs) = true;
  vol_lesions(:,:,:,i) = vol_tmp;
end

% Sort lesions by descending RSIrs_max
vec_rsirs_max = [relevant_lesions.rsirs_max]';
[vec_rsirs_max, I] = sort(vec_rsirs_max, 'descend');
vol_lesions = vol_lesions(:,:,:,I);

% Morphological closing to close holes in lesion masks
voxdim = sqrt(sum(ctx_rsirs.Mvxl2lph(1:3,1:3).^2,1))';
res_dwi = 2.5;
w = round(res_dwi/voxdim(1));
SE = strel('square', w);
for v = 1:size(vol_lesions, 4)
  for s = 1:size(vol_lesions, 3)
    vol_lesions(:,:,s,v) = imclose(vol_lesions(:,:,s,v), SE);
  end
end

% Only return the two lesions with highest RSIrs_max
N = 1:size(vol_lesions, 4);
I = N <= 2;
vol_lesions = vol_lesions(:,:,:,I);
vec_rsirs_max = vec_rsirs_max(I);

end
