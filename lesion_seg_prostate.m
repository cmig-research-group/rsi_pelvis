function [vol_lesions, vec_rsirs_max] = lesion_seg_prostate(ctx_rsirs, ctx_prostmask, smf)

vol_rsirs = ctx_rsirs.imgs;
mask = logical(ctx_prostmask.imgs);
size_prostate = sum(mask(:));

threshvol = vol_rsirs > 80;
CC = bwconncomp(threshvol);

relevant_lesions = struct;
lesion_cntr = 1;
for i = 1:length(CC.PixelIdxList)
  vol_blob = false(size(vol_rsirs));
  vol_blob(CC.PixelIdxList{i}) = true;
  size_blob = numel(CC.PixelIdxList{i});
  intersect_prostate = mask & vol_blob;
  touches_prostate = any(intersect_prostate(:));
  vec_rsirs_blob = vol_rsirs(vol_blob);
  at_least_100 = any(vec_rsirs_blob>=100);
  bigger_than_prostate = size_blob > size_prostate;
  too_high = any(vec_rsirs_blob>=1000);
  if touches_prostate && at_least_100 && ~bigger_than_prostate && ~too_high
    relevant_lesions(lesion_cntr,1).rsirs_max = max(vec_rsirs_blob);
    relevant_lesions(lesion_cntr,1).indxs = CC.PixelIdxList{i};
    lesion_cntr = lesion_cntr + 1;
  end
end

if ~isfield(relevant_lesions, 'indxs')
  vol_lesions = [];
  vec_rsirs_max = [];
  return
end

vol_lesions = false([size(vol_rsirs) length(relevant_lesions)]);
for i = 1:length(relevant_lesions)
  vol_tmp = vol_lesions(:,:,:,i);
  vol_tmp(relevant_lesions(i).indxs) = true;
  vol_lesions(:,:,:,i) = vol_tmp;
end

vec_rsirs_max = [relevant_lesions.rsirs_max]';
[vec_rsirs_max, I] = sort(vec_rsirs_max, 'descend');
vol_lesions = vol_lesions(:,:,:,I);

voxdim = sqrt(sum(ctx_rsirs.Mvxl2lph(1:3,1:3).^2,1))';
res_dwi = 2.5;
w = round(res_dwi/voxdim(1));
SE = strel('square', w);
for v = 1:size(vol_lesions, 4)
  for s = 1:size(vol_lesions, 3)
    vol_lesions(:,:,s,v) = imclose(vol_lesions(:,:,s,v), SE);
  end
end

vol_smooth = zeros(size(vol_lesions));
for v = 1:size(vol_lesions, 4)
  for s = 1:size(vol_lesions, 3)
    vol_smooth(:,:,s,v) = imgaussfilt(double(vol_lesions(:,:,s,v)), smf);
  end
end
vol_lesions = logical(vol_smooth);

if size(vol_lesions, 4) > 1
  vol_lesions_final = vol_lesions(:,:,:,1);
  vec_rsirs_max_final = vec_rsirs_max(1);
  for i = 2:size(vol_lesions, 4)
    A = vol_lesions_final(:,:,:,end);
    B = vol_lesions(:,:,:,i);
    intersect = A & B;
    if any(intersect(:))
      U = A | B;
      vol_lesions_final(:,:,:,end) = U;
    else
      vol_lesions_final = cat(4, vol_lesions_final, B);
      vec_rsirs_max_final = cat(1, vec_rsirs_max_final, vec_rsirs_max(i));
    end
  end
  vol_lesions = vol_lesions_final;
  vec_rsirs_max = vec_rsirs_max_final;
end

N = 1:size(vol_lesions, 4);
I = N <= 2;
vol_lesions = vol_lesions(:,:,:,I);
vec_rsirs_max = vec_rsirs_max(I);

end
