function [coordmat nvoxvec bbmat PixelIdxList] = RSI_find_centers(maskvol)

binvol = (maskvol>0);
conn = 26; CC = bwconncomp(binvol,conn);
[sv si] = sort(cellfun(@length,CC.PixelIdxList),'descend');
PixelIdxList = CC.PixelIdxList(si);
coordmat = NaN(length(PixelIdxList),3);
nvoxvec = NaN(length(PixelIdxList),1);
bbmat = NaN(length(PixelIdxList),6);
dims = size(binvol);
for i = 1:length(PixelIdxList)
  vol_tmp = false(dims);
  vol_tmp(PixelIdxList{i}) = true;
  nvoxvec(i) = length(PixelIdxList{i});
  bbmat(i,:) = getfield(regionprops(vol_tmp),'BoundingBox');
  coordmat(i,:) = mean(ind2sub_amd(dims,PixelIdxList{i}),1);
end

