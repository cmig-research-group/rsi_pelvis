function vol_lesions = lesion_id(ctx_prostmask, ctx_rsirs)

minvol  = 125; % Minimum lesion volume, in mm3
rsirsthresh = 150;

vol_prostmask = ctx_prostmask.imgs; 
vol_rsirs = ctx_rsirs.imgs;

vol_lesion = (vol_prostmask.*vol_rsirs) > rsirsthresh;

vol_lesions = [];
if sum(vol_lesion(:)) > 0
  voxdim = sqrt(sum(ctx_prostmask.Mvxl2lph(1:3,1:3).^2,1)); voxvol = prod(voxdim); % Change to use actual voxel size
  [coordmat, nvoxvec, bbmat, PixelIdxList] = RSI_find_centers(vol_lesion);
  volvec = voxvol*nvoxvec;

  lesioncnt = 0;
  for li = 1:min(4,size(coordmat,1))
    if volvec(li) > minvol % Minimum dimension
      vol_tmp = false(size(vol_lesion)); vol_tmp(PixelIdxList{li}) = true; vol_tmp0 = vol_tmp;
      for iter_plump = 1:3
        vol_tmp2 = imdilate(vol_tmp,strel('cube',3));
        vol_tmp = vol_tmp2&(vol_rsirs>120);
      end
      smf = 1;
      vol = vol_tmp;
      volsm = smooth3(vol,'gaussian',1+2*ceil(smf),smf);

      lesioncnt = lesioncnt+1;
      vol_lesions{lesioncnt} = volsm>0.25;
    end
  end

end

end
