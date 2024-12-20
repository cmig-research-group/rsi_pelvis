function bladder_mask = get_bladder_mask(C3_vol, threshold_percent)

% Finds largest connected component of C3 volume

thresh = threshold_percent .* prctile(C3_vol(:), 99);

init = C3_vol > thresh;

CC = bwconncomp(init, 6);
numPixels = cellfun(@numel, CC.PixelIdxList);
[biggest, idx] = max(numPixels);

bladder_mask = false(size(C3_vol));
bladder_mask(CC.PixelIdxList{idx}) = true;

end
