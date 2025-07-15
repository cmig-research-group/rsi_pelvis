function ctx_report = create_visual_report(ctx_T2, ctx_mask_prostate, ctx_RSIrs, ctx_lesions, prostate_detector_output, has_implant)

% Create overlay on which to print the report ----------------------------------------
if ~isempty(ctx_RSIrs)

  vol_prostate_contour = zeros(size(ctx_mask_prostate.imgs));
  for i = 1:size(vol_prostate_contour,3)
    vol_prostate_contour(:,:,i) = double(bwperim(logical(ctx_mask_prostate.imgs(:,:,i)), 8));
  end
  ctx_prostate_contour = ctx_mask_prostate;
  ctx_prostate_contour.imgs = vol_prostate_contour;
  ctx_overlay = create_color_overlay(ctx_prostate_contour, ctx_T2, [2 3], 'cool', ctx_prostate_contour, [], 1);

  for j = 1:size(ctx_lesions.imgs,4)
      vol_lesion_contour = zeros(size(ctx_mask_prostate.imgs));
      for i = 1:size(vol_lesion_contour,3)
	vol_lesion_contour(:,:,i) = double(bwperim(logical(ctx_lesions.imgs(:,:,i,j)), 8));
      end
      ctx_lesion_contour = ctx_mask_prostate;
      ctx_lesion_contour.imgs = vol_lesion_contour;
      ctx_overlay = create_color_overlay(ctx_lesion_contour, ctx_overlay, [0 0.5], 'cool', ctx_lesion_contour, [], 1);
  end

  ctx_overlay = create_color_overlay(ctx_RSIrs, ctx_overlay, [80 180]);

else
  ctx_overlay = ctx_T2;
  ctx_overlay.imgs = repmat(ctx_T2.imgs, 1, 1, 1, 3);
  ctx_overlay.imgs = ctx_overlay.imgs ./ max(ctx_T2.imgs(:));
end


% Create report -------------------------------------------------------------
report_title = '----- Visual Report -----';
box_color = 'c';
box_opacity = 0.9;
font_size = 12;

if prostate_detector_output
  text_prostate = 'Detected';
else
  text_prostate = 'Not detected';
  box_color = 'y';
end

if has_implant
  text_implant = 'Detected';
  box_color = 'y';
else
  text_implant = 'None detected';
end

header = sprintf('%s\nAI feature-detection results:\n    Prostate: %s\n    Hip implants: %s\n', report_title, text_prostate, text_implant);

vol_report = ctx_overlay.imgs;
[rows, cols, slices, ~] = size(vol_report);
row_mid = round(rows/2);

ctx_report = ctx_overlay;


% Alert user to missing prostate ---------------------------------------------
if isempty(ctx_mask_prostate)
  box_color = 'r';
  for i = 1:slices
    im = permute(squeeze(vol_report(:,:,i,:)), [2 1 3]);
    text_alert = sprintf('%sWARNING: AI contouring software could not identify a prostate\n(Beware that false negatives may occur)\nNo quantitative RSI values were returned', header);
    im_alert = insertText(im, [0 0], text_alert, 'AnchorPoint', 'LeftTop', 'BoxColor', box_color, 'BoxOpacity', box_opacity, 'FontSize', font_size);
    ctx_report.imgs(:,:,i,:) = permute(im_alert, [2 1 3]);
  end
  return
end


% Prostate overview ----------------------------------------------------------
mask_prostate = logical(ctx_mask_prostate.imgs);
vol_RSIrs = ctx_RSIrs.imgs;
max_RSIrs = max(vol_RSIrs(mask_prostate));

for i = 1:slices

  im = permute(squeeze(vol_report(:,:,i,:)), [2 1 3]);

  if max_RSIrs < 500
    text_RSIrs = 'Typical';
  elseif max_RSIrs >= 500 && max_RSIrs < 1000
    box_color = 'y';
    text_RSIrs = 'Abnormal; may indicate imaging artifact';
  elseif max_RSIrs >= 1000
    box_color = 'r';
    text_RSIrs = 'Highly abnormal; likely indicates image reconstruction failure or other severe artifact';
  end

  text_overview = sprintf('%sRSIrs range: %s', header, text_RSIrs);
  im_overview = insertText(im, [0 0], text_overview, 'AnchorPoint', 'LeftTop', 'BoxColor', box_color, 'BoxOpacity', box_opacity, 'FontSize', font_size);

  ctx_report.imgs(:,:,i,:) = permute(im_overview, [2 1 3]);

end

if isempty(ctx_lesions)
  return
end


% Per lesion ----------------------------------------------------------------
vol_lesions = logical(ctx_lesions.imgs);

for i = 1:size(vol_lesions,4)
  mask_lesion = vol_lesions(:,:,:,i);

  vol_masked = mask_lesion .* vol_RSIrs;
  [coord_row, coord_col, coord_slice] = ind2sub(size(vol_masked), find(vol_masked==max(vol_masked(:))));

  im = permute(squeeze(ctx_report.imgs(:,:,coord_slice,:)), [2 1 3]);
  max_RSIrs = max(vol_RSIrs(mask_lesion));

  text_lesion = sprintf('Max RSIrs: %0.2f', max_RSIrs);

  if coord_row < row_mid
    anchor = 'RightTop';
  else
    anchor = 'LeftTop';
  end

  im_lesion = insertText(im, [coord_row coord_col], text_lesion, 'AnchorPoint', anchor, 'BoxColor', 'w', 'BoxOpacity', 0.3, 'TextColor', 'black', 'FontSize', font_size);
  ctx_report.imgs(:,:,coord_slice,:) = permute(im_lesion, [2 1 3]);

end

end
