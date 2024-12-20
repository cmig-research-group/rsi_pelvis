function ADC_volume = compute_ADCs(DWI_volume, b_values)

% Assumes averaged DWI volume, including b=0
% b_values should correspond to the 4th dimension of DWI_volume

b0_volume = DWI_volume(:,:,:,1);
[rows, cols, slices] = size(b0_volume);

mask = b0_volume < 0;

log_b0_volume = log(b0_volume);
log_b0_volume(imag(log_b0_volume)~=0) = 0;
log_b0_volume(isinf(log_b0_volume)) = 0;
log_DWI_volume = log(DWI_volume);
log_DWI_volume(imag(log_DWI_volume)~=0) = 0;
log_DWI_volume(isinf(log_DWI_volume)) = 0;

Y = zeros( numel(b_values), numel(b0_volume) );
for b = 1:length(b_values)
    diff_vol = log_DWI_volume(:,:,:,b) - log_b0_volume;
    Y(b,:) = diff_vol(:)';
end

B = -b_values;

ADCs = B\Y;

ADC_volume = reshape(ADCs, rows, cols, slices);
ADC_volume(isnan(ADC_volume)) = 0;
ADC_volume(mask) = 0;

end
