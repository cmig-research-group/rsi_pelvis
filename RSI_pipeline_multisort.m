function RSI_pipeline_multisort(input_dir, output_dir, processing_params)

output_dir_dicoms = fullfile(output_dir, 'raw_dicoms');
output_dir_processed = fullfile(output_dir, 'processed_outputs');

disp('Unpacking data...');
unpack_data(input_dir, output_dir_dicoms);
input_dir = output_dir_dicoms;

patients = dir(input_dir);
exclude = {'.', '..', '.DS_Store'};

for p = 1:length(patients)
    if ~any(strcmp(patients(p).name, exclude))
       
       path_patient = fullfile(patients(p).folder, patients(p).name);
       RSI_pipeline(path_patient, output_dir_processed, processing_params);

    end
end

end
