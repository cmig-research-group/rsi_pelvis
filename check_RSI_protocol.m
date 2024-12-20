function problem_flag = check_RSI_protocol(output_dir, scans_dir, mandatory_series_list, ref_info_dwi_fwd, scan_info_dwi_fwd, ref_info_T2, scan_info_T2, ref_info_dwi_rev, scan_info_dwi_rev)

problem_flag = 0;

fname_report = fullfile(output_dir, 'protocol_compliance_report.txt');
fid = fopen(fname_report, 'w');


% Basic info --------------------------------------------------
institute = '';
if isfield(scan_info_dwi_fwd.dcminfo, 'InstitutionName')
  institute = scan_info_dwi_fwd.dcminfo.InstitutionName;
end

station = '';
if isfield(scan_info_dwi_fwd.dcminfo, 'StationName')
  station = scan_info_dwi_fwd.dcminfo.StationName;
end

protocol = '';
if isfield(scan_info_dwi_fwd.dcminfo, 'ProtocolName')
  protocol = scan_info_dwi_fwd.dcminfo.ProtocolName;
end

operator = '';
if isfield(scan_info_dwi_fwd.dcminfo, 'OperatorsName')
   names = fieldnames(scan_info_dwi_fwd.dcminfo.OperatorsName);
   for i = 1:length(names)
       name = scan_info_dwi_fwd.dcminfo.OperatorsName.(names{i});
       operator = [operator, name];
   end
end


% Begin compliance check --------------------------------------
fprintf(fid, '##############################\n');
fprintf(fid, '# Protocol Compliance Report #\n');
fprintf(fid, '##############################\n\n');
fprintf(fid, 'Institution: %s, Station: %s\n', institute, station);
fprintf(fid, 'Protocol: %s, Operator: %s\n\n', protocol, operator);

fprintf(fid, 'Generated %s\n\n', datetime('now'));


% Check that operator is listed
if isempty(operator) || ~isempty(regexpi(operator, 'none'))
   problem_flag = 1;
   fprintf(fid, 'No operator reported for this exam\n\n');
end


% Check for mandatory series
if ~isempty(mandatory_series_list)
  scans = dir(scans_dir);
  scans = {scans.name};
  for i = 1:length(mandatory_series_list)

      match = any(~cellfun(@isempty, regexpi(scans, mandatory_series_list{i})));
      
      if match == 0 && ~isempty(regexpi(mandatory_series_list{i}, 'DISCO'))
	match1 = any(~cellfun(@isempty, regexpi(scans, 'LAVA')));
	match2 = any(~cellfun(@isempty, regexpi(scans, 'DCE')));
	match = match1 | match2;
	if match == 0
          problem_flag = 1;
          fprintf(fid, 'Missing series ---------------------------------------\n');
          fprintf(fid, 'DCE series (either DISCO or LAVA) not found\n\n');
	end
	continue
      end

      if match == 0
	 problem_flag = 1;
	 fprintf(fid, 'Missing series ---------------------------------------\n');
	 fprintf(fid, '%s not found\n\n', mandatory_series_list{i});
      end

  end
end
   

% Check b-values 
ref_bvals_fwd = ref_info_dwi_fwd.bvals;
ref_ubvals_fwd = unique(ref_bvals_fwd);
scan_bvals_fwd = scan_info_dwi_fwd.bvals;
scan_ubvals_fwd = unique(scan_bvals_fwd);

if ~isequal(ref_ubvals_fwd, scan_ubvals_fwd)
   problem_flag = 1;
   fprintf(fid, 'b-value deviations: Forward RSI scan ---------------------------------------\n');
   fprintf(fid, 'Expected b-values: [');
   fprintf(fid, '%d, ', ref_ubvals_fwd(1:end-1));
   fprintf(fid, '%d]\n', ref_ubvals_fwd(end));
   fprintf(fid, 'Acquired b-values: [');
   fprintf(fid, '%d, ', scan_ubvals_fwd(1:end-1));
   fprintf(fid, '%d]\n\n', scan_ubvals_fwd(end));
end

if exist('scan_info_dwi_rev', 'var')
   ref_bvals_rev = ref_info_dwi_rev.bvals;
   ref_ubvals_rev = unique(ref_bvals_rev);
   scan_bvals_rev = scan_info_dwi_rev.bvals;
   scan_ubvals_rev = unique(scan_bvals_rev);
   if ~isequal(ref_ubvals_rev, scan_ubvals_rev)
      problem_flag = 1;
      fprintf(fid, 'b-value deviations: Reverse RSI scan ---------------------------------------\n');
      fprintf(fid, 'Expected b-values: [');
      fprintf(fid, '%d, ', ref_ubvals_rev(1:end-1));
      fprintf(fid, '%d]\n', ref_ubvals_rev(end));
      fprintf(fid, 'Acquired b-values: [');
      fprintf(fid, '%d, ', scan_ubvals_rev(1:end-1));
      fprintf(fid, '%d]\n\n', scan_ubvals_rev(end));
   end
end


% Check other scan parameters
ref_fields_fwd = fieldnames(ref_info_dwi_fwd.dcminfo);
for i = 1:length(ref_fields_fwd)
  try
    if ~isequal( ref_info_dwi_fwd.dcminfo.(ref_fields_fwd{i}), scan_info_dwi_fwd.dcminfo.(ref_fields_fwd{i}) )
      problem_flag = 1;
      fprintf(fid, 'Scan parameter deviation for RSI forward scan: %s ---------------------------------------\n', ref_fields_fwd{i});
      if isnumeric(ref_info_dwi_fwd.dcminfo.(ref_fields_fwd{i}))
	fprintf(fid, 'Expected value: %.2f, Observed value: %.2f\n', [ref_info_dwi_fwd.dcminfo.(ref_fields_fwd{i})'; scan_info_dwi_fwd.dcminfo.(ref_fields_fwd{i})']);
      else
	fprintf(fid, 'Expected value: %s, Observed value: %s\n', ref_info_dwi_fwd.dcminfo.(ref_fields_fwd{i}), scan_info_dwi_fwd.dcminfo.(ref_fields_fwd{i}));
      end
      fprintf(fid, '\n');
    end
  catch ME
    fprintf('%s\n', ME.message);
    continue
  end
end

ref_fields_T2 = fieldnames(ref_info_T2.dcminfo);
for i = 1:length(ref_fields_T2)
  try
    if ~isequal( ref_info_T2.dcminfo.(ref_fields_T2{i}), scan_info_T2.dcminfo.(ref_fields_T2{i}) )
      problem_flag = 1;
      fprintf(fid, 'Scan parameter deviation for axial T2 scan: %s ---------------------------------------\n', ref_fields_T2{i});
      if isnumeric(ref_info_T2.dcminfo.(ref_fields_T2{i}))
        fprintf(fid, 'Expected value: %.2f, Observed value: %.2f\n', [ref_info_T2.dcminfo.(ref_fields_T2{i})'; scan_info_T2.dcminfo.(ref_fields_T2{i})']);
      else
        fprintf(fid, 'Expected value: %s, Observed value: %s\n', ref_info_T2.dcminfo.(ref_fields_T2{i}), scan_info_T2.dcminfo.(ref_fields_T2{i}));
      end
      fprintf(fid, '\n');
    end
  catch ME
    fprintf('%s\n', ME.message);
    continue
  end
end

if exist('scan_info_dwi_rev', 'var')
  ref_fields_rev = fieldnames(ref_info_dwi_rev.dcminfo);
  for i = 1:length(ref_fields_rev)
    try
      if ~isequal( ref_info_dwi_rev.dcminfo.(ref_fields_rev{i}), scan_info_dwi_rev.dcminfo.(ref_fields_rev{i}) )
	problem_flag = 1;
	fprintf(fid, 'Scan parameter deviation for RSI reverse scan: %s ---------------------------------------\n', ref_fields_rev{i});
	if isnumeric(ref_info_dwi_rev.dcminfo.(ref_fields_rev{i}))
	  fprintf(fid, 'Expected value: %.2f, Observed value: %.2f\n', [ref_info_dwi_rev.dcminfo.(ref_fields_rev{i})'; scan_info_dwi_rev.dcminfo.(ref_fields_rev{i})']);
	else
	  fprintf(fid, 'Expected value: %s, Observed value: %s\n', ref_info_dwi_rev.dcminfo.(ref_fields_rev{i}), scan_info_dwi_rev.dcminfo.(ref_fields_rev{i}));
	end
	fprintf(fid, '\n');
      end
    catch
      fprintf('%s\n', ME.message);
      continue
    end
  end
end


% Check number of slices
slices_fwd = scan_info_dwi_fwd.slices;
slices_T2 = scan_info_T2.slices;

if slices_fwd ~= ref_info_dwi_fwd.slices;
   problem_flag = 1;
   fprintf(fid, 'Slice number deviation: Forward RSI scan ---------------------------------------\n');
   fprintf(fid, 'Expected number: %d, Observed number: %d\n\n', ref_info_dwi_fwd.slices, slices_fwd);
end

if slices_T2 ~= ref_info_T2.slices;
   problem_flag = 1;
   fprintf(fid, 'Slice number deviation: Axial T2 scan ---------------------------------------\n');
   fprintf(fid, 'Expected number: %d, Observed number: %d\n\n', ref_info_T2.slices, slices_T2);
end

if exist('scan_info_dwi_rev', 'var')
  slices_rev = scan_info_dwi_rev.slices;
  if slices_rev ~= ref_info_dwi_rev.slices;
    problem_flag = 1;
    fprintf(fid, 'Slice number deviation: Reverse RSI scan ---------------------------------------\n');
    fprintf(fid, 'Expected number: %d, Observed number: %d\n\n', ref_info_dwi_rev.slices, slices_rev);
  end
end


% Check scan coverage
M_fwd = scan_info_dwi_fwd.M;
M_T2 = scan_info_T2.M;
diff = abs(M_fwd(:,4) - M_T2(:,4));

if any(diff > 1)
   problem_flag = 1;
   fprintf(fid, 'Scan coverage deviation: Ax T2 and forward RSI scans do not match ---------------------------------------\n');
   fprintf(fid, 'T2 coordinates:          ');
   fprintf(fid, '%.4f ', M_T2(:,4));
   fprintf(fid, '\n');
   fprintf(fid, 'Forward RSI coordinates: ');
   fprintf(fid, '%.4f ', M_fwd(:,4));
   fprintf(fid, '\n\n');
end

if exist('scan_info_dwi_rev', 'var')
   M_rev = scan_info_dwi_rev.M;
   diff = abs(M_rev(:,4) - M_T2(:,4));
   if any(diff > 1)
      problem_flag = 1;
      fprintf(fid, 'Scan coverage deviation: Ax T2 and reverse RSI scans do not match ---------------------------------------\n');
      fprintf(fid, 'T2 coordinates:          ');
      fprintf(fid, '%.4f ', M_T2(:,4));
      fprintf(fid, '\n');
      fprintf(fid, 'Reverse RSI coordinates: ');
      fprintf(fid, '%.4f ', M_rev(:,4));
      fprintf(fid, '\n\n');
   end

   diff = abs(M_fwd(:,4) - M_rev(:,4));
   if any(diff > 1)
      problem_flag = 1;
      fprintf(fid, 'Scan coverage deviation: Forward and reverse RSI scans do not match ---------------------------------------\n');
      fprintf(fid, 'Forward RSI coordinates: ');
      fprintf(fid, '%.4f ', M_fwd(:,4));
      fprintf(fid, '\n');
      fprintf(fid, 'Reverse RSI coordinates: ');
      fprintf(fid, '%.4f ', M_rev(:,4));
      fprintf(fid, '\n\n');
   end
end


if problem_flag == 0
   fprintf(fid, 'No protocol deviations detected\n');
end

fclose(fid);

end
