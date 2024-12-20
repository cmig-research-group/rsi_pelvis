function dcminfo_fixed = fix_impax_dcm_tags(dcminfo)

load('private_tag_ref.mat');

tags2fix = {'Private_0043_1039', 'Private_0019_10bb', 'Private_0019_10bc', 'Private_0019_10bd'};

fields_in = fieldnames(dcminfo);
dcminfo_fixed = dcminfo;
for i = 1:length(fields_in)
  if isempty(strfind(fields_in{i},'Private')) == 0

    type_in = class( dcminfo.(fields_in{i}) );
    if isfield(private_tag_ref, fields_in{i}) == 1
      type_ref = private_tag_ref.(fields_in{i});

      if strcmp(type_ref, type_in) == 0 %any(strcmp(tags2fix, fields_in{i})) && strcmp(type_ref, type_in) == 0

	dcminfo_fixed.(fields_in{i}) = char( dcminfo_fixed.(fields_in{i})' );

	if strcmp(fields_in{i},'Private_0043_1039') == 1
	  stringtag = dcminfo_fixed.(fields_in{i});
	  [startIndex, endIndex] = regexp(stringtag, '[1-9]\d{1,3}\\');
	  bval_char = stringtag(startIndex:endIndex-1);
	  if isempty(bval_char) == 0
	    dcminfo_fixed.(fields_in{i}) = bval_char;
	  else
	    dcminfo_fixed.(fields_in{i}) = '0';
	  end
	
	elseif strcmp(type_ref, 'double')
	  dcminfo_fixed.(fields_in{i}) = str2double(dcminfo_fixed.(fields_in{i}));

	end

      end

    end

  end
end

end
