function nixy_string = nixify(string)

% Convert a string to something that will be easier to use in a *nix shell

string = strrep(string, ' ', '_');
string = strrep(string, '(', '');
string = strrep(string, ')', '');
string = strrep(string, '''', '');
string = strrep(string, '*', '');
string = strrep(string, '[', '');
string = strrep(string, ']', '');
string = strrep(string, '$', '');
string = strrep(string, ',', '');

% Convert string to basic ASCII 
unicode_vals = double(string);
unicode_vals(unicode_vals>127) = 63;

nixy_string = char(unicode_vals);

end
