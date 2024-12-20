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

nixy_string = string;

end
