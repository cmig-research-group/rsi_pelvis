function params = eval_file(path_file, var_name)

if ~exist('var_name', 'var')
  var_name = 'params';
end
  
fid = fopen(path_file);
line = fgetl(fid);
while ischar(line)
  fprintf('%s\n', line);
  eval(line);
  line = fgetl(fid);
end
fclose(fid);

params = eval(var_name);

end
