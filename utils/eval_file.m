function params = eval_file(path_file)

fid = fopen(path_file);
line = fgetl(fid);
while ischar(line)
  fprintf('%s\n', line);
  eval(line);
  line = fgetl(fid);
end
fclose(fid);

end
