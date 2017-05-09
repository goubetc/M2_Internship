function easy_test
  cd /users/goubet/Documents/Test/SIMPLE_M/
  'hello_world'
  file_id = fopen('OUTPUT/test.txt', 'w');

  fdisp(file_id, 'hello_world');

  fclose(file_id);

end
