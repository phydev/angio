echo "f2py2.7 compiling get_data_m.F90"
f2py2.7 -c get_data_m.F90 -m get_data_m --f90flags='-fdefault-real-8'
echo "get_data_m.so library was created"
