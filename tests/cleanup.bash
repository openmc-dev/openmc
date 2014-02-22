# This simple script ensures that all binary
# output files have been deleted in all the
# folders. This can occur if a previous error
# occurred and the test suite was rerun without
# deleting left over binary files. This will
# cause an assertion error in some of the 
# tests.
for i in ./*; do
  if [ -d "$i" ]; then
    cd $i
    rm -f *.binary *.h5
    cd ..
  fi
done
