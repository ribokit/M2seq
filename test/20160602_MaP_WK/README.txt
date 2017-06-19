# to run the tests, copy the example_dir to a new test directory,
# do the run with *master*, and save the directory into test_master
cp -rf example_files/ test/
git checkout master
cd test/
source README
cd ../
mv test test_master

# then rerun after checking out your favorite branch.
cp -rf example_files/ test/
git checkout <<my new branch>>
cd test/
source README
cd ../
diff -r test_master/ test/
