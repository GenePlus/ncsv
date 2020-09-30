set -e
scptpath=`dirname $0`
scptpath=`realpath $scptpath`
cd $scptpath

# build extension
rm -rf build dist lib/*.egg-info

cd extension/sv_helper
sh autogen.sh
./configure
make -j 10
make install

cd $scptpath
pip install 'cython>=0.28.2' --user # cython 0.28.2 is test, add quotes for min version requirement
rm -f lib/sv_helper.so && python setup.py build_ext --inplace

# test
cd $scptpath
python setup.py nosetests
if [ $? -ne 0 ]; then exit; fi;

# remove files, nosetests also generate *.pyc
rm -rf build dist lib/*.egg-info
find lib -type f -name '*.pyc'|xargs rm

# dos2unix
for f in `find lib -name "*.py" -type f ` ;do dos2unix $f ;done;

# zip variantCalling.zip -r \
    # lib/{dbwrapper,nccnv,ncsv,asmgraph,*.py,*.sh} setup* \
    # inst.sh README.md
