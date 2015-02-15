#!/bin/sh 

INCFLAGS=`../../FoX-config --fcflags --wxml`

rm -f test.xml
make test_xml_Openfile_1.exe
./test_xml_Openfile_1.exe > /dev/null 2>&1
if [ -f test.xml ]; then
   echo we passed
   PASSED=yes
else
   echo we failed
   PASSED=no
fi
./passed.sh $PASSED "file creation with default unit number"
rm -f test.xml

# check for creation with non-default file handle.
make test_xml_Openfile_2.exe
./test_xml_Openfile_2.exe> /dev/null 2>&1
if [ -f test.xml ]; then
   PASSED=yes
else
   PASSED=no
fi
./passed.sh $PASSED "file creation with non-default unit number"
rm -f test.xml

# check that we replace existing file by default
echo TESTINPUT > test.xml
make test_xml_Openfile_1.exe
./test_xml_Openfile_1.exe> /dev/null 2>&1
if grep -v TESTINPUT test.xml >/dev/null; then
   PASSED=yes
else
   PASSED=no
fi
./passed.sh $PASSED "file creation replacement by default"
rm -f test.xml

# check that we create a new file when we ask to.
# This should fail ...
echo TESTINPUT > test.xml
make test_xml_Openfile_3.exe
./test_xml_Openfile_3.exe >/dev/null 2>&1
if [ $? -ne 0 ]; then
   PASSED=yes
else
   PASSED=no
fi
./passed.sh $PASSED "file creation without replacement"
rm -f test.xml

echo TESTINPUT > test.xml
make test_xml_Openfile_1.exe
./test_xml_Openfile_1.exe> /dev/null 2>&1
if grep -v TESTINPUT test.xml >/dev/null; then
   PASSED=yes
else
   PASSED=no
fi
./passed.sh $PASSED "file creation specifying REPLACE"
rm -f test.xml

make test_xml_Openfile_4.exe
rm -f test.xml
./test_xml_Openfile_4.exe> test.out 2>&1
if diff test.out test_xml_Openfile_4.out | grep "^>" > /dev/null; then
   PASSED=no
else
   PASSED=yes
fi
./passed.sh $PASSED "Trying to create file twice"
rm -f test.xml
