#!/bin/sh -e

VN=$(grep dummy common/FoX_common.F90 | awk '{print $6}' | sed -e s/-dummy//)

if test x$VN != x\'$1\'; then
  echo Wrong version number
  exit 1 
fi

git archive --format=tar --prefix=FoX-$1/ HEAD | gzip -9 > ../FoX-$1-full.tar.gz
openssl dgst -md5 ../FoX-$1-full.tar.gz > ../FoX-$1-digest
openssl dgst -sha1 ../FoX-$1-full.tar.gz >> ../FoX-$1-digest

mkdir tmpFoX
cd tmpFoX
tar xzf ../../FoX-$1-full.tar.gz
(
  cd FoX-$1
  make cutdown
)
tar czf ../../FoX-$1.tar.gz FoX-$1
openssl dgst -md5 ../../FoX-$1.tar.gz >> ../../FoX-$1-digest
openssl dgst -sha1 ../../FoX-$1.tar.gz >> ../../FoX-$1-digest
rm -rf FoX-$1

for i in wxml wcml wkml sax dom
do
  tar xzf ../../FoX-$1-full.tar.gz
  (
    cd FoX-$1
    (
      cd config; \
       sed -e "s/CUTDOWN_TARGET=.*/CUTDOWN_TARGET=$i/" configure.ac > configure.ac.tmp ; \
       mv configure.ac.tmp configure.ac ; \
       make
    )
    make cutdown-$i
  )
  tar czf FoX-$1-$i.tar.gz FoX-$1
  openssl dgst -md5 FoX-$1-$i.tar.gz >> ../../FoX-$1-digest
  openssl dgst -sha1 FoX-$1-$i.tar.gz >> ../../FoX-$1-digest
  rm -rf FoX-$1
done

cp FoX-$1-* ../..

cd ..
rm -rf tmpFoX

