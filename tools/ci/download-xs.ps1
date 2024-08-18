
# Download HDF5 data
if (-not (Test-Path "$Env:USERPROFILE\nndc_hdf5")) {
  Invoke-WebRequest https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz -OutFile hdf5.xz
  tar -xvzf hdf5.xz
}

# Download ENDF/B-VII.1 distribution
$Env:ENDF = "$Env:USERPROFILE\endf-b-vii.1"
if (-not (Test-Path $Env:ENDF)) {
  Invoke-WebRequest https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz -OutFile endf.xs
  tar -xvzf endf.xz
}
