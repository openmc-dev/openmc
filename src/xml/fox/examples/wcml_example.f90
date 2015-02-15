 program wcml_example

   use FoX_wcml

  integer,  parameter ::  sp = selected_real_kind(6,30)
  integer,  parameter ::  dp = selected_real_kind(14,100)
!
! NB normally you will be writting to the xml file
! from mulitple fortran files/subroutines, therefore
! type(xmlf_t) :: myfile   (below)
! would normally need to be treated as a global
! variable, either in a module or a common block.
! 
  type(xmlf_t) :: myfile
!

  integer           :: num, na
  character(len=2)  :: elements(3)
  real(kind=dp)      :: coords(3,3)
  real(kind=dp)      :: adp, matrix(3,3)

  character(len=10) :: filename

  data coords(1:3,1)/0.0d0, 0.0d0, 0.0d0/
  data coords(1:3,2)/0.5d0, 0.5d0, 0.5d0/
  data coords(1:3,3)/0.4d0, 0.4d0, 0.4d0/

  matrix = 3.141_dp

  adp=1.234567890
  na=3
  elements(1) = 'Ca'
  elements(2) = 'Si'
  elements(3) = 'O'
  num = 20

  filename = 'mytest.cml'
  
  call cmlBeginFile(myfile, filename=filename, unit=-1)

  call cmlAddNamespace(myfile, prefix="myDict", URI="http://www.example.com/dict")

  call cmlStartCml(myfile)

  ! Add parameter
  call cmlStartParameterList(xf=myfile, title="Input parameters")
  call cmlAddParameter(xf=myfile, name='inputSize', value='True')
  call cmlAddParameter(xf=myfile, name='inputSize', value=.True.)
  call cmlAddParameter(xf=myfile, name='inputSize', value=3, units="cmlUnits:Angstrom")
  call cmlAddParameter(xf=myfile, name='inputSize', value=3.0, units="cmlUnits:Angstrom")
  call cmlAddParameter(xf=myfile, name='inputSize', value=3.0d0, units="cmlUnits:Angstrom")
  call cmlAddParameter(xf=myfile, name='inputSize', value=(3.0d0, 0.0d0), units="cmlUnits:Angstrom")
  call cmlEndParameterList(xf=myfile)

  call cmlStartPropertyList(xf=myfile, title="Scalars")
  call cmlAddProperty(xf=myfile, title='inputSize', value='True')
  call cmlAddProperty(xf=myfile, title='inputSize', value=.True.)
  call cmlAddProperty(xf=myfile, title='inputSize', value=3, units="cmlUnits:Angstrom")
  call cmlAddProperty(xf=myfile, title='inputSize', value=3.0, units="cmlUnits:Angstrom")
  call cmlAddProperty(xf=myfile, title='inputSize', value=3.0d0, units="cmlUnits:Angstrom")
  call cmlAddProperty(xf=myfile, title='inputSize', value=(3.0d0, 0.0d0), units="cmlUnits:Angstrom")
  call cmlEndPropertyList(xf=myfile)

  call cmlStartPropertyList(xf=myfile, title="Scalars")
  call cmlAddProperty(xf=myfile, title='inputArray', value=(/'one  ', 'two  ', 'three'/))
  call cmlAddProperty(xf=myfile, title='inputArray', value=(/.true., .false./))
  call cmlAddProperty(xf=myfile, title='inputArray', value=(/1, 2/), units="cmlUnits:Angstrom")
  call cmlAddProperty(xf=myfile, title='inputArray', value=(/1.0, 2.0/), units="cmlUnits:Angstrom")
  call cmlAddProperty(xf=myfile, title='inputArray', value=(/1.0d0, 2.0d0/), units="cmlUnits:Angstrom")
  call cmlAddProperty(xf=myfile, title='inputArray', value=(/(1.0d0,0.0d0), (2.0d0,0.0d0)/), units="cmlUnits:Angstrom")
  call cmlEndPropertyList(xf=myfile)

  call cmlStartPropertyList(xf=myfile, title="Scalars")

  call cmlStartPropertyList(myfile, dictref="castep:kpt_band")
  call cmlAddProperty(myfile, 1, dictref="castep:spin", &
             units="castepunits:spin")
  call cmlAddKpoint(myfile, (/1.0d0, 2.0d0, 3.0d0/), &
              weight=1.0d0, id="okpt:")
  call cmlAddProperty(myfile, (/1.0, 2.0/), title="Eigenvalues", &
              dictref="castep:eigenenergies", &
              units="castepunits:")
  call cmlEndPropertyList(myfile)

  call cmlAddProperty(xf=myfile, &
    value=matrix,&
    title='Elastic Constant Matrix',units='units:GPa')

  call cmlEndPropertyList(xf=myfile)

  ! Add molecule
  call cmlAddMolecule(xf=myfile, elements=elements,coords=coords)                               
  ! Add molecule output in short style
  call cmlAddMolecule(xf=myfile, elements=elements,coords=coords, style='xyz3')

  ! Add molecule output in short style in user supplied format
  call cmlAddMolecule(xf=myfile, elements=elements,coords=coords, style='xyz3', fmt='r6')

  call cmlAddCrystal(xf=myfile, a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0)

  ! End and Close
  call cmlEndCml(myfile)
  call cmlFinishFile(myfile)
 
end program wcml_example
