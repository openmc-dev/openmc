!(doc)TITLE{LIB\_VTK\_IO}
!(doc)SUBTITLE{VTK InPut/OutPut Fortran Library}
!(doc)VERSION{0.2}
!(doc)AUTHOR{Stefano Zaghi}
!(doc)COAUTHORS{Enrico Cavallini and Renato N. Elias}
!(doc)DATE{07-01-2008}

!!\newcommand{\LIBVTKIO}{\MaiuscolettoBS{LIB\_VTK\_IO }}

!(doc)header

!(doc)titlepage

!!\chapter{Acknowledgements}
!!\label{cap:Acknowledgements}
!!
!!I am very grateful to Renato N. Elias: whitout his support \LIBVTKIO would not born. As a matter of facts \LIBVTKIO is a
!!collection of his tips. Despite the fact that Renato does not write the code he is a \virgo{moral co-author} of the code.
!!
!!I thank Enrico Cavallini for his help in debugging the code. He also develop the MS Windows version
!!of \LIBVTKIO and he is the first co-author that I found.
!!
!!Finally I thank the ParaView mailing list for the great support of its members.
!!
!!\chapter{Introduction}
!!\label{cap:Introduction}
!!\begin{epigraphs}
!! \qitem{\emph{I have not failed. I've just found ten thousand ways that don't work.}}{{\sc Thomas Edison}}
!!\end{epigraphs}
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf L}}{IB\_VTK\_IO} is a Fortran library to write and read
!!(actually only to write) data conforming the VTK standard both binary and ascii. Even though there are many
!!wrappers/porting of the VTK source code (C++ code), there is not a fortran one. This library is not a porting
!!or a wrapper of the VTK code, but it only an exporter/importer of the VTK data format written in pure fortran
!!language (standard Fortran 95 with some extensions of non standard Fortran 2003) that can be used by fortran
!!coders (yes, there are still a lot of these brave coders...) without mixing fortran with C++ language.
!!
!!The library is still in developing and testing, this is first usable release, but there are not all the features
!!of the stable release (the importer is totaly absent and the exporter is not complete). Surely there are a lot of
!!bugs and the progamming style is not the best, but the exporter is usable for the 90\% of the VTK data format.
!!
!!The \LIBVTKIO is an open source project, it is distribuited under the GPL v3 (see appendix \ref{cap:GPL}). Anyone is
!!interest to use, to develop or contribuite to \LIBVTKIO is welcome.
!!
!!\section*{VTK Standard}
!!\label{sec:VTK Standard}
!!VTK, Visualization Toolkit, is an open source software that provides a powerful framework for the computer grafich, for
!!the images processing and for 3D rendering. It is widely used in the world and so it has a very large comunity of users;
!!besides the Kitware\footnote{The Kitware homepage can be found here: \href{http://public.kitware.com}{http://public.kitware.com}.}
!!company provides professional support. The toolkit is written in C++ and a lot of porting/wrappers for Tcl/Tk,
!!Java and Python are provided; unlucky there aren't wrappers for Fortran.
!!
!!Because of its good features the VTK toolkit has been used to develop a large set of open source programs. For my work
!!the most important family of programs is the scientific visualization programs. A lot of high-quality scientific visualization
!!tool are available on the web but for me the best is ParaView: I think that it is one of the best scintific visualization
!!program in the world and it is open source! Paraview is based on VTK.
!!
!!\section*{Paraview}
!!\label{sec:Paraview}
!!ParaView\footnote{The ParaView homepage can be found here: \href{http://www.paraview.org}{http://www.paraview.org}.}
!!is an open source software voted to scientific visualization and able to use the power of parallel architectures. It
!!has an architecture client-server in order to make easy the remote visualization of very large set of data. Because it is based
!!on VTK it inherits all VTK features. ParaView is very useful for Computational Fluid Dynamics visualizations because it provides
!!powerful post-processing tools; it provides a very large set of importers for the most used format like Plot3D and HDF (the list
!!is very large). It is easy to extend ParaView because it supports all the scripting language supported by VTK.
!!
!!\section*{LIB\_VTK\_IO}
!!\label{sec:LIB_VTK_IO}
!!Even though the VTK toolkit is written in C++ and so it is possible to use it in mixed fortran/c++ code this is not the easiest
!!way. Fortran is still the best language for high performance computing for scientific purpose, like CFD computing. It necessary a
!!tool to deal with VTK standard directly by fortran code. The library \LIBVTKIO was made to fill this empty: it is a simple
!!fortran module able to export native fortran data into VTK data format and to import VTK data into a fortran code (actually this
!!feature is missing), both in ascii and binary file format.
!!
!!The library provides an automatic way to deal with VTK data format: all the formatting processes is nested into the library and
!!the users comunicate with it by a simple API passing only native fortran data (native fortran scalar, vector and matrix).
!!
!!The library \LIBVTKIO is distribuited under the GNU GPL v3 license (see appendix \ref{cap:GPL}). Beyond to the source code there
!!are some precompiled binaries for GNU-Linux (amd x86, amd x86\_64, intel x86, intel x86\_64) and WindowsXP (amd x86, intel x86).
!!
!!Actually the library is still in developing/testing phase (a lot of features are missing); this is not a stable release, but the
!!exporter is quite complete and its API is quite stable. The exporter is usable and I use it for my work.
!!
!!\chapter{News and Changes}
!!\label{cap:NewsChanges}
!!
!!\section*{Version v0.2}
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} version v0.2 is the second testing release. From version v0.1 there are
!!only minor changes; this new version does not introduce new features and does not fix bugs: it is a simple code-cleaning. The
!!character variables are now case-insensitive; the names of some variables have been changed. The comments have been translated in
!!English (very poor translation...).
!!
!!\begin{boxred}{List of changes from v0.1}
!!\begin{enumerate1Red}
!!  \item variable {\color{Maroon}formato} is changed in {\color{Maroon}output\_format} and now appears only in VTK\_INI and
!!        in VTK\_INI\_XML.
!!  \item variable {\color{Maroon}nomefile} is changed in {\color{Maroon}filename}.
!!  \item variable {\color{Maroon}titolo} is changed in {\color{Maroon}title}.
!!  \item variable {\color{Maroon}topologia} is changed in {\color{Maroon}mesh\_topology} and now appears only in VTK\_INI and
!!        in VTK\_INI\_XML.
!!  \item variable {\color{Maroon}NCelle} is changed in {\color{Maroon}NC}.
!!  \item variable {\color{Maroon}Nnodi} is changed in {\color{Maroon}NN}.
!!  \item variable {\color{Maroon}tipo} in VTK\_CON and VTK\_CON\_XML is changed in {\color{Maroon}cell\_type}.
!!  \item variable {\color{Maroon}tipo} in VTK\_DAT and VTK\_DAT\_XML is changed in {\color{Maroon}var\_location}.
!!  \item variable {\color{Maroon}azione} in VTK\_DAT\_XML is changed in {\color{Maroon}var\_block\_action}.
!!  \item variable {\color{Maroon}tipo} in VTK\_VAR is changed in {\color{Maroon}vec\_type}.
!!  \item variable {\color{Maroon}nomevar} is changed in {\color{Maroon}varname}.
!!\end{enumerate1Red}
!!\end{boxred}
!!
!!The only relevant news in the v0.2 version is about this guide: now the guide is integrated in the code. The code has particular
!!comments: if the code is processed by the program FortranDOC\footnote{FortranDOC is an open-source fortran code available at:
!!\href{http://www.paraview.org}{http://www.paraview.org}. This code processing a free-format fortran code generates a corresponding
!!pretty-latex documentation file of the code structure.} a latex source of this guide will be made; compiling the latex file with
!!\virgo{pdflatex} you will obtain this guide in PDF.
!!
!!\section*{Version v0.1}
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} version v0.1 is the first testing release. There are not news and changes.
!!
!!\mainmatter
!!
!!\part{Compile and Install LIB\_VTK\_IO}
!!\label{part:Compile and Install}
!!
!!\chapter{Compile LIB\_VTK\_IO}
!!\label{cap:Compiling Library}
!!\minitoc
!!\vspace*{3mm}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} \LIBVTKIO is open source and so anyone is encouraged to use the source
!!code and to \virgo{patch} it.
!!
!!The code is written in Fortran: the standard adopted is the Fortran 95 standard that is a minor upgrade to the Fortran 90 standard
!!and that is widely supported by the almost all compilers actually available. Unluckily Fortran 95 does not allow the creation of
!!C-binary file (Fortran inserts some bytes before and after each records despite the C standard) that is the standard adopted by
!!VTK. Therefore in order to create binary files that are compatible whit VTK standard the only way is to use a non-standard 95
!!instructions. At today only Fortran 2003 can create C-binary file, but there are not any compilers that completely implement this
!!standard. In the next year (2008) maybe a new minor upgrade of Fortran standard (unofficial named Fortran 2008) will be born
!!and so the support to Fortran 2003/2008 probably will be improved. Luckily we need to use only some features of fortran 2003
!!that are supported by many compilers.
!!
!!The Fortran 2003 instructions are focused on the opening of the binary file, in particular in the functions
!!\MaiuscolettoBS{VTK\_INI} and \MaiuscolettoBS{VTK\_INI\_XML}. In these functions there are opening instructions like the following:
!!
!!\begin{boxred}{Fortran 2003 instructions}
!!\begin{verbatim}
!!open(unit       = ...,           &
!!     file       = ...,           &
!!     form       = ...,           &
!!     access     = ...,           &
!!     action     = ...,           &
!!     convert    = 'BIG_ENDIAN',  &
!!     recordtype = 'STREAM',      &
!!     buffered   = 'YES',         &
!!     iostat     = ...)
!!\end{verbatim}
!!\end{boxred}
!!
!!The specifiers \MaiuscolettoBS{convert}, \MaiuscolettoBS{recordtype} and \MaiuscolettoBS{buffered} are non standard for Fortran 95.
!!The \MaiuscolettoBS{buffered} specifier is not necessary and so can be commented or eliminated. The specifiers
!!\MaiuscolettoBS{convert} and \MaiuscolettoBS{recordtype} are instead necessary to write binary file but can be replaced by other
!!specifiers/instructions. In particular an alternative is opening the file with the specifier
!!\MaiuscolettoBS{form = BINARY}\footnote{Remember that also the value \MaiuscolettoBS{BINARY} for form specifier is non standard
!!for Fortran 95.} and using a compiler's option\footnote{Each compilers adopt differents option to achieve conversion of bytes
!!order (if it allows conversion). See the user guide of your compiler. Intel Fortran allows the conversion both by open specifier
!!and by compiling option.} to ensure the \MaiuscolettoBS{BIG\_ENDIAN} encoding. \MaiuscolettoBS{BIG\_ENDIAN} encoding is strictly
!!necessary only for legacy binary file; for XML binary file one can choice also the \MaiuscolettoBS{LITTLE\_ENDIAN} and so the
!!conversion is not necessary.
!!
!!Actually there is also another instruction that is non-standard for Fortran 95: the instruction \MaiuscolettoBS{sizeof}. This
!!instruction is used to comptuing the number of bytes of the saved data in the XML binary files. Maybe there are others
!!alternatives that are Fortran 95 compatible but at the moment I have not the time to implement them.
!!
!!Before you compile \LIBVTKIO ensure that your compiler allows these Fortran 2003 extensions. I use the Intel Fortran
!!Compiler\footnote{\href{http://www.intel.com}{http://www.intel.com}.} that is free for non-commercial use and it has a strong
!!support for Fortran 2003.
!!
!!\section{Compile under GNU/Linux}
!!\label{sec:CompileLinux}
!!
!!\LIBVTKIO can be compiled as a stand-alone library or it can be integrated directly in your code. It is a self-contained module
!!that can be safely included into others fortran codes. There are no any advices for compile \LIBVTKIO excluding the above non
!!standard instructions.
!!
!!For the GNU/Linux users there is available a makefile already set to compile \LIBVTKIO both as static and dynamic library with
!!Intel Fortran. The makefile has only one option: \MaiuscolettoBS{SHARED}. This variable (default set to \virgo{no}) can assume
!!two values:
!!\begin{enumerate1Blu}
!!\item {\color{RoyalBlue}\MaiuscolettoBS{no}}:  makefile creates a \MaiuscolettoBS{static} library
!!\item {\color{RoyalBlue}\MaiuscolettoBS{yes}}: makefile creates a \MaiuscolettoBS{dynamic} library
!!\end{enumerate1Blu}
!!
!!\section{Compile under MS Windows}
!!\label{sec:CompileWin}
!!
!!For MS Windows users there is not any support at the moment. As soon as I have the time I will make available a MS Visual Studio
!!Project to compile \LIBVTKIO with Intel Visual Fortran for Windows.
!!
!!\clearpage
!!
!!\chapter{Install and Link (Pre)Compiled LIB\_VTK\_IO}
!!\label{cap:Install and Linking}
!!\minitoc
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} \LIBVTKIO is distribuited in two different version (other than source code): the first is a static linking version (extensions are \emph{.a} and \emph{.lib}) and the second is dynamic linking version (extensions are \emph{.so} and \emph{.dll}). The use of these two version is different and it depends on the OS used. The library is been tested only on GNU/Linux (several different distro) and on MS Windows (Windows XP).
!!
!!The library is distribuited with two different archive: \MaiuscolettoBS{LIB\_VTK\_IO-bin-x.x}.tar for GNU/Linux systems and
!!\MaiuscolettoBS{LIB\_VTK\_IO-bin-x.x}.zip for MS Windows systems. Into the archives there is the source code of the library
!!(\MaiuscolettoBS{LIB\_VTK\_IO}.f90), there are both static and dynamic version of the librabry and there is also this guide
!!(\MaiuscolettoBS{LIB\_VTK\_IO\_Guide}.pdf).
!!
!!\section{GNU/Linux}
!!
!!\subsection{Static Library}
!!\label{sec:Linux Static}
!!The static version of the precompiled library (\MaiuscolettoBS{LIB\_VTK\_IO}.a) does not require any kind of installations. It is
!!enough to link against it in the linking phase. It is important to use the interface module \emph{lib\_vtk\_io.mod} distribuited
!!with the library: this is the interface of the subroutines and functions that constitute the library.
!!
!!To use the functions and subroutines of the library it is mandatory to \MaiuscolettoBS{USE} the module. Suppose one has a program
!!(or subprogram) named \emph{test} that use the library; the correct \MaiuscolettoBS{USE} is:
!!
!!\begin{boxred}{The \LIBVTKIO must to be loaded with the USE statement}
!!\begin{verbatim}
!!program test
!!USE LIB_VTK_IO
!!...
!!...
!!...
!!endprogram test
!!\end{verbatim}
!!\end{boxred}
!!
!!With the instruction \verb|USE LIB\_VTK\_IO| the program \emph{test} can use the functions and subroutines of the library. To
!!compile, without link, this code one must give the module interface \emph{lib\_vtk\_io.mod} to the compiler:
!!
!!\begin{boxred}{Static Compiling Phase}
!!\begin{verbatim}
!!ifort -c lib_vtk_io.mod test.f90 -o test.o
!!\end{verbatim}
!!\end{boxred}
!!
!!In this example \emph{ifort} is the Intel Fortran Compiler\footnote{Da aggiungere.} and the \verb|-c| flag compiles preventing
!! linking; the compiler must \virgo{see} the module interface: the file \emph{lib\_vtk\_io.mod} must be placed in a folder visible
!!by the compiler.
!!
!!In the linking phase one simply give the library to the compiler:
!!
!!\begin{boxred}{Static Linking Phase}
!!\begin{verbatim}
!!ifort test.o LIB_VTK_IO.a -o test.out
!!\end{verbatim}
!!\end{boxred}
!!
!!The library must be placed in a folder visible by the compiler.
!!
!!\subsection{Dynamic Library}
!!\label{sec:Linux Dynamic}
!!The dynamic version of the precompiled library must be installed. The operating system must know where is the library so it is
!!necessary to install the library in a folder where the OS search its shared objects. In the most of the GNU/Linux distro the
!!folder \emph{/usr/lib/} is scanned to find shared objects. After you have copied the \MaiuscolettoBS{LIB\_VTK\_IO}.so file in
!!this folder, update the list of the shared objects with the command \verb|ldconfig -v| and the OS is ready to use the library.
!!
!!After you set your OS the compiling and linking phase is identical to the previous (remember to you the module interface at
!!the compiling phase). The only difference is to use the dynamic library at the linking phase:
!!
!!\begin{boxred}{Dynamic Linking Phase}
!!\begin{verbatim}
!!ifort test.o LIB_VTK_IO.so -o test.out
!!\end{verbatim}
!!\end{boxred}
!!
!!\section{MS Windows}
!!
!!Unluckily for MS Windows there is not any support at the moment. As soon as I have the time, I make some instructions on how
!!use \LIBVTKIO with MS Visual Studio and Intel Visual Fortran for MS Windows.
!!
!!\chapter{LIB\_VTK\_IO Programming Style}
!!\label{cap:Programming Style}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf A}}{ll} the \LIBVTKIO functions are \MaiuscolettoBS{4-byte integer functions}:
!!the output of these functions is an integer that is $0$ if the function calling has been done right while it is $> 0$  if some
!!errors occur (the error handling is only at its embryonal phase). Therefore the functions calling must be done in the following
!!way:
!!
!!\begin{boxred}{Functions Calling}
!!\begin{verbatim}
!!...
!!integer(4):: E_IO
!!...
!!E_IO = VTK_INI(....
!!...
!!\end{verbatim}
!!\end{boxred}
!!
!!The \LIBVTKIO programming style is based on two main principles: \MaiuscolettoBS{portable kind-precision} of reals and integers
!!variables and \MaiuscolettoBS{dynamic dispatching}. In the appendix \ref{cap:kind precision} and \ref{cap:Dynamic Dispatching}
!!there are more details about these choices. I just remark some consequences of these choices. Using \MaiuscolettoBS{dynamic
!!dispatching} the \LIBVTKIO has a simple API. The user calls a generic procedure (VTK\_INI, VTK\_GEO,...) and the library,
!!depending on the type of the inputs passed, calls the correct internal function (i.e. VTK\_GEO for 8-byte real type if the input
!!passed is 8-byte real type). By this interface only few functions are used whitout the necessity of calling a different function
!!for every different inputs type. \MaiuscolettoBS{Dynamic dispatching} is valid also for the different kind of topology and
!!variables-data-dimensions; the function VTK\_GEO is the same for all topologies, just the inputs passed to the functions change
!!as the topology changes. Also the dimensions of variables-data use the \MaiuscolettoBS{dynamic dispatching}: the function
!!(VTK\_VAR) used to save vectorial data is identical to the one used for scalar data, depending on the dimensions of the data
!!\LIBVTKIO calls the correct internal function. \MaiuscolettoBS{Dynamic dispatching} is based on the internal kind-precision
!!selecting convention: Fortran 90/95 standard has some useful functions to achive the portability of reals and integers precision
!!and \LIBVTKIO uses these functions to define portable kind-precision; because it is important to make portable the code on
!!different architectures I suggest to use this programming style.
!!
!!The data handled by \LIBVTKIO can be classified into two main categories:
!!
!!\begin{enumerate1Red}
!!\item Geometric Data. These are the geometric informations of the mesh and they can be of different kind and different number
!!      depending on the topology choiced. The mesh points coordinates type must be of 4-byte real type or 8-byte real type.
!!\item Variable Data. These are the scalar or vectorial variables appended to the mesh points (both at the cell-nodes and the
!!      cell-centers of the mesh). The type of these data can be of 8-byte real type, 4-byte real type and 4-byte integer type
!!      (for the XML output there are also the 8-byte integer type, 2-byte integer type and 1-byte integer type).
!!\end{enumerate1Red}
!!
!!In the following chapters theare the details of \LIBVTKIO API.
!!
!!\part{LIB\_VTK\_IO API}
!!\label{part:LIBVTKIO API}
module vtk_writer 
!----------------------------------------------------------------------------------------------------------------------------------
!!\LIBVTKIO is a library of functions for Input and Output pure fortran data (both ascii and binary) in VTK format.
!!
!!The VTK standard can be separated into two main catagories: the \MaiuscolettoBS{VTK Legacy Standard} and the
!!\MaiuscolettoBS{VTK XML Standard}. The latter is more powerful and will has a stronger support from VTk comunity than legacy
!!standard; XML file format would to be preferred despite the legacy one.
!!
!!At the present only a few functions of the final library have been implemented. The InPut functions are totaly absent, but the
!!OutPut functions are almost complete (the \virgo{polydata} functions are the only missing).
!!
!!The functions actually present are:
!!
!!\begin{boxred}{Functions for Legacy VTK file format}
!!\begin{enumerate1Red}
!! \item \MaiuscolettoS{VTK\_INI}
!! \item \MaiuscolettoS{VTK\_GEO}
!! \item \MaiuscolettoS{VTK\_CON}
!! \item \MaiuscolettoS{VTK\_DAT}
!! \item \MaiuscolettoS{VTK\_VAR}
!! \item \MaiuscolettoS{VTK\_END}
!!\end{enumerate1Red}
!!\end{boxred}
!!
!!\begin{boxred}{Functions for XML VTK file format}
!!\begin{enumerate1Red}
!! \item \MaiuscolettoS{VTK\_INI\_XML}
!! \item \MaiuscolettoS{VTK\_GEO\_XML}
!! \item \MaiuscolettoS{VTK\_CON\_XML}
!! \item \MaiuscolettoS{VTK\_DAT\_XML}
!! \item \MaiuscolettoS{VTK\_VAR\_XML}
!! \item \MaiuscolettoS{VTK\_END\_XML}
!!\end{enumerate1Red}
!!\end{boxred}
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
! functions for VTK LEGACY
public:: VTK_INI
public:: VTK_GEO
public:: VTK_CON
public:: VTK_DAT
public:: VTK_VAR
public:: VTK_END
! functions for VTK XML
public:: VTK_INI_XML
public:: VTK_GEO_XML
public:: VTK_CON_XML
public:: VTK_DAT_XML
public:: VTK_VAR_XML
public:: VTK_END_XML
! portable kind-precision
public:: R16P, FR16P
public:: R8P,  FR8P
public:: R4P,  FR4P
public:: R_P,  FR_P
public:: I8P,  FI8P
public:: I4P,  FI4P
public:: I2P,  FI2P
public:: I1P,  FI1P
public:: I_P,  FI_P
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! overloading of VTK_GEO
interface VTK_GEO
  module procedure VTK_GEO_UNST_R8, & ! real(R8P) UNSTRUCTURED_GRID
                   VTK_GEO_UNST_R4, & ! real(R4P) UNSTRUCTURED_GRID
                   VTK_GEO_STRP_R8, & ! real(R8P) STRUCTURED_POINTS
                   VTK_GEO_STRP_R4, & ! real(R4P) STRUCTURED_POINTS
                   VTK_GEO_STRG_R8, & ! real(R8P) STRUCTURED_GRID
                   VTK_GEO_STRG_R4, & ! real(R4P) STRUCTURED_GRID
                   VTK_GEO_RECT_R8, & ! real(R8P) RECTILINEAR_GRID
                   VTK_GEO_RECT_R4    ! real(R4P) RECTILINEAR_GRID
endinterface
! overloading of VTK_VAR
interface VTK_VAR
  module procedure VTK_VAR_SCAL_R8, & ! real(R8P)    scalar
                   VTK_VAR_SCAL_R4, & ! real(R4P)    scalar
                   VTK_VAR_SCAL_I4, & ! integer(I4P) scalar
                   VTK_VAR_VECT_R8, & ! real(R8P)    vectorial
                   VTK_VAR_VECT_R4, & ! real(R4P)    vectorial
                   VTK_VAR_VECT_I4, & ! integer(I4P) vectorial
                   VTK_VAR_TEXT_R8, & ! real(R8P)    vectorial (texture)
                   VTK_VAR_TEXT_R4    ! real(R4P)    vectorial (texture)
endinterface
! overloading of VTK_GEO_XML
interface VTK_GEO_XML
  module procedure VTK_GEO_XML_STRG_R4, & ! real(R4P) StructuredGrid
                   VTK_GEO_XML_STRG_R8, & ! real(R8P) StructuredGrid
                   VTK_GEO_XML_RECT_R8, & ! real(R8P) RectilinearGrid
                   VTK_GEO_XML_RECT_R4, & ! real(R4P) RectilinearGrid
                   VTK_GEO_XML_UNST_R8, & ! real(R8P) UnstructuredGrid
                   VTK_GEO_XML_UNST_R4, & ! real(R4P) UnstructuredGrid
                   VTK_GEO_XML_CLOSEP     ! closing tag "Piece" function
endinterface
! overloading of VTK_VAR_XML
interface VTK_VAR_XML
  module procedure VTK_VAR_XML_SCAL_R8, & ! real(R8P)    scalar
                   VTK_VAR_XML_SCAL_R4, & ! real(R4P)    scalar
                   VTK_VAR_XML_SCAL_I8, & ! integer(I8P) scalar
                   VTK_VAR_XML_SCAL_I4, & ! integer(I4P) scalar
                   VTK_VAR_XML_SCAL_I2, & ! integer(I2P) scalar
                   VTK_VAR_XML_SCAL_I1, & ! integer(I1P) scalar
                   VTK_VAR_XML_VECT_R8, & ! real(R8P)    vectorial
                   VTK_VAR_XML_VECT_R4, & ! real(R4P)    vectorial
                   VTK_VAR_XML_VECT_I8, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I4, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I2, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I1    ! integer(I4P) vectorial
endinterface
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!!\LIBVTKIO has a small set of internal variables and parameters some of which have public visibility.
!!
!!The \LIBVTKIO uses a partable kind parameters for real and integer variables. The following are the kind parameters used: these
!!parameters are public and their use is strong encouraged.
!!
!!Real precision definitions:
!!
integer, parameter:: R16P = selected_real_kind(33,4931) ! 33  digits, range $[\pm 10^{-4931}  ,\pm 10^{+4931}   -1]$
integer, parameter:: R8P  = selected_real_kind(15,307)  ! 15  digits, range $[\pm 10^{-307}~~ ,\pm 10^{+307}~~  -1]$
integer, parameter:: R4P  = selected_real_kind(6,37)    ! 6~~~digits, range $[\pm 10^{-37}~~~~,\pm 10^{+37}~~~~ -1]$
integer, parameter:: R_P  = R8P                         ! default real precision
!!Integer precision definitions:
!!
integer, parameter:: I8P  = selected_int_kind(18)       ! range $[-2^{63} ,+2^{63}  -1]$
integer, parameter:: I4P  = selected_int_kind(9)        ! range $[-2^{31} ,+2^{31}  -1]$
integer, parameter:: I2P  = selected_int_kind(4)        ! range $[-2^{15} ,+2^{15}  -1]$
integer, parameter:: I1P  = selected_int_kind(2)        ! range $[-2^{7}~~,+2^{7}~~ -1]$
integer, parameter:: I_P  = I4P                         ! default integer precision
!!
!!Besides the kind parameters there are also the format parameters useful for writing in a well-ascii-format numeric variables.
!!Also these parameters are public.
!!
!! Real output formats:
!!
character(10), parameter:: FR16P = '(E41.33E4)'         ! R16P  output format
character(10), parameter:: FR8P  = '(E23.15E3)'         ! R8P   output format
character(9),  parameter:: FR4P  = '(E14.6E2)'          ! R4P   output format
character(10), parameter:: FR_P  = '(E23.15E3)'         ! R\_P  output format
!! Integer output formats:
!!
character(5), parameter:: FI8P  = '(I21)'               ! I8P  output format
character(5), parameter:: FI4P  = '(I12)'               ! I4P  output format
character(4), parameter:: FI2P  = '(I7)'                ! I2P  output format
character(4), parameter:: FI1P  = '(I5)'                ! I1P  output format
character(5), parameter:: FI_P  = '(I12)'               ! I\_P output format
!!
!!\LIBVTKIO uses a small set of internal variables that are private (not accessible from the outside). The following are
!! private variables:
!!
integer(I4P), parameter:: maxlen       = 500         ! max number of characters os static string
character(1), parameter:: end_rec      = char(10)    ! end-character for binary-record finalize
integer(I4P), parameter:: f_out_ascii  = 0           ! ascii-output-format parameter identifier
integer(I4P), parameter:: f_out_binary = 1           ! binary-output-format parameter identifier
integer(I4P)::            f_out        = f_out_ascii ! current output-format (initialized to ascii format)
character(len=maxlen)::   topology                   ! mesh topology
integer(I4P)::            Unit_VTK                   ! internal logical unit
integer(I4P)::            Unit_VTK_Append            ! internal logical unit for raw binary XML append file
integer(I4P)::            N_Byte                     ! number of byte to be written/read
real(R8P)::               tipo_R8                    ! prototype of R8P real
real(R4P)::               tipo_R4                    ! prototype of R4P real
integer(I8P)::            tipo_I8                    ! prototype of I8P integer
integer(I4P)::            tipo_I4                    ! prototype of I4P integer
integer(I2P)::            tipo_I2                    ! prototype of I2P integer
integer(I1P)::            tipo_I1                    ! prototype of I1P integer
integer(I4P)::            ioffset                    ! offset pointer
integer(I4P)::            indent                     ! indent pointer
!----------------------------------------------------------------------------------------------------------------------------------

!!In the following chapters there is the API reference of all functions of \LIBVTKIO.
contains
  !!\chapter{Auxiliary functions}
  !!\minitoc
  !!\vspace*{8mm}
  !!
  !!\LIBVTKIO uses two auxiliary functions that are not connected with the VTK standard. These functions are private and so they
  !!cannot be called outside the library.
  function GetUnit() result(Free_Unit)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The GetUnit function is used for getting a free logic unit. The users of \LIBVTKIO does not know which is
  !!the logical unit: \LIBVTKIO handels this information without boring the users. The logical unit used is safe-free: if the
  !!program calling \LIBVTKIO has others logical units used \LIBVTKIO will never use these units, but will choice one that is free.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P):: Free_Unit ! free logic unit
  integer(I4P):: n1        ! counter
  integer(I4P):: ios       ! inquiring flag
  logical(4)::   lopen     ! inquiring flag
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  !!The following is the code snippet of GetUnit function: the units 0, 5, 6, 9 and all non-free units are discarded.
  !!
  !(\doc)codesnippet
  Free_Unit = -1_I4P                                      ! initializing free logic unit
  n1=1_I4P                                                ! initializing counter
  do
    if ((n1/=5_I4P).AND.(n1/=6_I4P).AND.(n1/=9_I4P)) then
      inquire (unit=n1,opened=lopen,iostat=ios)           ! verify logic units
      if (ios==0_I4P) then
        if (.NOT.lopen) then
          Free_Unit = n1                                  ! assignment of free logic
          return
        endif
      endif
    endif
    n1=n1+1_I4P                                           ! updating counter
  enddo
  return
  !(doc/)codesnippet
  !!GetUnit function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope to public.
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction GetUnit

  function Upper_Case(string)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The Upper\_Case function converts the lower case characters of a string to upper case one. \LIBVTKIO uses this function in
  !!order to achieve case-insensitive: all character variables used within \LIBVTKIO functions are pre-processed by
  !!Uppper\_Case function before these variables are used. So the users can call \LIBVTKIO functions whitout pay attention of the
  !!case of the kwywords passed to the functions: calling the function VTK\_INI with the string \code{E_IO = VTK_INI('Ascii',...)}
  !!or with the string  \code{E_IO = VTK_INI('AscII',...)} is equivalent.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(len=*), intent(IN):: string     ! string to be converted
  character(len=len(string))::   Upper_Case ! converted string
  integer::                      n1         ! characters counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  !!The following is the code snippet of Upper\_Case function.
  !!
  !(\doc)codesnippet
  Upper_Case = string
  do n1=1,len(string)
    select case(ichar(string(n1:n1)))
    case(97:122)
      Upper_Case(n1:n1)=char(ichar(string(n1:n1))-32) ! Upper case conversion
    endselect
  enddo
  return
  !(doc/)codesnippet
  !!Upper\_Case function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope
  !!to public.
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction Upper_Case

  !!\chapter{VTK LEGACY functions}
  !!\minitoc
  !!\vspace*{8mm}
  !!
  function VTK_INI(output_format,filename,title,mesh_topology) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The VTK\_INI function is used for initializing file. This function must be the first to be called.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: output_format ! output format: ASCII or BINARY
  character(*), intent(IN):: filename      ! name of file
  character(*), intent(IN):: title         ! title
  character(*), intent(IN):: mesh_topology ! mesh topology
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!The VTK\_INI variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}output\_format}] indicates the \virgo{format} of output file. It can assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{ascii} (it is case insensitive) $\rightarrow$ creating an ascii output file.
  !!  \item \emph{binary} (it is case insensitive) $\rightarrow$ creating a binary (big\_endian encoding) output file.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}filename}] contains the name (with its path) of the output file.
  !! \item[{\color{RoyalBlue}title}] contains the title of the VTK dataset.
  !! \item[{\color{RoyalBlue}topology}] indicates the topology of the mesh and can assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{STRUCTURED\_POINTS}.
  !!  \item \emph{STRUCTURED\_GRID}.
  !!  \item \emph{UNSTRUCTURED\_GRID}.
  !!  \item \emph{RECTILINEAR\_GRID}.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_INI calling:
  !!
  !!\begin{boxred}{VTK\_INI Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_INI('Binary','example.vtk','VTK legacy file','UNSTRUCTURED_GRID')
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!\noindent Note that the \virgo{.vtk} extension is necessary in the file name.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  topology = trim(mesh_topology)
  Unit_VTK=GetUnit()
  select case(trim(Upper_Case(output_format)))
  case('ASCII')
    f_out = f_out_ascii
    open(unit     = Unit_VTK,       &
         file     = trim(filename), &
         form     = 'FORMATTED',    &
         access   = 'SEQUENTIAL',   &
!        action   = 'WRITE')        &
         iostat   = E_IO)

    ! writing header of file
!   write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'# vtk DataFile Version 3.0'
!   write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)trim(title)
!   write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)trim(Upper_Case(output_format))
!   write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'DATASET '//trim(topology)
  case('BINARY')
    f_out = f_out_binary
    open(unit       = Unit_VTK,       &
         file       = trim(filename), &
         form       = 'UNFORMATTED',  &
         access     = 'SEQUENTIAL',   &
         action     = 'WRITE',        &
         iostat     = E_IO)
    ! writing header of file
    write(unit=Unit_VTK,iostat=E_IO)'# vtk DataFile Version 3.0'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)trim(title)//end_rec
    write(unit=Unit_VTK,iostat=E_IO)trim(Upper_Case(output_format))//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'DATASET '//trim(topology)//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_INI

  !!\section{VTK\_GEO}
  !!
  !!VTK\_GEO is an interface to 8 different functions; there are 2 functions for each 4 different topologies actually supported:
  !!one function for mesh coordinates with R8P precision and one for mesh coordinates with R4P precision.
  !!This function must be called after VTK\_INI. It saves the mesh geometry. The inputs that must be passed change depending on
  !!the topologies choiced. Not all VTK topologies have been implemented (\virgo{polydata} topologies are absent). The signatures
  !!for all implemented topologies are now reported.
  !!
  !!\subsection{VTK\_GEO STRUCTURED POINTS}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Structured Points Signature}]
  !! function VTK_GEO(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{structured points} is useful for structured grid with uniform discretization steps.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Structured Points Variables}]
  !!integer(I4P),     intent(IN):: Nx   ! number of nodes in x direction
  !!integer(I4P),     intent(IN):: Ny   ! number of nodes in y direction
  !!integer(I4P),     intent(IN):: Nz   ! number of nodes in z direction
  !!real(R8P or R4P), intent(IN):: X0   ! x coordinate of origin
  !!real(R8P or R4P), intent(IN):: Y0   ! y coordinate of origin
  !!real(R8P or R4P), intent(IN):: Z0   ! z coordinate of origin
  !!real(R8P or R4P), intent(IN):: Dx   ! space step in x
  !!real(R8P or R4P), intent(IN):: Dy   ! space step in y
  !!real(R8P or R4P), intent(IN):: Dz   ! space step in z
  !!integer(I4P)::                 E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X0,Y0,Z0,Dx,Dy,Dz} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be sure
  !!that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO structured point variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}Nx}] indicates the number of nodes in $X$ direction.
  !! \item[{\color{RoyalBlue}Ny}] indicates the number of nodes in $Y$ direction.
  !! \item[{\color{RoyalBlue}NZ}] indicates the number of nodes in $Z$ direction.
  !! \item[{\color{RoyalBlue}X0}] indicates the $X$ value of coordinates system origin. It is a scalar.
  !! \item[{\color{RoyalBlue}Y0}] indicates the $Y$ value of coordinates system origin. It is a scalar.
  !! \item[{\color{RoyalBlue}Z0}] indicates the $Z$ value of coordinates system origin. It is a scalar.
  !! \item[{\color{RoyalBlue}Dx}] indicates the uniform grid step discretization in $X$ direction. It is a scalar.
  !! \item[{\color{RoyalBlue}Dy}] indicates the uniform grid step discretization in $Y$ direction. It is a scalar.
  !! \item[{\color{RoyalBlue}DZ}] indicates the uniform grid step discretization in $Z$ direction. It is a scalar.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO structured point calling:
  !!
  !!\begin{boxred}{VTK\_GEO Structured Points Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4):: Nx,Ny,Nz
  !!real(8):: X0,Y0,Z0
  !!real(8):: Dx,Dy,Dz
  !!...
  !!E_IO = VTK_GEO(Nx,Ny,Nz, &
  !!               X0,Y0,Z0,Dx,Dy,Dz)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_GEO STRUCTURED GRID}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Structured Grid Signature}]
  !!function VTK_GEO(Nx,Ny,Nz,NN,X,Y,Z) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{structured grid} is useful for structured grid with non-uniform discretization steps.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Structured Grid Variables}]
  !!integer(I4P),     intent(IN):: Nx      ! number of nodes in x direction
  !!integer(I4P),     intent(IN):: Ny      ! number of nodes in y direction
  !!integer(I4P),     intent(IN):: Nz      ! number of nodes in z direction
  !!integer(I4P),     intent(IN):: NN      ! number of all nodes
  !!real(R8P or R4P), intent(IN):: X(1:NN) ! x coordinates
  !!real(R8P or R4P), intent(IN):: Y(1:NN) ! y coordinates
  !!real(R8P or R4P), intent(IN):: Z(1:NN) ! z coordinates
  !!integer(I4P)::                 E_IO    ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X,Y,Z} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be
  !!sure that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO structured grid variables have the following meaning:
  !!
  !!\begin{description}
  !!  \item[{\color{RoyalBlue}Nx}] indicates the number of nodes in $X$ direction.
  !!  \item[{\color{RoyalBlue}Ny}] indicates the number of nodes in $Y$ direction.
  !!  \item[{\color{RoyalBlue}NZ}] indicates the number of nodes in $Z$ direction.
  !!  \item[{\color{RoyalBlue}NN}] indicates the number of all nodes, $NN= Nx\cdot Ny\cdot Nz$.
  !!  \item[{\color{RoyalBlue}X}] contains the $X$ coordinates values of all nodes. It is a vector of $[1:NN]$.
  !!  \item[{\color{RoyalBlue}Y}] contains the $Y$ coordinates values of all nodes. It is a vector of $[1:NN]$.
  !!  \item[{\color{RoyalBlue}Z}] contains the $Z$ coordinates values of all nodes. It is a vector of $[1:NN]$.
  !!  \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO structured grid calling:
  !!
  !!\begin{boxred}{VTK\_GEO Structured Grid Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: Nx=10,Ny=10,Nz=10
  !!integer(4), parameter:: Nnodi=Nx*Ny*Nz
  !!real(8):: X(1:Nnodi),Y(1:Nnodi),Z(1:Nnodi)
  !!...
  !!E_IO = VTK_GEO(Nx,Ny,Nz,Nnodi,X,Y,Z)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_GEO RECTILINEAR GRID}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Rectilinear Grid Signature}]
  !!function VTK_GEO(Nx,Ny,Nz,X,Y,Z) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{rectilinear grid} is useful for structured grid with non-uniform discretization steps even
  !!in generalized coordinates.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Rectilinear Grid Signature}]
  !!integer(I4P),     intent(IN):: Nx      ! number of nodes in x direction
  !!integer(I4P),     intent(IN):: Ny      ! number of nodes in y direction
  !!integer(I4P),     intent(IN):: Nz      ! number of nodes in z direction
  !!real(R8P or R4P), intent(IN):: X(1:Nx) ! x coordinates
  !!real(R8P or R4P), intent(IN):: Y(1:Ny) ! y coordinates
  !!real(R8P or R4P), intent(IN):: Z(1:Nz) ! z coordinates
  !!integer(I4P)::                 E_IO    ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X,Y,Z} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be
  !!sure that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO rectilinear grid variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}Nx}] indicates the number of nodes in $X$ direction.
  !! \item[{\color{RoyalBlue}Ny}] indicates the number of nodes in $Y$ direction.
  !! \item[{\color{RoyalBlue}Nz}] indicates the number of nodes in $Z$ direction.
  !! \item[{\color{RoyalBlue}X}] contains the $X$ coordinates values of nodes. It is a vector of $[1:Nx]$.
  !! \item[{\color{RoyalBlue}Y}] contains the $Y$ coordinates values of nodes. It is a vector of $[1:Ny]$.
  !! \item[{\color{RoyalBlue}Z}] contains the $Z$ coordinates values of nodes. It is a vector of $[1:Nz]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO rectilinear grid calling:
  !!
  !!\begin{boxred}{VTK\_GEO Rectilinear Grid Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: Nx=10,Ny=20,Nz=30
  !!real(4):: X(1:Nx),Y(1:Ny),Z(1:Nz)
  !!...
  !!E_IO = VTK_GEO(Nx,Ny,Nz,X,Y,Z)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_GEO UNSTRUCTURED GRID}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Unstructured Grid Signature}]
  !!function VTK_GEO(Nnodi,X,Y,Z) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{unstructured grid} is necessary for unstructured grid, the most general mesh format. This
  !!topology is also useful for scructured mesh in order to save only a non-structured clip of mesh.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Unstructured Grid Variables}]
  !!integer(I4P),     intent(IN):: NN      ! number of nodes
  !!real(R8P or R4P), intent(IN):: X(1:NN) ! x coordinates of all nodes
  !!real(R8P or R4P), intent(IN):: Y(1:NN) ! y coordinates of all nodes
  !!real(R8P or R4P), intent(IN):: Z(1:NN) ! z coordinates of all nodes
  !!integer(I4P)::                 E_IO    ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X,Y,Z} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be
  !!sure that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO unstructured grid variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NN}] indicates the number of all nodes.
  !! \item[{\color{RoyalBlue}X}] contains the $X$ coordinates values of nodes. It is a vector of $[1:NN]$.
  !! \item[{\color{RoyalBlue}Y}] contains the $Y$ coordinates values of nodes. It is a vector of $[1:NN]$.
  !! \item[{\color{RoyalBlue}Z}] contains the $Z$ coordinates values of nodes. It is a vector of $[1:NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO unstructured grid calling:
  !!
  !!\begin{boxred}{VTK\_GEO Unstructured Grid Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NN=100
  !!real(4):: X(1:NN),Y(1:NN),Z(1:NN)
  !!...
  !!E_IO = VTK_GEO(NN,X,Y,Z)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!In order to use the \virgo{unstructured grid} it is necessary to save also the \virgo{connectivity} of the grid. The
  !!connectivity must be saved with the function \MaiuscolettoBS{VTK\_CON}.
  !!
  !(\doc)skippedblock
  function VTK_GEO_STRP_R8(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = STRUCTURED\_POINTS (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
  integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
  integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
  real(R8P),    intent(IN):: X0        ! x coordinate of origin
  real(R8P),    intent(IN):: Y0        ! y coordinate of origin
  real(R8P),    intent(IN):: Z0        ! z coordinate of origin
  real(R8P),    intent(IN):: Dx        ! space step in x direction
  real(R8P),    intent(IN):: Dy        ! space step in y direction
  real(R8P),    intent(IN):: Dz        ! space step in z direction
  integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer  ! buffer string
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,fmt='(A,3'//FR8P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=Unit_VTK,fmt='(A,3'//FR8P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
  case(f_out_binary)
    write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,3'//FR8P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,3'//FR8P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRP_R8

  function VTK_GEO_STRP_R4(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = STRUCTURED\_POINTS (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
  integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
  integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
  real(R4P),    intent(IN):: X0        ! x coordinate of origin
  real(R4P),    intent(IN):: Y0        ! y coordinate of origin
  real(R4P),    intent(IN):: Z0        ! z coordinate of origin
  real(R4P),    intent(IN):: Dx        ! space step in x direction
  real(R4P),    intent(IN):: Dy        ! space step in y direction
  real(R4P),    intent(IN):: Dz        ! space step in z direction
  integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer  ! buffer string
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,fmt='(A,3'//FR4P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=Unit_VTK,fmt='(A,3'//FR4P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
  case(f_out_binary)
    write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,3'//FR4P//')', iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,3'//FR4P//')', iostat=E_IO)'SPACING ',Dx,Dy,Dz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRP_R4

  function VTK_GEO_STRG_R8(Nx,Ny,Nz,NN,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = STRUCTURED\_GRID (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: Nx       ! number of nodes in x direction
  integer(I4P), intent(IN):: Ny       ! number of nodes in y direction
  integer(I4P), intent(IN):: Nz       ! number of nodes in z direction
  integer(I4P), intent(IN):: NN       ! number of all nodes
  real(R8P),    intent(IN):: X(1:NN)  ! x coordinates
  real(R8P),    intent(IN):: Y(1:NN)  ! y coordinates
  real(R8P),    intent(IN):: Z(1:NN)  ! z coordinates
  integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer ! buffer string
  integer(I4P)::             n1       ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=Unit_VTK,fmt='(3'//FR8P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRG_R8

  function VTK_GEO_STRG_R4(Nx,Ny,Nz,NN,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = STRUCTURED\_GRID (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: Nx       ! number of nodes in x direction
  integer(I4P), intent(IN):: Ny       ! number of nodes in y direction
  integer(I4P), intent(IN):: Nz       ! number of nodes in z direction
  integer(I4P), intent(IN):: NN       ! number of all nodes
  real(R4P),    intent(IN):: X(1:NN)  ! x coordinates
  real(R4P),    intent(IN):: Y(1:NN)  ! y coordinates
  real(R4P),    intent(IN):: Z(1:NN)  ! z coordinates
  integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer ! buffer string
  integer(I4P)::             n1       ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=Unit_VTK,fmt='(3'//FR4P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRG_R4

  function VTK_GEO_RECT_R8(Nx,Ny,Nz,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = RECTILINEAR\_GRID (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
  integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
  integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
  real(R8P),    intent(IN):: X(1:Nx)   ! x coordinates
  real(R8P),    intent(IN):: Y(1:Ny)   ! y coordinates
  real(R8P),    intent(IN):: Z(1:Nz)   ! z coordinates
  integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer  ! buffer string
  integer(I4P)::             n1        ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' double'
    write(unit=Unit_VTK,fmt=FR8P,              iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' double'
    write(unit=Unit_VTK,fmt=FR8P,              iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' double'
    write(unit=Unit_VTK,fmt=FR8P,              iostat=E_IO)(Z(n1),n1=1,Nz)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' double'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' double'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' double'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(Z(n1),n1=1,Nz)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_RECT_R8

  function VTK_GEO_RECT_R4(Nx,Ny,Nz,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = RECTILINEAR\_GRID (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: Nx        ! number of nodes in x direction
  integer(I4P), intent(IN):: Ny        ! number of nodes in y direction
  integer(I4P), intent(IN):: Nz        ! number of nodes in z direction
  real(R4P),    intent(IN):: X(1:Nx)   ! x coordinates
  real(R4P),    intent(IN):: Y(1:Ny)   ! y coordinates
  real(R4P),    intent(IN):: Z(1:Nz)   ! z coordinates
  integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer  ! buffer string
  integer(I4P)::             n1        ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' float'
    write(unit=Unit_VTK,fmt=FR4P,              iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' float'
    write(unit=Unit_VTK,fmt=FR4P,              iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' float'
    write(unit=Unit_VTK,fmt=FR4P,              iostat=E_IO)(Z(n1),n1=1,Nz)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' float'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' float'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' float'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(Z(n1),n1=1,Nz)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_RECT_R4

  function VTK_GEO_UNST_R8(NN,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = UNSTRUCTURED\_GRID (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NN        ! number of nodes
  real(R8P),    intent(IN):: X(1:NN)   ! x coordinates of all nodes
  real(R8P),    intent(IN):: Y(1:NN)   ! y coordinates of all nodes
  real(R8P),    intent(IN):: Z(1:NN)   ! z coordinates of all nodes
  integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer  ! buffer string
  integer(I4P)::             n1        ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=Unit_VTK,fmt='(3'//FR8P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_UNST_R8

  function VTK_GEO_UNST_R4(NN,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = UNSTRUCTURED\_GRID (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NN        ! number of nodes
  real(R4P),    intent(IN):: X(1:NN)   ! x coordinates of all nodes
  real(R4P),    intent(IN):: Y(1:NN)   ! y coordinates of all nodes
  real(R4P),    intent(IN):: Z(1:NN)   ! z coordinates of all nodes
  integer(I4P)::             E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer  ! buffer string
  integer(I4P)::             n1        ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=Unit_VTK,fmt='(3'//FR4P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=Unit_VTK,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                       iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_UNST_R4
  !(doc/)skippedblock

  function VTK_CON(NC,connect,cell_type) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!This function \MaiuscolettoBS{must} be used when unstructured grid is used. It saves the connectivity of the unstructured
  !!mesh.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC              ! number of cells
  integer(I4P), intent(IN):: connect(:)      ! mesh connectivity
  integer(I4P), intent(IN):: cell_type(1:NC) ! VTK cell type
  integer(I4P)::             E_IO            ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer        ! buffer string
  integer(I4P)::             ncon            ! dimension of connectivity vector 
  !!The VTK\_CON variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC}] indicates the number of all cells.
  !! \item[{\color{RoyalBlue}connect}] contains the connectivity of the mesh. It is a vector.
  !! \item[{\color{RoyalBlue}cell\_type}] contains the type of every cells. It is a vector of $[1:NC]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The vector \MaiuscolettoBS{connect} must follow the VTK legacy standard. It is passed as \MaiuscolettoBS{assumed-shape} array
  !!because its dimensions is related to the mesh dimensions in a complex way. Its dimensions can be calculated by the following
  !!equation:
  !!
  !!\begin{equation}
  !!dc = NC + \sum\limits_{i = 1}^{NC} {nvertex_i }
  !!\label{eq:connectivity dimensions}
  !!\end{equation}
  !!
  !!\noindent where $dc$ is connectivity vector dimension and $nvertex_i$ is the number of vertices of $i^{th}$ cell. The VTK
  !!legacy standard for the mesh connectivity is quite obscure at least at first sight. It is more simple analizing an example.
  !!Suppose we have a mesh composed by 2 cells, one hexahedron (8 vertices) and one pyramid with square basis (5 vertices); suppose
  !!that the basis of pyramid is constitute by a face of the hexahedron and so the two cells share 4 vertices. The equation
  !!\ref{eq:connectivity dimensions} gives $dc=2+8+5=15$; the connectivity vector for this mesh can be:
  !!
  !!\begin{boxred}{Connectivity vector example for VTK legacy standard}
  !!\begin{verbatim}
  !!! first cell
  !!connect(1)  = 8  => number of vertices of 1° cell
  !!connect(2)  = 0  => identification flag of 1° vertex of 1° cell
  !!connect(3)  = 1  => identification flag of 2° vertex of 1° cell
  !!connect(4)  = 2  => identification flag of 3° vertex of 1° cell
  !!connect(5)  = 3  => identification flag of 4° vertex of 1° cell
  !!connect(6)  = 4  => identification flag of 5° vertex of 1° cell
  !!connect(7)  = 5  => identification flag of 6° vertex of 1° cell
  !!connect(8)  = 6  => identification flag of 7° vertex of 1° cell
  !!connect(9)  = 7  => identification flag of 8° vertex of 1° cell
  !!! second cell
  !!connect(10) = 5  => number of vertices of 2° cell
  !!connect(11) = 0  => identification flag of 1° vertex of 2° cell
  !!connect(12) = 1  => identification flag of 2° vertex of 2° cell
  !!connect(13) = 2  => identification flag of 3° vertex of 2° cell
  !!connect(14) = 3  => identification flag of 4° vertex of 2° cell
  !!connect(15) = 8  => identification flag of 5° vertex of 2° cell
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\noindent Note that the first 4 identification flags of pyramid vertices as the same of the first 4 identification flags of
  !!the hexahedron because the two cells share this face. It is also important to note that the identification flags start
  !!form $0$ value: this is impose to the VTK standard. The function VTK\_CON does not calculate the connectivity vector: it
  !!writes the connectivity vector conforming the VTK standard, but does not calculate it. In the future release of \LIBVTKIO will
  !!be included a function to calculate the connectivity vector.
  !!
  !!The vector variable \MaiuscolettoBS{tipo} must conform the VTK standard \footnote{See the file VTK-Standard at the Kitware
  !!homepage.}. It contains the \emph{type} of each cells. For the above example this vector is:
  !!
  !!\begin{boxred}{Cell-Type vector example for VTK legacy standard}
  !!\begin{verbatim}
  !!tipo(1) = 12  => VTK hexahedron type of 1° cell
  !!tipo(2) = 14  => VTK pyramid type of 2° cell
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!The following is an example of VTK\_CON calling:
  !!
  !!\begin{boxred}{VTK\_CON Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC=2
  !!integer(4), parameter:: Nvertex1=8
  !!integer(4), parameter:: Nvertex2=5
  !!integer(4), parameter:: dc=NC+Nvertex1+Nvertex2
  !!integer(4)::            connect(1:dc)
  !!integer(4)::            cell_type(1:NC)
  !!...
  !!E_IO = VTK_CON(NC,connect,cell_type)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  ncon = size(connect,1)
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,2'//FI4P//')',iostat=E_IO)'CELLS ',NC,ncon
    write(unit=Unit_VTK,fmt=FI4P,             iostat=E_IO)connect
    write(unit=Unit_VTK,fmt='(A,'//FI4P//')', iostat=E_IO)'CELL_TYPES ',NC
    write(unit=Unit_VTK,fmt=FI4P,             iostat=E_IO)cell_type
  case(f_out_binary)
    write(s_buffer,     fmt='(A,2'//FI4P//')',iostat=E_IO)'CELLS ',NC,ncon
    write(unit=Unit_VTK,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                      iostat=E_IO)connect
    write(unit=Unit_VTK,                      iostat=E_IO)end_rec
    write(s_buffer,     fmt='(A,'//FI4P//')', iostat=E_IO)'CELL_TYPES ',NC
    write(unit=Unit_VTK,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                      iostat=E_IO)cell_type
    write(unit=Unit_VTK,                      iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_CON

  function VTK_DAT(NC_NN,var_location) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!This function \MaiuscolettoBS{must} be called before saving the data related to geometric mesh. This function initializes the
  !!saving of data variables indicating the \emph{type} of variables that will be saved.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes of field
  character(*), intent(IN):: var_location ! location of saving variables: cell for cell-centered, node for node-centered
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  !!The VTK\_DAT variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of {\color{RoyalBlue}tipo}.
  !! \item[{\color{RoyalBlue}var\_location}] contains the location-type of variables that will be saved after VTK\_DAT. It is a scalar and cab assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{cell} (it is case insensitive) $\rightarrow$ variables will be cell-centered.
  !!  \item \emph{node} (it is case insensitive) $\rightarrow$ variables will be node-centered.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Of course a single file can contain both cell and node centered variables; in this case the VTK\_DAT function must be called two times, before saving cell-centered variables and before saving node-centered variables.
  !!
  !!The following is an example of VTK\_DAT calling:
  !!
  !!\begin{boxred}{VTK\_DAT Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_DAT(50,'node')
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      write(unit=Unit_VTK,fmt='(A,'//FI4P//')',iostat=E_IO)'CELL_DATA ',NC_NN
    case('NODE')
      write(unit=Unit_VTK,fmt='(A,'//FI4P//')',iostat=E_IO)'POINT_DATA ',NC_NN
    endselect
  case(f_out_binary)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      write(s_buffer,     fmt='(A,'//FI4P//')',iostat=E_IO)'CELL_DATA ',NC_NN
      write(unit=Unit_VTK,                     iostat=E_IO)trim(s_buffer)//end_rec
    case('NODE')
      write(s_buffer,     fmt='(A,'//FI4P//')',iostat=E_IO)'POINT_DATA ',NC_NN
      write(unit=Unit_VTK,                     iostat=E_IO)trim(s_buffer)//end_rec
    endselect
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_DAT

  !!\section{VTK\_VAR}
  !!
  !!VTK\_VAR is an interface to 8 different functions; there are 3 functions for scalar variables, 3 functions for vectorial
  !!variables and 2 function texture variables.
  !!This function saves the data variables related to geometric mesh. The inputs that must be passed change depending on the data
  !!variables type.
  !!
  !!\subsection{VTK\_VAR SCALAR DATA}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Scalar Data Signature}]
  !!function VTK_VAR(formato,NC_NN,varname,var) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!This kind of call is used to save scalar data.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Scalar Data Variables}]
  !!integer(I4P),                     intent(IN):: NC_NN        ! number of nodes or cells
  !!character(*),                     intent(IN):: varname      ! variable name
  !!real(R8P or R4P) or integer(I4P), intent(IN):: var(1:NC_NN) ! variable to be saved
  !!integer(I4P)::                                 E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_VAR variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of
  !!                                  {\color{RoyalBlue}tipo} passed to VTK\_DAT.
  !! \item[{\color{RoyalBlue}varname}] contains the name attribuited the variable saved.
  !! \item[{\color{RoyalBlue}var}] contains the values of variables in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Note that the variables \texttt{var} can be passed both as 8-byte real kind, 4-byte real kind and 4-byte integer; the
  !!dynamic displacement interface will call the correct function.
  !!
  !!The following is an example of VTK\_VAR scalar data calling:
  !!
  !!\begin{boxred}{VTK\_VAR Scalar Data Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC_NN=100
  !!real(4)::               var(1:NC_NN)
  !!...
  !!E_IO = VTK_VAR(NC_NN,'Scalar Data',var)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_VAR REAL VECTORIAL DATA}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Real Vectorial Data Signature}]
  !!function VTK_VAR(tipo,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!This kind of call is used to save real vectorial data.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Real Vectorial Data Variables}]
  !!character(*),     intent(IN):: vec_type      ! vector type: vect = generic vector , norm = normal vector
  !!integer(I4P),     intent(IN):: NC_NN         ! number of nodes or cells
  !!character(*),     intent(IN):: varname       ! variable name
  !!real(R8P or R4P), intent(IN):: varX(1:NC_NN) ! x component of vector
  !!real(R8P or R4P), intent(IN):: varY(1:NC_NN) ! y component of vector
  !!real(R8P or R4P), intent(IN):: varZ(1:NC_NN) ! z component of vector
  !!integer(I4P)::                 E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_VAR variables have the following meaning:
  !!
  !!\begin{description}
  !! \item [{\color{RoyalBlue}tipo}] indicates the type of vector. It can assume the following value:
  !! \begin{enumerateABlu}
  !!  \item \emph{vect} $\rightarrow$ generic vector.
  !!  \item \emph{norm} $\rightarrow$ normal vector of face.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of
  !!                                  {\color{RoyalBlue}tipo} passed to VTK\_DAT.
  !! \item[{\color{RoyalBlue}varname}] contains the name attribuited the variable saved.
  !! \item[{\color{RoyalBlue}varX}] contains the values of $X$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}varY}] contains the values of $Y$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}varZ}] contains the values of $Z$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Note that the variables \texttt{varX,varY,varZ} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function.
  !!
  !!The following is an example of VTK\_VAR real vectorial data calling:
  !!
  !!\begin{boxred}{VTK\_VAR Real Vectorial Data Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC_NN=100
  !!real(4)::               varX(1:NC_NN)
  !!real(4)::               varZ(1:NC_NN)
  !!real(4)::               varZ(1:NC_NN)
  !!...
  !!E_IO = VTK_VAR('vect',NC_NN,'Real Vectorial Data',...
  !!            ...varX,varY,varZ)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_VAR INTEGER VECTORIAL DATA}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Integer Vectorial Data Signature}]
  !!function VTK_VAR(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!This kind of call is used to save integer vectorial data.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Integer Vectorial Data Variables}]
  !!integer(R4P),   intent(IN):: NC_NN         ! number of nodes or cells
  !!character(*),   intent(IN):: varname       ! variable name
  !!integer(R4P),   intent(IN):: varX(1:NC_NN) ! x component of vector
  !!integer(R4P),   intent(IN):: varY(1:NC_NN) ! y component of vector
  !!integer(R4P),   intent(IN):: varZ(1:NC_NN) ! z component of vector
  !!integer(R4P)::               E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_VAR variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of
  !!                                  {\color{RoyalBlue}tipo} passed to VTK\_DAT.
  !! \item[{\color{RoyalBlue}varname}] contains the name attribuited the variable saved.
  !! \item[{\color{RoyalBlue}varX}] contains the values of $X$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}varY}] contains the values of $Y$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}varZ}] contains the values of $Z$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_VAR real vectorial data calling:
  !!
  !!\begin{boxred}{VTK\_VAR Integer Vectorial Data Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC_NN=100
  !!integer(4)::            varX(1:NC_NN)
  !!integer(4)::            varZ(1:NC_NN)
  !!integer(4)::            varZ(1:NC_NN)
  !!...
  !!E_IO = VTK_VAR(NC_NN,'Integer Vectorial Data', &
  !!               varX,varY,varZ)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_VAR TEXTURE DATA}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Texture Data Signature}]
  !!function VTK_VAR(NC_NN,,dimm,varname,textCoo) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!This kind of call is used to save texture data.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR Texture Data Variables}]
  !!integer(R4P),     intent(IN):: NC_NN                   ! number of nodes or cells
  !!integer(R4P),     intent(IN):: dimm                    ! texture dimensions
  !!character(*),     intent(IN):: varname                 ! variable name
  !!real(R8P or R4P), intent(IN):: textCoo(1:NC_NN,1:dimm) ! texture
  !!integer(R4P)::                 E_IO                    ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_VAR variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of
  !!                                  {\color{RoyalBlue}tipo} passed to VTK\_DAT.
  !! \item[{\color{RoyalBlue}dimm}] indicates the dimensions of the texture coordinates. It can assume the value:
  !! \begin{enumerateABlu}
  !!  \item \emph{1} $\rightarrow$ scalar texture.
  !!  \item \emph{2} $\rightarrow$ twodimensional texture.
  !!  \item \emph{3} $\rightarrow$ threedimensional texture.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}varname}] contains the name attribuited the variable saved.
  !! \item[{\color{RoyalBlue}textCoo}] contains the coordinates of texture in each nodes or cells. It is a vector of
  !!                                   $[1:NC\_NN,1:dimm]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Note that the variable \texttt{textCoo} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function.
  !!
  !!The following is an example of VTK\_VAR texture data calling:
  !!
  !!\begin{boxred}{VTK\_VAR Texture Data Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC_NN=100
  !!integer(4), parameter:: dimm=2
  !!real(4)::               textCoo(1:NC_NN,1:dimm)
  !!...
  !!E_IO = VTK_VAR(NC_NN,dimm,'Texture Data',textCoo)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !(\doc)skippedblock
  function VTK_VAR_SCAL_R8(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving field of scalar variable (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of nodes or cells
  character(*), intent(IN):: varname      ! variable name
  real(R8P),    intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'SCALARS '//trim(varname)//' double 1'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
    write(unit=Unit_VTK,fmt=FR8P, iostat=E_IO)var
  case(f_out_binary)
    write(unit=Unit_VTK,iostat=E_IO)'SCALARS '//trim(varname)//' double 1'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'LOOKUP_TABLE default'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)var
    write(unit=Unit_VTK,iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_SCAL_R8

  function VTK_VAR_SCAL_R4(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving field of scalar variable (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of nodes or cells
  character(*), intent(IN):: varname      ! variable name
  real(R4P),    intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'SCALARS '//trim(varname)//' float 1'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
    write(unit=Unit_VTK,fmt=FR4P, iostat=E_IO)var
  case(f_out_binary)
    write(unit=Unit_VTK,iostat=E_IO)'SCALARS '//trim(varname)//' float 1'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'LOOKUP_TABLE default'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)var
    write(unit=Unit_VTK,iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_SCAL_R4

  function VTK_VAR_SCAL_I4(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving field of scalar variable (I4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of nodes or cells
  character(*), intent(IN):: varname      ! variable name
  integer(I4P), intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'SCALARS '//trim(varname)//' int 1'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
    write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)var
  case(f_out_binary)
    write(unit=Unit_VTK,iostat=E_IO)'SCALARS '//trim(varname)//' int 1'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'LOOKUP_TABLE default'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)var
    write(unit=Unit_VTK,iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_SCAL_I4

  function VTK_VAR_VECT_R8(vec_type,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving field of vectorial variable (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: vec_type      ! vector type: vect = generic vector , norm = normal vector
  integer(I4P), intent(IN):: NC_NN         ! number of nodes or cells
  character(*), intent(IN):: varname       ! variable name
  real(R8P),    intent(IN):: varX(1:NC_NN) ! x component of vector
  real(R8P),    intent(IN):: varY(1:NC_NN) ! y component of vector
  real(R8P),    intent(IN):: varZ(1:NC_NN) ! z component of vector
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  integer(I8P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    select case(Upper_Case(trim(vec_type)))
    case('VECT')
      write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'VECTORS '//trim(varname)//' double'
    case('NORM')
      write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'NORMALS '//trim(varname)//' double'
    endselect
    write(unit=Unit_VTK,fmt='(3'//FR8P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  case(f_out_binary)
    select case(Upper_Case(trim(vec_type)))
    case('VECT')
      write(unit=Unit_VTK,iostat=E_IO)'VECTORS '//trim(varname)//' double'//end_rec
    case('NORM')
      write(unit=Unit_VTK,iostat=E_IO)'NORMALS '//trim(varname)//' double'//end_rec
    endselect
    write(unit=Unit_VTK,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_VECT_R8

  function VTK_VAR_VECT_R4(vec_type,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving field of vectorial variable (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: vec_type      ! vector type: vect = generic vector , norm = normal vector
  integer(I4P), intent(IN):: NC_NN         ! number of nodes or cells
  character(*), intent(IN):: varname       ! variable name
  real(R4P),    intent(IN):: varX(1:NC_NN) ! x component of vector
  real(R4P),    intent(IN):: varY(1:NC_NN) ! y component of vector
  real(R4P),    intent(IN):: varZ(1:NC_NN) ! z component of vector
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  integer(I8P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    select case(Upper_Case(trim(vec_type)))
    case('vect')
      write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'VECTORS '//trim(varname)//' float'
    case('norm')
      write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'NORMALS '//trim(varname)//' float'
    endselect
    write(unit=Unit_VTK,fmt='(3'//FR4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  case(f_out_binary)
    select case(Upper_Case(trim(vec_type)))
    case('vect')
      write(unit=Unit_VTK,iostat=E_IO)'VECTORS '//trim(varname)//' float'//end_rec
    case('norm')
      write(unit=Unit_VTK,iostat=E_IO)'NORMALS '//trim(varname)//' float'//end_rec
    endselect
    write(unit=Unit_VTK,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_VECT_R4

  function VTK_VAR_VECT_I4(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving field of vectorial variable (I4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of nodes or cells
  character(*), intent(IN):: varname       ! variable name
  integer(I4P), intent(IN):: varX(1:NC_NN) ! x component of vector
  integer(I4P), intent(IN):: varY(1:NC_NN) ! y component of vector
  integer(I4P), intent(IN):: varZ(1:NC_NN) ! z component of vector
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  integer(I8P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)'VECTORS '//trim(varname)//' int'
    write(unit=Unit_VTK,fmt='(3'//FI4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  case(f_out_binary)
    write(unit=Unit_VTK,iostat=E_IO)'VECTORS '//trim(varname)//' int'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_VECT_I4

  function VTK_VAR_TEXT_R8(NC_NN,dimm,varname,textCoo) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving texture variable (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN                   ! number of nodes or cells
  integer(I4P), intent(IN):: dimm                    ! texture dimensions
  character(*), intent(IN):: varname                 ! variable name
  real(R8P),    intent(IN):: textCoo(1:NC_NN,1:dimm) ! texture
  integer(I4P)::             E_IO                    ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer                ! buffer string
  integer(I8P)::             n1,n2                   ! counters
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,1X,'//FI4P//'1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' double'
    write(s_buffer,     fmt='(I1)',                 iostat=E_IO)dimm
    s_buffer='('//trim(s_buffer)//FR4P//')'
    write(unit=Unit_VTK,fmt=trim(s_buffer),         iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,1X,'//FI4P//'1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' double'
    write(unit=Unit_VTK,                            iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                            iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
    write(unit=Unit_VTK,                            iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_TEXT_R8

  function VTK_VAR_TEXT_R4(NC_NN,dimm,varname,textCoo) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving texture variable (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN                   ! number of nodes or cells
  integer(I4P), intent(IN):: dimm                    ! texture dimensions
  character(*), intent(IN):: varname                 ! variable name
  real(R4P),    intent(IN):: textCoo(1:NC_NN,1:dimm) ! texture
  integer(I4P)::             E_IO                    ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer                ! buffer string
  integer(I8P)::             n1,n2                   ! counters
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,1X,'//FI4P//'1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' float'
    write(s_buffer,     fmt='(I1)',                 iostat=E_IO)dimm
    s_buffer='('//trim(s_buffer)//FR4P//')'
    write(unit=Unit_VTK,fmt=trim(s_buffer),         iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
  case(f_out_binary)
    write(s_buffer,     fmt='(A,1X,'//FI4P//'1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' float'
    write(unit=Unit_VTK,                            iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                            iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
    write(unit=Unit_VTK,                            iostat=E_IO)end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_TEXT_R4
  !(doc/)skippedblock

  function VTK_END() result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!This function is used to finalize the file opened and it has not inputs. The \LIBVTKIO manages the file unit without the
  !!user's action.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P):: E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!The VTK\_END variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_END calling:
  !!
  !!\begin{boxred}{VTK\_END Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_END()
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  close(unit=Unit_VTK,iostat=E_IO)
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_END

  !!\chapter{VTK XML functions}
  !!\minitoc
  !!\vspace*{8mm}
  !!
  !!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} XML standard is more powerful than legacy one. It is more flexible
  !!and free but on the other hand is more (but not so more using a library like \LIBVTKIO...) complex than legacy standard. The
  !!output of XML functions is a well-formated XML file at least for the ascii format (in the binary format \LIBVTKIO use
  !!raw-data format that does not produce a well formated XML file).
  !!
  !!The XML functions follow the same calling-convention of the legacy functions; all the \LIBVTKIO XML functions are
  !!\MaiuscolettoBS{4-byte integer function}: the output of these functions is an integer that is $0$ if the function calling
  !!has been done right while it is $> 0$  if some errors occur. The functions calling is the same as legacy functions:
  !!
  !!\begin{boxred}{Functions Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4):: E_IO
  !!...
  !!E_IO = VTK_INI_XML(....
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\noindent Note that the XML functions have the same name of legacy functions with the suffix \virgo{\_XML}.
  !!
  function VTK_INI_XML(output_format,filename,mesh_topology,nx1,nx2,ny1,ny2,nz1,nz2) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The VTK\_INI\_XML function is used for initializing file. This function must be the first to be called.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN)::           output_format ! output format: ASCII or BINARY
  character(*), intent(IN)::           filename      ! file name
  character(*), intent(IN)::           mesh_topology ! mesh topology
  integer(I4P), intent(IN), optional:: nx1,nx2       ! initial and final nodes of x axis
  integer(I4P), intent(IN), optional:: ny1,ny2       ! initial and final nodes of y axis
  integer(I4P), intent(IN), optional:: nz1,nz2       ! initial and final nodes of z axis
  integer(I4P)::                       E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::              s_buffer      ! buffer string
  !!The VTK\_INI\_XML variables have the following meaning:
  !!
  !!\begin{description}
  !!\item[{\color{RoyalBlue}output\_format}] indicates the \virgo{format} of output file. It can assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{ascii} (it is case insensitive) $\rightarrow$ creating an ascii output file.
  !!  \item \emph{binary} (it is case insensitive) $\rightarrow$ creating a binary (big\_endian encoding) output file.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}filename}] contains the name (with its path) of the output file.
  !! \item[{\color{RoyalBlue}topology}] indicates the topology of the mesh and can assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{StructuredGrid}.
  !!  \item \emph{RectilinearGrid}.
  !!  \item \emph{UnstructuredGrid}.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}nx1,nx2}] contains the extent of X axis; $nx1$ is the initial node and $nx2$ is the final.
  !! \item[{\color{RoyalBlue}ny1,ny2}] contains the extent of Y axis; $ny1$ is the initial node and $ny2$ is the final.
  !! \item[{\color{RoyalBlue}nz1,nz2}] contains the extent of Z axis; $nz1$ is the initial node and $nz2$ is the final.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!This function is quite more complex than the rispective legacy function; it needs more inputs: the XML standard needs more
  !!informations to initialize the file.
  !!
  !!The following is an example of VTK\_INI\_XML calling:
  !!
  !!\begin{boxred}{VTK\_INI\_XML Calling}
  !!\begin{verbatim}
  !!...
  !!...
  !!E_IO = VTK_INI_XML('BINARY','XML_RECT_BINARY.vtr', &
  !!                   'RectilinearGrid',              &
  !!                   nx1=nx1,nx2=nx2,                &
  !!                   ny1=ny1,ny2=ny2,                &
  !!                   nz1=nz1,nz2=nz2)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\noindent Note that the file extension is necessary in the file name. The XML standard has different extensions for each
  !!different topologies (i.e. \MaiuscolettoBS{.vtr} for rectilinear topology). See the VTK-standard file for more information.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  topology = trim(mesh_topology)
  Unit_VTK=GetUnit()
  select case(trim(Upper_Case(output_format)))
  case('ASCII')
    f_out = f_out_ascii
    open(unit   = Unit_VTK,       &
         file   = trim(filename), &
         form   = 'FORMATTED',    &
         access = 'SEQUENTIAL',   &
!        action = 'WRITE',        &
!        buffered   = 'YES',      &
         iostat = E_IO)
    ! writing header of file
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'<?xml version="1.0"?>'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'<VTKFile type="'//trim(topology)//'" version="0.1" byte_order="BigEndian">'
    indent = 2
    select case(trim(topology))
    case('RectilinearGrid','StructuredGrid')
      write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent) &
              //'<'//trim(topology)//' WholeExtent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    case('UnstructuredGrid')
      write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//'>'
    endselect
    indent = indent + 2
  case('BINARY')
    f_out = f_out_binary
    open(unit       = Unit_VTK,       &
         file       = trim(filename), &
         form       = 'UNFORMATTED',  &
         access     = 'SEQUENTIAL',   &
         action     = 'WRITE',        &
!        convert    = 'BIG_ENDIAN',   &
!        recordtype = 'STREAM',       &
!        buffered   = 'YES',          &
         iostat     = E_IO)
    ! writing header of file
    write(unit=Unit_VTK,iostat=E_IO)'<?xml version="1.0"?>'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'<VTKFile type="'//trim(topology)//'" version="0.1" byte_order="BigEndian">'//end_rec
    indent = 2
    select case(trim(topology))
    case('RectilinearGrid','StructuredGrid')
      write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//&
                 ' WholeExtent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    case('UnstructuredGrid')
      write(s_buffer,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<'//trim(topology)//'>'
    endselect
    write(unit=Unit_VTK,iostat=E_IO)trim(s_buffer)//end_rec
    indent = indent + 2
    Unit_VTK_Append=GetUnit()
    ! opening the SCRATCH file used for appending raw binary data
    open(unit       = Unit_VTK_Append, &
         form       = 'UNFORMATTED',   &
         access     = 'SEQUENTIAL',    &
         action     = 'WRITE',         &
!        convert    = 'BIG_ENDIAN',    &
!        recordtype = 'STREAM',        &
!        buffered   = 'YES',           &
         status     = 'SCRATCH',       &
         iostat     = E_IO)
    ioffset = 0 ! initializing offset puntator
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_INI_XML

  !!\section{VTK\_GEO\_XML}
  !!
  !!VTK\_GEO\_XML is an interface to 6 different functions; there are 2 functions for each 3 topologies supported.
  !!This function must be called after VTK\_INI\_XML. It saves the mesh geometry. The inputs that must be passed change
  !!depending on the topologies choiced. Not all VTK topologies have been implemented (\virgo{polydata} topologies are absent).
  !!The signatures for all implemented topologies are now reported.
  !!
  !!\subsection{VTK\_GEO\_XML STRUCTURED GRID}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Structured Grid Signature}]
  !!function VTK_GEO_XML(nx1,nx2,ny1,ny2,nz1,nz2,NN, &
  !!                     X,Y,Z) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{structured grid} is useful for structured grid with non-uniform discretization steps.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Structured Grid Variables}]
  !!integer(I4P),     intent(IN):: nx1,nx2  ! initial and final nodes of x axis
  !!integer(I4P),     intent(IN):: ny1,ny2  ! initial and final nodes of y axis
  !!integer(I4P),     intent(IN):: nz1,nz2  ! initial and final nodes of z axis
  !!integer(I4P),     intent(IN):: NN       ! number of all nodes
  !!real(R8P or R4P), intent(IN):: X(1:NN)  ! x coordinates
  !!real(R8P or R4P), intent(IN):: Y(1:NN)  ! y coordinates
  !!real(R8P or R4P), intent(IN):: Z(1:NN)  ! z coordinates
  !!integer(I4P)::                 E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X,Y,Z} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic displacement
  !!interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be sure that all
  !!variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO\_XML structured grid variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}nx1,nx2}] contains the extent of X axis; $nx1$ is the initial node and $nx2$ is the final.
  !! \item[{\color{RoyalBlue}ny1,ny2}] contains the extent of Y axis; $ny1$ is the initial node and $ny2$ is the final.
  !! \item[{\color{RoyalBlue}nz1,nz2}] contains the extent of Z axis; $nz1$ is the initial node and $nz2$ is the final.
  !! \item[{\color{RoyalBlue}NN}] contains the global number of nodes $NN=(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)$.
  !! \item[{\color{RoyalBlue}X}] contains the $X$ coordinates values of all nodes. It is a vector of $[1:NN]$.
  !! \item[{\color{RoyalBlue}Y}] contains the $Y$ coordinates values of all nodes. It is a vector of $[1:NN]$.
  !! \item[{\color{RoyalBlue}Z}] contains the $Z$ coordinates values of all nodes. It is a vector of $[1:NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO\_XML structured grid calling:
  !!
  !!\begin{boxred}{VTK\_GEO\_XML Structured Grid Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4):: nx1,nx2
  !!integer(4):: ny1,ny2
  !!integer(4):: nz1,nz2
  !!integer(4):: NN
  !!real(4):: X(1:NN),Y(1:NN),Z(1:NN)
  !!...
  !!E_IO = VTK_GEO_XML(nx1,nx2,ny1,ny2,nz1,nz2, &
  !!                   NN,                      &
  !!                   X,Y,Z)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_GEO\_XML RECTILINEAR GRID}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Rectilinear Grid Signature}]
  !!function VTK_GEO_XML(nx1,nx2,ny1,ny2,nz1,nz2, &
  !!                     X,Y,Z) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{rectilinear grid} is useful for structured grid with non-uniform discretization steps even in
  !!generalized coordinates.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Rectilinear Grid Variables}]
  !!integer(I4P),     intent(IN):: nx1,nx2    ! initial and final nodes of x axis
  !!integer(I4P),     intent(IN):: ny1,ny2    ! initial and final nodes of y axis
  !!integer(I4P),     intent(IN):: nz1,nz2    ! initial and final nodes of z axis
  !!real(R8P or R4P), intent(IN):: X(nx1:nx2) ! x coordinates
  !!real(R8P or R4P), intent(IN):: Y(ny1:ny2) ! y coordinates
  !!real(R8P or R4P), intent(IN):: Z(nz1:nz2) ! z coordinates
  !!integer(I4P)::                 E_IO       ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X,Y,Z} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be sure
  !!that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO\_XML rectilinear grid variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}nx1,nx2}] contains the extent of X axis; $nx1$ is the initial node and $nx2$ is the final.
  !! \item[{\color{RoyalBlue}ny1,ny2}] contains the extent of Y axis; $ny1$ is the initial node and $ny2$ is the final.
  !! \item[{\color{RoyalBlue}nz1,nz2}] contains the extent of Z axis; $nz1$ is the initial node and $nz2$ is the final.
  !! \item[{\color{RoyalBlue}X}] contains the $X$ coordinates values of X nodes. It is a vector of $[nx1:nx2]$.
  !! \item[{\color{RoyalBlue}Y}] contains the $Y$ coordinates values of Y nodes. It is a vector of $[ny1:ny2]$.
  !! \item[{\color{RoyalBlue}Z}] contains the $Z$ coordinates values of Z nodes. It is a vector of $[nz1:nz2]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO\_XML rectilinear grid calling:
  !!
  !!\begin{boxred}{VTK\_GEO\_XML Structured Grid Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4):: nx1,nx2
  !!integer(4):: ny1,ny2
  !!integer(4):: nz1,nz2
  !!real(4):: X(nx1:nx2),Y(ny1:ny2),Z(nz1:nz2)
  !!...
  !!E_IO = VTK_GEO_XML(nx1,nx2,ny1,ny2,nz1,nz2, &
  !!                   X,Y,Z)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_GEO\_XML UNSTRUCTURED GRID}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Unstructured Grid Signature}]
  !!function VTK_GEO_XML(Nnodi,NCelle,X,Y,Z) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{unstructured grid} is necessary for unstructured grid, the most general mesh format. This topology
  !!is also useful for scructured mesh in order to save only a non-structured clip of mesh.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Unstructured Grid Variables}]
  !!integer(I4P),     intent(IN):: NN       ! number of nodes
  !!integer(I4P),     intent(IN):: NC       ! number of cells
  !!real(R8P or R4P), intent(IN):: X(1:NN)  ! x coordinates
  !!real(R8P or R4P), intent(IN):: Y(1:NN)  ! y coordinates
  !!real(R8P or R4P), intent(IN):: Z(1:NN)  ! z coordinates
  !!integer(I4P)::                 E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X,Y,Z} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be
  !!sure that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO\_XML unstructured grid variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}Nnodi}] indicates the number of all nodes.
  !! \item[{\color{RoyalBlue}NCelle}] indicates the number of all cells.
  !! \item[{\color{RoyalBlue}X}] contains the $X$ coordinates values of nodes. It is a vector of $[1:Nnodi]$.
  !! \item[{\color{RoyalBlue}Y}] contains the $Y$ coordinates values of nodes. It is a vector of $[1:Nnodi]$.
  !! \item[{\color{RoyalBlue}Z}] contains the $Z$ coordinates values of nodes. It is a vector of $[1:Nnodi]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO\_XML unstructured grid calling:
  !!
  !!\begin{boxred}{VTK\_GEO\_XML Unstructured Grid Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: Nnodi=100
  !!integer(4), parameter:: NCelle=50
  !!real(4):: X(1:Nnodi),Y(1:Nnodi),Z(1:Nnodi)
  !!...
  !!E_IO = VTK_GEO_XML('ascii',Nnodi,NCelle,X,Y,Z)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!In order to use the \virgo{unstructured grid} it is necessary to save also the \virgo{connectivity} of the grid.
  !!The connectivity must be saved with the function \MaiuscolettoBS{VTK\_CON\_XML}.
  !!
  !!\subsection{VTK\_GEO\_XML CLOSE PIECE}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Close Piece Signature}]
  !!function VTK_GEO_XML() result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!As we said before the XML standard is more powerful than legacy. XML file can contain more than 1 mesh with its
  !!associated variables. Thus there is the necessity to close each \virgo{pieces} that compose the data-set saved in the
  !!XML file. The \MaiuscolettoBS{VTK\_GEO\_XML} called in the \virgo{close piece} format is used just to close the
  !!current piece before saving another piece or closing the file.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO\_XML Close Piece Variables}]
  !!integer(I4P):: E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_GEO\_XML close piece variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_GEO\_XML close piece calling:
  !!
  !!\begin{boxred}{VTK\_GEO\_XML Unstructured Grid Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_GEO_XML()
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !(\doc)skippedblock
  function VTK_GEO_XML_STRG_R8(nx1,nx2,ny1,ny2,nz1,nz2,NN,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = StructuredGrid (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: nx1,nx2  ! initial and final nodes of x axis
  integer(I4P), intent(IN):: ny1,ny2  ! initial and final nodes of y axis
  integer(I4P), intent(IN):: nz1,nz2  ! initial and final nodes of z axis
  integer(I4P), intent(IN):: NN       ! number of all nodes
  real(R8P),    intent(IN):: X(1:NN)  ! x coordinates
  real(R8P),    intent(IN):: Y(1:NN)  ! y coordinates
  real(R8P),    intent(IN):: Z(1:NN)  ! z coordinates
  integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer ! buffer string
  integer(I4P)::             n1       ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<Points>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)// &
                          '<DataArray type="Float64" NumberOfComponents="3" Name="Point" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FR8P//')',iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</Points>'
  case(f_out_binary)
    write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'<Points>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)',iostat=E_IO)ioffset
    write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" NumberOfComponents="3" Name="Point"i &
   &                        format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NN*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R8',3*NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</Points>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_STRG_R8

  function VTK_GEO_XML_STRG_R4(nx1,nx2,ny1,ny2,nz1,nz2,NN,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = StructuredGrid (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: nx1,nx2  ! initial and final nodes of x axis
  integer(I4P), intent(IN):: ny1,ny2  ! initial and final nodes of y axis
  integer(I4P), intent(IN):: nz1,nz2  ! initial and final nodes of z axis
  integer(I4P), intent(IN):: NN       ! number of all nodes
  real(R4P),    intent(IN):: X(1:NN)  ! x coordinates
  real(R4P),    intent(IN):: Y(1:NN)  ! y coordinates
  real(R4P),    intent(IN):: Z(1:NN)  ! z coordinates
  integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer ! buffer string
  integer(I4P)::             n1       ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<Points>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)// &
                              '<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FR4P//')',    iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</Points>'
  case(f_out_binary)
    write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,                   iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<Points>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)// &
      '<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NN*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',3*NN
    write(unit=Unit_VTK_Append,            iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</Points>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_STRG_R4

  function VTK_GEO_XML_RECT_R8(nx1,nx2,ny1,ny2,nz1,nz2,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = RectilinearGrid (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: nx1,nx2    ! initial and final nodes of x axis
  integer(I4P), intent(IN):: ny1,ny2    ! initial and final nodes of y axis
  integer(I4P), intent(IN):: nz1,nz2    ! initial and final nodes of z axis
  real(R8P),    intent(IN):: X(nx1:nx2) ! x coordinates
  real(R8P),    intent(IN):: Y(ny1:ny2) ! y coordinates
  real(R8P),    intent(IN):: Z(nz1:nz2) ! z coordinates
  integer(I4P)::             E_IO       ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer   ! buffer string
  integer(I4P)::             n1         ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<Coordinates>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" Name="X" format="ascii">'
    write(unit=Unit_VTK,fmt=FR8P,               iostat=E_IO)(X(n1),n1=nx1,nx2)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" Name="Y" format="ascii">'
    write(unit=Unit_VTK,fmt=FR8P,               iostat=E_IO)(Y(n1),n1=ny1,ny2)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" Name="Z" format="ascii">'
    write(unit=Unit_VTK,fmt=FR8P,               iostat=E_IO)(Z(n1),n1=nz1,nz2)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</Coordinates>'
  case(f_out_binary)
    write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,                   iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<Coordinates>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)// &
                  '<DataArray type="Float64" Name="X" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = (nx2-nx1+1)*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R8',nx2-nx1+1
    write(unit=Unit_VTK_Append,            iostat=E_IO)(X(n1),n1=nx1,nx2)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)// &
                      '<DataArray type="Float64" Name="Y" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = (ny2-ny1+1)*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R8',ny2-ny1+1
    write(unit=Unit_VTK_Append,            iostat=E_IO)(Y(n1),n1=ny1,ny2)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)// &
                     '<DataArray type="Float64" Name="Z" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = (nz2-nz1+1)*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R8',nz2-nz1+1
    write(unit=Unit_VTK_Append,            iostat=E_IO)(Z(n1),n1=nz1,nz2)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</Coordinates>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_RECT_R8

  function VTK_GEO_XML_RECT_R4(nx1,nx2,ny1,ny2,nz1,nz2,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = RectilinearGrid (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: nx1,nx2    ! initial and final nodes of x axis
  integer(I4P), intent(IN):: ny1,ny2    ! initial and final nodes of y axis
  integer(I4P), intent(IN):: nz1,nz2    ! initial and final nodes of z axis
  real(R4P),    intent(IN):: X(nx1:nx2) ! x coordinates
  real(R4P),    intent(IN):: Y(ny1:ny2) ! y coordinates
  real(R4P),    intent(IN):: Z(nz1:nz2) ! z coordinates
  integer(I4P)::             E_IO       ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer   ! buffer string
  integer(I4P)::             n1         ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<Coordinates>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="X" format="ascii">'
    write(unit=Unit_VTK,fmt=FR4P,               iostat=E_IO)(X(n1),n1=nx1,nx2)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="Y" format="ascii">'
    write(unit=Unit_VTK,fmt=FR4P,               iostat=E_IO)(Y(n1),n1=ny1,ny2)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="Z" format="ascii">'
    write(unit=Unit_VTK,fmt=FR4P,               iostat=E_IO)(Z(n1),n1=nz1,nz2)
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',              iostat=E_IO)repeat(' ',indent)//'</Coordinates>'
  case(f_out_binary)
    write(s_buffer,fmt='(A,6'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//'<Piece Extent="',nx1,nx2,ny1,ny2,nz1,nz2,'">'
    indent = indent + 2
    write(unit=Unit_VTK,                   iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'<Coordinates>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//&
                '<DataArray type="Float32" Name="X" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = (nx2-nx1+1)*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',nx2-nx1+1
    write(unit=Unit_VTK_Append,            iostat=E_IO)(X(n1),n1=nx1,nx2)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)// &
                      '<DataArray type="Float32" Name="Y" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = (ny2-ny1+1)*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',ny2-ny1+1
    write(unit=Unit_VTK_Append,            iostat=E_IO)(Y(n1),n1=ny1,ny2)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    write(s_buffer,fmt='(I8)',             iostat=E_IO)ioffset
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)// &
                     '<DataArray type="Float32" Name="Z" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = (nz2-nz1+1)*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,            iostat=E_IO)N_Byte,'R4',nz2-nz1+1
    write(unit=Unit_VTK_Append,            iostat=E_IO)(Z(n1),n1=nz1,nz2)
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,                   iostat=E_IO)repeat(' ',indent)//'</Coordinates>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_RECT_R4

  function VTK_GEO_XML_UNST_R8(NN,NC,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = UnstructuredGrid (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NN       ! number of nodes
  integer(I4P), intent(IN):: NC       ! number of cells
  real(R8P),    intent(IN):: X(1:NN)  ! x coordinates
  real(R8P),    intent(IN):: Y(1:NN)  ! y coordinates
  real(R8P),    intent(IN):: Z(1:NN)  ! z coordinates
  integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer ! buffer string
  integer(I4P)::             n1       ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A,'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)// &
                              '<Piece NumberOfPoints="',NN,'" NumberOfCells="',NC,'">'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'<Points>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//& 
                   '<DataArray type="Float64" NumberOfComponents="3" Name="Point" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FR8P//')',                iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'</Points>'
  case(f_out_binary)
    write(s_buffer,fmt='(A,'//FI4P//',A,'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//&
                    '<Piece NumberOfPoints="',NN,'" NumberOfCells="',NC,'">'
    indent = indent + 2
    write(unit=Unit_VTK,                               iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'<Points>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)',                         iostat=E_IO)ioffset
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//&
             '<DataArray type="Float64" NumberOfComponents="3" Name="Point" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NN*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,                        iostat=E_IO)N_Byte,'R8',3*NN
    write(unit=Unit_VTK_Append,                        iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'</Points>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_UNST_R8

  function VTK_GEO_XML_UNST_R4(NN,NC,X,Y,Z) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving mesh; topology = UnstructuredGrid (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NN       ! number of nodes
  integer(I4P), intent(IN):: NC       ! number of cells
  real(R4P),    intent(IN):: X(1:NN)  ! x coordinates
  real(R4P),    intent(IN):: Y(1:NN)  ! y coordinates
  real(R4P),    intent(IN):: Z(1:NN)  ! z coordinates
  integer(I4P)::             E_IO     ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer ! buffer string
  integer(I4P)::             n1       ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A,'//FI4P//',A,'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)// &
                    '<Piece NumberOfPoints="',NN,'" NumberOfCells="',NC,'">'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'<Points>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//&
                     '<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FR4P//')',                iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',                          iostat=E_IO)repeat(' ',indent)//'</Points>'
  case(f_out_binary)
    write(s_buffer,fmt='(A,'//FI4P//',A,'//FI4P//',A)',iostat=E_IO)repeat(' ',indent)//&
                     '<Piece NumberOfPoints="',NN,'" NumberOfCells="',NC,'">'
    indent = indent + 2
    write(unit=Unit_VTK,                               iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'<Points>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)',                         iostat=E_IO)ioffset
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//& 
        '<DataArray type="Float32" NumberOfComponents="3" Name="Point" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NN*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,                        iostat=E_IO)N_Byte,'R4',3*NN
    write(unit=Unit_VTK_Append,                        iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,                               iostat=E_IO)repeat(' ',indent)//'</Points>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_UNST_R4

  function VTK_GEO_XML_CLOSEP() result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for closing mesh block data.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P):: E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  indent = indent - 2
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</Piece>'
  case(f_out_binary)
    write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</Piece>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_CLOSEP
  !(doc/)skippedblock

  function VTK_CON_XML(NC,connect,offset,cell_type) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!This function \MaiuscolettoBS{must} be used when unstructured grid is used. It saves the connectivity of the unstructured mesh.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC              ! number of cells
  integer(I4P), intent(IN):: connect(:)      ! mesh connectivity
  integer(I4P), intent(IN):: offset(1:NC)    ! cell offset
  integer(I1P), intent(IN):: cell_type(1:NC) ! VTK cell type
  integer(I4P)::             E_IO            ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer        ! buffer string
  integer(I4P)::             n1              ! counter
  !!The VTK\_CON\_XML variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NCelle}] indicates the number of all cells.
  !! \item[{\color{RoyalBlue}connect}] contains the connectivity of the mesh. It is a vector.
  !! \item[{\color{RoyalBlue}offset}] contains the offset\footnote{The summ of nodes of all previous cells included the
  !!                                  current cell.} of every cells. It is a vector of $[1:NCelle]$.
  !! \item[{\color{RoyalBlue}tipo}] contains the type of every cells. It is a vector of $[1:NCelle]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The vector \MaiuscolettoBS{connect} must follow the VTK XML standard. It is passed as \MaiuscolettoBS{assumed-shape}
  !!array because its dimensions is related to the mesh dimensions in a complex way. Its dimensions can be calculated by
  !!the following equation:
  !!
  !!\begin{equation}
  !!dc = \sum\limits_{i = 1}^{NCelle} {nvertex_i }
  !!\label{eq:xml connectivity dimensions}
  !!\end{equation}
  !!
  !!\noindent where $dc$ is connectivity vector dimension and $nvertex_i$ is the number of vertices of $i^{th}$ cell.
  !!Note that this equation is different from the legacy one (eq. \ref{eq:connectivity dimensions}). The XML connectivity
  !!convention is quite different from the legacy standard. As an example considering the same mesh of section \ref{sec:VTKCON}:
  !!suppose we have a mesh composed by 2 cells, one hexahedron (8 vertices) and one pyramid with square basis (5 vertices);
  !!suppose that the basis of pyramid is constitute by a face of the hexahedron and so the two cells share 4 vertices. The
  !!equation \ref{eq:xml connectivity dimensions} gives $dc=8+5=13$; the connectivity vector for this mesh can be:
  !!
  !!\begin{boxred}{Connectivity vector example for VTK XML standard}
  !!\begin{verbatim}
  !!! first cell
  !!connect(1)  = 0  => identification flag of 1° vertex of 1° cell
  !!connect(2)  = 1  => identification flag of 2° vertex of 1° cell
  !!connect(3)  = 2  => identification flag of 3° vertex of 1° cell
  !!connect(4)  = 3  => identification flag of 4° vertex of 1° cell
  !!connect(5)  = 4  => identification flag of 5° vertex of 1° cell
  !!connect(6)  = 5  => identification flag of 6° vertex of 1° cell
  !!connect(7)  = 6  => identification flag of 7° vertex of 1° cell
  !!connect(8)  = 7  => identification flag of 8° vertex of 1° cell
  !!! second cell
  !!connect(9)  = 0  => identification flag of 1° vertex of 2° cell
  !!connect(10) = 1  => identification flag of 2° vertex of 2° cell
  !!connect(11) = 2  => identification flag of 3° vertex of 2° cell
  !!connect(12) = 3  => identification flag of 4° vertex of 2° cell
  !!connect(13) = 8  => identification flag of 5° vertex of 2° cell
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!Therefore this connectivity vector convention is more simple than the legacy convention, now we must create also the
  !!\MaiuscolettoBS{offset} vector that contains the data now missing in the \MaiuscolettoBS{connect} vector. The offset
  !!vector for this mesh can be:
  !!
  !!\begin{boxred}{Offset vector example for VTK XML standard}
  !!\begin{verbatim}
  !!! first cell
  !!offset(1) = 8  => summ of nodes of 1° cell
  !!! second cell
  !!offset(2) = 13 => summ of nodes of 1° and 2° cells
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\noindent The value of every cell-offset can be calculated by the following equation:
  !!
  !!\begin{equation}
  !!offset_c = \sum\limits_{i = 1}^{c} {nvertex_i }
  !!\label{eq:xml offset vlue}
  !!\end{equation}
  !!
  !!\noindent where $offset_c$ is the value of $c^{th}$ cell and $nvertex_i$ is the number of vertices of $i^{th}$ cell.
  !!
  !!The function VTK\_CON\_XML does not calculate the connectivity and offset vectors: it writes the connectivity and offset
  !!vectors conforming the VTK XML standard, but does not calculate them. In the future release of \LIBVTKIO will be included
  !!a function to calculate the connectivity and offset vector.
  !!
  !!The vector variable \MaiuscolettoBS{tipo} must conform the VTK XML standard \footnote{See the file VTK-Standard at the
  !!Kitware homepage.} that is the same of the legacy standard presented previous (sec. \ref{sec:VTKCON}). It contains the
  !!\emph{type} of each cells. For the above example this vector is:
  !!
  !!\begin{boxred}{Cell-Type vector example for VTK legacy standard}
  !!\begin{verbatim}
  !!tipo(1) = 12  => VTK hexahedron type of 1° cell
  !!tipo(2) = 14  => VTK pyramid type of 2° cell
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!The following is an example of VTK\_CON\_XML calling:
  !!
  !!\begin{boxred}{VTK\_CON\_XML Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NCelle=2
  !!integer(4), parameter:: Nvertex1=8
  !!integer(4), parameter:: Nvertex2=5
  !!integer(4), parameter:: dc=Nvertex1+Nvertex2
  !!integer(4)::            connect(1:dc)
  !!integer(4)::            offset(1:NCelle)
  !!integer(4)::            tipo(1:NCelle)
  !!...
  !!E_IO = VTK_CON_XML(NCelle,connect,offset,tipo)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<Cells>'
    indent = indent + 2
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="connectivity" format="ascii">'
    write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)(connect(n1),n1=1,size(connect))
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="offsets" format="ascii">'
    write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)(offset(n1),n1=1,NC)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="types" format="ascii">'
    write(unit=Unit_VTK,fmt=FI1P, iostat=E_IO)(cell_type(n1),n1=1,NC)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)', iostat=E_IO)repeat(' ',indent)//'</Cells>'
  case(f_out_binary)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<Cells>'//end_rec
    indent = indent + 2
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)// &
             '<DataArray type="Int32" Name="connectivity" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = size(connect)*sizeof(Tipo_I4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',size(connect)
    write(unit=Unit_VTK_Append,iostat=E_IO)(connect(n1),n1=1,size(connect))
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)// &
                '<DataArray type="Int32" Name="offsets" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC*sizeof(Tipo_I4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',NC
    write(unit=Unit_VTK_Append,iostat=E_IO)(offset(n1),n1=1,NC)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//&
                 '<DataArray type="Int8" Name="types" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC*sizeof(Tipo_I1)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I1',NC
    write(unit=Unit_VTK_Append,iostat=E_IO)(cell_type(n1),n1=1,NC)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
    indent = indent - 2
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</Cells>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_CON_XML

  function VTK_DAT_XML(var_location,var_block_action) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!This function \MaiuscolettoBS{must} be called before saving the data related to geometric mesh. This function initializes
  !!the saving of data variables indicating the \emph{type} of variables that will be saved.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: var_location     ! location of saving variables: CELL for cell-centered, NODE for node-centered
  character(*), intent(IN):: var_block_action ! variables block action: OPEN or CLOSE block
  integer(I4P)::             E_IO             ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!The VTK\_DAT\_XML variables have the following meaning:
  !!
  !!\begin{description}
  !!\item[{\color{RoyalBlue}var\_location}] contains the location-type of variables that will be saved after VTK\_DAT.
  !!It is a scalar and cab assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{cell} (it is case insensitive) $\rightarrow$ variables will be cell-centered.
  !!  \item \emph{node} (it is case insensitive) $\rightarrow$ variables will be node-centered.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}var\_block\_action}] indicates if the block-data-variables is being opened or closed; it can
  !!                                              assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{open}  (it is case insensitive) $\rightarrow$ block-data is being opened.
  !!  \item \emph{close} (it is case insensitive) $\rightarrow$ block-data is being closed.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Of course a single file can contain both cell and node centered variables. The \MaiuscolettoBS{VTK\_DAT\_XML} must be
  !!called two times, before saving a block-data-variables in order to open the block, and after the block-data-variables
  !!has been saved in order to close the block. XML file can contains as many blocks as you want.
  !!
  !!The following is an example of VTK\_DAT\_XML calling:
  !!
  !!\begin{boxred}{VTK\_DAT\_XML Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_DAT_XML('node','OPEN')
  !!...
  !!SAVE YOUR DATA WITH VTK_VAR_XML
  !!...
  !!E_IO = VTK_DAT_XML('node','CLOSE')
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<CellData>'
        indent = indent + 2
      case('CLOSE')
        indent = indent - 2
        write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</CellData>'
      endselect
    case('NODE')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<PointData>'
        indent = indent + 2
      case('CLOSE')
        indent = indent - 2
        write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</PointData>'
      endselect
    endselect
  case(f_out_binary)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'<CellData>'//end_rec
        indent = indent + 2
      case('CLOSE')
        indent = indent - 2
        write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</CellData>'//end_rec
      endselect
    case('NODE')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'<PointData>'//end_rec
        indent = indent + 2
      case('CLOSE')
        indent = indent - 2
        write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</PointData>'//end_rec
      endselect
    endselect
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_DAT_XML

  !!\section{VTK\_VAR\_XML}
  !!
  !!VTK\_VAR\_XML is an interface to 12 different functions; there are 6 functions for scalar variables (1 for each supported
  !!precision: R8P, R4P, I8P, I4P, I2P and I1P) and 6 for vectorial variables (1 for each supported precision: R8P, R4P, I8P,
  !!I4P, I2P and I1P)
  !!This function saves the data variables related to geometric mesh. The inputs that must be passed change depending on the
  !!data variables type.
  !!
  !!\subsection{VTK\_VAR\_XML SCALAR DATA}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR\_XML Scalar Data Signature}]
  !!function VTK_VAR_XML(NC_NN,varname,var) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!This kind of call is used to save scalar data.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR\_XML Scalar Data Variables}]
  !!integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  !!character(*), intent(IN):: varname      ! variable name
  !!real(R8P or...
  !!     R4P) or...
  !!integer(I8P or...
  !!        I4P or...
  !!        I2P or...
  !!        I1P), intent(IN):: var(1:NC_NN) ! variable to be saved
  !!integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_VAR\_XML variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of
  !!                                  {\color{RoyalBlue}var\_location} passed to VTK\_DAT\_XML.
  !! \item[{\color{RoyalBlue}varname}] contains the name attribuited the variable saved.
  !! \item[{\color{RoyalBlue}var}] contains the values of variables in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Note that the variables \texttt{var} can be passed both 8-byte real kind, 4-byte real kind, 8-byte integer, 4-byte integer, 2-byte integer and 1-byte integer; XML is very flexible; the dynamic displacement interface will call the correct function.
  !!
  !!The following is an example of VTK\_VAR\_XML scalar data calling:
  !!
  !!\begin{boxred}{VTK\_VAR\_XML Scalar Data Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC_NN=100
  !!integer(2)::            var(1:NC_NN)
  !!...
  !!E_IO = VTK_VAR_XML(NC_NN,'Scalar Data',var)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !!\subsection{VTK\_VAR\_XML VECTORIAL DATA}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR\_XML Vectorial Data Signature}]
  !!function VTK_VAR_XML(NC_NN,varname, &
  !!                     varX,varY,varZ) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!This kind of call is used to save vectorial data.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_VAR\_XML Vectorial Data Variables}]
  !!integer(I4P),                      intent(IN):: NC_NN         ! number of cells or nodes
  !!character(*),                      intent(IN):: varname       ! variable name
  !!real(R8P or R4P) or...
  !!integer(I8P or I4P or I2P or I1P), intent(IN):: varX(1:NC_NN) ! x component
  !!real(R8P or R4P) or...
  !!integer(I8P or I4P or I2P or I1P), intent(IN):: varY(1:NC_NN) ! y component
  !!real(R8P or R4P) or...
  !!integer(I8P or I4P or I2P or I1P), intent(IN):: varZ(1:NC_NN) ! z component
  !!integer(I4P)::                                  E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The VTK\_VAR\_XML variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}NC\_NN}] indicates the number of all cells or all nodes according to the value of
  !!                                  {\color{RoyalBlue}var\_location} passed to VTK\_DAT\_XML.
  !! \item[{\color{RoyalBlue}varname}] contains the name attribuited the variable saved.
  !! \item[{\color{RoyalBlue}varX}] contains the values of $X$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}varY}] contains the values of $Y$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}varZ}] contains the values of $Z$ component in each nodes or cells. It is a vector of $[1:NC\_NN]$.
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!Note that the variables \texttt{varX,varY,varZ} can be passed both 8-byte real kind, 4-byte real kind, 8-byte integer,
  !!4-byte integer, 2-byte integer and 1-byte integer; XML is very flexible; the dynamic displacement interface will call
  !!the correct function.
  !!
  !!The following is an example of VTK\_VAR\_XML vectorial data calling:
  !!
  !!\begin{boxred}{VTK\_VAR\_XML Vectorial Data Calling}
  !!\begin{verbatim}
  !!...
  !!integer(4), parameter:: NC_NN=100
  !!integer(4)::            varX(1:NC_NN)
  !!integer(4)::            varZ(1:NC_NN)
  !!integer(4)::            varZ(1:NC_NN)
  !!...
  !!E_IO = VTK_VAR_XML(NC_NN,'Vectorial Data', &
  !!                   varX,varY,varZ)
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!
  !(\doc)skippedblock
  function VTK_VAR_XML_SCAL_R8(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving scalar variable (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  character(*), intent(IN):: varname      ! variable name
  real(R8P),    intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  integer(I4P)::             n1           ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//&
       '<DataArray type="Float64" Name="'//trim(varname)//'" NumberOfComponents="1" format="ascii">'
    write(unit=Unit_VTK,fmt=FR8P, iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" Name="'//trim(varname)&
                //'" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC_NN*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R8',NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_R8

  function VTK_VAR_XML_SCAL_R4(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving scalar variable (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  character(*), intent(IN):: varname      ! variable name
  real(R4P),    intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  integer(I4P)::             n1           ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//&
                 '" NumberOfComponents="1" format="ascii">'
    write(unit=Unit_VTK,fmt=FR4P, iostat=E_IO)var
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//&
              '" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC_NN*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R4',NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_R4

  function VTK_VAR_XML_SCAL_I8(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving scalar variable (I8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  character(*), intent(IN):: varname      ! variable name
  integer(I8P), intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  integer(I4P)::             n1           ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int64" Name="'&
                //trim(varname)//'" NumberOfComponents="1" format="ascii">'
    write(unit=Unit_VTK,fmt=FI8P, iostat=E_IO)var
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int64" Name="'//trim(varname)//&
             '" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC_NN*sizeof(Tipo_I8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I8',NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I8

  function VTK_VAR_XML_SCAL_I4(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving scalar variable (I4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  character(*), intent(IN):: varname      ! variable name
  integer(I4P), intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  integer(I4P)::             n1           ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'&
              //trim(varname)//'" NumberOfComponents="1" format="ascii">'
    write(unit=Unit_VTK,fmt=FI4P, iostat=E_IO)var
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'&
           //trim(varname)//'" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC_NN*sizeof(Tipo_I4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I4

  function VTK_VAR_XML_SCAL_I2(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving scalar variable (I2P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  character(*), intent(IN):: varname      ! variable name
  integer(I2P), intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  integer(I4P)::             n1           ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int16" Name="'&
               //trim(varname)//'" NumberOfComponents="1" format="ascii">'
    write(unit=Unit_VTK,fmt=FI2P, iostat=E_IO)var
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int16" Name="'&
               //trim(varname)//'" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC_NN*sizeof(Tipo_I2)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I2',NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I2

  function VTK_VAR_XML_SCAL_I1(NC_NN,varname,var) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving scalar variable (I1P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN        ! number of cells or nodes
  character(*), intent(IN):: varname      ! variable name
  integer(I1P), intent(IN):: var(1:NC_NN) ! variable to be saved
  integer(I4P)::             E_IO         ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer     ! buffer string
  integer(I4P)::             n1           ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="'&
               //trim(varname)//'" NumberOfComponents="1" format="ascii">'
    write(unit=Unit_VTK,fmt=FI1P, iostat=E_IO)var
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="'&
              //trim(varname)//'" NumberOfComponents="1" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = NC_NN*sizeof(Tipo_I1)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I1',NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I1

  function VTK_VAR_XML_VECT_R8(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving vectorial variable (R8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of cells or nodes
  character(*), intent(IN):: varname       ! variable name
  real(R8P),    intent(IN):: varX(1:NC_NN) ! x component
  real(R8P),    intent(IN):: varY(1:NC_NN) ! y component
  real(R8P),    intent(IN):: varZ(1:NC_NN) ! z component
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer      ! buffer string
  integer(I4P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" Name="'//&
                       trim(varname)//'" NumberOfComponents="3" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FR8P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float64" Name="'//&
           trim(varname)//'" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NC_NN*sizeof(Tipo_R8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R8',3*NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_R8

  function VTK_VAR_XML_VECT_R4(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving vectorial variable (R4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of cells or nodes
  character(*), intent(IN):: varname       ! variable name
  real(R4P),    intent(IN):: varX(1:NC_NN) ! x component
  real(R4P),    intent(IN):: varY(1:NC_NN) ! y component
  real(R4P),    intent(IN):: varZ(1:NC_NN) ! z component
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer      ! buffer string
  integer(I4P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'&
            //trim(varname)//'" NumberOfComponents="3" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FR4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Float32" Name="'//trim(varname)//&
                  '" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NC_NN*sizeof(Tipo_R4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'R4',3*NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_R4

  function VTK_VAR_XML_VECT_I8(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving vectorial variable (I8P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of cells or nodes
  character(*), intent(IN):: varname       ! variable name
  integer(I8P), intent(IN):: varX(1:NC_NN) ! x component
  integer(I8P), intent(IN):: varY(1:NC_NN) ! y component
  integer(I8P), intent(IN):: varZ(1:NC_NN) ! z component
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer      ! buffer string
  integer(I4P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//&
                '<DataArray type="Int64" Name="'//trim(varname)//'" NumberOfComponents="3" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FI8P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int64" Name="'//trim(varname)//&
              '" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NC_NN*sizeof(Tipo_I8)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I8',3*NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I8

  function VTK_VAR_XML_VECT_I4(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving vectorial variable (I4P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of cells or nodes
  character(*), intent(IN):: varname       ! variable name
  integer(I4P), intent(IN):: varX(1:NC_NN) ! x component
  integer(I4P), intent(IN):: varY(1:NC_NN) ! y component
  integer(I4P), intent(IN):: varZ(1:NC_NN) ! z component
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer      ! buffer string
  integer(I4P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//&
              '<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="3" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FI4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int32" Name="'//&
                  trim(varname)//'" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NC_NN*sizeof(Tipo_I4)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I4',3*NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I4

  function VTK_VAR_XML_VECT_I2(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving vectorial variable (I2P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of cells or nodes
  character(*), intent(IN):: varname       ! variable name
  integer(I2P), intent(IN):: varX(1:NC_NN) ! x component
  integer(I2P), intent(IN):: varY(1:NC_NN) ! y component
  integer(I2P), intent(IN):: varZ(1:NC_NN) ! z component
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer      ! buffer string
  integer(I4P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//&
           '<DataArray type="Int16" Name="'//trim(varname)//'" NumberOfComponents="3" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FI2P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int16" Name="'//&
      trim(varname)//'" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NC_NN*sizeof(Tipo_I2)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I2',3*NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I2

  function VTK_VAR_XML_VECT_I1(NC_NN,varname,varX,varY,varZ) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving vectorial variable (I1P).
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P), intent(IN):: NC_NN         ! number of cells or nodes
  character(*), intent(IN):: varname       ! variable name
  integer(I1P), intent(IN):: varX(1:NC_NN) ! x component
  integer(I1P), intent(IN):: varY(1:NC_NN) ! y component
  integer(I1P), intent(IN):: varZ(1:NC_NN) ! z component
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(len=maxlen)::    s_buffer      ! buffer string
  integer(I4P)::             n1            ! counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="'//trim(varname)//&
    '" NumberOfComponents="3" format="ascii">'
    write(unit=Unit_VTK,fmt='(3'//FI1P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,fmt='(A)',          iostat=E_IO)repeat(' ',indent)//'</DataArray>'
  case(f_out_binary)
    write(s_buffer,fmt='(I8)', iostat=E_IO)ioffset
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<DataArray type="Int8" Name="'//trim(varname)//&
     '" NumberOfComponents="3" format="appended" offset="',trim(s_buffer),'">'//end_rec
    N_Byte  = 3*NC_NN*sizeof(Tipo_I1)
    ioffset = ioffset + sizeof(Tipo_I4) + N_Byte
    write(unit=Unit_VTK_Append,iostat=E_IO)N_Byte,'I1',3*NC_NN
    write(unit=Unit_VTK_Append,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</DataArray>'//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I1
  !(doc/)skippedblock

  function VTK_END_XML() result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!This function is used to finalize the file opened. The \LIBVTKIO manages the file unit without the user's action.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P)::              E_IO      ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  character(2)::              var_type  ! var\_type = R8,R4,I8,I4,I2,I1
  real(R8P),    allocatable:: v_R8(:)   ! R8 vector for IO in AppendData
  real(R4P),    allocatable:: v_R4(:)   ! R4 vector for IO in AppendData
  integer(I8P), allocatable:: v_I8(:)   ! I8 vector for IO in AppendData
  integer(I4P), allocatable:: v_I4(:)   ! I4 vector for IO in AppendData
  integer(I2P), allocatable:: v_I2(:)   ! I2 vector for IO in AppendData
  integer(I1P), allocatable:: v_I1(:)   ! I1 vector for IO in AppendData
  integer(I4P)::              N_v       ! vector dimension
  integer(I4P)::              n1        ! counter
  !!The following is an example of VTK\_END\_XML calling:
  !!
  !!\begin{boxred}{VTK\_END\_XML Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_END_XML()
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  select case(f_out)
  case(f_out_ascii)
    indent = indent - 2
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)repeat(' ',indent)//'</'//trim(topology)//'>'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'</VTKFile>'
  case(f_out_binary)
    indent = indent - 2
    write(unit  =Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec
    write(unit  =Unit_VTK,       iostat=E_IO)repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
    write(unit  =Unit_VTK,       iostat=E_IO)'_'
    endfile(unit=Unit_VTK_Append,iostat=E_IO)
    rewind(unit =Unit_VTK_Append,iostat=E_IO)
    do
      read(unit=Unit_VTK_Append,iostat=E_IO,end=100)N_Byte,var_type,N_v
      select case(var_type)
      case('R8')
        allocate(v_R8(1:N_v))
        read(unit =Unit_VTK_Append,iostat=E_IO)(v_R8(n1),n1=1,N_v)
        write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_R8(n1),n1=1,N_v)
        deallocate(v_R8)
      case('R4')
        allocate(v_R4(1:N_v))
        read(unit =Unit_VTK_Append,iostat=E_IO)(v_R4(n1),n1=1,N_v)
        write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_R4(n1),n1=1,N_v)
        deallocate(v_R4)
      case('I8')
        allocate(v_I8(1:N_v))
        read(unit =Unit_VTK_Append,iostat=E_IO)(v_I8(n1),n1=1,N_v)
        write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_I8(n1),n1=1,N_v)
        deallocate(v_I8)
      case('I4')
        allocate(v_I4(1:N_v))
        read(unit =Unit_VTK_Append,iostat=E_IO)(v_I4(n1),n1=1,N_v)
        write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_I4(n1),n1=1,N_v)
        deallocate(v_I4)
      case('I2')
        allocate(v_I2(1:N_v))
        read(unit =Unit_VTK_Append,iostat=E_IO)(v_I2(n1),n1=1,N_v)
        write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_I2(n1),n1=1,N_v)
        deallocate(v_I2)
      case('I1')
        allocate(v_I1(1:N_v))
        read(unit =Unit_VTK_Append,iostat=E_IO)(v_I1(n1),n1=1,N_v)
        write(unit=Unit_VTK,       iostat=E_IO)N_Byte,(v_I1(n1),n1=1,N_v)
        deallocate(v_I1)
      endselect
    enddo
    100 continue
    write(unit=Unit_VTK,iostat=E_IO)end_rec
    write(unit=Unit_VTK,iostat=E_IO)repeat(' ',indent)//'</AppendedData>'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'</VTKFile>'//end_rec
    ! closing AppendData file
    close(unit=Unit_VTK_Append,iostat=E_IO)
  endselect
  close(unit=Unit_VTK,iostat=E_IO)
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_END_XML
endmodule vtk_writer 
!!
!!\appendix
!!
!!\chapter{LIB\_VTK\_IO Usage Example}
!!\label{cap:example}
!!\minitoc
!!
!!\vspace*{8mm}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} usage of \LIBVTKIO is quite simple. In this chapter there are some
!!example of \LIBVTKIO usage. Some of the following examples are present also in the file \MaiuscolettoBS{Test\_LIB\_VTK\_IO.f90}
!!distributed within the \LIBVTKIO.
!!
!!\section{Legacy Rectilinear Grid}
!!\label{sec:example LRECTG}
!!
!!\begin{boxred}{Legacy Rectilinear Grid}
!!\begin{verbatim}
!!...
!!integer(4), intent(IN)::   Nx
!!real(8),    intent(IN)::           p(1:Nx)
!!real(8),    intent(IN)::         rho(1:Nx)
!!real(8),    intent(IN)::           u(1:Nx)
!!real(8),    intent(IN)::       gamma(1:Nx)
!!character(*), intent(IN):: filename
!!real(8)::                  x(1:Nx)
!!integer(4)::               i
!!...
!!x=(/(i, i=1, Nx, 1)/)
!!E_IO = VTK_INI(output_format = 'ascii',                &
!!               filene        = trim(filename)//'.vtk', &
!!               title         = 'Field',                &
!!               mesh_topology = 'RECTILINEAR_GRID')
!!E_IO = VTK_GEO(Nx        = Nx,        &
!!               Ny        = 1,         &
!!               Nz        = 1,         &
!!               X         = x,         &
!!               Y         = (/0.0_8/), &
!!               Z         = (/0.0_8/))
!!E_IO = VTK_DAT(NC_NN   = Nx,      &
!!               tipo    = 'node')
!!E_IO = VTK_VAR(NC_NN   = Nx,      &
!!               varname = 'p',     &
!!               var     = p)
!!E_IO = VTK_VAR(NC_NN   = Nx,      &
!!               varname = 'rho',   &
!!               var     = rho)
!!E_IO = VTK_VAR(NC_NN   = Nx,      &
!!               varname = 'u',     &
!!               var     = u)
!!E_IO = VTK_VAR(NC_NN   = Nx,      &
!!               varname = 'gamma', &
!!               var     = gamma)
!!E_IO = VTK_VAR(NC_NN   = Nx,      &
!!               varname = 'a',     &
!!               var     = sqrt(gamma*p/rho))
!!E_IO = VTK_END()
!!...
!!\end{verbatim}
!!\end{boxred}
!!
!!\section{XML Rectilinear Grid}
!!\label{sec:example XRECTG}
!!
!!\begin{boxred}{XML Rectilinear Grid}
!!\begin{verbatim}
!!...
!!integer(4),   intent(IN):: n
!!integer(4),   intent(IN):: Nx
!!real(8),      intent(IN)::     p(1:Nx)
!!real(8),      intent(IN)::   rho(1:Nx)
!!real(8),      intent(IN)::     u(1:Nx)
!!real(8),      intent(IN):: gamma(1:Nx)
!!character(*), intent(IN):: filename
!!real(8)::                  x(1:Nx)
!!integer(4)::               i
!!...
!!x=(/(i, i=1, Nx, 1)/)
!!E_IO = VTK_INI_XML(output_format = 'ascii',                &
!!                   filename      = trim(filename)//'.vtr', &
!!                   mesh_topology = 'RectilinearGrid',      &
!!                   nx1=1,nx2=Nx,ny1=1,ny2=1,nz1=1,nz2=1)
!!E_IO = VTK_GEO_XML(nx1=1,nx2=Nx,ny1=1,ny2=1,nz1=1,nz2=1, &
!!                   X=x,Y=(/0.0_8/),Z=(/0.0_8/))
!!E_IO = VTK_DAT_XML(tipo    = 'node',   &
!!                   azione  = 'OPEN')
!!E_IO = VTK_VAR_XML(NC_NN   = Nx,      &
!!                   varname = 'p',     &
!!                   var     = p)
!!E_IO = VTK_VAR_XML(NC_NN   = Nx,      &
!!                   varname = 'rho',   &
!!                   var     = rho)
!!E_IO = VTK_VAR_XML(NC_NN   = Nx,      &
!!                   varname = 'u',     &
!!                   var     = u)
!!E_IO = VTK_VAR_XML(NC_NN   = Nx,      &
!!                   varname = 'gamma', &
!!                   var     = gamma)
!!E_IO = VTK_VAR_XML(NC_NN   = Nx,      &
!!                   varname = 'a',     &
!!                   var     = sqrt(gamma*p/rho))
!!E_IO = VTK_DAT_XML(tipo    = 'node',  &
!!                   azione  = 'CLOSE')
!!E_IO = VTK_GEO_XML()
!!E_IO = VTK_END_XML()
!!...
!!\end{verbatim}
!!\end{boxred}
!!
!!\section{Legacy Unstructured Grid}
!!\label{sec:example LUNSTG}
!!
!!\begin{boxred}{Legacy Unstructured Grid}
!!\begin{verbatim}
!!...
!!integer(4), parameter::       Nn   = 27
!!integer(4), parameter::       Ne   = 11
!!real(4),    dimension(1:Nn):: x_uns
!!real(4),    dimension(1:Nn):: y_uns
!!real(4),    dimension(1:Nn):: z_uns
!!integer(4), dimension(1:Ne):: tipo
!!integer(4), dimension(1:60):: connect
!!real(8),    dimension(1:Nn):: var_uns_grid
!!integer(4), dimension(1:Nn):: var_uns_grid_X
!!integer(4), dimension(1:Nn):: var_uns_grid_Y
!!integer(4), dimension(1:Nn):: var_uns_grid_Z
!!...
!!E_IO = VTK_INI(output_format  = 'BINARY',                   &
!!               filename       = 'UNST_GRID_BIN.vtk',        &
!!               title          = 'Unstructured Grid Example' &
!!               mesh_topology  = 'UNSTRUCTURED_GRID')
!!
!!x_uns=(/0,1,2,0,1,2, &
!!        0,1,2,0,1,2, &
!!        0,1,2,0,1,2, &
!!        0,1,2,0,1,2, &
!!        0,1,2/)
!!y_uns=(/0,0,0,1,1,1, &
!!        0,0,0,1,1,1, &
!!        1,1,1,1,1,1, &
!!        1,1,1,1,1,1, &
!!        1,1,1/)
!!z_uns=(/0,0,0,0,0,0, &
!!        1,1,1,1,1,1, &
!!        2,2,2,3,3,3, &
!!        4,4,4,5,5,5, &
!!        6,6,6/)
!!
!!E_IO = VTK_GEO(Nnodi = Nn, &
!!               X=x_uns,Y=y_uns,Z=z_uns)
!!
!!connect = (/ 8, 0, 1, 4, 3, 6, 7,10, 9, &
!!             8, 1, 2, 5, 4, 7, 8,11,10, &
!!             4, 6,10, 9,12,             &
!!             4, 5,11,10,14,             &
!!             6,15,16,17,14,13,12,       &
!!             6,18,15,19,16,20,17,       &
!!             4,22,23,20,19,             &
!!             3,21,22,18,                &
!!             3,22,19,18,                &
!!             2,26,25,                   &
!!             1,24/)
!!tipo = (/12, &
!!         12, &
!!         10, &
!!         10, &
!!          7, &
!!          6, &
!!          9, &
!!          5, &
!!          5, &
!!          3, &
!!          1/)
!!E_IO = VTK_CON(NCelle  = Ne,       &
!!               connect = connect,  &
!!               tipo    = tipo)
!!E_IO = VTK_DAT(NC_NN   = Nn,       &
!!               tipo    = 'node')
!!
!!var_uns_grid =(/ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, &
!!                 6.0, 7.0, 8.0, 9.0,10.0,11.0, &
!!                12.0,13.0,14.0,15.0,16.0,17.0, &
!!                18.0,19.0,20.0,21.0,22.0,23.0, &
!!                24.0,25.0,26.0/)
!!
!!E_IO = VTK_VAR(NC_NN   = Nn,        &
!!               varname = 'scalars', &
!!               var     = var_uns_grid)
!!
!!var_uns_grid_X=(/1,1,0,1,1,0, &
!!                 1,1,0,1,1,0, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0/)
!!var_uns_grid_Y=(/0,1,2,0,1,2, &
!!                 0,1,2,0,1,2, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0/)
!!var_uns_grid_Z=(/0,0,0,0,0,0, &
!!                 0,0,0,0,0,0, &
!!                 1,1,1,1,1,1, &
!!                 1,1,1,1,1,1, &
!!                 1,1,1/)
!!E_IO = VTK_VAR(NC_NN   = Nn,             &
!!               varname = 'vectors',      &
!!               varX    = var_uns_grid_X, &
!!               varY    = var_uns_grid_Y, &
!!               varZ    = var_uns_grid_Z)
!!E_IO = VTK_END()
!!...
!!\end{verbatim}
!!\end{boxred}
!!
!!\section{XML Unstructured Grid}
!!\label{sec:example XUNSTG}
!!
!!\begin{boxred}{XML Unstructured Grid}
!!\begin{verbatim}
!!...
!!integer(4), parameter::       Nn   = 27
!!integer(4), parameter::       Ne   = 11
!!real(4),    dimension(1:Nn):: x_uns
!!real(4),    dimension(1:Nn):: y_uns
!!real(4),    dimension(1:Nn):: z_uns
!!integer(4), dimension(1:Ne):: tipo
!!integer(4), dimension(1:49):: connect_xml
!!integer(4), dimension(1:Ne):: offset_xml
!!real(8),    dimension(1:Nn):: var_uns_grid
!!integer(4), dimension(1:Nn):: var_uns_grid_X
!!integer(4), dimension(1:Nn):: var_uns_grid_Y
!!integer(4), dimension(1:Nn):: var_uns_grid_Z
!!...
!!E_IO = VTK_INI_XML(output_format = 'BINARY',              &
!!                   filename      = 'XML_UNST_BINARY.vtu', &
!!                   mesh_topology = 'UnstructuredGrid')
!!
!!x_uns=(/0,1,2,0,1,2, &
!!        0,1,2,0,1,2, &
!!        0,1,2,0,1,2, &
!!        0,1,2,0,1,2, &
!!        0,1,2/)
!!y_uns=(/0,0,0,1,1,1, &
!!        0,0,0,1,1,1, &
!!        1,1,1,1,1,1, &
!!        1,1,1,1,1,1, &
!!        1,1,1/)
!!z_uns=(/0,0,0,0,0,0, &
!!        1,1,1,1,1,1, &
!!        2,2,2,3,3,3, &
!!        4,4,4,5,5,5, &
!!        6,6,6/)
!!
!!E_IO = VTK_GEO_XML(Nnodi     = Nn, &
!!                   NCelle    = Ne, &
!!                   X=x_uns,Y=y_uns,Z=z_uns)
!!
!!connect_xml = (/ 0, 1, 4, 3, 6, 7,10, 9, &
!!                 1, 2, 5, 4, 7, 8,11,10, &
!!                 6,10, 9,12,             &
!!                 5,11,10,14,             &
!!                15,16,17,14,13,12,       &
!!                18,15,19,16,20,17,       &
!!                22,23,20,19,             &
!!                21,22,18,                &
!!                22,19,18,                &
!!                26,25,                   &
!!                24/)
!!offset_xml = (/ 8, &
!!               16, &
!!               20, &
!!               24, &
!!               30, &
!!               36, &
!!               40, &
!!               43, &
!!               46, &
!!               48, &
!!               49/)
!!
!!E_IO = VTK_CON_XML(NCelle  = Ne,          &
!!                   connect = connect_xml, &
!!                   offset  = offset_xml,  &
!!                   tipo    = (/12_1, &
!!                               12_1, &
!!                               10_1, &
!!                               10_1, &
!!                                7_1, &
!!                                6_1, &
!!                                9_1, &
!!                                5_1, &
!!                                5_1, &
!!                                3_1, &
!!                                1_1/))
!!
!!E_IO = VTK_DAT_XML(tipo    = 'node', &
!!                   azione  = 'OPEN')
!!
!!var_uns_grid =(/ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, &
!!                 6.0, 7.0, 8.0, 9.0,10.0,11.0, &
!!                12.0,13.0,14.0,15.0,16.0,17.0, &
!!                18.0,19.0,20.0,21.0,22.0,23.0, &
!!                24.0,25.0,26.0/)
!!
!!E_IO = VTK_VAR_XML(NC_NN   = Nn,        &
!!                   varname = 'scalars', &
!!                   var     = var_uns_grid)
!!
!!var_uns_grid_X=(/1,1,0,1,1,0, &
!!                 1,1,0,1,1,0, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0/)
!!var_uns_grid_Y=(/0,1,2,0,1,2, &
!!                 0,1,2,0,1,2, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0,0,0,0, &
!!                 0,0,0/)
!!var_uns_grid_Z=(/0,0,0,0,0,0, &
!!                 0,0,0,0,0,0, &
!!                 1,1,1,1,1,1, &
!!                 1,1,1,1,1,1, &
!!                 1,1,1/)
!!
!!E_IO = VTK_VAR_XML(NC_NN   = Nn,             &
!!                   varname = 'vector',       &
!!                   varX    = var_uns_grid_X, &
!!                   varY    = var_uns_grid_Y, &
!!                   varZ    = var_uns_grid_Z)
!!E_IO = VTK_DAT_XML(tipo    = 'node',   &
!!                   azione  = 'CLOSE')
!!E_IO = VTK_GEO_XML()
!!E_IO = VTK_END_XML()
!!...
!!\end{verbatim}
!!\end{boxred}
!!
!!\chapter{Fortran \& Portable-Kind-Precision Selection}
!!\label{cap:kind precision}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf F}}{ortran} is the most popular programming language for scientific computing.
!!With fortran it is quite simple obtain fast code and manage large multidimensional array. Because fortran permits the achivment
!!of high performance it is also used on great range of different computer-architettures, and often on the fastest supercomputer
!!in the world. Therefore fortran programs must be \MaiuscolettoBS{portable}: portability means that the code will give the same
!!results on every different computer-architettures. One of the most important goal of the numeric code is to control the
!!\MaiuscolettoBS{the numeric error} due to finite precision of numerical operations. Fortran uses the \MaiuscolettoBS{IEEE
!!rappresentations}; integers and reals (floating point) are represented with a finite precision. So when the code computes an
!!operation it has a \MaiuscolettoBS{trunction error} due to the truncation of the numerical finite rappresentaions. For numerical
!!and more in general scientific applications this source of errors must be controlled. The programmer must know which is the
!!precision associated to the code variables. Before the standard fortran 90/95 there are not any way to select the precision of
!!the numerical variables in a portable fashion. With the possibility to specify a kind parameter for variables, the standard
!!fortran 90/95 makes avaible two useful functions to select the kind precision of integers and reals:
!!
!!\begin{boxred}{selected\_real\_kind \& selected\_int\_kind}
!!\begin{verbatim}
!!function selected_real_kind(p,r) result(kind_id)
!!integer, intent(IN), optional:: p
!!integer, intent(IN), optional:: r
!!integer::                       kind_id
!!
!!The result, kind_id, is a scalar of type default integer.
!!If both arguments are absent, the result is zero.
!!Otherwise, the result has a value equal to a value of
!!the kind parameter of a real data type with decimal
!!precision, as returned by the function PRECISION, of at
!!least p digits and a decimal exponent range, as returned
!!by the function RANGE, of at least r.
!!
!!function selected_int_kind(p) result(kind_id)
!!integer, intent(IN), optional:: p
!!integer::                       kind_id
!!
!!The result, kind_id, is a scalar of type default integer.
!!The result has a value equal to the value of the kind
!!parameter of the integer data type that represents all
!!values n in the range of about values n with
!!-10^p < n < 10^p.
!!\end{verbatim}
!!\end{boxred}
!!
!!Using these two functions the programmer can accurately control the precision of its own variables in a portable manner.
!!Note that specifing the kind precision without using these two functions is not portable: $real(8)$ means different
!!precisions on different architettures. Parametrizing the kind of all numerical variables using these two functions makes
!!the portable. The \LIBVTKIO uses this principle to achive portable-kind-precision selection; in the library are defined
!!some parameters by which all variables kind-precisions are parametrized:
!!
!!\begin{boxblu}{\LIBVTKIO Kind-Precision Parameters}
!!{\color{RoyalBlue}\MaiuscolettoS{Real Precision Definitions}}
!!\begin{description}
!! \item [{\color{RoyalBlue}R16P}] real with $33$ digits, range $[+-10^{-4931},+-10^{+4931}-1]$
!! \item [{\color{RoyalBlue}R8P}]  real with $15$ digits, range $[+-10^{-307} ,+-10^{+307}-1 ]$
!! \item [{\color{RoyalBlue}R4P}]  real with $6$  digits, range $[+-10^{-37}  ,+-10^+{37}-1  ]$
!!\end{description}
!!{\color{RoyalBlue}\MaiuscolettoS{Integer Precision Definitions}}
!!\begin{description}
!! \item [{\color{RoyalBlue}I8P}] range $[-2^{63},+2^{63}-1]$
!! \item [{\color{RoyalBlue}I4P}] range $[-2^{31},+2^{31}-1]$
!! \item [{\color{RoyalBlue}I2P}] range $[-2^{15},+2^{15}-1]$
!! \item [{\color{RoyalBlue}I1P}] range $[-2^{7} ,+2^{7} -1]$
!!\end{description}
!!\end{boxblu}
!!
!!In order to avoid strange results porting your code the use of parametrized-kind-precision is very useful. The \LIBVTKIO
!!makes avaible to the external its own kind-parameters that can be used to parametrize the code.
!!
!!\chapter{Dynamic Dispatching}
!!\label{cap:Dynamic Dispatching}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf F}}{ortran} is not an \MaiuscolettoBS{object oriented} (OOp) programming
!!language. It is a procedural language with some of the the goals (ineritance, user-definited data type, polimorphism...)
!!of OOp. Fortran most important aim is to ensure the performance of the code not its \virgo{friendliness}... Despite its
!!nature, fortran 90/95 makes avaible some interesting features: it permits the dynamic dispatching of functions and
!!subroutine ensuring the best performance. This goal is achived with use of $interface$ construct. In the \LIBVTKIO there are,
!!at today, 4 interface blocks:
!!
!!\begin{boxred}{\LIBVTKIO Interface Blocks}
!!\begin{verbatim}
!!interface VTK_GEO
!!  module procedure VTK_GEO_UNST_R8, &
!!                   VTK_GEO_UNST_R4, &
!!                   VTK_GEO_STRP_R8, &
!!                   VTK_GEO_STRP_R4, &
!!                   VTK_GEO_STRG_R8, &
!!                   VTK_GEO_STRG_R4, &
!!                   VTK_GEO_RECT_R8, &
!!                   VTK_GEO_RECT_R4
!!endinterface
!!
!!interface VTK_VAR
!!  module procedure VTK_VAR_SCAL_R8, &
!!                   VTK_VAR_SCAL_R4, &
!!                   VTK_VAR_SCAL_I4, &
!!                   VTK_VAR_VECT_R8, &
!!                   VTK_VAR_VECT_R4, &
!!                   VTK_VAR_VECT_I4, &
!!                   VTK_VAR_TEXT_R8, &
!!                   VTK_VAR_TEXT_R4
!!endinterface
!!
!!interface VTK_GEO_XML
!!  module procedure VTK_GEO_XML_STRG_R4, &
!!                   VTK_GEO_XML_STRG_R8, &
!!                   VTK_GEO_XML_RECT_R8, &
!!                   VTK_GEO_XML_RECT_R4, &
!!                   VTK_GEO_XML_UNST_R8, &
!!                   VTK_GEO_XML_UNST_R4, &
!!                   VTK_GEO_XML_CLOSEP
!!endinterface
!!
!!interface VTK_VAR_XML
!!  module procedure VTK_VAR_XML_SCAL_R8, &
!!                   VTK_VAR_XML_SCAL_R4, &
!!                   VTK_VAR_XML_SCAL_I8, &
!!                   VTK_VAR_XML_SCAL_I4, &
!!                   VTK_VAR_XML_SCAL_I2, &
!!                   VTK_VAR_XML_SCAL_I1, &
!!                   VTK_VAR_XML_VECT_R8, &
!!                   VTK_VAR_XML_VECT_R4, &
!!                   VTK_VAR_XML_VECT_I8, &
!!                   VTK_VAR_XML_VECT_I4, &
!!                   VTK_VAR_XML_VECT_I2, &
!!                   VTK_VAR_XML_VECT_I1
!!endinterface
!!\end{verbatim}
!!\end{boxred}
!!
!!By the interface construct \LIBVTKIO has a more simple API. The user deals with a few functions without non-sense-long-name...
!!Dynamic dispatching is not the magic wand to solve all problems but it is an useful tool to simplify the code API. It is
!!not powerful as the C++ template, but it is a \MaiuscolettoBS{quantum-leap} for fortran programmers.
!!
!!\chapter{Known Bugs}
!!\label{cap:BUG}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} \LIBVTKIO is a very young project and it is a good example of wrong
!!programming style... It is unstable and not tested. It is used by only one user (... me of course!) and there are a lot of
!!bugs that are still hidden. At the moment several features are missing (the input functions and the poly-data topology...),
!!but it is useful to export fortran data to VTK standard, and this goal was the most important for me.
!!
!!At today only one main bug was found. Fortran allows the automatic reshape of arrays: as an example 2D array can be
!!automatically (in the function calling) transformed  to a 1D array with the same number of element of 2D array. The use of
!!dynamic dispatching had disable this feature: dynamic dispatching use the array-shape information to dectet, at compile-time,
!!the correct function to be called. So reshape arrays at calling phase is not allowed. In the next release I will fix this bug
!!introducing the function to reshape arrays between 1D, 2D and 3D arrays.
!!
!!A possible, not already found, bug is the non correct kind detection. It is possible that a code uses kind-precision parameter
!!that does not match the \LIBVTKIO parameters. I never observe this bug but it is possible. To avoid it the simple way is to use
!!always the \LIBVTKIO kind-precision parameters; if the parameters actually present do not match your necessities, define new
!!parameters in \LIBVTKIO and redistribuite \LIBVTKIO with your pacth!
!!
!!Finally there is a strong inefficiency when saving XML binary file. To write XML binary \LIBVTKIO uses a temporary scratch file
!!to save binary data while saving all formatting data to the final XML file; only when all XML formatting data have been written
!!the scratch file is rewinded and the binary data is saved in the final tag of XML file as \MaiuscolettoBS{raw} data. This
!!algorithm is obviously inefficient. Any tip is welcome!
!!
!!\chapter{GNU GENERAL PUBLIC LICENSE}
!!\label{cap:GPL}
!!
!!\begin{center}
!!\MaiuscolettoS{Version 3, 29 June 2007}
!!
!!{\parindent 0in
!!
!!Copyright \copyright\ 2007 Free Software Foundation, Inc. \texttt{http://fsf.org/}
!!
!!\bigskip
!!Everyone is permitted to copy and distribute verbatim copies of this
!!
!!license document, but changing it is not allowed.}
!!
!!\end{center}
!!
!!
!!\section*{Preamble}
!!The GNU General Public License is a free, copyleft license for
!!software and other kinds of works.
!!
!!The licenses for most software and other practical works are designed
!!to take away your freedom to share and change the works.  By contrast,
!!the GNU General Public License is intended to guarantee your freedom to
!!share and change all versions of a program--to make sure it remains free
!!software for all its users.  We, the Free Software Foundation, use the
!!GNU General Public License for most of our software; it applies also to
!!any other work released this way by its authors.  You can apply it to
!!your programs, too.
!!
!!When we speak of free software, we are referring to freedom, not
!!price.  Our General Public Licenses are designed to make sure that you
!!have the freedom to distribute copies of free software (and charge for
!!them if you wish), that you receive source code or can get it if you
!!want it, that you can change the software or use pieces of it in new
!!free programs, and that you know you can do these things.
!!
!!To protect your rights, we need to prevent others from denying you
!!these rights or asking you to surrender the rights.  Therefore, you have
!!certain responsibilities if you distribute copies of the software, or if
!!you modify it: responsibilities to respect the freedom of others.
!!
!!For example, if you distribute copies of such a program, whether
!!gratis or for a fee, you must pass on to the recipients the same
!!freedoms that you received.  You must make sure that they, too, receive
!!or can get the source code.  And you must show them these terms so they
!!know their rights.
!!
!!Developers that use the GNU GPL protect your rights with two steps:
!!(1) assert copyright on the software, and (2) offer you this License
!!giving you legal permission to copy, distribute and/or modify it.
!!
!!For the developers' and authors' protection, the GPL clearly explains
!!that there is no warranty for this free software.  For both users' and
!!authors' sake, the GPL requires that modified versions be marked as
!!changed, so that their problems will not be attributed erroneously to
!!authors of previous versions.
!!
!!Some devices are designed to deny users access to install or run
!!modified versions of the software inside them, although the manufacturer
!!can do so.  This is fundamentally incompatible with the aim of
!!protecting users' freedom to change the software.  The systematic
!!pattern of such abuse occurs in the area of products for individuals to
!!use, which is precisely where it is most unacceptable.  Therefore, we
!!have designed this version of the GPL to prohibit the practice for those
!!products.  If such problems arise substantially in other domains, we
!!stand ready to extend this provision to those domains in future versions
!!of the GPL, as needed to protect the freedom of users.
!!
!!Finally, every program is threatened constantly by software patents.
!!States should not allow patents to restrict development and use of
!!software on general-purpose computers, but in those that do, we wish to
!!avoid the special danger that patents applied to a free program could
!!make it effectively proprietary.  To prevent this, the GPL assures that
!!patents cannot be used to render the program non-free.
!!
!!The precise terms and conditions for copying, distribution and
!!modification follow.
!!
!!
!!\begin{center}
!!{\Large \sc Terms and Conditions}
!!\end{center}
!!
!!\begin{enumerate}
!!
!!\addtocounter{enumi}{-1}
!!
!!\item Definitions.
!!
!!``This License'' refers to version 3 of the GNU General Public License.
!!
!!``Copyright'' also means copyright-like laws that apply to other kinds of
!!works, such as semiconductor masks.
!!
!!``The Program'' refers to any copyrightable work licensed under this
!!License.  Each licensee is addressed as ``you''.  ``Licensees'' and
!!``recipients'' may be individuals or organizations.
!!
!!To ``modify'' a work means to copy from or adapt all or part of the work
!!in a fashion requiring copyright permission, other than the making of an
!!exact copy.  The resulting work is called a ``modified version'' of the
!!earlier work or a work ``based on'' the earlier work.
!!
!!A ``covered work'' means either the unmodified Program or a work based
!!on the Program.
!!
!!To ``propagate'' a work means to do anything with it that, without
!!permission, would make you directly or secondarily liable for
!!infringement under applicable copyright law, except executing it on a
!!computer or modifying a private copy.  Propagation includes copying,
!!distribution (with or without modification), making available to the
!!public, and in some countries other activities as well.
!!
!!To ``convey'' a work means any kind of propagation that enables other
!!parties to make or receive copies.  Mere interaction with a user through
!!a computer network, with no transfer of a copy, is not conveying.
!!
!!An interactive user interface displays ``Appropriate Legal Notices''
!!to the extent that it includes a convenient and prominently visible
!!feature that (1) displays an appropriate copyright notice, and (2)
!!tells the user that there is no warranty for the work (except to the
!!extent that warranties are provided), that licensees may convey the
!!work under this License, and how to view a copy of this License.  If
!!the interface presents a list of user commands or options, such as a
!!menu, a prominent item in the list meets this criterion.
!!
!!\item Source Code.
!!
!!The ``source code'' for a work means the preferred form of the work
!!for making modifications to it.  ``Object code'' means any non-source
!!form of a work.
!!
!!A ``Standard Interface'' means an interface that either is an official
!!standard defined by a recognized standards body, or, in the case of
!!interfaces specified for a particular programming language, one that
!!is widely used among developers working in that language.
!!
!!The ``System Libraries'' of an executable work include anything, other
!!than the work as a whole, that (a) is included in the normal form of
!!packaging a Major Component, but which is not part of that Major
!!Component, and (b) serves only to enable use of the work with that
!!Major Component, or to implement a Standard Interface for which an
!!implementation is available to the public in source code form.  A
!!``Major Component'', in this context, means a major essential component
!!(kernel, window system, and so on) of the specific operating system
!!(if any) on which the executable work runs, or a compiler used to
!!produce the work, or an object code interpreter used to run it.
!!
!!The ``Corresponding Source'' for a work in object code form means all
!!the source code needed to generate, install, and (for an executable
!!work) run the object code and to modify the work, including scripts to
!!control those activities.  However, it does not include the work's
!!System Libraries, or general-purpose tools or generally available free
!!programs which are used unmodified in performing those activities but
!!which are not part of the work.  For example, Corresponding Source
!!includes interface definition files associated with source files for
!!the work, and the source code for shared libraries and dynamically
!!linked subprograms that the work is specifically designed to require,
!!such as by intimate data communication or control flow between those
!!subprograms and other parts of the work.
!!
!!The Corresponding Source need not include anything that users
!!can regenerate automatically from other parts of the Corresponding
!!Source.
!!
!!The Corresponding Source for a work in source code form is that
!!same work.
!!
!!\item Basic Permissions.
!!
!!All rights granted under this License are granted for the term of
!!copyright on the Program, and are irrevocable provided the stated
!!conditions are met.  This License explicitly affirms your unlimited
!!permission to run the unmodified Program.  The output from running a
!!covered work is covered by this License only if the output, given its
!!content, constitutes a covered work.  This License acknowledges your
!!rights of fair use or other equivalent, as provided by copyright law.
!!
!!You may make, run and propagate covered works that you do not
!!convey, without conditions so long as your license otherwise remains
!!in force.  You may convey covered works to others for the sole purpose
!!of having them make modifications exclusively for you, or provide you
!!with facilities for running those works, provided that you comply with
!!the terms of this License in conveying all material for which you do
!!not control copyright.  Those thus making or running the covered works
!!for you must do so exclusively on your behalf, under your direction
!!and control, on terms that prohibit them from making any copies of
!!your copyrighted material outside their relationship with you.
!!
!!Conveying under any other circumstances is permitted solely under
!!the conditions stated below.  Sublicensing is not allowed; section 10
!!makes it unnecessary.
!!
!!\item Protecting Users' Legal Rights From Anti-Circumvention Law.
!!
!!No covered work shall be deemed part of an effective technological
!!measure under any applicable law fulfilling obligations under article
!!11 of the WIPO copyright treaty adopted on 20 December 1996, or
!!similar laws prohibiting or restricting circumvention of such
!!measures.
!!
!!When you convey a covered work, you waive any legal power to forbid
!!circumvention of technological measures to the extent such circumvention
!!is effected by exercising rights under this License with respect to
!!the covered work, and you disclaim any intention to limit operation or
!!modification of the work as a means of enforcing, against the work's
!!users, your or third parties' legal rights to forbid circumvention of
!!technological measures.
!!
!!\item Conveying Verbatim Copies.
!!
!!You may convey verbatim copies of the Program's source code as you
!!receive it, in any medium, provided that you conspicuously and
!!appropriately publish on each copy an appropriate copyright notice;
!!keep intact all notices stating that this License and any
!!non-permissive terms added in accord with section 7 apply to the code;
!!keep intact all notices of the absence of any warranty; and give all
!!recipients a copy of this License along with the Program.
!!
!!You may charge any price or no price for each copy that you convey,
!!and you may offer support or warranty protection for a fee.
!!
!!\item Conveying Modified Source Versions.
!!
!!You may convey a work based on the Program, or the modifications to
!!produce it from the Program, in the form of source code under the
!!terms of section 4, provided that you also meet all of these conditions:
!!  \begin{enumerate}
!!  \item The work must carry prominent notices stating that you modified
!!  it, and giving a relevant date.
!!
!!  \item The work must carry prominent notices stating that it is
!!  released under this License and any conditions added under section
!!  7.  This requirement modifies the requirement in section 4 to
!!  ``keep intact all notices''.
!!
!!  \item You must license the entire work, as a whole, under this
!!  License to anyone who comes into possession of a copy.  This
!!  License will therefore apply, along with any applicable section 7
!!  additional terms, to the whole of the work, and all its parts,
!!  regardless of how they are packaged.  This License gives no
!!  permission to license the work in any other way, but it does not
!!  invalidate such permission if you have separately received it.
!!
!!  \item If the work has interactive user interfaces, each must display
!!  Appropriate Legal Notices; however, if the Program has interactive
!!  interfaces that do not display Appropriate Legal Notices, your
!!  work need not make them do so.
!!\end{enumerate}
!!A compilation of a covered work with other separate and independent
!!works, which are not by their nature extensions of the covered work,
!!and which are not combined with it such as to form a larger program,
!!in or on a volume of a storage or distribution medium, is called an
!!``aggregate'' if the compilation and its resulting copyright are not
!!used to limit the access or legal rights of the compilation's users
!!beyond what the individual works permit.  Inclusion of a covered work
!!in an aggregate does not cause this License to apply to the other
!!parts of the aggregate.
!!
!!\item Conveying Non-Source Forms.
!!
!!You may convey a covered work in object code form under the terms
!!of sections 4 and 5, provided that you also convey the
!!machine-readable Corresponding Source under the terms of this License,
!!in one of these ways:
!!  \begin{enumerate}
!!  \item Convey the object code in, or embodied in, a physical product
!!  (including a physical distribution medium), accompanied by the
!!  Corresponding Source fixed on a durable physical medium
!!  customarily used for software interchange.
!!
!!  \item Convey the object code in, or embodied in, a physical product
!!  (including a physical distribution medium), accompanied by a
!!  written offer, valid for at least three years and valid for as
!!  long as you offer spare parts or customer support for that product
!!  model, to give anyone who possesses the object code either (1) a
!!  copy of the Corresponding Source for all the software in the
!!  product that is covered by this License, on a durable physical
!!  medium customarily used for software interchange, for a price no
!!  more than your reasonable cost of physically performing this
!!  conveying of source, or (2) access to copy the
!!  Corresponding Source from a network server at no charge.
!!
!!  \item Convey individual copies of the object code with a copy of the
!!  written offer to provide the Corresponding Source.  This
!!  alternative is allowed only occasionally and noncommercially, and
!!  only if you received the object code with such an offer, in accord
!!  with subsection 6b.
!!
!!  \item Convey the object code by offering access from a designated
!!  place (gratis or for a charge), and offer equivalent access to the
!!  Corresponding Source in the same way through the same place at no
!!  further charge.  You need not require recipients to copy the
!!  Corresponding Source along with the object code.  If the place to
!!  copy the object code is a network server, the Corresponding Source
!!  may be on a different server (operated by you or a third party)
!!  that supports equivalent copying facilities, provided you maintain
!!  clear directions next to the object code saying where to find the
!!  Corresponding Source.  Regardless of what server hosts the
!!  Corresponding Source, you remain obligated to ensure that it is
!!  available for as long as needed to satisfy these requirements.
!!
!!  \item Convey the object code using peer-to-peer transmission, provided
!!  you inform other peers where the object code and Corresponding
!!  Source of the work are being offered to the general public at no
!!  charge under subsection 6d.
!!  \end{enumerate}
!!
!!A separable portion of the object code, whose source code is excluded
!!from the Corresponding Source as a System Library, need not be
!!included in conveying the object code work.
!!
!!A ``User Product'' is either (1) a ``consumer product'', which means any
!!tangible personal property which is normally used for personal, family,
!!or household purposes, or (2) anything designed or sold for incorporation
!!into a dwelling.  In determining whether a product is a consumer product,
!!doubtful cases shall be resolved in favor of coverage.  For a particular
!!product received by a particular user, ``normally used'' refers to a
!!typical or common use of that class of product, regardless of the status
!!of the particular user or of the way in which the particular user
!!actually uses, or expects or is expected to use, the product.  A product
!!is a consumer product regardless of whether the product has substantial
!!commercial, industrial or non-consumer uses, unless such uses represent
!!the only significant mode of use of the product.
!!
!!``Installation Information'' for a User Product means any methods,
!!procedures, authorization keys, or other information required to install
!!and execute modified versions of a covered work in that User Product from
!!a modified version of its Corresponding Source.  The information must
!!suffice to ensure that the continued functioning of the modified object
!!code is in no case prevented or interfered with solely because
!!modification has been made.
!!
!!If you convey an object code work under this section in, or with, or
!!specifically for use in, a User Product, and the conveying occurs as
!!part of a transaction in which the right of possession and use of the
!!User Product is transferred to the recipient in perpetuity or for a
!!fixed term (regardless of how the transaction is characterized), the
!!Corresponding Source conveyed under this section must be accompanied
!!by the Installation Information.  But this requirement does not apply
!!if neither you nor any third party retains the ability to install
!!modified object code on the User Product (for example, the work has
!!been installed in ROM).
!!
!!The requirement to provide Installation Information does not include a
!!requirement to continue to provide support service, warranty, or updates
!!for a work that has been modified or installed by the recipient, or for
!!the User Product in which it has been modified or installed.  Access to a
!!network may be denied when the modification itself materially and
!!adversely affects the operation of the network or violates the rules and
!!protocols for communication across the network.
!!
!!Corresponding Source conveyed, and Installation Information provided,
!!in accord with this section must be in a format that is publicly
!!documented (and with an implementation available to the public in
!!source code form), and must require no special password or key for
!!unpacking, reading or copying.
!!
!!\item Additional Terms.
!!
!!``Additional permissions'' are terms that supplement the terms of this
!!License by making exceptions from one or more of its conditions.
!!Additional permissions that are applicable to the entire Program shall
!!be treated as though they were included in this License, to the extent
!!that they are valid under applicable law.  If additional permissions
!!apply only to part of the Program, that part may be used separately
!!under those permissions, but the entire Program remains governed by
!!this License without regard to the additional permissions.
!!
!!When you convey a copy of a covered work, you may at your option
!!remove any additional permissions from that copy, or from any part of
!!it.  (Additional permissions may be written to require their own
!!removal in certain cases when you modify the work.)  You may place
!!additional permissions on material, added by you to a covered work,
!!for which you have or can give appropriate copyright permission.
!!
!!Notwithstanding any other provision of this License, for material you
!!add to a covered work, you may (if authorized by the copyright holders of
!!that material) supplement the terms of this License with terms:
!!  \begin{enumerate}
!!  \item Disclaiming warranty or limiting liability differently from the
!!  terms of sections 15 and 16 of this License; or
!!
!!  \item Requiring preservation of specified reasonable legal notices or
!!  author attributions in that material or in the Appropriate Legal
!!  Notices displayed by works containing it; or
!!
!!  \item Prohibiting misrepresentation of the origin of that material, or
!!  requiring that modified versions of such material be marked in
!!  reasonable ways as different from the original version; or
!!
!!  \item Limiting the use for publicity purposes of names of licensors or
!!  authors of the material; or
!!
!!  \item Declining to grant rights under trademark law for use of some
!!  trade names, trademarks, or service marks; or
!!
!!  \item Requiring indemnification of licensors and authors of that
!!  material by anyone who conveys the material (or modified versions of
!!  it) with contractual assumptions of liability to the recipient, for
!!  any liability that these contractual assumptions directly impose on
!!  those licensors and authors.
!!  \end{enumerate}
!!
!!All other non-permissive additional terms are considered ``further
!!restrictions'' within the meaning of section 10.  If the Program as you
!!received it, or any part of it, contains a notice stating that it is
!!governed by this License along with a term that is a further
!!restriction, you may remove that term.  If a license document contains
!!a further restriction but permits relicensing or conveying under this
!!License, you may add to a covered work material governed by the terms
!!of that license document, provided that the further restriction does
!!not survive such relicensing or conveying.
!!
!!If you add terms to a covered work in accord with this section, you
!!must place, in the relevant source files, a statement of the
!!additional terms that apply to those files, or a notice indicating
!!where to find the applicable terms.
!!
!!Additional terms, permissive or non-permissive, may be stated in the
!!form of a separately written license, or stated as exceptions;
!!the above requirements apply either way.
!!
!!\item Termination.
!!
!!You may not propagate or modify a covered work except as expressly
!!provided under this License.  Any attempt otherwise to propagate or
!!modify it is void, and will automatically terminate your rights under
!!this License (including any patent licenses granted under the third
!!paragraph of section 11).
!!
!!However, if you cease all violation of this License, then your
!!license from a particular copyright holder is reinstated (a)
!!provisionally, unless and until the copyright holder explicitly and
!!finally terminates your license, and (b) permanently, if the copyright
!!holder fails to notify you of the violation by some reasonable means
!!prior to 60 days after the cessation.
!!
!!Moreover, your license from a particular copyright holder is
!!reinstated permanently if the copyright holder notifies you of the
!!violation by some reasonable means, this is the first time you have
!!received notice of violation of this License (for any work) from that
!!copyright holder, and you cure the violation prior to 30 days after
!!your receipt of the notice.
!!
!!Termination of your rights under this section does not terminate the
!!licenses of parties who have received copies or rights from you under
!!this License.  If your rights have been terminated and not permanently
!!reinstated, you do not qualify to receive new licenses for the same
!!material under section 10.
!!
!!\item Acceptance Not Required for Having Copies.
!!
!!You are not required to accept this License in order to receive or
!!run a copy of the Program.  Ancillary propagation of a covered work
!!occurring solely as a consequence of using peer-to-peer transmission
!!to receive a copy likewise does not require acceptance.  However,
!!nothing other than this License grants you permission to propagate or
!!modify any covered work.  These actions infringe copyright if you do
!!not accept this License.  Therefore, by modifying or propagating a
!!covered work, you indicate your acceptance of this License to do so.
!!
!!\item Automatic Licensing of Downstream Recipients.
!!
!!Each time you convey a covered work, the recipient automatically
!!receives a license from the original licensors, to run, modify and
!!propagate that work, subject to this License.  You are not responsible
!!for enforcing compliance by third parties with this License.
!!
!!An ``entity transaction'' is a transaction transferring control of an
!!organization, or substantially all assets of one, or subdividing an
!!organization, or merging organizations.  If propagation of a covered
!!work results from an entity transaction, each party to that
!!transaction who receives a copy of the work also receives whatever
!!licenses to the work the party's predecessor in interest had or could
!!give under the previous paragraph, plus a right to possession of the
!!Corresponding Source of the work from the predecessor in interest, if
!!the predecessor has it or can get it with reasonable efforts.
!!
!!You may not impose any further restrictions on the exercise of the
!!rights granted or affirmed under this License.  For example, you may
!!not impose a license fee, royalty, or other charge for exercise of
!!rights granted under this License, and you may not initiate litigation
!!(including a cross-claim or counterclaim in a lawsuit) alleging that
!!any patent claim is infringed by making, using, selling, offering for
!!sale, or importing the Program or any portion of it.
!!
!!\item Patents.
!!
!!A ``contributor'' is a copyright holder who authorizes use under this
!!License of the Program or a work on which the Program is based.  The
!!work thus licensed is called the contributor's ``contributor version''.
!!
!!A contributor's ``essential patent claims'' are all patent claims
!!owned or controlled by the contributor, whether already acquired or
!!hereafter acquired, that would be infringed by some manner, permitted
!!by this License, of making, using, or selling its contributor version,
!!but do not include claims that would be infringed only as a
!!consequence of further modification of the contributor version.  For
!!purposes of this definition, ``control'' includes the right to grant
!!patent sublicenses in a manner consistent with the requirements of
!!this License.
!!
!!Each contributor grants you a non-exclusive, worldwide, royalty-free
!!patent license under the contributor's essential patent claims, to
!!make, use, sell, offer for sale, import and otherwise run, modify and
!!propagate the contents of its contributor version.
!!
!!In the following three paragraphs, a ``patent license'' is any express
!!agreement or commitment, however denominated, not to enforce a patent
!!(such as an express permission to practice a patent or covenant not to
!!sue for patent infringement).  To ``grant'' such a patent license to a
!!party means to make such an agreement or commitment not to enforce a
!!patent against the party.
!!
!!If you convey a covered work, knowingly relying on a patent license,
!!and the Corresponding Source of the work is not available for anyone
!!to copy, free of charge and under the terms of this License, through a
!!publicly available network server or other readily accessible means,
!!then you must either (1) cause the Corresponding Source to be so
!!available, or (2) arrange to deprive yourself of the benefit of the
!!patent license for this particular work, or (3) arrange, in a manner
!!consistent with the requirements of this License, to extend the patent
!!license to downstream recipients.  ``Knowingly relying'' means you have
!!actual knowledge that, but for the patent license, your conveying the
!!covered work in a country, or your recipient's use of the covered work
!!in a country, would infringe one or more identifiable patents in that
!!country that you have reason to believe are valid.
!!
!!If, pursuant to or in connection with a single transaction or
!!arrangement, you convey, or propagate by procuring conveyance of, a
!!covered work, and grant a patent license to some of the parties
!!receiving the covered work authorizing them to use, propagate, modify
!!or convey a specific copy of the covered work, then the patent license
!!you grant is automatically extended to all recipients of the covered
!!work and works based on it.
!!
!!A patent license is ``discriminatory'' if it does not include within
!!the scope of its coverage, prohibits the exercise of, or is
!!conditioned on the non-exercise of one or more of the rights that are
!!specifically granted under this License.  You may not convey a covered
!!work if you are a party to an arrangement with a third party that is
!!in the business of distributing software, under which you make payment
!!to the third party based on the extent of your activity of conveying
!!the work, and under which the third party grants, to any of the
!!parties who would receive the covered work from you, a discriminatory
!!patent license (a) in connection with copies of the covered work
!!conveyed by you (or copies made from those copies), or (b) primarily
!!for and in connection with specific products or compilations that
!!contain the covered work, unless you entered into that arrangement,
!!or that patent license was granted, prior to 28 March 2007.
!!
!!Nothing in this License shall be construed as excluding or limiting
!!any implied license or other defenses to infringement that may
!!otherwise be available to you under applicable patent law.
!!
!!\item No Surrender of Others' Freedom.
!!
!!If conditions are imposed on you (whether by court order, agreement or
!!otherwise) that contradict the conditions of this License, they do not
!!excuse you from the conditions of this License.  If you cannot convey a
!!covered work so as to satisfy simultaneously your obligations under this
!!License and any other pertinent obligations, then as a consequence you may
!!not convey it at all.  For example, if you agree to terms that obligate you
!!to collect a royalty for further conveying from those to whom you convey
!!the Program, the only way you could satisfy both those terms and this
!!License would be to refrain entirely from conveying the Program.
!!
!!\item Use with the GNU Affero General Public License.
!!
!!Notwithstanding any other provision of this License, you have
!!permission to link or combine any covered work with a work licensed
!!under version 3 of the GNU Affero General Public License into a single
!!combined work, and to convey the resulting work.  The terms of this
!!License will continue to apply to the part which is the covered work,
!!but the special requirements of the GNU Affero General Public License,
!!section 13, concerning interaction through a network will apply to the
!!combination as such.
!!
!!\item Revised Versions of this License.
!!
!!The Free Software Foundation may publish revised and/or new versions of
!!the GNU General Public License from time to time.  Such new versions will
!!be similar in spirit to the present version, but may differ in detail to
!!address new problems or concerns.
!!
!!Each version is given a distinguishing version number.  If the
!!Program specifies that a certain numbered version of the GNU General
!!Public License ``or any later version'' applies to it, you have the
!!option of following the terms and conditions either of that numbered
!!version or of any later version published by the Free Software
!!Foundation.  If the Program does not specify a version number of the
!!GNU General Public License, you may choose any version ever published
!!by the Free Software Foundation.
!!
!!If the Program specifies that a proxy can decide which future
!!versions of the GNU General Public License can be used, that proxy's
!!public statement of acceptance of a version permanently authorizes you
!!to choose that version for the Program.
!!
!!Later license versions may give you additional or different
!!permissions.  However, no additional obligations are imposed on any
!!author or copyright holder as a result of your choosing to follow a
!!later version.
!!
!!\item Disclaimer of Warranty.
!!
!!\begin{sloppypar}
!! THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
!! APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE
!! COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM ``AS IS''
!! WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
!! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
!! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE
!! RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
!! SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
!! NECESSARY SERVICING, REPAIR OR CORRECTION.
!!\end{sloppypar}
!!
!!\item Limitation of Liability.
!!
!! IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN
!! WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES
!! AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR
!! DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL
!! DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM
!! (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED
!! INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE
!! OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH
!! HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
!! DAMAGES.
!!
!!\item Interpretation of Sections 15 and 16.
!!
!!If the disclaimer of warranty and limitation of liability provided
!!above cannot be given local legal effect according to their terms,
!!reviewing courts shall apply local law that most closely approximates
!!an absolute waiver of all civil liability in connection with the
!!Program, unless a warranty or assumption of liability accompanies a
!!copy of the Program in return for a fee.
!!
!!\begin{center}
!!{\Large\sc End of Terms and Conditions}
!!
!!\bigskip
!!How to Apply These Terms to Your New Programs
!!\end{center}
!!
!!If you develop a new program, and you want it to be of the greatest
!!possible use to the public, the best way to achieve this is to make it
!!free software which everyone can redistribute and change under these terms.
!!
!!To do so, attach the following notices to the program.  It is safest
!!to attach them to the start of each source file to most effectively
!!state the exclusion of warranty; and each file should have at least
!!the ``copyright'' line and a pointer to where the full notice is found.
!!
!!{\footnotesize
!!\begin{verbatim}
!!<one line to give the program's name and a brief idea of what it does.>
!!
!!Copyright (C) <textyear>  <name of author>
!!
!!This program is free software: you can redistribute it and/or modify
!!it under the terms of the GNU General Public License as published by
!!the Free Software Foundation, either version 3 of the License, or
!!(at your option) any later version.
!!
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details.
!!
!!You should have received a copy of the GNU General Public License
!!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!\end{verbatim}
!!}
!!
!!Also add information on how to contact you by electronic and paper mail.
!!
!!If the program does terminal interaction, make it output a short
!!notice like this when it starts in an interactive mode:
!!
!!{\footnotesize
!!\begin{verbatim}
!!<program>  Copyright (C) <year>  <name of author>
!!
!!This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
!!This is free software, and you are welcome to redistribute it
!!under certain conditions; type `show c' for details.
!!\end{verbatim}
!!}
!!
!!The hypothetical commands {\tt show w} and {\tt show c} should show
!!the appropriate
!!parts of the General Public License.  Of course, your program's commands
!!might be different; for a GUI interface, you would use an ``about box''.
!!
!!You should also get your employer (if you work as a programmer) or
!!school, if any, to sign a ``copyright disclaimer'' for the program, if
!!necessary.  For more information on this, and how to apply and follow
!!the GNU GPL, see \texttt{http://www.gnu.org/licenses/}.
!!
!!The GNU General Public License does not permit incorporating your
!!program into proprietary programs.  If your program is a subroutine
!!library, you may consider it more useful to permit linking proprietary
!!applications with the library.  If this is what you want to do, use
!!the GNU Lesser General Public License instead of this License.  But
!!first, please read \newline\texttt{http://www.gnu.org/philosophy/why-not-lgpl.html}.
!!
!!\end{enumerate}

!(doc)footer
