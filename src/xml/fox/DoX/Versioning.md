# FoX versioning

This documentation describes version 4.1 of the FoX library.

This version includes output modules for general XML, and for CML; and a fully validating XML parser, exposed through a Fortran version of the SAX2 input parser and a Fortran mapping of the W3C DOM interface.

This is a stable branch, which will be maintained with important bugfixes.

<a name="Changes"/>

## FoX Changes

As of FoX-3.0, there is one user-visible change that should be noted.

### Configuration/compilation

In previous versions of FoX, the configure script was accessible as `config/configure`. Version 3.0 now follows common practice by placing the script in the main directory, so it is now called as `./configure`.

Previous versions of FoX made it quite hard to compile only portions of the library (eg only the CML output portion; or just the SAX input). This is now possible by specifying arguments to the configuration script. For example,

`./configure --enable-wcml`

will cause the generated Makefile to only compile the CML writing module and its dependencies.

See [Compilation](|Compilation|) for further details.


