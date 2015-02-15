# FoX_common

FoX_common is a module exporting interfaces to a set of convenience functions common to all of the FoX modules, which are of more general use.

Currently, there are three publically available functions and four subroutines:

*  The subroutine `str` converts primitive datatypes into strings in a consistent fashion, conformant with the expectations of XML processors.

It is fully described in [StringFormatting](|StringFormatting|)

* The subroutine `rts` performs the reverse function, taking a string (obtained from an XML document) and converts it into a primitive Fortran datatype.

* The function `countrts` examinies a string and determines the size of array requiered to hold all its data, once converted to a primitive Fortran datatype.

It is fully described in [StringConversion](|StringConversion|)

The final four procedures change the way that errors and warnings are handled when encounterd by any FoX modules. Using these procedures it is possible to convert non-fatal warnings and fatal errors to calls to the internal about routine. This generally has the effect of generating a stack trace or core dump of the program before temination. This is a global setting for all XML documents being manipulated. Two subroutines take a single logical argument to turn on (true) and off (false) the feature for warnings and errors respectivly:
    
* `FoX_set_fatal_warnings` for warnings 

* `FoX_set_fatal_errors` for errors
    
and two functions (without arguments) allow the state to be checked:

* `FoX_get_fatal_warnings` for warnings
    
* `FoX_get_fatal_errors` for errors
 
Both fatal warnings and errors are off by default. This corresponds to the previous behaviour. 
