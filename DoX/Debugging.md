#Debugging with FoX.

Following experience integrating `FoX_wxml` into several codes, here are a few tips for debugging any problems you may encounter.


## Compilation problems

You may encounter problems at the compiling or linking stage, with error messages along the lines of:
     'No Specific Function can be found for this Generic Function'
(exact phrasing depending on compiler, of course.)

If this is the case, it is possible that you have accidentally got the arguments to the offending out of order. If so, then use the keyword form of the argument to ensure correctness; that is, instead of doing:

    call cmlAddProperty(file, name, value)

do:

    call cmlAddProperty(xf=file, name=name, value=value)

This will prevent argument mismatches, and is recommended practise in any case.

## Runtime problems

You may encounter run-time issues. FoX performs many run-time checks to ensure the validity of the resultant XML code. In so far as it is possible, FoX will either issue warnings about potential problems, or try and safely handle any errors it encounters. In both cases, warning will be output on stderr, which will hopefully help diagnose the problem.

Sometimes, however, FoX will encounter a problem it can do nothing about, and must stop. In all cases, it will try and write out an error message highlighting the reason, and generate a backtrace pointing to the offending line. Occasionally though, the compiler will not generate this information, and the error message will be lost.

If this is the case, you can either investigate the coredump to find the problem, or (if you are on a Mac) look in ~/Library/Logs/CrashReporter to find a human-readable log.

If this is not enlightening, or you cannot find the problem, then some of the most common issues we have encountered are listed below. Many of them are general Fortran problems, but sometimes are not easily spotted in the context of FoX.

### Incorrect formatting.

Make sure, whenever you are writing out a real number through one of FoX's routines, and specifying a format, that the format is correct according to [StringFormatting](|StringFormatting|). Fortran-style formats are **not** permitted, and will cause crashes at runtime.

### Array overruns

If you are outputting arrays or matrices, and are doing so in the traditional Fortran style - by passing both the array and its length to the routine, like so:

     call xml_AddAttribute(xf=file, name=name, value=array, nvalue=n)

then if `n` is wrong, you may end up with an array overrun, and cause a crash.

We highly recommend wherever possible using the Fortran-90 style, like so:

     call xml_AddAttribute(xf=file, name=name, value=array)

where the array length will be passed automatically.

### Uninitialized variables

If you are passing variables to FoX which have not been initialized, you may well cause a crash. This is especially true, and easy to cause if you are passing in an array which (due to a bug elsewhere) has been partly but not entirely initialized. To diagnose this, try printing out suspect variables just before passing them to FoX, and look for suspiciously wrong values.

### Invalid floating point numbers.

If during the course of your calculation you accidentally generate Infinities, or NaNs, then passing them to any Fortran subroutine can result in a crash - therefore trying to pass them to FoX for output may result in a crash.

If you suspect this is happening, try printing out suspect variables before calling FoX. 
