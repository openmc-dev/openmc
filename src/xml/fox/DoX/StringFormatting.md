# String handling in FoX

Many of the routines in wxml, and indeed in wcml which is built on top of wxml, are overloaded so that data may be passed to the same routine as string, integer, logical, real, or complex data.

In such cases, a few notes on the conversion of non-textual data to text is in order. The
standard Fortran I/O formatting routines do not offer the control required for useful XML output, so FoX performs all its own formatting.

This formatting is done internally through a function which is also available publically to the user, `str`.

To use this in your program, import it via:

    use FoX_common, only; str

and use it like so:

     print*, str(data)

In addition, for ease of use, the `//` concatenation operator is overloaded, such that strings can easily be formed by concatenation of strings to other datatypes. To use this you must import it via:

     use FoX_common, only: operator(//)

and use it like so:

     integer :: data
     print*, "This is a number "//data

This will work for all native Fortran data types - but no floating point formatting is available as described below with concatenation, only with str()

You may pass data of the following primitive types to `str`:

## Scalar data

### Character (default kind)

Character data is returned unchanged.

### Logical (default kind)

Logical data is output such that True values are converted to the string 'true', and False to the string 'false'.

### Integer (default kind)

Integer data is converted to the standard decimal representation.

### Real numbers (single and double precision)

Real numbers, both single and double precision, are converted to strings in one of two ways, with some control offered to the user. The output will conform to the real number formats specified by XML Schema Datatypes.

This may be done in one of two ways:

 1. Exponential notation, with variable number of significant figures. Format strings of the form "`s`**n**"  are accepted, where **n** is the number of significant figures.

 Thus the number `111`, when output with various formats, will produce the following output:

<table class="format">
<tr>
  <td class="format"> s1 </td><td> 1e2 </td>
</tr><tr>
  <td> s2 </td><td> 1.1e2 </td>
</tr><tr>
  <td> s3 </td><td> 1.11e2 </td>
</tr><tr>
  <td> s4 </td><td> 1.110e2 </td>
</tr>
</table>

 The number of significant figures should lie between 1 and the number of digits precision provided by the real kind. If a larger or smaller number is specified, output will be truncated accordingly. If unspecified, then a sensible default will be chosen.

  This format is not permitted by XML Schema Datatypes 1.0, though it is in 2.0

 2. Non-exponential notation, with variable number of digits after the decimal point. Format strings of the form "`r`**n**", where **n** is the number of digits after the decimal point.

 Thus the number `3.14159`, when output with various formats, will produce the following output:

<table class="format">
<tr>
  <td> r0 </td><td> 3 </td>
</tr><tr>
  <td> r1 </td><td> 3.1</td>
</tr><tr>
  <td> r2 </td><td> 3.14</td>
</tr><tr>
  <td> r3 </td><td> 3.142 </td>
</tr>
</table>

 The number of decimal places must lie between 0 and whatever would output the maximum digits precision for that real kind.  If a larger or smaller number is specified, output will be truncated accorsingly. If unspecified, then a sensible default will be chosen.

 This format is the only one permitted by XML Schema Datatypes 1.0

 If no format is specified, then a default of exponential notation will be used.

 If a format is specified not conforming to either of the two forms above, a run-time error will be generated.

**NB** Since by using FoX or str, you are passing real numbers through various functions, this means that
       they must be valid real numbers. A corollary of this is that if you pass in +/-Infinity, or NaN, then
       the behaviour of FoX is unpredictable, and may well result in a crash. This is a consequence of the
       Fortran standard, which strictly disallows doing anything at all with such numbers, including even
       just passing them to a subroutine.

## Complex numbers (single and double precision)

Complex numbers will be output as pairs of real numbers, in the following way:

`(1.0e0)+i(1.0e0)`

where the two halves can be formatted in the way described for 'Real numbers' above; only one format may be specified, and it will apply to both.

All the caveats described above apply for complex number as well; that is, output of complex numbers either of whose components are infinite or NaN is illegal in Fortran, and more than likely will cause a crash in FoX.

## Arrays and matrices

All of the above types of data may be passed in as arrays and matrices as well. In this case, a string containing all the individual elements will be returned, ordered as they would be in memory, each element separated by a single space.

If the data is character data, then there is an additional option to str, `delimiter` which may be any single-character string, and will replace a space as the delimiter.

## wxml/wcml wrappers.

All functions in wxml which can accept arbitrary data (roughly, wherever you put anything that is not an XML name; attribute values, pseudo-attribute values, character data) will take scalars, arrays, and matrices of any of the above data types, with `fmt=` and `delimiter=` optional arguments where appropriate.

Similarly, wcml functions which can accept varied data will behave similarly.
