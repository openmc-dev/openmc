# String conversion

Two procedures are provided to simplify reading data retreved from XML documents into Fortran variables. The subroutine `rts` performs the data conversion step and the function `countrts` can be used to allocate an array of the correct size for the incomming data. 

## `rts` subroutine

The `rts` subroutine can be imported from `FoX_common`. In its simplest form, it is called in this fashion:

    call rts(string, data)

`string` is a simple Fortran string (probably retrieved from an XML file.)

`data` is any native Fortran datatype: `logical`, `character`, `integer`, `real`, `double precision`, `complex`, `double complex`, and may be a scalar, 1D or 2D array.

`rts` will attempt to parse the contents of `string` into the appropriate datatype, and return the value in `data`.

Additional information or error handling is accomplished with the following optional arguments:

### `num`

`num` is an integer; on returning from the function it indicates the number of data items read before either:

* an error occurred
* the string was exhausted of data items  
* `data` was filled.

### `iostat`

`iostat` is an integer, which on return from the function has the values:

* `0` for no problems
* `-1` if too few elements were found in `string` to fill up `data`
* `1` if `data` was filled, but there were still data items left in `string`
* `2` if the characters found in `string` could not be converted to the appropriate type for `data`.

NB if `iostat` is not specified, and a non-zero value is returned, then the program will stop with an error message.

## String formatting

When `string` is expected to be an array of strings, the following options are used to break `string` into its constituent elements:

* By default it is assumed that the elements are separated by whitespace, and that multiple whitespace characters are not significant. No zero-length elements are possible, nor are elements containing whitespace.

* An optional argument, `separator` may be specified, which is a single character. In this case, each element consists of all characters between subsequent occurences of the `separator`. Zero-length elements are possible, but no escaping mechanism is possible.

* Alternatively, an optional logical argument `csv` may be specified. In this case, the value of `delimiter` is ignored, and the string is parsed as a Comma-Separated-Value string, according to [RFC 4180](http://tools.ietf.org/html/rfc4180).

## Numerical formatting.

Numbers are expected to be formatted according to the usual conventions for Fortran input.

## Complex number formatting.

Complex numbers may be formatted according to either normal Fortran conventions (comma-separated pairs) or [CMLComp conventions](http://cmlcomp.org/t/wiki/FpxStandard)

## Logical variable formatting.

Logical variables must be encoded according to the conventions of [XML Schema Datatypes](http://www.w3.org/TR/xmlschema-2/#boolean)  - that is, True may be written as "true" or "1", and False may be written as "false" or "0".

## `countrts` function

The `countrts` function can also be imported from `FoX_common`. In its simplest form, it is called in this fashion:

    countrts(string, datatype)

`string` is a simple Fortran string (probably retrived from an XML file)

`datatype` is a scalar argument of any native Fortran datatype (`logical`, `character`, `integer`, `real`, `double precision`, `complex` or `double complex`).

The function returns a default integer equal to the number of elements that rts would
return if called with a sufficently large array of the same type as `datatype`. `countrts` returns 0 to indicate that characters were found in the string that could not be converted. If datatype is a character, the optional arguments `seperator` and `csv` are avalable as described in "string formatting" above. The `countrts` function is pure and can be used as a specification function.

