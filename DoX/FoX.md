#FoX documentation.

[All in one page](FoX_DoX.html)

[Separate pages](FoX.html)

## Introduction

This document is the primary documentation for FoX, the Fortan/XML library. See below for [other sources of documentation](#otherdoc). It consists of:

Reference information on [versions](|Versioning|), 
[standards compliance](|Standards|),
and [licensing](|Licensing|).

Information about how to [get up and running with FoX](|Compilation|)
and how to [use FoX in an existing project](|Embedding|).

Finally, there is [full API reference documentation](#apidoc).

## Other documentation

<a name="otherdoc"/>

This documentation is largely reference in nature. For new users it is best to start elsewhere:

### iFaX workshops

Two workshops, entitled iFaX (Integrating Fortran and XML) have been run teaching the use of FoX, [one in January 2007](http://www.niees.ac.uk/events/ifax/index.shtml), and [one in January 2008](http://www.nesc.ac.uk/esi/events/841/). The full documentation and lectures from these may be found at:

* [iFaX I](http://buffalo.niees.group.cam.ac.uk/archive2.php?event_details=ifax)
* [iFaX II](http://www.nesc.ac.uk/action/esi/contribution.cfm?Title=841)

### Tutorials

Out of the above workshops, some [tutorial material](http://www1.gly.bris.ac.uk/~walker/FoX/iFaX) has been written, focussing on different use cases. Currently two are available:

* [SAX input](http://www1.gly.bris.ac.uk/~walker/FoX/iFaX/iFaX.4/iFaX.4.html)
* [DOM input](http://www1.gly.bris.ac.uk/~walker/FoX/iFaX/iFaX.5/iFaX.5.html)

There is also tutorial information on the use of WKML [here](http://web.me.com/dove_family/xml/kml.html).

## API documentation

<a name="apidoc"/>

* FoX has seven sets of publically exported interfaces. These are documented here:

### COMMON interfaces

* [FoX_common](|FoX_common|)
* [FoX_utils](|FoX_utils|)

### OUTPUT interfaces

* [FoX_wxml](|FoX_wxml|)
* [FoX_wcml](|FoX_wcml|)
* [FoX_wkml](|FoX_wkml|)

### INPUT interface

* [FoX_sax](|FoX_sax|)
* [FoX_dom](|FoX_dom|)

These documents describe all publically usable APIs.

Worked examples of the use of some of these APIs may be found in the `examples/` subdirectory, and tutorial-style documentaion is available from the links [above](#otherdoc).

## Other things

* [Hints for debugging](|Debugging|)

* [Further information](|Information|)

