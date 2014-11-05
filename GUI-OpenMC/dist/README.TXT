
GUI-OpenMC is an open source application written by jaafar el bakkali (bahmedj@gmail.com) which offers a GUI for OpenMC Monte Carlo particle transport simulation code, allowing creation/edition separate mandatory XML files describing the geometry, the materials involved and the settings of a given OpenMC simulation. The creation/edition facility is also allowed to optional xml files such tallies (results of given simulation), plotting and CMFD (Coarse-Mesh Finite Difference) acceleration. The GUI-OpenMC will help the OpenMC users to avoid hardcoding full content of those XML files when creating, editing these files.
In order to displaying the user application xml files in a more readable way, the syntax highlighter approach has been considered where the open source java library jsyntaxpane was used as highlighter for those xml files.
The java-based open source application was developed using NetBeans IDE 7.4 under Ubuntu 12 operating system. The GUI application calls three external programs: openmc which used as neutron calculations engine, xterm used as application console and eog (Eyes of GENOME) used as PPM image viewer.

Requirements:
Netbeans IDE 7.4
JDK 7

