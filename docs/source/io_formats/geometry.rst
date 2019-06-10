.. _io_geometry:

======================================
Geometry Specification -- geometry.xml
======================================

.. _surface_element:

---------------------
``<surface>`` Element
---------------------

Each ``<surface>`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the surface.

    *Default*: None

  :name:
    An optional string name to identify the surface in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :type:
    The type of the surfaces. This can be "x-plane", "y-plane", "z-plane",
    "plane", "x-cylinder", "y-cylinder", "z-cylinder", "sphere", "x-cone",
    "y-cone", "z-cone", or "quadric".

    *Default*: None

  :coeffs:
    The corresponding coefficients for the given type of surface. See below for
    a list a what coefficients to specify for a given surface

    *Default*: None

  :boundary:
     The boundary condition for the surface. This can be "transmission",
     "vacuum", "reflective", or "periodic". Periodic boundary conditions can
     only be applied to x-, y-, and z-planes. Only axis-aligned periodicity is
     supported, i.e., x-planes can only be paired with x-planes. Specify which
     planes are periodic and the code will automatically identify which planes
     are paired together.

    *Default*: "transmission"

  :periodic_surface_id:
     If a periodic boundary condition is applied, this attribute identifies the
     ``id`` of the corresponding periodic sufrace.

The following quadratic surfaces can be modeled:

  :x-plane:
    A plane perpendicular to the x axis, i.e. a surface of the form :math:`x -
    x_0 = 0`. The coefficients specified are ":math:`x_0`".

  :y-plane:
    A plane perpendicular to the y axis, i.e. a surface of the form :math:`y -
    y_0 = 0`. The coefficients specified are ":math:`y_0`".

  :z-plane:
    A plane perpendicular to the z axis, i.e. a surface of the form :math:`z -
    z_0 = 0`. The coefficients specified are ":math:`z_0`".

  :plane:
    An arbitrary plane of the form :math:`Ax + By + Cz = D`. The coefficients
    specified are ":math:`A \: B \: C \: D`".

  :x-cylinder:
    An infinite cylinder whose length is parallel to the x-axis. This is a
    quadratic surface of the form :math:`(y - y_0)^2 + (z - z_0)^2 = R^2`. The
    coefficients specified are ":math:`y_0 \: z_0 \: R`".

  :y-cylinder:
    An infinite cylinder whose length is parallel to the y-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (z - z_0)^2 = R^2`. The
    coefficients specified are ":math:`x_0 \: z_0 \: R`".

  :z-cylinder:
    An infinite cylinder whose length is parallel to the z-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (y - y_0)^2 = R^2`. The
    coefficients specified are ":math:`x_0 \: y_0 \: R`".

  :sphere:
    A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 =
    R^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0 \: R`".

  :x-cone:
    A cone parallel to the x-axis of the form :math:`(y - y_0)^2 + (z - z_0)^2 =
    R^2 (x - x_0)^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0
    \: R^2`".

  :y-cone:
    A cone parallel to the y-axis of the form :math:`(x - x_0)^2 + (z - z_0)^2 =
    R^2 (y - y_0)^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0
    \: R^2`".

  :z-cone:
    A cone parallel to the x-axis of the form :math:`(x - x_0)^2 + (y - y_0)^2 =
    R^2 (z - z_0)^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0
    \: R^2`".

  :quadric:
     A general quadric surface of the form :math:`Ax^2 + By^2 + Cz^2 + Dxy +
     Eyz + Fxz + Gx + Hy + Jz + K = 0` The coefficients specified are ":math:`A
     \: B \: C \: D \: E \: F \: G \: H \: J \: K`".

.. _cell_element:

------------------
``<cell>`` Element
------------------

Each ``<cell>`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the cell.

    *Default*: None

  :name:
    An optional string name to identify the cell in summary output files.
    This string is limmited to 52 characters for formatting purposes.

    *Default*: ""

  :universe:
    The ``id`` of the universe that this cell is contained in.

    *Default*: 0

  :fill:
    The ``id`` of the universe that fills this cell.

    .. note:: If a fill is specified, no material should be given.

    *Default*: None

  :material:
    The ``id`` of the material that this cell contains. If the cell should
    contain no material, this can also be set to "void". A list of materials
    can be specified for the "distributed material" feature. This will give each
    unique instance of the cell its own material.

    .. note:: If a material is specified, no fill should be given.

    *Default*: None

  :region:
    A Boolean expression of half-spaces that defines the spatial region which
    the cell occupies. Each half-space is identified by the unique ID of the
    surface prefixed by `-` or `+` to indicate that it is the negative or
    positive half-space, respectively. The `+` sign for a positive half-space
    can be omitted. Valid Boolean operators are parentheses, union `|`,
    complement `~`, and intersection. Intersection is implicit and indicated by
    the presence of whitespace. The order of operator precedence is parentheses,
    complement, intersection, and then union.

    As an example, the following code gives a cell that is the union of the
    negative half-space of surface 3 and the complement of the intersection of
    the positive half-space of surface 5 and the negative half-space of surface
    2:

    .. code-block:: xml

        <cell id="1" material="1" region="-3 | ~(5 -2)" />

    .. note:: The ``region`` attribute/element can be omitted to make a cell
              fill its entire universe.

    *Default*: A region filling all space.

  :temperature:
    The temperature of the cell in Kelvin. If windowed-multipole data is
    avalable, this temperature will be used to Doppler broaden some cross
    sections in the resolved resonance region. A list of temperatures can be
    specified for the "distributed temperature" feature. This will give each
    unique instance of the cell its own temperature.

    *Default*: If a material default temperature is supplied, it is used. In the
    absence of a material default temperature, the :ref:`global default
    temperature <temperature_default>` is used.

  :rotation:
    If the cell is filled with a universe, this element specifies the angles in
    degrees about the x, y, and z axes that the filled universe should be
    rotated. Should be given as three real numbers. For example, if you wanted
    to rotate the filled universe by 90 degrees about the z-axis, the cell
    element would look something like:

    .. code-block:: xml

        <cell fill="..." rotation="0 0 90" />

    The rotation applied is an intrinsic rotation whose Tait-Bryan angles are
    given as those specified about the x, y, and z axes respectively. That is to
    say, if the angles are :math:`(\phi, \theta, \psi)`, then the rotation
    matrix applied is :math:`R_z(\psi) R_y(\theta) R_x(\phi)` or

    .. math::

       \left [ \begin{array}{ccc} \cos\theta \cos\psi & -\cos\phi \sin\psi +
       \sin\phi \sin\theta \cos\psi & \sin\phi \sin\psi + \cos\phi \sin\theta
       \cos\psi \\ \cos\theta \sin\psi & \cos\phi \cos\psi + \sin\phi \sin\theta
       \sin\psi & -\sin\phi \cos\psi + \cos\phi \sin\theta \sin\psi \\
       -\sin\theta & \sin\phi \cos\theta & \cos\phi \cos\theta \end{array}
       \right ]

    *Default*: None

  :translation:
    If the cell is filled with a universe, this element specifies a vector that
    is used to translate (shift) the universe. Should be given as three real
    numbers.

    .. note:: Any translation operation is applied after a rotation, if also
              specified.

    *Default*: None


---------------------
``<lattice>`` Element
---------------------

The ``<lattice>`` can be used to represent repeating structures (e.g. fuel pins
in an assembly) or other geometry which fits onto a rectilinear grid. Each cell
within the lattice is filled with a specified universe. A ``<lattice>`` accepts
the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the lattice.

  :name:
    An optional string name to identify the lattice in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :dimension:
    Two or three integers representing the number of lattice cells in the x- and
    y- (and z-) directions, respectively.

    *Default*: None

  :lower_left:
    The coordinates of the lower-left corner of the lattice. If the lattice is
    two-dimensional, only the x- and y-coordinates are specified.

    *Default*: None

  :pitch:
    If the lattice is 3D, then three real numbers that express the distance
    between the centers of lattice cells in the x-, y-, and z- directions.  If
    the lattice is 2D, then omit the third value.

    *Default*: None

  :outer:
    The unique integer identifier of a universe that will be used to fill all
    space outside of the lattice.  The universe will be tiled repeatedly as if
    it were placed in a lattice of infinite size.  This element is optional.

    *Default*: An error will be raised if a particle leaves a lattice with no
    outer universe.

  :universes:
    A list of the universe numbers that fill each cell of the lattice.

    *Default*: None

Here is an example of a properly defined 2d rectangular lattice:

.. code-block:: xml

    <lattice id="10" dimension="3 3" outer="1">
        <lower_left> -1.5 -1.5 </lower_left>
        <pitch> 1.0 1.0 </pitch>
        <universes>
          2 2 2
          2 1 2
          2 2 2
        </universes>
    </lattice>

-------------------------
``<hex_lattice>`` Element
-------------------------

The ``<hex_lattice>`` can be used to represent repeating structures (e.g. fuel
pins in an assembly) or other geometry which naturally fits onto a hexagonal
grid or hexagonal prism grid. Each cell within the lattice is filled with a
specified universe. This lattice uses the "flat-topped hexagon" scheme where two
of the six edges are perpendicular to the y-axis.  A ``<hex_lattice>`` accepts
the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the lattice.

  :name:
    An optional string name to identify the hex_lattice in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :n_rings:
    An integer representing the number of radial ring positions in the xy-plane.
    Note that this number includes the degenerate center ring which only has one
    element.

    *Default*: None

  :n_axial:
    An integer representing the number of positions along the z-axis.  This
    element is optional.

    *Default*: None

  :orientation:
    The orientation of the hexagonal lattice. The string "x" indicates that two
    sides of the lattice are parallel to the x-axis, whereas the string "y"
    indicates that two sides are parallel to the y-axis.

    *Default*: "y"

  :center:
    The coordinates of the center of the lattice. If the lattice does not have
    axial sections then only the x- and y-coordinates are specified.

    *Default*: None

  :pitch:
    If the lattice is 3D, then two real numbers that express the distance
    between the centers of lattice cells in the xy-plane and along the z-axis,
    respectively.  If the lattice is 2D, then omit the second value.

    *Default*: None

  :outer:
    The unique integer identifier of a universe that will be used to fill all
    space outside of the lattice.  The universe will be tiled repeatedly as if
    it were placed in a lattice of infinite size.  This element is optional.

    *Default*: An error will be raised if a particle leaves a lattice with no
    outer universe.

  :universes:
    A list of the universe numbers that fill each cell of the lattice.

    *Default*: None

Here is an example of a properly defined 2d hexagonal lattice:

.. code-block:: xml

    <hex_lattice id="10" n_rings="3" outer="1">
        <center> 0.0 0.0 </center>
        <pitch> 1.0 </pitch>
        <universes>
                  202
               202   202
            202   202   202
               202   202
            202   101   202
               202   202
            202   202   202
               202   202
                  202
        </universes>
    </hex_lattice>
