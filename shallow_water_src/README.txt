Shallow Water Demo by Stephen Thompson
======================================

This archive contains the source code for my Shallow Water Demo.

For more information on the demo itself, and pre-built binaries,
please refer to http://www.solarflare.org.uk/shallow_water_demo/index.html.


Compiling
=========

A Visual Studio 2008 solution file is provided ("shallow_water.sln" in
the "msvc" directory). Compiling should be straightforward, although
you will need to have the latest DirectX SDK correctly installed on
your machine.

To run the compiled exe, make sure you set Working Directory to 
"$(SolutionDir)\..", otherwise the data files will not be found.


Roadmap
=======

The source consists of three projects:

1) Coercri -- This is a simple wrapper library around operating system
functions (graphics, sounds and so on) which I use in several of my
projects. Here I have included only the parts relevant to the shallow
water demo (which is mostly the graphics functionality and Guichan
interfacing in this case).

2) Guichan -- This is an open source GUI library which was used to
create the buttons and sliders at the right-hand side of the screen.

3) Shallow_Water -- This is the main project containing the shallow
water simulations. Briefly, the main files are as follows:

 * engine.cpp -- Main "engine" for the simulation, contains all the
   code that drives the GPU. The bulk of the code is found here.

 * kp07.hlsl -- Contains shaders for doing the numerical simulation of
   the shallow water equations on the GPU.

 * shallow_water.fx -- Contains shaders for creating the graphical
   appearance of the land and water surfaces.

 * gui_manager.cpp -- Creates the GUI at the right-hand side of the
   screen.

 * settings.cpp -- Stores and manages simulation settings.

 * terrain_heightfield.cpp -- Contains formulas for determining the
   shape of the terrain.


Copyright / Legal Information
=============================

The Shallow Water Demo source code consists of three components each
of which has their own copyright and licence conditions. The three
components are:

1) The Shallow Water Demo itself
2) Coercri
3) Guichan

Each source code file identifies clearly (at the top of the file)
which licence applies to that file.

The licences for the three components are as follows:

1) The Shallow Water Demo is Copyright (C) 2012 Stephen Thompson.

The Shallow Water Demo is free software: you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see GPL.txt). If not, see
<http://www.gnu.org/licenses/>.

2) Coercri is Copyright (C) 2011 Stephen Thompson.

Coercri is licensed under the Boost Software Licence
(http://www.boost.org/LICENSE_1_0.txt), the text of which can be found
in the relevant source code files.

3) Guichan is Copyright (C) 2004 - 2008 Olof Naessén and Per Larsson.

Guichan is licensed under a BSD licence; the exact licence terms can
be found in the relevant source code files.


Contact
=======

I can be contacted by email at stephen@solarflare.org.uk.
