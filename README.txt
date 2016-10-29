    **************************************************************
                                POSEST
                              version 1.1
                 By Manolis Lourakis & Xenophon Zabulis

                     Institute of Computer Science
            Foundation for Research and Technology - Hellas
                       Heraklion, Crete, Greece
    **************************************************************


GENERAL
This is posest, a copylefted C/C++ library for 6DoF pose estimation. The
computation is based on a set of 3D points and their corresponding 2D
projections on an imaging sensor. The library estimates the relative motion
between the 3D points and the camera. See posest_demo.c for examples of use.

The theory behind posest is described in the paper entitled "Model-Based Pose
Estimation for Rigid Objects" by M. Lourakis and X. Zabulis, Springer LNCS 7963,
pp. 83-92, 2013.

Please cite this paper if your use posest in your published work.

LICENSE
posest is released under the GNU Public License (GPL), which can be found in
the included LICENSE file. Note that under the terms of GPL, commercial use
is allowed only if a software employing posest is also published in source
under the GPL. However, if you are interested in using posest in a proprietary
commercial application, a commercial license for posest can be obtained by
contacting the author using the email address at the end of this file.

COMPILATION
posest can be compiled with the CMake cross-platform build system. The included
CMakeLists.txt file can be used to generate makefiles for Unix systems or project
files for Windows systems.  CMakeLists.txt defines some configuration variables
that control certain aspects of posest and can be modified from CMake's user
interface.
More information on how to use CMake can be found at http://www.cmake.org

Prior to building posest, the levmar (http://www.ics.forth.gr/~lourakis/levmar)
optimization library should be installed.

The demo program posest_demo and binocposest_demo provide working examples
of using posest. Please refer to examples/README.txt for sample runs.

MATLAB INTERFACE
The 'matlab' subdirectory in posest's distribution includes a matlab mex
interface. See the accompanying 'README.txt' for more information and examples
of use.

Send your comments/bug reports to lourakis (at) ics (dot) forth (dot) gr
