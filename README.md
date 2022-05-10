# Curves #

### This repository contains a set of classes to create N-dimensional curves in C++ ###

<a name="toc"></a>
#### Table of contents: ####
- [Curves](#curves)
- [Instructions](#instructions)
- [Table of contents](#toc)
- [Code Overview](#overview)
- [TODOs](#todos)

<a name="curves"></a>
#### Curves: ####
1. B-spline
2. Bezier Curves
3. Catmull-Rom Splines
4. Hermite Splines

A special class for cubic hermite spline has also been created that might be deprecated in the future, as the generalized class of hermite spline can be used to create the cubic hermite spline.

<a name="instructions"></a>
#### To build and run: ####
    git clone https://github.com/vss2sn/curves.git  
    cd curves
    g++-10 main.cpp -std=c++2a  #(Might need to specify g++-10 or higher)
    ./a.out

<a name="overview"></a>
#### Code Overview: ####
* Currently this is a header-only repository containing the classes for the curves listed above, implemented using `std::array`
* This allows for compile time checking of degrees, number of points, etc., and ensures that most mismatches are caught at compile time, and should be easier to understand and debug
* The constructors for each of the classes are `constexpr` allowing them to be created and used at compile time (this has been tested using the `std::is_constant_evaluated`)
* The main file has multiple examples of how to use each of the classes and the util and curve conversion functions.

<a name="todos"></a>
#### TO DOs: ####
* Add documentation
* Add references
* Add implementations using vectors
* Add tests
* Add knot insertion algorithm for converting `bspline` to group of `bezier` curves
* Add functions to transform one type of curve into another
* Add (python?) plotting utility
* Move to standard C++ project file arrangement and format
