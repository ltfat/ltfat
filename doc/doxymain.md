\mainpage libltfat - Large Time-Frequency Alalysis Toolbox backend library

General conventions
-------------------

\anchor stft
 \f[
(\mathcal{V}_g f) (\sfreq, \stime)
    = \int_{\mathbb{R}}\! f(\stimearg+\stime)
    \overline{g(\stimearg)} \me^{-\mi 2\pi \sfreq \stimearg  } \,
    \mathrm{d}\stimearg,
    \ \ \sfreq, \stime \in\mathbb{R}\\
\f]

Papers which use this toolbox
\cite ltfatnote014

Function naming convention
--------------------------

The function names are in the following format:

\c ltfat_<function_name>[_<d|s|dc|sc>](<parameters>)

The \c ltfat_ prefix is present in all function names while the suffix
is optional and identifies the data type the function is working with.

<table>
<caption id="multi_row">Data type suffix</caption>
<tr><th>Suffix</th><th>Data type</th></tr>
<tr><td>d</td><td>double</td></tr>
<tr><td>s</td><td>float</td></tr>
<tr><td>dc</td><td>complex double</td></tr>
<tr><td>sc</td><td>complex float</td></tr>
</table>

\note In the documentation the prefix and the suffix will be omitted when
introducing a non-unique function and when referring to the function group.


Arrays, matrices and naming conventions
---------------------------------------

The multidimensional arrays are contiguous in memory and therefore, they
are reffered to by a single pointer and the individual dimensions are
accessed trough an offset.

When an array represents a matrix, it is assumed that the columns are the
first dimension and therefore they are stored continuously in memory.

In the function headers, the data arrays are denoted with array brackets []
and a pointer is used whenever reffering to a single object. This distinction
is just cosmetics as the arrays decay to pointers anyway.

Further, the following naming conventions are used consistently:
<table>
<caption id="multi_row">Argument naming</caption>
<tr><th>Arg. name</th><th>Definition</th></tr>
<tr><td>f</td><td>Time domain signal</td></tr>
<tr><td>g,gd,gt</td><td>Window, canotical dual window, canonical tight window </td></tr>
<tr><td>c</td><td>Coefficients</td></tr>
<tr><td>a</td><td>Integer. Time hop size.</td></tr>
<tr><td>M</td><td>Integer. Number of frequency channels (FFT length).</td></tr>
<tr><td>M2</td><td>Integer. Number of unique frequency channels for real signals:
M2=M/2+1.</td></tr>
<tr><td>L</td><td>Integer. Length of the signal. </td></tr>
<tr><td>N</td><td>Integer. Nuber of time shifts: N=L/a.</td></tr>
<tr><td>W</td><td>Integer. Number of signal channels.</td></tr>
</table>


Error handling
--------------

Every function which returns a status code should be checked by the user.
Additionally, the error message is printed to the standard error stream.
This behavior can be turned off or a custom error handler can be registered.
For details see \ref error

Plans
-----

When repeated computations with the same settings are desired, it is
convenient to create a __plan__ using the appropriate *_init function,
call the *_execute function multiple times and destroy the plan by
calling the *_done function.
The plan usually contains some precomputed read-only data,
working arrays and FFTW plans.
The plan is represented as a pointer to an opaque structure and here
is an example how to use it: 
~~~~~~~~~~~~~~~{.c}
dgt_long_plan_d* plan = NULL;

dgt_long_init_d(f, g, L, W, a, M, cout, ptype, FFTW_ESTIMATE, &plan)
// Fill in c and f after calling init. FFTW plan might overwrite it.

dgt_long_execute_d(plan);
// Refresh data in f and call execute again  

dgt_long_done_d(&plan);
~~~~~~~~~~~~~~~

\note Please note that due to the
<a href="https://github.com/FFTW/fftw3/issues/16">limitation of FFTW</a>
the *_init routines are not re-entrant because of the FFTW planning happening in them.
Therefore, the *_init functions cannot be called simultaneously on different threads even 
when creating completely separate plans.

\note Further, the *_execute function is reentrant and thread-safe, but not when executed 
on a single plan concurrently on separate threads. 
This limitation comes from the fact that the plan contains some working buffers.

States
------

A __state__ is a plan which additionally holds some data which persists
between the execute calls.

Compatibility
-------------

The library internally uses complex numbers from ISO C99
[complex.h](http://en.cppreference.com/w/c/numeric/complex).
Their memory layout is interleved and they are
binary compatible with the complex class from C++.
One can cast pointers back and forth in the following way:

~~~~~~~~~~~~~~~{.cpp}
double ccomp[][2] = {{1.0,2.0},{3.0,4.0},{5.0,6.0}};
std::complex<double>* ccpp = reinterpret_cast<std::complex<double>*>(ccomp);
double (*ccomp2)[2] = reinterpret_cast<double(*)[2]>(ccpp);
~~~~~~~~~~~~~~~






