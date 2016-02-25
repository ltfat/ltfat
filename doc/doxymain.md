\mainpage

libltfat - Large Time-Frequency Alalysis Toolbox backend library
================================================================

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


Compatibility
-------------

The library internally uses complex numbers from ISO C99
[complex.h](http://en.cppreference.com/w/c/numeric/complex) and
[aligned_malloc](http://en.cppreference.com/w/c/memory/aligned_alloc)
 defined in ISO C11, therefore it must be compiled with a
fairly modern compiler.

When linking the library however, none of the above is required. It is not
required to supply memory aligned output arrays and the complex data type
in the headers will become a length 2 array.

Such format is binary compatible with the complex class from C++.
One can cast pointers back and forth in the following way:
~~~~~~~~~~~~~~~{.cpp}
double ccomp[][2] = {{1.0,2.0},{3.0,4.0},{5.0,6.0}};
std::complex<double>* ccpp = reinterpret_cast<std::complex<double>*>(ccomp);
double (*ccomp2)[2] = reinterpret_cast<double(*)[2]>(ccpp);
~~~~~~~~~~~~~~~

A complete example:


> Blockquote
> 2nd line

- More text for this item.
- Further
  + nested list item.
  + another nested item.
- Prd

1. First item.
2. Second item.

This a normal paragraph

    This is a code block

We continue with a normal paragraph again.


Code inline `int prd`

Cross link to \ref stft "(formula)"



