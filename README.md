# INHDF5
A mini API for using (programming) the parallel HDF5 library in Fortran

SUMMARY

1. Quick ways to create and use files/groups/datasets in a distributed parallel
program written in Fortran.

2. Writes and reads of either a portion or the entire dataset can be done, and
both can be throttled by having collective write/read operations be wrapped
in sweeps (this can be used to tune efficiency of the underlying parallel
file system). Segmentation or throttling of writes/reads can be performed by
proper use of the "offset" and "dimension" specification inputs.

3. Writing and reading is done with the underlying Fortran2003 API of HDF5.
This is with good reason. The HDF5 library's Fortran API has "Fortran
interfaces" to call the appropriate function that specifies a particular
atomic datatype (integer, double, float, etc). However, HDF5 has Fortran 2003
bindings, which require the use of "ISO C bindings" in the Fortran code. Those
are useful because there is a single Fortran interface that wraps the HDF5 C
calls that do the actual wrok. (This is not as complicated as it sounds.) The
Fortran user needs only to specify the C-pointer of the begining of the buffer
of data that is to be written (or read), like so:
```
      Real*8 :: buffer(10)
      Type(C_PTR) :: ptr
      ptr = C_LOC( buffer)

      Integer*4 :: int_buffer(10)
      Type(C_PTR) :: ptr
      ptr = C_LOC( int_buffer)
```

This makes the API (and this library) agnostic to the data type that is to
be written or read.

4. This API is agnostic to datatypes of the data that is to be written or read.
This is related to (3) above. Basically, we can use the same API to write/read
any kind of data, including "Fortran derived types" (like C structs). One word
of caution of derived types: I have encountered errors during parallel writes
of chunks of derived type data to a LUSTRE filesystem.


USING THE LIBRARY

- Compiling

Compile just as you would a distributed parallel program with MPI and HDF5.


- Using this API

Look at what the code does (it is obvious) and just build on top of it or use
it as is. There is nothing to it... All these are shortcuts for the simple
operations that would otherwise need error-trapping to be repeated over and
over. I put the error-trapping there for you.


MOTIVATION

The things that can be kept generic, should be kept as generic as possible.


TO BE CONTINUED...


IN 2019/04/29


Added an example and a Makefile.

IN 2020/06/16
