# interactive-IB-NS-code
interactive immersed boundary incompressible Navier-Stokes solver 

![porsche](http://dutw1479.wbmt.tudelft.nl/~renep/images/porsche.png)

Required libraries: 

1.) openMP (if homebrew installed on MAC OS, simply use: brew install clang-omp)

2.) FFTW, http://www.fftw.org/download.html (use: ./configure --enable-float --enable-threads)

3.) OpenGL utility toolkit (glut) (use: brew install freeglut)

Adapt locations of libs in the Makefile (for OSX) and run "make".  

To change number of threads for openMP, change the setting in "defines.h". 
There are two defines: NUM_THREADS for the loops and NUM_THREADS_FFT for FFTW. 



Known issues / things to do: 

1.) IB not completely implemented for QUICK scheme

2.) for the sake of computational speed, pressure only solved after RK sub-steps



have fun...


