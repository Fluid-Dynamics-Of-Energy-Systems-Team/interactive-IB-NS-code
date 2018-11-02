# interactive-IB-NS-code
interactive immersed boundary incmpressible Navier-Stokes code 

Required libraries: 

1.) openMP

2.) FFTW (use ./configure --enable-float --enable-threads)

3.) OpenGL utility toolkit (glut). If homebrew installed on MAC OS, simply use: "brew install freeglut"

Adapt locations of libs in the Makefile (for OSX) and make the code. 

To change number of threads for openMP, change the setting in the "defines.h". 
There are two defines: NUM_THREADS for the loops and NUM_THREADS_FFT for FFTW. 



Known issues / things to do: 

1.) IB not completely implemented for QUICK scheme

2.) for the sake of computational speed, pressure only solved after RK sub-steps



have fun...

Rene

