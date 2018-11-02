# interactive-IB-NS-code
interactive immersed boundary Navier-Stokes code 

In order to install and run the code, install: 
1.) openMP
2.) FFTW
3.) OpenGL utility toolkit (glut). For MAC os simply use homebrew: "brew install freeglut"

adapt the Makefile (for OSX) accordinlgy and install. 

To change number of threads for openMP change the setting in the file "defines.h". 
There are two defines: NUM_THREADS for the loops and NUM_THREADS_FFT for FFTW. 

have fun...


