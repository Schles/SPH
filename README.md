# Startparameter
./FluidSim [Szene] [ParticlePerAxis]

Bei mir läuft es mit 12 Partikeln pro Achse noch in realtime

# Szenen
2 - Breaking Dam
3 - Double Breaking Dam
4 - Moving Wall


# Tasten

p - Start / Stop Simulation
r - Restart Simulation
g - Gravitation
t - Adaptive Stepsize
SPC - Step once

# Probleme
Adaptive Stepsize kann manchmal funktionieren. Bei sehr dichten Stellen, tendiert es allerdings instabil zu werden


# Hinweise

## Ubuntu dependencies
Ich vermute, folgende Pakete müssen installiert werden, damit das Projekt auf Ubuntu (16.04) kompiliert.

    sudo apt-get install freeglut3-dev dh-autoreconf libasound2-dev libgl1-mesa-dev libpulse-dev libudev-dev libdbus-1-dev libx11-dev libxcursor-dev libxext-dev libxi-dev libxinerama-dev libxrandr-dev libxss-dev libxt-dev libxxf86vm-dev 	

## CMake build

Das Projekt muss out-of-source (?!) kompiliert werden. Falls ihr nicht die cmake-gui verwendet wollt, habe ich ein cmake.sh Skript hinterlegt. Das setzt einfach nur einen (externen) Build Ordner setzt (out-of-source anstatt in-source ?!).



# GLviz

GLviz is a small collection of C++ classes and GLSL shaders to facilitate the development of OpenGL demos. It is built on top of [CMake](http://www.cmake.org/), [SDL](http://libsdl.org), [GLEW](http://glew.sourceforge.net), [Eigen](http://eigen.tuxfamily.org/), and [AntTweakBar](http://anttweakbar.sourceforge.net/) and requires at least OpenGL 3.3. GLviz has been tested on a NVIDIA GTX 680 GPU using driver version 347.52 on Windows (compiled with MSVC 2013) and 340.65 on Linux (compiled with GCC 4.9). It includes all external dependencies to make compilation on either Windows or Linux as simple and convenient as possible.

**Author**: [Sebastian Lipponer](http://sebastianlipponer.de), **License**: MIT

## Features

* Camera with trackball navigation.
* Shader management.
* Triangle mesh shader (supports flat or Phong shading and an optional high-quality wireframe<sup>1</sup> overlay).
* Sphere shader (for the sake of performance not perspectively correct spheres).
* Supports embedding of shaders in the executable.

## Screenshots

[![](http://sebastianlipponer.github.io/glviz/dragon_thumbnail.png)](http://sebastianlipponer.github.io/glviz/dragon.png)[![](http://sebastianlipponer.github.io/glviz/dragon_wireframe_thumbnail.png)](http://sebastianlipponer.github.io/glviz/dragon_wireframe.png)[![](http://sebastianlipponer.github.io/glviz/dragon_spheres_thumbnail.png)](http://sebastianlipponer.github.io/glviz/dragon_spheres.png)

## Documentation

Currently there is no documentation, but the code includes an example. From the source code it should become apparent how GLviz is intended to be used.

## References

[1] Bærentzen, J. A., Munk-Lund, S., Gjøl, M., Larsen, B. D.: **Two Methods for Antialiased Wireframe Drawing with Hidden Line Removal**. In Proceedings of the Spring Conference in Computer Graphics, 2008.
