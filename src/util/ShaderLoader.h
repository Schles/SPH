#ifndef LOADSHADER_H
#define LOADSHADER_H

#include "../../extern/glad/glad.h"
#include <GLFW/glfw3.h>

#include <string>
#include <fstream>


class ShaderLoader {
public:
    int static loadShader(std::string vs, std::string fs);

private:
    GLuint loadShaderFromFile(std::string path, GLenum shaderType);
    void printShaderLog( GLuint shader );
    void printProgramLog( GLuint program );
};

#endif
