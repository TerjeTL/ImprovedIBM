#pragma once

#include "Shader.h"

#include "glm.hpp"
#include "Eigen/Eigen"

#include <vector>
#include <filesystem>

template<typename T>
std::vector<T> EigenMatrixToArray(const Eigen::MatrixX<T>& eigen_mat)
{
    std::vector<T> mat_vector(eigen_mat.rows() * eigen_mat.cols());

    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
    RowMajMat::Map(mat_vector.data(), eigen_mat.rows(), eigen_mat.cols()) = eigen_mat;

    return mat_vector;
}

class MatrixObject
{
public:
    MatrixObject(std::vector<int> grid_flags);
    ~MatrixObject()
    {
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
        glDeleteBuffers(1, &FBO);
        glDeleteBuffers(1, &texturebuffer);
    }
    void SetMVP(glm::mat4 mvp)
    {
        MVP = mvp;
    }

    void DrawToTexture();
    unsigned int GetTextureBuffer() const
    {
        return texturebuffer;
    }

    std::vector<int> grid_flag_mat;
    uint32_t size = 0;

private:
    unsigned int texture1, texture2;
    unsigned int VBO, VAO, EBO, FBO;
    unsigned int texturebuffer;

    static constexpr float vertices[] = {
        //pos                //texcoords
         0.5f,  0.5f, 0.0f,  1.0f, 1.0f, //tr
         0.5f, -0.5f, 0.0f,  1.0f, 0.0f, //br
        -0.5f, -0.5f, 0.0f,  0.0f, 0.0f, //bl
        -0.5f,  0.5f, 0.0f,  0.0f, 1.0f  //tl
    };
    static constexpr unsigned int indices[] = {
        0, 1, 3, // first triangle
        1, 2, 3  // second triangle
    };

    glm::mat4 MVP;
    glm::vec3 color;

    Eigen::Matrix4f test_tex = (Eigen::Matrix4f() <<
        0.0f, 0.05f, 0.1f, 0.15f,
        0.2f, 0.25f, 0.3f, 0.35f,
        0.4f, 0.45f, 0.5f, 0.55f,
        0.6f, 0.65f, 0.7f, 1.0f).finished();

    std::filesystem::path vert_path = std::filesystem::current_path().parent_path().parent_path() / "ibm_application/src/shaders/MatrixViewer.vert";
    std::filesystem::path frag_path = std::filesystem::current_path().parent_path().parent_path() / "ibm_application/src/shaders/MatrixViewer.frag";

    Shader shaders{ vert_path.string().c_str(), frag_path.string().c_str() };
};