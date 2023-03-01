#pragma once

#include "Mesh.h"
#include "Shader.h"
#include "ibm_application/DataStructs.h"

// holds
// - mesh to draw
// - shader to use
// - solution and logic controlling what to draw

class SolutionModel
{
public:
	SolutionModel() {}

	void SetSolution(std::shared_ptr<Solution> solution)
	{
		m_solution = solution;
	}

	void SetMVP(const std::tuple<glm::mat4, glm::mat4, glm::mat4>& mvp_in)
	{
		m_mesh.SetMVP(mvp_in);
	}

	void DrawSolutionModelToTexture()
	{
		shaders.use();
		shaders.SetMat4fv("mvp", m_mesh.mvp);
		m_mesh.DrawToTexture(shaders);
	}

	void AddIntegerDataTexture(const Eigen::MatrixXi& eigen_mat)
	{
		m_mesh.textures.push_back( Texture{ TextureFromMatrixXi(eigen_mat), "data_texture" } );
	}

	unsigned int& GetTextureBuffer()
	{
		return m_mesh.out_texture_buffer;
	}

	Mesh& GetMesh()
	{
		return m_mesh;
	}

private:
	template<typename T>
	std::pair<std::vector<T>, int> EigenMatToVector(const Eigen::MatrixX<T>& eigen_mat)
	{
		std::vector<T> mat_vector(eigen_mat.rows() * eigen_mat.cols());

		typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajMat;
		RowMajMat::Map(mat_vector.data(), eigen_mat.rows(), eigen_mat.cols()) = eigen_mat;

		return { mat_vector,  eigen_mat.rows() };
	}
	
	unsigned int TextureFromMatrixXi(const Eigen::MatrixX<int>& eigen_mat)
	{
		auto mat = EigenMatToVector<int>(eigen_mat);

		unsigned int texture_id;
		glGenTextures(1, &texture_id);

		int* data = mat.first.data();
		if (data)
		{
			glBindTexture(GL_TEXTURE_2D, texture_id);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R16I, mat.second, mat.second, 0, GL_RED_INTEGER, GL_INT, data);
			glGenerateMipmap(GL_TEXTURE_2D);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);	// set texture wrapping to GL_REPEAT (default wrapping method)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
			float borderColor[] = { 1.0f, 1.0f, 0.0f, 1.0f };
			glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);
			// set texture filtering parameters
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		}
		else
		{
			std::cout << "Could not load texture" << std::endl;
		}

		return texture_id;
	}

	// define the square
	std::vector<Vertex> verts{
		{ { 0.5f,  0.5f, 0.0f}, {1.0f, 1.0f} },
		{ { 0.5f, -0.5f, 0.0f}, {1.0f, 0.0f} },
		{ {-0.5f, -0.5f, 0.0f}, {0.0f, 0.0f} },
		{ {-0.5f,  0.5f, 0.0f}, {0.0f, 1.0f} }
	};

	std::vector<unsigned int> indices{
		0, 1, 3, // first triangle
		1, 2, 3  // second triangle
	};

	Mesh m_mesh{verts, indices, {}};
	Shader shaders{ "D:/Development/C++/ibm_application/ibm_application/src/shaders/MatrixViewer.vert",
	"D:/Development/C++/ibm_application/ibm_application/src/shaders/MatrixViewer.frag" };

	std::shared_ptr<Solution> m_solution = nullptr;
};