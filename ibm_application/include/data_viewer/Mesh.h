#pragma once

#include <glad/glad.h> // holds all OpenGL type declarations

#include <glm.hpp>
#include <gtc/matrix_transform.hpp>

#include <data_viewer/Shader.h>

#include <string>
#include <vector>

struct Vertex {
    // position
    glm::vec3 Position;
    // normal
    //glm::vec3 Normal;
    // texCoords
    glm::vec2 TexCoords;
    // tangent
    //glm::vec3 Tangent;
    // bitangent
    //glm::vec3 Bitangent;
    //bone indexes which will influence this vertex
    //int m_BoneIDs[MAX_BONE_INFLUENCE];
    //weights from each bone
    //float m_Weights[MAX_BONE_INFLUENCE];
};

struct Texture {
    unsigned int id;
    std::string type;
};

class Mesh {
public:
    // mesh Data
    std::vector<Vertex>       vertices;
    std::vector<unsigned int> indices;
    std::vector<Texture>      textures;
    unsigned int VAO;

    unsigned int out_texture_buffer;
    std::pair<int, int> viewport_width_height{ 800, 800 };

    glm::vec3 pos{ 0.0 };

    glm::mat4 mvp = glm::mat4(1.0f);

    void SetMVP(glm::mat4 model, glm::mat4 view, glm::mat4 proj)
    {
        auto t_model = glm::translate(model, pos);
        mvp = proj * t_model * view;
    }
    void SetMVP(const std::tuple<glm::mat4, glm::mat4, glm::mat4>& mvp_in)
    {
        auto t_model = glm::translate(std::get<0>(mvp_in), pos);
        mvp = std::get<2>(mvp_in) * t_model * std::get<1>(mvp_in);
    }

    // constructor
    Mesh(std::vector<Vertex> vertices, std::vector<unsigned int> indices, std::vector<Texture> textures)
    {
        this->vertices = vertices;
        this->indices = indices;
        this->textures = textures;

        // now that we have all the required data, set the vertex buffers and its attribute pointers.
        setupMesh();
    }

    // render the mesh
    void Draw(Shader& shader)
    {
        // bind uniforms


        // bind appropriate textures
        unsigned int diffuseNr = 1;
        unsigned int specularNr = 1;
        unsigned int normalNr = 1;
        unsigned int heightNr = 1;
        for (unsigned int i = 0; i < textures.size(); i++)
        {
            glActiveTexture(GL_TEXTURE0 + i); // active proper texture unit before binding
            // retrieve texture number (the N in diffuse_textureN)
            std::string number;
            std::string name = textures[i].type;
            if (name == "texture_diffuse")
                number = std::to_string(diffuseNr++);
            else if (name == "texture_specular")
                number = std::to_string(specularNr++); // transfer unsigned int to string
            else if (name == "texture_normal")
                number = std::to_string(normalNr++); // transfer unsigned int to string
            else if (name == "texture_height")
                number = std::to_string(heightNr++); // transfer unsigned int to string

            // now set the sampler to the correct texture unit
            glUniform1i(glGetUniformLocation(shader.ID, (name + number).c_str()), i);
            // and finally bind the texture
            glBindTexture(GL_TEXTURE_2D, textures[i].id);
        }

        // draw mesh
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(indices.size()), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        // always good practice to set everything back to defaults once configured.
        glActiveTexture(GL_TEXTURE0);
    }

    void DrawToTexture(Shader& shader)
    {
        // store previous values
        int prev_viewport[4];
        glGetIntegerv(GL_VIEWPORT, prev_viewport);

        // bind to framebuffer and draw scene as we normally would to color texture 
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        //glEnable(GL_DEPTH_TEST); // enable depth testing (is disabled for rendering screen-space quad)

        // make sure we clear the framebuffer's content and set viewport size
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT); // DEPTH BUFFER BIT
        glViewport(0, 0, viewport_width_height.first, viewport_width_height.second);

        // bind appropriate textures
        unsigned int data_nr = 1;
        unsigned int specularNr = 1;
        for (unsigned int i = 0; i < textures.size(); i++)
        {
            glActiveTexture(GL_TEXTURE0 + i); // active proper texture unit before binding
            // retrieve texture number (the N in diffuse_textureN)
            std::string number;
            std::string name = textures[i].type;
            if (name == "data_texture")
            {
                number = std::to_string(data_nr++);
            }
            else if (name == "texture_specular")
            {
                number = std::to_string(specularNr++); // transfer unsigned int to string   
            }

            // now set the sampler to the correct texture unit
            //glUniform1i(glGetUniformLocation(shader.ID, (name + number).c_str()), i);
            // and finally bind the texture
            glBindTexture(GL_TEXTURE_2D, textures[i].id);
        }

        // draw mesh
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, static_cast<unsigned int>(indices.size()), GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        // always good practice to set everything back to defaults once configured.
        glActiveTexture(GL_TEXTURE0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(prev_viewport[0], prev_viewport[1], prev_viewport[2], prev_viewport[3]);
    }

private:
    // render data 
    unsigned int VBO, EBO, FBO;

    // initializes all the buffer objects/arrays
    void setupMesh()
    {
        // create buffers/arrays
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);
        // load data into vertex buffers
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        // A great thing about structs is that their memory layout is sequential for all its items.
        // The effect is that we can simply pass a pointer to the struct and it translates perfectly to a glm::vec3/2 array which
        // again translates to 3/2 floats which translates to a byte array.
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

        // set the vertex attribute pointers
        // vertex Positions
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
        // vertex normals
        //glEnableVertexAttribArray(1);
        //glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, Normal));
        // vertex texture coords
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, TexCoords));
        // vertex tangent
        //glEnableVertexAttribArray(3);
        //glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, Tangent));
        // vertex bitangent
        //glEnableVertexAttribArray(4);
        //glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, Bitangent));
        // ids
        //glEnableVertexAttribArray(5);
        //glVertexAttribIPointer(5, 4, GL_INT, sizeof(Vertex), (void*)offsetof(Vertex, m_BoneIDs));

        // weights
        //glEnableVertexAttribArray(6);
        //glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, m_Weights));

        // framebuffer
        glGenFramebuffers(1, &FBO);
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);

        glGenTextures(1, &out_texture_buffer);
        glBindTexture(GL_TEXTURE_2D, out_texture_buffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, viewport_width_height.first, viewport_width_height.second, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, out_texture_buffer, 0);

    	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // finish up
    	glBindVertexArray(0);
    	glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
};