#pragma once

#include <tuple>

#include "ext.hpp"

class Camera
{
public:
	Camera()
	{
		view = glm::translate(view, loc);
		model = glm::scale(model, scale);
	}

	std::tuple<glm::mat4, glm::mat4, glm::mat4> GetMVP() const
	{
		return { model, view, projection };
	}
private:
	glm::vec3 loc = glm::vec3(0.0f, 0.0f, -1.0f);

	glm::mat4 model = glm::mat4(1.0f);
	glm::vec3 scale = glm::vec3(1.f);

	glm::mat4 view = glm::mat4(1.0f);
	glm::mat4 projection = glm::ortho(-0.5f, 0.5f, -0.5f, 0.5f, 0.1f, 100.0f);
};