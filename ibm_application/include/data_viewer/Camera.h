#pragma once

#include <tuple>

#include "SDL.h"
#include "ext.hpp"

class Camera
{
public:
	Camera()
	{
		view = glm::translate(view, loc);
		model = glm::scale(model, scale);
	}

	void MouseTranslate(glm::vec2 delta)
	{
		glm::vec3 input = glm::vec3{ -delta, 0.0 };
		input.x /= model[0][0];
		input.y /= model[1][1];
		view = glm::translate(view, input);
	}

	void ScrollZoom(float scroll)
	{
		const float speed = 0.08f;
		glm::vec3 input{ glm::vec3{speed, speed, 0.0} * scroll };
		model = glm::scale(model, glm::vec3{ 1.0 } + input);
	}

	void ProcessKeyboardCommands(SDL_Keycode key)
	{
		KeyboardTranslate(key);
		KeyboardZoom(key);
	}

	void KeyboardTranslate(SDL_Keycode key)
	{
		glm::vec3 translate{ 0.0 };

		if (key == SDLK_UP || key == SDLK_w)
		{
			translate.y = -1.0;
		}
		if (key == SDLK_DOWN || key == SDLK_s)
		{
			translate.y = 1.0;
		}
		if (key == SDLK_LEFT || key == SDLK_a)
		{
			translate.x = 1.0;
		}
		if (key == SDLK_RIGHT || key == SDLK_d)
		{
			translate.x = -1.0;
		}

		translate *= 0.02;

		view = glm::translate(view, translate);
	}

	void KeyboardZoom(SDL_Keycode key)
	{
		glm::vec3 input{ 0.0 };
		if (key == SDLK_z)
		{
			input.x = 1.0;
			input.y = 1.0;
		}
		if (key == SDLK_x)
		{
			input.x = -1.0;
			input.y = -1.0;
		}

		input *= 0.04;
		model = glm::scale(model, glm::vec3{1.0} + input);
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