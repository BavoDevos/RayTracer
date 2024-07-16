#pragma once
#include <cassert>
#include <SDL_keyboard.h>
#include <SDL_mouse.h>

#include "Math.h"
#include "Timer.h"
#include <iostream>
#include <algorithm>

namespace dae
{
	struct Camera
	{
		Camera() = default;

		Camera(const Vector3& _origin, float _fovAngle) :
			origin{ _origin },
			fovAngle{ _fovAngle }
		{
		}


		Vector3 origin{};
		float fovAngle{ 45.f };
		const float movementSpeed{ 7.0f };
		const float rotationSpeed{ 20.0f };
		const float keyboardRotationSpeed{ 80.0f };

		Vector3 forward{ Vector3::UnitZ };
		Vector3 up{ Vector3::UnitY };
		Vector3 right{ Vector3::UnitX };

		float totalPitch{ 0.f };
		float totalYaw{ 0.f };

		Matrix cameraToWorld{};


		Matrix CalculateCameraToWorld()
		{
			right = Vector3::Cross(Vector3::UnitY, forward).Normalized();
			up = Vector3::Cross(forward, right).Normalized();

			cameraToWorld = Matrix{
				{right.x  , right.y   , right.z  , 0},
				{up.x     , up.y      , up.z     , 0},
				{forward.x, forward.y , forward.z, 0},
				{origin.x , origin.y  , origin.z , 1}
			};

			return cameraToWorld;
		}

		void SetFovAngle(float fov)
		{

		}

		void Update(Timer* pTimer)
		{
			const float deltaTime = pTimer->GetElapsed();

			//Keyboard Input
			const uint8_t* pKeyboardState = SDL_GetKeyboardState(nullptr);

			// Keyboard movement of the camera
			if (pKeyboardState[SDL_SCANCODE_W]|| pKeyboardState[SDL_SCANCODE_UP])
			{
				origin += forward * movementSpeed * deltaTime;
			}
			if (pKeyboardState[SDL_SCANCODE_S] || pKeyboardState[SDL_SCANCODE_DOWN])
			{
				origin -= forward * movementSpeed * deltaTime;
			}
			if (pKeyboardState[SDL_SCANCODE_D] || pKeyboardState[SDL_SCANCODE_RIGHT])
			{
				origin += right * movementSpeed * deltaTime;
			}
			if (pKeyboardState[SDL_SCANCODE_A] || pKeyboardState[SDL_SCANCODE_LEFT])
			{
				origin -= right * movementSpeed * deltaTime;
			}
			if (pKeyboardState[SDL_SCANCODE_SPACE])
			{
				origin += up * movementSpeed * deltaTime;
			}
			if (pKeyboardState[SDL_SCANCODE_LSHIFT])
			{
				origin -= up * movementSpeed * deltaTime;
			}
			//Mouse Input
			int mouseX{}, mouseY{};
			const uint32_t mouseState = SDL_GetRelativeMouseState(&mouseX, &mouseY);

			// Mouse movements / rotation of the camera
			if ((mouseState & SDL_BUTTON(SDL_BUTTON_LEFT)) && mouseState & SDL_BUTTON(SDL_BUTTON_RIGHT))
			{
				// mouseX yaw left & right, mouse Y moves forwards & backwards
				const float upwards = -mouseY * movementSpeed * deltaTime;
				origin += up * upwards;
			}
			else if (mouseState & SDL_BUTTON(SDL_BUTTON_LEFT))
			{
				// mouseX yaw left & right, mouse Y moves forwards & backwards
				const float forwards = -mouseY * deltaTime;
				const float yaw = mouseX * deltaTime;

				origin += forward * forwards;
				totalYaw += yaw;
			}
			else if (mouseState & SDL_BUTTON(SDL_BUTTON_RIGHT))
			{
				// Look around the current origin
				const float pitch = -mouseY * rotationSpeed * deltaTime;
				const float yaw = mouseX * rotationSpeed * deltaTime;

				totalPitch += pitch;
				totalYaw += yaw;
			}


			totalPitch = std::clamp(totalPitch, -88.0f, 88.0f);
			if (totalYaw > 360.0f)
				totalYaw -= 360.0f;
			else if (totalYaw < 0.0f)
				totalYaw += 360.0f;

			const Matrix finalRotation = Matrix::CreateRotationX(totalPitch * TO_RADIANS) * Matrix::CreateRotationY(totalYaw * TO_RADIANS);
			forward = finalRotation.TransformVector(Vector3::UnitZ);
			forward.Normalize();
		}
	};
}
