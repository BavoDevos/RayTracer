//External includes
#include "SDL.h"
#include "SDL_surface.h"

//Project includes
#include "Renderer.h"
#include "Math.h"
#include "Matrix.h"
#include "Material.h"
#include "Scene.h"
#include "Utils.h"
#include <thread>
#include "camera.h"
#include <future>
#include <ppl.h>
#include "SDL_pixels.h"
#include <vector>
#include <cmath>
#include <numeric>

using namespace dae;

// Both in comment -> synchronous execution
//#define ASYNC
#define PARALLEL_FOR


Renderer::Renderer(SDL_Window* pWindow) :
	m_pWindow(pWindow),
	m_pBuffer(SDL_GetWindowSurface(pWindow))
{
	//Initialize
	SDL_GetWindowSize(pWindow, &m_Width, &m_Height);
	m_pBufferPixels = static_cast<uint32_t*>(m_pBuffer->pixels);
	RunTests();
}

void Renderer::Render(Scene* pScene) const
{
	Camera& camera = pScene->GetCamera();
	auto& materials = pScene->GetMaterials();
	auto& lights = pScene->GetLights();

	const float aspectRatio{ m_Width / float(m_Height) };

	const float fovRatio{ tan(camera.fovAngle * TO_RADIANS / 2.0f) };
	camera.CalculateCameraToWorld();

	const uint32_t numPixels = m_Width * m_Height;



#if defined(ASYNC)
	// ASYNC EXECUTION WITH THREADS
	// Get the number of cores of the system
	const uint32_t numCores{ std::thread::hardware_concurrency() };

	// Vector to keep track of all the async futures
	std::vector<std::future<void>> async_futures{};

	// Calculate how many pixels per task
	const uint32_t pixelsPerTask{ numPixels / numCores };
	uint32_t unassignedPixels{ numPixels % numCores }; // Pixels that are not assigned to a task
	uint32_t currPixelIndex{ 0 };

	// Create a task for each core
	for (uint32_t coreId{ 0 }; coreId < numCores; ++coreId)
	{
		uint32_t taskSize = pixelsPerTask;
		if (unassignedPixels > 0)
		{
			++taskSize;
			--unassignedPixels;
		}

		async_futures.push_back(
			std::async(std::launch::async, [=, this]
				{
					const uint32_t endPixel = currPixelIndex + taskSize;
		for (uint32_t pixelIndex{ currPixelIndex }; pixelIndex < endPixel; ++pixelIndex)
		{
			RenderPixel(pScene, pixelIndex, fovRatio, aspectRatio, camera, lights, materials);
		}
				}
			)
		);

		currPixelIndex += taskSize;
	}

	// Wait for all tasks to be finished
	for (const std::future<void>& future : async_futures)
	{
		future.wait();
	}
#elif defined(PARALLEL_FOR)
	// PARALLEL FOR EXECUTION
	//concurrency::parallel_for()
	concurrency::parallel_for((uint32_t)0, numPixels, [=, this](int pixelIndex)
		{
			RenderPixel(pScene, pixelIndex, fovRatio, aspectRatio, camera, lights, materials);
		});

#else
	// SYNCHRONOUS EXECUTION
	for (uint32_t pixelIndex{}; pixelIndex < numPixels; ++pixelIndex)
	{
		RenderPixel(pScene, pixelIndex, fovRatio, aspectRatio, camera, lights, materials);
	}

#endif
	//@END

	// Update SDL Surface
	SDL_UpdateWindowSurface(m_pWindow);
}

void Renderer::RenderPixel(Scene* pScene, uint32_t pixelIndex, float fov, float aspectRatio, const Camera& camera, const std::vector<Light>& lights, const std::vector<Material*>& materials) const
{
	uint32_t px{ pixelIndex % m_Width };
	uint32_t py{ pixelIndex / m_Width };

	const float cx{ ((2.0f * (px + 0.5f) / float(m_Width)) - 1) * aspectRatio * fov };
	const float cy{ (1.0f - ((2.0f * (py + 0.5f)) / float(m_Height))) * fov };

	float multiplier = 1.0f;

	const Vector3 rayDirection{ camera.cameraToWorld.TransformVector(Vector3{cx, cy, 1}).Normalized() };
	Ray viewRay{ camera.origin,  rayDirection };

	ColorRGB finalColor{};
	for (int i{}; i < m_Bounces; i++)
	{
		HitRecord closestHit{};
		pScene->GetClosestHit(viewRay, closestHit);  // Checks EVERY object in the scene and returns the closest one hit.
		if (closestHit.didHit)
		{
			for (const Light& light : lights)
			{
				finalColor += LightUtils::GetFinalColorFromLight(light, closestHit, m_ShadowsEnabled, m_CurrentLightingMode, pScene, materials, rayDirection, multiplier);
			}

			multiplier *= 0.7f;
			viewRay.origin = closestHit.origin + closestHit.normal * 0.0001f;
			viewRay.direction = Vector3::Reflect(viewRay.direction, closestHit.normal);
		}
		else
		{
			ColorRGB skyColor{ colors::White };
			finalColor += skyColor;
		}
		//Update Color in Buffer
		finalColor.MaxToOne();
		m_pBufferPixels[px + (py * m_Width)] = SDL_MapRGB(m_pBuffer->format,
			static_cast<uint8_t>(finalColor.r * 255),
			static_cast<uint8_t>(finalColor.g * 255),
			static_cast<uint8_t>(finalColor.b * 255));
	}

}


bool Renderer::SaveBufferToImage() const
{
	return SDL_SaveBMP(m_pBuffer, "RayTracing_Buffer.bmp");
}

void dae::Renderer::CycleLightingMode()
{
	switch (m_CurrentLightingMode)
	{
	case dae::Renderer::LightingMode::ObservedArea:
		m_CurrentLightingMode = LightingMode::Radiance;
		std::cout << "LightingMode: Radiance\n";
		break;
	case dae::Renderer::LightingMode::Radiance:
		m_CurrentLightingMode = LightingMode::BRDF;
		std::cout << "LightingMode: BRDF\n";
		break;
	case dae::Renderer::LightingMode::BRDF:
		m_CurrentLightingMode = LightingMode::Combined;
		std::cout << "LightingMode: Combined\n";
		break;
	case dae::Renderer::LightingMode::Combined:
		m_CurrentLightingMode = LightingMode::ObservedArea;
		std::cout << "LightingMode: ObservedArea\n";
		break;
	}
}

int Renderer::RunTests()
{
	// Test dot & cross product for vector3 & vector4
	assert(int(roundf(Vector3::Dot(Vector3::UnitX, Vector3::UnitX))) == 1);  // Should be 1 -> same direction
	assert(int(roundf(Vector3::Dot(Vector3::UnitX, -Vector3::UnitX))) == -1);  // Should be -1 -> opposite direction
	assert(int(roundf(Vector3::Dot(Vector3::UnitX, Vector3::UnitY))) == 0);  // Should be 0 -> perpendicular direction


	return 0;
}
