//External includes
#include "vld.h"
#include "SDL.h"
#include "SDL_surface.h"
#undef main

//Standard includes
#include <iostream>

//Project includes
#include "Timer.h"
#include "Renderer.h"
#include "Scene.h"

using namespace dae;

void ShutDown(SDL_Window* pWindow)
{
	SDL_DestroyWindow(pWindow);
	SDL_Quit();
}

int main(int argc, char* args[])
{
	//Unreferenced parameters
	(void)argc;
	(void)args;

	//Create window + surfaces
	SDL_Init(SDL_INIT_VIDEO);

	const uint32_t width = 640;
	const uint32_t height = 480;

	SDL_Window* pWindow = SDL_CreateWindow(
		"RayTracer - Bavo Devos/DAE07",
		SDL_WINDOWPOS_UNDEFINED,
		SDL_WINDOWPOS_UNDEFINED,
		width, height, 0);

	if (!pWindow)
		return 1;
	
	//Initialize "framework"
	const auto pTimer = new Timer();
	const auto pRenderer = new Renderer(pWindow);

	const auto pScene = new Scene_W4_BunnyScene;
	pScene->Initialize();

	//Start loop
	pTimer->Start();
	float printTimer = 0.f;
	bool isLooping = true;
	bool takeScreenshot = false;
	while (isLooping)
	{
		//--------- Get input events ---------
		SDL_Event e;
		while (SDL_PollEvent(&e))
		{
			switch (e.type)
			{
			case SDL_QUIT:
				isLooping = false;
				break;
			case SDL_KEYUP:
				switch (e.key.keysym.scancode)
				{
					case SDL_SCANCODE_X:
						takeScreenshot = true;
						break;
					case SDL_SCANCODE_F2:
						if (not e.key.repeat) pRenderer->ToggleShadows();
						break;
					case SDL_SCANCODE_F3:
						if (not e.key.repeat) pRenderer->CycleLightingMode();
						break;
					case SDL_SCANCODE_F6:
						if (not e.key.repeat) pTimer->StartBenchmark(10);
						break;
				}
			}
			
		}

		//--------- Update ---------
		pScene->Update(pTimer);

		//--------- Render ---------
		pRenderer->Render(pScene);

		//--------- Timer ---------
		pTimer->Update();
		printTimer += pTimer->GetElapsed();
		if (printTimer >= 1.f)
		{
			printTimer = 0.f;
			std::cout << "dFPS: " << pTimer->GetdFPS()
				<< "\tms/frame: " << pTimer->GetElapsed() * 1000.0f<<'\n';
			/*	<< "\tavg(" << pTimer->GetHistorySize() << ") ms/frame: " << pTimer->GetAverageElapsed() * 1000.0f
				<< "\n";*/
		}

		//Save screenshot after full render
		if (takeScreenshot)
		{
			if (!pRenderer->SaveBufferToImage())
				std::cout << "Screenshot saved!" << "\n";
			else
				std::cout << "Something went wrong. Screenshot not saved!" << "\n";
			takeScreenshot = false;
		}
	}
	pTimer->Stop();

	//Shutdown "framework"
	delete pScene;
	delete pRenderer;
	delete pTimer;

	ShutDown(pWindow);
	return 0;
}