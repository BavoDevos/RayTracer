#include "Scene.h"
#include "Utils.h"
#include "Material.h"

namespace dae
{

#pragma region Base Scene
	//Initialize Scene with Default Solid Color Material (RED)
	Scene::Scene() :
		m_Materials({ new Material_SolidColor({1,0,0}) })
	{
		m_SphereGeometries.reserve(32);
		m_PlaneGeometries.reserve(32);
		m_TriangleMeshGeometries.reserve(32);
		m_Lights.reserve(32);
	}

	Scene::~Scene()
	{
		for (auto& pMaterial : m_Materials)
		{
			delete pMaterial;
			pMaterial = nullptr;
		}

		m_Materials.clear();
	}

	void dae::Scene::GetClosestHit(const Ray& ray, HitRecord& closestHit) const
	{
		//todo W1

		// Check the planes
		const size_t planeGeometriesSize{ m_PlaneGeometries.size() };
		for (size_t i{}; i < planeGeometriesSize; ++i)
		{
			HitRecord hitInfo{};
			GeometryUtils::HitTest_Plane(m_PlaneGeometries[i], ray, hitInfo);
			if (hitInfo.t > 0.0f && hitInfo.t < closestHit.t)
			{
				closestHit = hitInfo;
			}
		}

		// Check the spheres
		const size_t sphereGeometriesSize{ m_SphereGeometries.size() };
		for (size_t i{}; i < sphereGeometriesSize; ++i)
		{
			HitRecord hitInfo{};
			GeometryUtils::HitTest_Sphere(m_SphereGeometries[i], ray, hitInfo);
			if (hitInfo.t > 0.0f && hitInfo.t < closestHit.t)
			{
				closestHit = hitInfo;
			}
		}
		// Triangles
		const size_t triangleMeshGeometriesSize{ m_TriangleMeshGeometries.size() };
		for (size_t i{}; i < triangleMeshGeometriesSize; ++i)
		{
			HitRecord hitInfo{};
			GeometryUtils::HitTest_TriangleMesh(m_TriangleMeshGeometries[i], ray, hitInfo, false, true);
			if (hitInfo.t > 0.0f && hitInfo.t < closestHit.t)
			{
				closestHit = hitInfo;
			}
		}
	}

	bool Scene::DoesHit(const Ray& ray) const
	{
		for (const Plane& plane : m_PlaneGeometries)
		{
			if (GeometryUtils::HitTest_Plane(plane, ray))
				return true;
		}

		for (const Sphere& sphere : m_SphereGeometries)
		{
			if (GeometryUtils::HitTest_Sphere(sphere, ray))
				return true;
		}

		for (const Triangle& triangle : m_Triangles)
		{
			if (GeometryUtils::HitTest_Triangle(triangle, ray))
				return true;
		}

		for (const TriangleMesh& triangleMesh : m_TriangleMeshGeometries)
		{
			if (GeometryUtils::HitTest_TriangleMesh(triangleMesh, ray))
				return true;
		}

		return false;
	}

#pragma region Scene Helpers
	Sphere* Scene::AddSphere(const Vector3& origin, float radius, unsigned char materialIndex)
	{
		Sphere s;
		s.origin = origin;
		s.radius = radius;
		s.materialIndex = materialIndex;

		m_SphereGeometries.emplace_back(s);
		return &m_SphereGeometries.back();
	}

	Plane* Scene::AddPlane(const Vector3& origin, const Vector3& normal, unsigned char materialIndex)
	{
		Plane p;
		p.origin = origin;
		p.normal = normal;
		p.materialIndex = materialIndex;

		m_PlaneGeometries.emplace_back(p);
		return &m_PlaneGeometries.back();
	}

	TriangleMesh* Scene::AddTriangleMesh(TriangleCullMode cullMode, unsigned char materialIndex)
	{
		TriangleMesh m{};
		m.cullMode = cullMode;
		m.materialIndex = materialIndex;

		m_TriangleMeshGeometries.emplace_back(m);
		return &m_TriangleMeshGeometries.back();
	}

	Light* Scene::AddPointLight(const Vector3& origin, float intensity, const ColorRGB& color)
	{
		Light l;
		l.origin = origin;
		l.intensity = intensity;
		l.color = color;
		l.type = LightType::Point;
		l.samplePoints.push_back(origin);
		m_Lights.emplace_back(l);
		return &m_Lights.back();
	}

	Light* Scene::AddDirectionalLight(const Vector3& direction, float intensity, const ColorRGB& color)
	{
		Light l;
		l.direction = direction;
		l.intensity = intensity;
		l.color = color;
		l.type = LightType::Directional;

		m_Lights.emplace_back(l);
		return &m_Lights.back();
	}

	Light* Scene::AddSphereLight(const Vector3& origin, float intensity, const ColorRGB& color, float radius, int samplePoints)
	{
		Light l;
		l.origin = origin;
		l.intensity = intensity;
		l.color = color;
		l.radius = radius;
		l.type = LightType::Sphere;

		//this ensures an even distribution
		int samplesPerDimension = static_cast<int>(std::sqrt(samplePoints));

		for (int i = 0; i < samplesPerDimension; i++) {
			for (int j = 0; j < samplesPerDimension; j++) {
				float theta = 2.0f * static_cast<float>(M_PI) * ((i + Random()) / samplesPerDimension);
				float phi = static_cast<float>(M_PI) * ((j + Random()) / samplesPerDimension);
				float r = radius * sqrt(Random());
				float x = r * sin(phi) * cos(theta);
				float y = r * sin(phi) * sin(theta);
				float z = r * cos(phi);
				auto sampledPointOnLight = origin + Vector3{ x, y, z };
				l.samplePoints.push_back(sampledPointOnLight);
			}
		}

		m_Lights.emplace_back(l);
		return &m_Lights.back();
	}

	Light* Scene::AddPlaneLight(const Vector3& origin, const Vector3& normal, float intensity, const ColorRGB& color, float width, float height, int samplePoints)
	{
		Light l;
		l.origin = origin;
		l.intensity = intensity;
		l.color = color;
		l.type = LightType::Plane;

		int samplesPerDimension = static_cast<int>(std::sqrt(samplePoints));
		Vector3 xAxis, yAxis;

		// Ensure the normal is normalized
		Vector3 planeNormal = normal.Normalized();

		// Generate orthogonal axes for the plane
		if (std::abs(Vector3::Dot(planeNormal, Vector3::Up())) < 0.9f) {
			xAxis = Vector3::Cross(Vector3::Up(), planeNormal).Normalized();
		}
		else {
			xAxis = Vector3::Cross(Vector3::Right(), planeNormal).Normalized();
		}
		yAxis = Vector3::Cross(planeNormal, xAxis).Normalized();


		//ensures even distribution
		for (int i = 0; i < samplesPerDimension; i++) {
			for (int j = 0; j < samplesPerDimension; j++) {
				float u = (i + Random()) / samplesPerDimension;
				float v = (j + Random()) / samplesPerDimension;
				Vector3 sampledPointOnLight = origin + (u - 0.5f) * width * xAxis + (v - 0.5f) * height * yAxis;
				l.samplePoints.push_back(sampledPointOnLight);
			}
		}

		m_Lights.emplace_back(l);
		return &m_Lights.back();
	}


	unsigned char Scene::AddMaterial(Material* pMaterial)
	{
		m_Materials.push_back(pMaterial);
		return static_cast<unsigned char>(m_Materials.size() - 1);
	}
#pragma endregion

	void Scene_W4_ReferenceScene::Initialize()
	{
		sceneName = "Reference Scene";
		m_Camera.origin = { 0.f, 3.f, -9.f };
		m_Camera.fovAngle = 45.f;

		// Materials
		const auto matCt_GrayRoughMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 1.f));
		const auto matCt_GrayMediumMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 0.6f));
		const auto matCt_GraySmoothMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 0.1f));

		const auto matCt_GrayRoughPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 1.f));
		const auto matCt_GrayMediumPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 0.6f));
		const auto matCt_GraySmoothPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 0.1f));

		const auto matLambert_GrayBlue = AddMaterial(new Material_Lambert({ 0.49f, 0.57f, 0.57f }, 1.f));
		const auto matLambert_White = AddMaterial(new Material_Lambert(colors::White, 1.f));

		// Planes
		AddPlane({ 0.f, 0.f, 10.f }, { 0.f, 0.f, -1.f }, matLambert_GrayBlue);	// BACK
		AddPlane({ 0.f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, matLambert_GrayBlue);	// BOTTOM
		AddPlane({ 0.f, 10.f, 0.f }, { 0.f, -1.f, 0.f }, matLambert_GrayBlue);  // TOP
		AddPlane({ 5.f, 0.f, 0.f }, { -1.f, 0.f, 0.f }, matLambert_GrayBlue);	// RIGHT
		AddPlane({ -5.f, 0.f, 0.f }, { 1.f, 0.f, 0.f }, matLambert_GrayBlue);	// LEFT

		// Spheres
		AddSphere({ -1.75f, 1.0f, 0.0f }, 0.75f, matCt_GrayRoughMetal);
		AddSphere({ 0.0f, 1.0f, 0.0f }, 0.75f, matCt_GrayMediumMetal);
		AddSphere({ 1.75f, 1.0f, 0.0f }, 0.75f, matCt_GraySmoothMetal);

		AddSphere({ -1.75f, 3.0f, 0.0f }, 0.75f, matCt_GrayRoughPlastic);
		AddSphere({ 0.0f, 3.0f, 0.0f }, 0.75f, matCt_GrayMediumPlastic);
		AddSphere({ 1.75f, 3.0f, 0.0f }, 0.75f, matCt_GraySmoothPlastic);

		// Triangles
		const Triangle baseTriangle = { { -.75f, 1.5f, 0.f }, { .75f, 0.f, 0.f }, { -.75f, 0.f, 0.f } };

		m_Meshes[0] = AddTriangleMesh(TriangleCullMode::BackFaceCulling, matLambert_White);
		m_Meshes[0]->AppendTriangle(baseTriangle, true);
		m_Meshes[0]->Translate({ -1.75f, 4.5f, 0.f });
		m_Meshes[0]->CalculateNormals();
		m_Meshes[0]->UpdateAABB();
		m_Meshes[0]->UpdateTransforms();

		m_Meshes[1] = AddTriangleMesh(TriangleCullMode::FrontFaceCulling, matLambert_White);
		m_Meshes[1]->AppendTriangle(baseTriangle, true);
		m_Meshes[1]->Translate({ 0.f, 4.5f, 0.f });
		m_Meshes[1]->CalculateNormals();
		m_Meshes[1]->UpdateAABB();
		m_Meshes[1]->UpdateTransforms();

		m_Meshes[2] = AddTriangleMesh(TriangleCullMode::NoCulling, matLambert_White);
		m_Meshes[2]->AppendTriangle(baseTriangle, true);
		m_Meshes[2]->Translate({ 1.75f, 4.5f, 0.f });
		m_Meshes[2]->CalculateNormals();
		m_Meshes[2]->UpdateAABB();
		m_Meshes[2]->UpdateTransforms();


		//// Lights
		AddPointLight({ 0.f, 5.f, 5.f }, 50.f, { 1.f, .61f, .45f }); // BACKLIGHT
		AddPointLight({ -2.5f, 5.f, -5.f }, 70.f, { 1.f, .8f, .45f }); // FRONT LIGHT LEFT
		AddPointLight({ 2.5f, 2.5f, -5.f }, 50.f, { 0.34f, .47f, .68f });
		
		//// Lights
		//AddSphereLight({ 0.f, 5.f, 5.f }, 50.f, { 1.f, .61f, .45f },.75f,50); // BACKLIGHT
		//AddSphereLight({ -2.5f, 5.f, -5.f }, 70.f, { 1.f, .8f, .45f },.75f,50); // FRONT LIGHT LEFT
		//AddSphereLight({ 2.5f, 2.5f, -5.f }, 50.f, { 0.34f, .47f, .68f },0.75f,50);

		//AddPlaneLight({ 0,2,-5 }, { 0,0,1 }, 70, { 1,.61f,.45f }, 5.1, 5.2, 50);
	}
	void Scene_W4_ReferenceScene::Update(dae::Timer* pTimer)
	{
		Scene::Update(pTimer);  // run base class function

		const float yawAngle = (cos(pTimer->GetTotal()) + 1.f) / 2.f * PI_2;
		for (TriangleMesh* mesh : m_Meshes)
		{
			mesh->RotateY(yawAngle);
			mesh->UpdateTransforms();
		}
		
	}
	void Scene_W4_BunnyScene::Initialize()
	{
		sceneName = "Reference Scene";
		m_Camera.origin = { 0.f, 3.f, -9.f };
		m_Camera.fovAngle = 45.f;

		// Materials
		const auto matCt_GrayRoughMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 1.f));
		const auto matCt_GrayMediumMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 0.6f));
		const auto matCt_GraySmoothMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 0.1f));

		const auto matCt_GrayRoughPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 1.f));
		const auto matCt_GrayMediumPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 0.6f));
		const auto matCt_GraySmoothPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 0.1f));

		const auto matLambert_GrayBlue = AddMaterial(new Material_Lambert({ 0.49f, 0.57f, 0.57f }, 1.f));
		const auto matLambert_White = AddMaterial(new Material_Lambert(colors::White, 1.f));

		// Planes
		AddPlane({ 0.f, 0.f, 10.f }, { 0.f, 0.f, -1.f }, matLambert_GrayBlue);	// BACK
		AddPlane({ 0.f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, matLambert_GrayBlue);	// BOTTOM
		AddPlane({ 0.f, 10.f, 0.f }, { 0.f, -1.f, 0.f }, matLambert_GrayBlue);  // TOP
		AddPlane({ 5.f, 0.f, 0.f }, { -1.f, 0.f, 0.f }, matLambert_GrayBlue);	// RIGHT
		AddPlane({ -5.f, 0.f, 0.f }, { 1.f, 0.f, 0.f }, matLambert_GrayBlue);	// LEFT

		// Bunny
		pMesh = AddTriangleMesh(TriangleCullMode::BackFaceCulling, matLambert_White);
		Utils::ParseOBJ("Resources/lowpoly_bunny.obj", pMesh->positions, pMesh->normals, pMesh->indices);

		//pMesh->CalculateNormals();
		pMesh->Scale({ 2.f, 2.f, 2.f });

		pMesh->UpdateAABB();
		pMesh->UpdateTransforms();


		// Lights
		AddPointLight({ 0.f, 5.f, 5.f }, 50.f, { 1.f, .61f, .45f }); // BACKLIGHT
		AddPointLight({ -2.5f, 5.f, -5.f }, 70.f, { 1.f, .8f, .45f }); // FRONT LIGHT LEFT
		AddPointLight({ 2.5f, 2.5f, -5.f }, 50.f, { 0.34f, .47f, .68f });
	}
	void Scene_W4_BunnyScene::Update(dae::Timer* pTimer)
	{
		Scene::Update(pTimer);
		pMesh->RotateY(PI_DIV_4 * pTimer->GetTotal() / 4.0f);
		pMesh->UpdateTransforms();
	}

	void Scene_AreaLightSphere::Initialize()
	{
		sceneName = "Reference Scene";
		m_Camera.origin = { 0.f, 3.f, -9.f };
		m_Camera.fovAngle = 45.f;

		// Materials
		const auto matCt_GrayRoughMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 1.f));
		const auto matCt_GrayMediumMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 0.6f));
		const auto matCt_GraySmoothMetal = AddMaterial(new Material_CookTorrence({ 0.972f, 0.96f, 0.915f }, 1.f, 0.1f));

		const auto matCt_GrayRoughPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 1.f));
		const auto matCt_GrayMediumPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 0.6f));
		const auto matCt_GraySmoothPlastic = AddMaterial(new Material_CookTorrence({ 0.75f, 0.75f, 0.75f }, 0.f, 0.1f));

		const auto matLambert_GrayBlue = AddMaterial(new Material_Lambert({ 0.49f, 0.57f, 0.57f }, 1.f));
		const auto matLambert_White = AddMaterial(new Material_Lambert(colors::White, 1.f));

		// Planes
		AddPlane({ 0.f, 0.f, 10.f }, { 0.f, 0.f, -1.f }, matLambert_GrayBlue);	// BACK
		AddPlane({ 0.f, 0.f, 0.f }, { 0.f, 1.f, 0.f }, matLambert_GrayBlue);	// BOTTOM
		AddPlane({ 0.f, 10.f, 0.f }, { 0.f, -1.f, 0.f }, matLambert_GrayBlue);  // TOP
		AddPlane({ 5.f, 0.f, 0.f }, { -1.f, 0.f, 0.f }, matLambert_GrayBlue);	// RIGHT
		AddPlane({ -5.f, 0.f, 0.f }, { 1.f, 0.f, 0.f }, matLambert_GrayBlue);	// LEFT

		// Spheres
		AddSphere({ -1.75f, 1.0f, 0.0f }, 0.75f, matCt_GrayRoughMetal);
		AddSphere({ 0.0f, 1.0f, 0.0f }, 0.75f, matCt_GrayMediumMetal);
		AddSphere({ 1.75f, 1.0f, 0.0f }, 0.75f, matCt_GraySmoothMetal);

		AddSphere({ -1.75f, 3.0f, 0.0f }, 0.75f, matCt_GrayRoughPlastic);
		AddSphere({ 0.0f, 3.0f, 0.0f }, 0.75f, matCt_GrayMediumPlastic);
		AddSphere({ 1.75f, 3.0f, 0.0f }, 0.75f, matCt_GraySmoothPlastic);

		// Triangles
		const Triangle baseTriangle = { { -.75f, 1.5f, 0.f }, { .75f, 0.f, 0.f }, { -.75f, 0.f, 0.f } };

		m_Meshes[0] = AddTriangleMesh(TriangleCullMode::BackFaceCulling, matLambert_White);
		m_Meshes[0]->AppendTriangle(baseTriangle, true);
		m_Meshes[0]->Translate({ -1.75f, 4.5f, 0.f });
		m_Meshes[0]->CalculateNormals();
		m_Meshes[0]->UpdateAABB();
		m_Meshes[0]->UpdateTransforms();

		m_Meshes[1] = AddTriangleMesh(TriangleCullMode::FrontFaceCulling, matLambert_White);
		m_Meshes[1]->AppendTriangle(baseTriangle, true);
		m_Meshes[1]->Translate({ 0.f, 4.5f, 0.f });
		m_Meshes[1]->CalculateNormals();
		m_Meshes[1]->UpdateAABB();
		m_Meshes[1]->UpdateTransforms();

		m_Meshes[2] = AddTriangleMesh(TriangleCullMode::NoCulling, matLambert_White);
		m_Meshes[2]->AppendTriangle(baseTriangle, true);
		m_Meshes[2]->Translate({ 1.75f, 4.5f, 0.f });
		m_Meshes[2]->CalculateNormals();
		m_Meshes[2]->UpdateAABB();
		m_Meshes[2]->UpdateTransforms();

		//// Lights
		AddSphereLight({ 0.f, 5.f,-5.f }, 75.f, { 1.f, .61f, .45f },1.75f,20); // BACKLIGHT
		//AddPlaneLight({ 0,2,-5 }, { 0,0,1 }, 70, { 1,.61f,.45f }, 2.1, 2.2, 50);
	}
	void Scene_AreaLightSphere::Update(dae::Timer* pTimer)
	{
		Scene::Update(pTimer);  // run base class function

		const float yawAngle = (cos(pTimer->GetTotal()) + 1.f) / 2.f * PI_2;
		for (TriangleMesh* mesh : m_Meshes)
		{
			mesh->RotateY(yawAngle);
			mesh->UpdateTransforms();
		}

	}
}
