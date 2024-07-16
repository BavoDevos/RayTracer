#pragma once
#include <string>
#include <vector>

#include "Math.h"
#include "DataTypes.h"
#include "Camera.h"

namespace dae
{
	//Forward Declarations
	class Timer;
	class Material;
	struct Plane;
	struct Sphere;
	struct Light;

	//Scene Base Class
	class Scene
	{
	public:
		Scene();
		virtual ~Scene();

		Scene(const Scene&) = delete;
		Scene(Scene&&) noexcept = delete;
		Scene& operator=(const Scene&) = delete;
		Scene& operator=(Scene&&) noexcept = delete;

		virtual void Initialize() = 0;
		virtual void Update(dae::Timer* pTimer)
		{
			m_Camera.Update(pTimer);
		}

		Camera& GetCamera() { return m_Camera; }
		void GetClosestHit(const Ray& ray, HitRecord& closestHit) const;
		bool DoesHit(const Ray& ray) const;

		const std::vector<Plane>& GetPlaneGeometries() const { return m_PlaneGeometries; }
		const std::vector<Sphere>& GetSphereGeometries() const { return m_SphereGeometries; }
		const std::vector<Light>& GetLights() const { return m_Lights; }
		const std::vector<Material*> GetMaterials() const { return m_Materials; }

	protected:
		std::string	sceneName;

		std::vector<Plane> m_PlaneGeometries{};
		std::vector<Sphere> m_SphereGeometries{};
		std::vector<TriangleMesh> m_TriangleMeshGeometries{};
		std::vector<Light> m_Lights{};
		std::vector<Material*> m_Materials{};
		
		// Temp for triangles
		std::vector<Triangle> m_Triangles{};

		Camera m_Camera{};

		Sphere* AddSphere(const Vector3& origin, float radius, unsigned char materialIndex = 0);
		Plane* AddPlane(const Vector3& origin, const Vector3& normal, unsigned char materialIndex = 0);
		TriangleMesh* AddTriangleMesh(TriangleCullMode cullMode, unsigned char materialIndex = 0);

		Light* AddPointLight(const Vector3& origin, float intensity, const ColorRGB& color);
		Light* AddDirectionalLight(const Vector3& direction, float intensity, const ColorRGB& color);
		Light* AddSphereLight(const Vector3& origin, float intensity,const ColorRGB& color,float radius,int samplePoints);
		Light* AddPlaneLight(const Vector3& origin, const Vector3& normal, float intensity, const ColorRGB& color, float width, float height, int samplePoints);
		unsigned char AddMaterial(Material* pMaterial);
	};


	//i removed previous scenes to perserve space as they were not used here anyway, also to have to scroll less.

	//+++++++++++++++++++++++++++++++++++++++++
	//Reference Scene
	class Scene_W4_ReferenceScene final : public Scene
	{
	public:
		Scene_W4_ReferenceScene() = default;
		~Scene_W4_ReferenceScene() override = default;

		Scene_W4_ReferenceScene(const Scene_W4_ReferenceScene&) = delete;
		Scene_W4_ReferenceScene(Scene_W4_ReferenceScene&&) noexcept = delete;
		Scene_W4_ReferenceScene& operator=(const Scene_W4_ReferenceScene&) = delete;
		Scene_W4_ReferenceScene& operator=(Scene_W4_ReferenceScene&&) noexcept = delete;

		void Initialize() override;
		void Update(dae::Timer* pTimer) override;

	private:
		TriangleMesh* m_Meshes[3]{};

	};

	//+++++++++++++++++++++++++++++++++++++++++
	//Bunny Scene
	class Scene_W4_BunnyScene final : public Scene
	{
	public:
		Scene_W4_BunnyScene() = default;
		~Scene_W4_BunnyScene() override = default;

		Scene_W4_BunnyScene(const Scene_W4_BunnyScene&) = delete;
		Scene_W4_BunnyScene(Scene_W4_BunnyScene&&) noexcept = delete;
		Scene_W4_BunnyScene& operator=(const Scene_W4_BunnyScene&) = delete;
		Scene_W4_BunnyScene& operator=(Scene_W4_BunnyScene&&) noexcept = delete;

		void Initialize() override;
		void Update(dae::Timer* pTimer) override;

	private:
		TriangleMesh* pMesh{ nullptr };

	};

	//+++++++++++++++++++++++++++++++++++++++++
	//AREALIGHT SPHERE Scene
	class Scene_AreaLightSphere final : public Scene
	{
	public:
		Scene_AreaLightSphere() = default;
		~Scene_AreaLightSphere() override = default;

		Scene_AreaLightSphere(const Scene_AreaLightSphere&) = delete;
		Scene_AreaLightSphere(Scene_AreaLightSphere&&) noexcept = delete;
		Scene_AreaLightSphere& operator=(const Scene_AreaLightSphere&) = delete;
		Scene_AreaLightSphere& operator=(Scene_AreaLightSphere&&) noexcept = delete;

		void Initialize() override;
		void Update(dae::Timer* pTimer) override;

	private:
		TriangleMesh* m_Meshes[3]{};

	};

}
