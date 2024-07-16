#pragma once
#include <cassert>

#include "Math.h"
#include "vector"
#include <iostream>

namespace dae
{
#pragma region GEOMETRY
	struct Sphere
	{
		Vector3 origin{};
		float radius{};

		unsigned char materialIndex{ 0 };
	};

	struct Plane
	{
		Vector3 origin{};
		Vector3 normal{};

		unsigned char materialIndex{ 0 };
	};

	enum class TriangleCullMode
	{
		FrontFaceCulling,
		BackFaceCulling,
		NoCulling
	};

	struct Triangle
	{
		Triangle() = default;
		Triangle(const Vector3& _v0, const Vector3& _v1, const Vector3& _v2, const Vector3& _normal) :
			v0{ _v0 }, v1{ _v1 }, v2{ _v2 }, normal{ _normal.Normalized() }
		{
		}

		Triangle(const Vector3& _v0, const Vector3& _v1, const Vector3& _v2) :
			v0{ _v0 }, v1{ _v1 }, v2{ _v2 }
		{
			const Vector3 edgeV0V1 = v1 - v0;
			const Vector3 edgeV0V2 = v2 - v0;
			normal = Vector3::Cross(edgeV0V1, edgeV0V2).Normalized();
		}

		Vector3 v0{};
		Vector3 v1{};
		Vector3 v2{};

		Vector3 normal{};

		TriangleCullMode cullMode{};
		unsigned char materialIndex{};
	};

	struct TriangleMesh
	{
		TriangleMesh() = default;
		TriangleMesh(const std::vector<Vector3>& _positions, const std::vector<int>& _indices, TriangleCullMode _cullMode) :
			positions(_positions), indices(_indices), cullMode(_cullMode)
		{
			//Calculate Normals
			CalculateNormals();

			//Update Transforms
			UpdateTransforms();
		}

		TriangleMesh(const std::vector<Vector3>& _positions, const std::vector<int>& _indices, const std::vector<Vector3>& _normals, TriangleCullMode _cullMode) :
			positions(_positions), indices(_indices), normals(_normals), cullMode(_cullMode)
		{
			UpdateTransforms();
		}

		std::vector<Vector3> positions{};
		std::vector<Vector3> normals{};
		std::vector<int> indices{};
		unsigned char materialIndex{};

		TriangleCullMode cullMode{ TriangleCullMode::BackFaceCulling };

		Matrix rotationTransform{};
		Matrix translationTransform{};
		Matrix scaleTransform{};

		// Min/Max Axis-Aligned-Bounding-Box
		Vector3 minAABB{};
		Vector3 maxAABB{};

		// Transformed Min/Max Axis-Aligned-Bounding-Box (after applying the transform)
		Vector3 transformedMinAABB{};
		Vector3 transformedMaxAABB{};


		std::vector<Vector3> transformedPositions{};
		std::vector<Vector3> transformedNormals{};


		void Translate(const Vector3& translation)
		{
			translationTransform = Matrix::CreateTranslation(translation);
		}

		void RotateY(float yaw)
		{
			rotationTransform = Matrix::CreateRotationY(yaw);
		}

		void Scale(float uniformScale)
		{
			Scale({ uniformScale, uniformScale, uniformScale });
		}

		void Scale(const Vector3& scale)
		{
			scaleTransform = Matrix::CreateScale(scale);
		}

		void AppendTriangle(const Triangle& triangle, bool ignoreTransformUpdate = false)
		{
			int startIndex = static_cast<int>(positions.size());

			positions.push_back(triangle.v0);
			positions.push_back(triangle.v1);
			positions.push_back(triangle.v2);

			indices.push_back(startIndex);
			indices.push_back(++startIndex);
			indices.push_back(++startIndex);

			//Not ideal, but making sure all vertices are updated
			if (!ignoreTransformUpdate)
				UpdateTransforms();
		}

		void CalculateNormals()
		{
			// Cross of 2 edges of the triangle (clockwise order of vertices because left handed system used
			// Each indice specifies the vertices of the triangle

			normals.clear();
			normals.reserve(indices.size() / 3);

			for (size_t i{}; i < indices.size(); i += 3)
			{
				const Vector3 v0{ positions[indices[i]] };
				const Vector3 v1{ positions[indices[i + 1]] };
				const Vector3 v2{ positions[indices[i + 2]] };

				const Vector3 edgeA = v1 - v0;
				const Vector3 edgeB = v2 - v1;
				// No need for 3rd edge

				normals.emplace_back(Vector3::Cross(edgeA, edgeB).Normalized());
			}
		}

		void UpdateTransforms()
		{
			//Calculate Final Transform 
			// First scale, then rotate, then translate
			const auto finalTransform = scaleTransform * rotationTransform * translationTransform;

			// Loop over every position & apply the transformation
			transformedPositions.clear();
			transformedPositions.reserve(positions.size());
			for (const Vector3& p : positions)
			{
				transformedPositions.emplace_back(finalTransform.TransformPoint(p));
			}


			//Transform Normals (normals > transformedNormals)
			//...
			transformedNormals.clear();
			transformedNormals.reserve(normals.size());
			for (const Vector3& n : normals)
			{
				transformedNormals.emplace_back(rotationTransform.TransformVector(n));
			}

			UpdateTransformedAABB(finalTransform);
		}

		void UpdateAABB()
		{
			//Update Min/Max Axis-Aligned-Bounding-Box
			// Init with smallest/biggest possible values
			minAABB = { FLT_MAX, FLT_MAX, FLT_MAX };
			maxAABB = { -FLT_MAX, -FLT_MAX, -FLT_MAX };

			// Loop over every position, compare with the current min & max, and update it with the new min / max
			for (const Vector3& p : positions)
			{
				minAABB = Vector3::Min(minAABB, p);
				maxAABB = Vector3::Max(maxAABB, p);
			}

		}

		void UpdateTransformedAABB(const Matrix& finalTransform)
		{
			// Instead of transforming every position, we can use the min/max AABB to calculate the transformed AAB
			// We loop over all points of the cube (8) and check which one is the new min & max after transforming			
			
			for (int x = 0; x < 2; x++)
			{
				for (int y = 0; y < 2; y++)
				{
					for (int z = 0; z < 2; z++)
					{
						const Vector3 p = finalTransform.TransformPoint(Vector3{ x ? minAABB.x : maxAABB.x, y ? minAABB.y : maxAABB.y, z ? minAABB.z : maxAABB.z });
						transformedMinAABB = Vector3::Min(transformedMinAABB, p);
						transformedMaxAABB = Vector3::Max(transformedMaxAABB, p);
					}
				}
			}	
		}
	};
#pragma endregion
#pragma region LIGHT
	enum class LightType
	{
		Point,
		Directional,
		Sphere,
		Plane
	};

	struct Light
	{
		Vector3 origin{};
		Vector3 direction{};
		Vector3 RectDimensions{};
		ColorRGB color{};
		float intensity{};
		float radius{};
		LightType type{};
		std::vector<Vector3>samplePoints{};
	};
#pragma endregion
#pragma region MISC
	struct Ray
	{
		Vector3 origin{};
		Vector3 direction{};

		float min{ 0.0001f };
		float max{ FLT_MAX };
	};

	struct HitRecord
	{
		Vector3 origin{};
		Vector3 normal{};
		float t = FLT_MAX;

		bool didHit{ false };
		unsigned char materialIndex{ 0 };
	};
#pragma endregion
}