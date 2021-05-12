#include <iostream>
#include <functional>

#include <GLFW/glfw3.h>
#include <boost/heap/priority_queue.hpp>
#include <stack>
#include <algorithm>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "imgui.h"

#include <iomanip>


using namespace geometrycentral;
using namespace geometrycentral::surface;

//Mesh and Geometry
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;



//Geodesic Properties (Might be stored in a class, or struct as a better and safer design)
//For now making these global to easilty access them from IMGUI
std::vector<std::set<std::pair<size_t, double>>> graph;
std::vector<int> parent;
std::list<Vector3> path;
int source;
int target;
double pathfindingRuntime;
double pathLength;
std::vector<Vector3> visitedVertices;
int lastVertex;



Vector3 projectToTangentPlane(const Vertex& v, const Vector3& p)
{
	const Vector3& vp = p - geometry->inputVertexPositions[v];
	Vector3 normal = Vector3::zero();
	geometry->requireVertexNormals();

	//for (const Face& f : v.adjacentFaces())
	//{
	//	normal += geometry->faceNormal(f);
	//}

	//normal = normal.normalize();

	normal = geometry->vertexNormals[v];

	Vector3 offset = vp - (dot(normal, vp) * normal);
	return geometry->inputVertexPositions[v] + offset;

}



//Heuristics for A*
double EucledianHeuristic(const Vector3& v1, const Vector3& v2)
{
	return (v1 - v2).norm();
}


double ManhattanHeuristic(const Vector3& vN, const Vector3& vT)
{
	double dx = fabs(vN.x - vT.x);
	double dy = fabs(vN.y - vT.y);
	double dz = fabs(vN.z - vT.z);

	const Vector3& vToT = vT - geometry->inputVertexPositions[lastVertex];
	const Vector3& vToN = vN - geometry->inputVertexPositions[lastVertex];
	double c = cross(vToT, vToN).norm();

	return c * (vN - vT).norm();
	//return (dx+dy+dz) + sqrt(5) * std::max({ dx, dy, dz });
}



class VertexComparator
{
public:
	//Comparator, that will form min heap wrt eucledian distance
	bool operator() (const std::pair<size_t, double>& lhs, const std::pair<size_t,
		double>& rhs) const
	{
		return lhs.second > rhs.second;
	}
};


void createGraph(const std::unique_ptr<ManifoldSurfaceMesh>& mesh, const std::unique_ptr<VertexPositionGeometry>& geometry)
{
	graph.resize(mesh->nVertices());
	//Traverse all the faces
	for (const Vertex& v : mesh->vertices())
	{
		for (const Vertex& n : v.adjacentVertices())
		{
			graph[v.getIndex()].insert({ n.getIndex(), (geometry->inputVertexPositions[v] - geometry->inputVertexPositions[n]).norm() });
		}
	}
}

//Dijkstra stops upon reaching to the target now. Since no need to find all the distances from source to all other vertices
void dijkstraPQ(int source, int target)
{
	//Could be initalized outside
	visitedVertices.clear();
	std::vector<double> dist(graph.size(), std::numeric_limits<double>::max());
	parent.resize(graph.size(), -1);
	boost::heap::priority_queue<std::pair<size_t, double>, boost::heap::compare<VertexComparator>> pq;
	//FOR BENCHMARK (initializations are skipped for benchmark)
	double t1 = glfwGetTime();


	pq.push({ source, 0.0 });
	dist[source] = 0.0;
	parent[source] = source;

	while (!pq.empty())
	{
		std::pair<size_t, double> vertex = pq.top();
		pq.pop();

		//Get the neighbours
		const std::set<std::pair<size_t, double>>& neighbours = graph[vertex.first];

		visitedVertices.push_back(geometry->inputVertexPositions[vertex.first]);

		//If the target is being explored that means the shortest distance is already found just quit.
		if (vertex.first == target)
		{
			break;
		}

		//Explore the neighbours
		for (const auto& n : neighbours)
		{
			//Relax condition
			if (dist[vertex.first] + n.second < dist[n.first])
			{
				dist[n.first] = dist[vertex.first] + n.second;
				parent[n.first] = vertex.first;
				pq.push({ n.first, dist[n.first] });
			}
		}
	}

	pathfindingRuntime = glfwGetTime() - t1;

}


void Astar(int source, int target, std::function<double(const Vector3&, const Vector3&)> heuristic)
{
	//Could be initalized outside
	visitedVertices.clear();
	std::vector<double> dist(graph.size(), std::numeric_limits<double>::max());
	parent.resize(graph.size(), -1);
	boost::heap::priority_queue<std::pair<size_t, double>, boost::heap::compare<VertexComparator>> pq;
	//FOR BENCHMARK (initializations are skipped for benchmark)
	double t1 = glfwGetTime();


	pq.push({ source, 0.0 });
	dist[source] = 0.0;
	parent[source] = source;

	while (!pq.empty())
	{
		std::pair<size_t, double> vertex = pq.top();
		pq.pop();
		//Get the neighbours
		const std::set<std::pair<size_t, double>>& neighbours = graph[vertex.first];

		visitedVertices.push_back(geometry->inputVertexPositions[vertex.first]);

		lastVertex = vertex.first;

		//If the target is being explored that means the shortest distance is already found just quit.
		if (vertex.first == target)
		{
			break;
		}

		//Explore the neighbours
		for (const auto& n : neighbours)
		{
			double newCost = dist[vertex.first] + n.second;
			//Relax condition
			if (newCost < dist[n.first])
			{
				dist[n.first] = newCost;
				parent[n.first] = vertex.first;
				double h = heuristic(geometry->inputVertexPositions[n.first], geometry->inputVertexPositions[target]);
				pq.push({ n.first, (dist[n.first] + h) });
			}
		}
	}

	pathfindingRuntime = glfwGetTime() - t1;

}




void retrievePath()
{
	pathLength = 0;
	std::stack<int> st;
	int u = target;
	while (u != source)
	{
		st.push(u);
		u = parent[u];
	}
	st.push(source);
	
	Vector3 lastVertex = geometry->inputVertexPositions[st.top()];
	while (!st.empty())
	{
		path.push_back({ geometry->inputVertexPositions[st.top()] });
		pathLength += (geometry->inputVertexPositions[st.top()] - lastVertex).norm();
		lastVertex = geometry->inputVertexPositions[st.top()];
		st.pop();
	}
}


//NOTE FOR MYSELF ABOUT GEOMETRY CENTRAL API:
//A mesh consists of two parts in Geometry Central. 
//SurfaceMesh: Topology  and  Geometry: Properties like positions, lengths, areas, curvatures etc.
//Whenever you add a new vertex or form a new connection via SurfaceMesh you only update its topology no 
//Position data is specified. So whenever you update the topology you should also update its geometry
//And always a SurfaceMesh and Geometry is tie'd together via std::tie. Therefore, geometry will be 
//aware of any update on Topology. Hence, after updating the topology without worrying about internal
//containers you can direclty update mesh.



//SUBDIVISION PART
//Returns the newly created Halfedge whose tail vertex is connected to the vertices of the incident face(s).
//This function only updates topology caller must set its position. (by default vertex is at origin)
Halfedge splitEdge(const Edge& e)
{
	//insertVertexAlongEdge just adds a vertex alligned with the direction of the e.halfedge(), nothing fantastic
	//at he.tailVertex() our newly created vertex lies, but topology is not adequate yet.
	const Halfedge& he = mesh->insertVertexAlongEdge(e);
	//Connect the vertex to the prime face's opposite vertex.
	const Halfedge& primeOpposite = he.next().next();
	mesh->connectVertices(he, primeOpposite);

	//If not boundary, there must be a secondary incident face
	if (he.twin().isInterior())
	{
		//Even though the tailvertex is the same to preserve orientation we need to use twin of the he.
		const Halfedge& heReverse = he.twin().next();
		const Halfedge& secondaryOpposite = heReverse.next().next();
		mesh->connectVertices(heReverse, secondaryOpposite);
	}

	return he;
}




void ftvSubdivision()
{
	//Keep track of original vertices and edges. It will be beneficial in two ways:
	//1) Since while we are looping we are updating the mesh (adding new edges vertices etc.) and newly created geometry should not be processed
	//2) These will help us to keep track of which edges should be flipped.

	//Geometry Central Data Containers are adaptive to mutations, meaning that even the topology changes these data containers
	//Adapt themselves and stays valid. (Resizes themselves etc.)
	VertexData<bool> isOriginalVertex(*mesh, true);
	EdgeData<bool> isOriginalEdge(*mesh, true);
	std::vector<Edge> edgesToFlip;

	for (const Edge& e : mesh->edges())
	{
		if (!isOriginalEdge[e]) continue; //Do not process newly created edges

		const Vertex& oldU = e.halfedge().tailVertex();
		const Vertex& oldV = e.halfedge().tipVertex();
		const Vector3& oldUPos = geometry->inputVertexPositions[oldU];
		const Vector3& oldVPos = geometry->inputVertexPositions[oldV];

		//Split the edge
		const Vertex& newVertex = splitEdge(e).vertex();
		//Mark the vertex as new
		isOriginalVertex[newVertex] = false;

		//Splitting only changes topology geometry is not initalized. Initalize the geometry
		const Vector3& newPos = 0.5 * (oldUPos + oldVPos);
		geometry->inputVertexPositions[newVertex] = newPos;

		//Traverse through adjacent edges of the new vertex
		for (const Edge& nEdge : newVertex.adjacentEdges())
		{
			//When split an edge the corresponding edge splitted into two, and also it is connected to two opposite vertices.
			//Thus resulting new edges (if interior it is 4, if boundary it is 3), are new edges.
			//Mark Them
			isOriginalEdge[nEdge] = false;
			const Vertex& otherVertex = nEdge.otherVertex(newVertex);

			//Mark the edges that forms crosses inside the triangles since they deform the shape
			//We need to flip them. No need to mark already existing edges.
			if (isOriginalVertex[otherVertex] && (otherVertex != oldU) && (otherVertex != oldV))
			{
				edgesToFlip.push_back(nEdge);
			}

		}
	}

	//After splitting as I stated before, there are crosses (direct lines between opposite vertices)
	//Get rid of them by flipping those edges
	for (const Edge& e : edgesToFlip)
	{
		mesh->flip(e);
	}


	polyscope::registerSurfaceMesh("4 to 1 Subdivision", geometry->inputVertexPositions, mesh->getFaceVertexList());


}


//Subject to change, not sure this is a nice "midpoint"
Vector3 faceMidpoint(const Face& f)
{
	std::vector<std::reference_wrapper<Vector3>> vertices;
	for (const Vertex& v : f.adjacentVertices())
	{
		vertices.push_back(geometry->inputVertexPositions[v]);
	}
	
	Vector3 mid; //This is not mid but centroid.
	//mid = vertices[0].get() + 0.5 * (vertices[1].get() - vertices[0].get());
	//mid += (1.0 / 3.0) * (vertices[2].get() - mid);
	mid = (1.0 / 3.0) * (vertices[0].get() + vertices[1].get() + vertices[2].get());
	return mid;
}

void stSubdivision()
{
	EdgeData<bool> isOriginalEdge(*mesh, true);
	VertexData<bool> isOriginalVertex(*mesh, true);
	//STEP1
	//Form new points inside the faces.
	for (const Face& f : mesh->faces())
	{
		const Vector3& midPoint = faceMidpoint(f);;
		//It is important to gather the inner vertex first since when you form a new vertex
		//inside the face topology changes and the resulting point will be deformed.
		const Vertex& newVertex = mesh->insertVertex(f);
		isOriginalVertex[newVertex] = false;
		geometry->inputVertexPositions[newVertex] = midPoint;
		for (const Edge& e : newVertex.adjacentEdges())
		{
			isOriginalEdge[e] = false;
		}
	}
	
	//STEP2
	//Shift the original vertices.
	for (const Vertex& vert : mesh->vertices())
	{
		if (!isOriginalVertex[vert]) continue;

		double n = vert.degree();
		double an = (4 - 2 * cos((2.0 * PI) / n)) / 9.0;

		Vector3 newPos = (1.0 - an) * geometry->inputVertexPositions[vert];
		Vector3 offset = Vector3::zero();
		for (const Vertex& v : vert.adjacentVertices())
		{
			offset += geometry->inputVertexPositions[v];
		}

		offset *= an / n;
		newPos += offset;

		geometry->inputVertexPositions[vert] = newPos;

	}

	for (const Edge& e : mesh->edges())
	{
		if (!isOriginalEdge[e]) continue;

		mesh->flip(e);
	}


	polyscope::registerSurfaceMesh("Sqrt3 Subdivision", geometry->inputVertexPositions, mesh->getFaceVertexList());

	
}


Vector3 computePStar(const Face& f, double shapeFactor)
{
	std::vector<std::reference_wrapper<Vector3>> vertices;
	std::vector<Vector3> tangentPlaneVertices;
	for (const Vertex& v : f.adjacentVertices())
	{
		vertices.push_back(geometry->inputVertexPositions[v]);
	}
	//Obtain p by barycentric interpolation (For now I could not figure out what should be u and v) 
	Vector3 p = (1.0 / 3.0) * (vertices[0].get() + vertices[1].get() + vertices[2].get());

	//Project every point to the tangent planes then interpolate with same u and v
	for (const Vertex& v : f.adjacentVertices())
	{
		tangentPlaneVertices.push_back(projectToTangentPlane(v, p));
	}
	
	//Interpolate back 
	Vector3 primePStar = (1.0 / 3.0) * (tangentPlaneVertices[0] + tangentPlaneVertices[1] + tangentPlaneVertices[2]);

	//return primePStar;

	return (1 - shapeFactor) * p + shapeFactor * primePStar;


}

void phongTesellation(double shapeFactor)
{
	EdgeData<bool> isOriginalEdge(*mesh, true);

	for (const Face& f : mesh->faces())
	{
		const Vector3& pStar = computePStar(f, shapeFactor);
		const Vertex& newVertex = mesh->insertVertex(f);
		geometry->inputVertexPositions[newVertex] = pStar;
		for (const Edge& e : newVertex.adjacentEdges())
		{
			isOriginalEdge[e] = false;
		}
	}

	for (const Edge& e : mesh->edges())
	{
		if (!isOriginalEdge[e]) continue;

		mesh->flip(e);
	}

	polyscope::registerSurfaceMesh("Phong Subdivision", geometry->inputVertexPositions, mesh->getFaceVertexList());

}




void callback()
{
	//Graph Construction Runtime
	if (ImGui::Button("Create Graph"))
	{
		graph.clear();
		createGraph(mesh, geometry);
	}

	//Dijkstra Callback
	ImGui::InputInt("Source", &source);
	if (ImGui::Button("Dijkstra"))
	{
		dijkstraPQ(source, target);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		std::vector<std::array<size_t, 2>> ind;
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Visited Vertices", visitedVertices, ind);
	}
	
	//A* Callback
	if (ImGui::Button("A*Eucledian"))
	{
		Astar(source, target, EucledianHeuristic);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		std::vector<std::array<size_t, 2>> ind;
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Visited Vertices", visitedVertices, ind);
	}

	if (ImGui::Button("A*Manhattan"))
	{
		Astar(source, target, ManhattanHeuristic);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		std::vector<std::array<size_t, 2>> ind;
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Visited Vertices", visitedVertices, ind);
	}

	//Show Dijkstra Runtime 
	ImGui::SameLine(200);
	ImGui::Text("Total Time Taken: %f", pathfindingRuntime);

	//Retrieve Path to Target
	ImGui::InputInt("Target", &target);
	if (ImGui::Button("Retrieve Path"))
	{
		path.clear();
		retrievePath();
		std::vector<std::array<size_t, 2>> ind;
		for (size_t i = 0; i < path.size() - 1; ++i)
		{
			ind.push_back({ i, i + 1 });
		}
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Path", path, ind);
	}
	//Show PathLength
	ImGui::SameLine(200);
	ImGui::Text("Path Length is: %f", pathLength);

}


int main()
{
	std::tie(mesh, geometry) = readManifoldSurfaceMesh("Models/horse0m.off");
	
	
	
	polyscope::init();
	polyscope::registerSurfaceMesh("mesh", geometry->inputVertexPositions, mesh->getFaceVertexList());

	//ftvSubdivision();

	//for (size_t i = 0; i < 1; ++i)
	//{
	//	stSubdivision();
	//}

	for (size_t i = 0; i < 1; ++i)
	{
		phongTesellation(0.75);
	}


	polyscope::state::userCallback = callback;


	polyscope::show(); // pass control to the gui until the user exits


	return 0;
}