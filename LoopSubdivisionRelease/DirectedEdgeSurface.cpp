///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

using namespace  std;
// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face " << face << " ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << " " << otherHalf[dirEdge] << std::endl;
    } // WriteObjectStream()


void DirectedEdgeSurface::LoopSubdivision()
{
    unsigned int i,j;
    map<unsigned int, unsigned int> edge;
    map<unsigned int, unsigned int> new_edge;
    map<unsigned int, unsigned int> edge_match_vertex;

    //volume before operation
    unsigned int old_v_size = vertices.size();
    unsigned int old_face_size = faceVertices.size()/3;

    //create new points & normals
    for (i = 0; i < otherHalf.size(); i++)
    {
        unsigned int mod = i % 3;
        unsigned int belong = i / 3.0f;

        edge[i] = otherHalf[i];
        if (i != edge[edge[i]] || i == 0)
        {
            if (mod == 0)
            {
                unsigned int otherhalf_index = otherHalf[i];
                unsigned int otherhalf_mod = otherhalf_index % 3;
                unsigned int otherhalf_belong = otherhalf_index / 3;

                unsigned int otherVertex;
                if (otherhalf_mod == 0) otherVertex = faceVertices[otherhalf_belong * 3 + 1];
                else if (otherhalf_mod == 1) otherVertex = faceVertices[otherhalf_belong * 3 + 2];
                else if (otherhalf_mod == 2) otherVertex = faceVertices[otherhalf_belong * 3];


                Cartesian3 newPoint_first;
                newPoint_first = (3.0f / 8.0f) * (vertices[faceVertices[belong * 3]] + vertices[faceVertices[belong * 3 + 2]])
                    + (1.0f / 8.0f) * (vertices[faceVertices[belong * 3 + 1]] + vertices[otherVertex]);

                Cartesian3 newNormal_first;
                newNormal_first = (3.0f / 8.0f) * (normals[faceVertices[belong * 3]] + normals[faceVertices[belong * 3 + 2]])
                    + (1 / 8) * (normals[faceVertices[belong * 3 + 1]] + normals[otherVertex]);

                vertices.push_back(newPoint_first);
                normals.push_back(newNormal_first);

                edge_match_vertex[i] = vertices.size() - 1;
            }
            else if (mod == 1)
            {
                unsigned int otherhalf_index = otherHalf[i];
                unsigned int otherhalf_mod = otherhalf_index % 3;
                unsigned int otherhalf_belong = otherhalf_index / 3;

                unsigned int otherVertex;
                if (otherhalf_mod == 0) otherVertex = faceVertices[otherhalf_belong * 3 + 1];
                else if (otherhalf_mod == 1) otherVertex = faceVertices[otherhalf_belong * 3 + 2];
                else if (otherhalf_mod == 2) otherVertex = faceVertices[otherhalf_belong * 3];


                Cartesian3 newPoint_second;
                newPoint_second = (3.0f / 8.0f) * (vertices[faceVertices[belong * 3]] + vertices[faceVertices[belong * 3 + 1]])
                    + (1.0f / 8.0f) * (vertices[faceVertices[belong * 3 + 2]] + vertices[otherVertex]);

                Cartesian3 newNormal_second;
                newNormal_second = (3.0f / 8.0f) * (normals[faceVertices[belong * 3]] + normals[faceVertices[belong * 3 + 1]])
                    + (1 / 8) * (normals[faceVertices[belong * 3 + 2]] + normals[otherVertex]);

                vertices.push_back(newPoint_second);
                normals.push_back(newNormal_second);

                edge_match_vertex[i] = vertices.size() - 1;
            }
            else if (mod == 2)
            {
                unsigned int otherhalf_index = otherHalf[i];
                unsigned int otherhalf_mod = otherhalf_index % 3;
                unsigned int otherhalf_belong = otherhalf_index / 3;

                unsigned int otherVertex;
                if (otherhalf_mod == 0) otherVertex = faceVertices[otherhalf_belong * 3 + 1];
                else if (otherhalf_mod == 1) otherVertex = faceVertices[otherhalf_belong * 3 + 2];
                else if (otherhalf_mod == 2) otherVertex = faceVertices[otherhalf_belong * 3];


                Cartesian3 newPoint_third;
                newPoint_third = (3.0f / 8.0f) * (vertices[faceVertices[belong * 3 + 1]] + vertices[faceVertices[belong * 3 + 2]])
                    + (1.0f / 8.0f) * (vertices[faceVertices[belong * 3]] + vertices[otherVertex]);

                Cartesian3 newNormal_third;
                newNormal_third = (3.0f / 8.0f) * (normals[faceVertices[belong * 3 + 1]] + normals[faceVertices[belong * 3 + 2]])
                    + (1 / 8) * (normals[faceVertices[belong * 3]] + normals[otherVertex]);

                vertices.push_back(newPoint_third);
                normals.push_back(newNormal_third);

                edge_match_vertex[i] = vertices.size() - 1;
            }
        }
        else if(i == edge[edge[i]])
        {
            edge_match_vertex[i] = edge_match_vertex[edge[i]];
        }
    }

    //create new faces
    for (i = 0 ; i < old_face_size ; i++)
    {
        map<unsigned int, unsigned int>::iterator it0 = edge_match_vertex.find(3 * i);
        map<unsigned int, unsigned int>::iterator it1 = edge_match_vertex.find(3 * i + 1);
        map<unsigned int, unsigned int>::iterator it2 = edge_match_vertex.find(3 * i + 2);

        faceVertices.push_back(faceVertices[3 * i]);
        faceVertices.push_back(it0->second);
        faceVertices.push_back(it1->second);

        faceVertices.push_back(faceVertices[3 * i + 1]);
        faceVertices.push_back(it1->second);
        faceVertices.push_back(it2->second);

        faceVertices.push_back(faceVertices[3 * i + 2]);
        faceVertices.push_back(it0->second);
        faceVertices.push_back(it2->second);

        faceVertices[3 * i] = it0->second;
        faceVertices[3 * i + 1] = it1->second;
        faceVertices[3 * i + 2] = it2->second;
    }

    //caculate first dirededge
    vector<vector<int>> firstedges;
    firstDirectedEdge.resize(vertices.size());
    for (i = 0; i < vertices.size(); i++)
    {
        vector<int> temp;
        for (j = 0; j < faceVertices.size()/3; j++)
        {
            if (faceVertices[3 * j] == i)
            {
                temp.push_back(faceVertices[3 * j + 1]);
            }
            else if (faceVertices[3 * j + 1] == i)
            {
                temp.push_back(faceVertices[3 * j + 2]);
            }
            else if (faceVertices[3 * j + 2] == i)
            {
                temp.push_back(faceVertices[3 * j]);
            }
        }
        firstedges.push_back(temp);
    }
    //save firstedges
    for (i = 0; i < firstedges.size(); i++)
    {
        firstDirectedEdge[i] = firstedges[i][0];
    }

    //sum otherHalf
    otherHalf.resize(faceVertices.size());
    vector<Edge> e;
    for (i = 0; i < faceVertices.size()/3; i++)
    {
        Edge temp_e, temp_e1, temp_e2;
        temp_e.first_index = faceVertices[3 * i + 2];
        temp_e.second_index = faceVertices[3 * i];

        temp_e1.first_index = faceVertices[3 * i];
        temp_e1.second_index = faceVertices[3 * i + 1];

        temp_e2.first_index = faceVertices[3 * i + 1];
        temp_e2.second_index = faceVertices[3 * i + 2];

        e.push_back(temp_e);
        e.push_back(temp_e1);
        e.push_back(temp_e2);
    }

    for (i = 0; i < e.size(); i++)
    {
        for (j = 0; j < e.size(); j++)
        {
            if (i != j && ((e[i].second_index == e[j].first_index && e[i].first_index == e[j].second_index) || (e[i].second_index == e[j].second_index && e[i].first_index == e[j].first_index)))
            {
                otherHalf[i] = j;
            }
        }
    }

    //find each vertex's degree
    vector<int> degree;
    vector<vector<int>> neighbors;
    for (i = 0; i < old_v_size; i++)
    {
        int time = 0;
        vector<int> temp_nei;
        for (j = 0; j < faceVertices.size() / 3; j++)
        {
            if (faceVertices[3 * j] == i)
            {
                temp_nei.push_back(faceVertices[3 * j + 1]);
                temp_nei.push_back(faceVertices[3 * j + 2]);
                time++;
            }
            else if (faceVertices[3 * j + 1] == i)
            {
                temp_nei.push_back(faceVertices[3 * j + 2]);
                temp_nei.push_back(faceVertices[3 * j]);
                time++;
            }
            else if (faceVertices[3 * j + 2] == i)
            {
                temp_nei.push_back(faceVertices[3 * j]);
                temp_nei.push_back(faceVertices[3 * j + 1]);
                time++;
            }
        }

        degree.push_back(time);
        neighbors.push_back(temp_nei);
    }

    //remove same points
    for (int k = 0; k < neighbors.size(); k++)
    {
        for (i = 0; i < neighbors[k].size() - 1; i++)
        {
            for (j = i + 1; j < neighbors[k].size(); j++)
            {
                if (neighbors[k][j] == neighbors[k][i])
                {
                    int temp = neighbors[k][neighbors[k].size() - 1];
                    neighbors[k][neighbors[k].size() - 1] = neighbors[k][j];
                    neighbors[k][j] = temp;
                    neighbors[k].pop_back();
                }
            }
        }
    }
    for (int k = 0; k < neighbors.size(); k++)
    {
        for (i = 0; i < neighbors[k].size() - 1; i++)
        {
            for (j = i + 1; j < neighbors[k].size(); j++)
            {
                if (neighbors[k][j] == neighbors[k][i])
                {
                    int temp = neighbors[k][neighbors[k].size() - 1];
                    neighbors[k][neighbors[k].size() - 1] = neighbors[k][j];
                    neighbors[k][j] = temp;
                    neighbors[k].pop_back();
                }
            }
        }
    }

    //modify old vertices
    for (i = 0; i < old_v_size; i++)
    {
        //caculate u
        float u;
        if (degree[i] == 3) u = 3.0f / 16.0f;
        else u = 3.0f / (8.0f * degree[i]);

        Cartesian3 sum_neighbour,sum_normal;
        for (j = 0; j < neighbors[i].size(); j++)
        {
            sum_neighbour = sum_neighbour + vertices[neighbors[i][j]];

            sum_normal = sum_normal + normals[neighbors[i][j]];
        }
        Cartesian3 modified_pos = (1 - degree[i] * u) * vertices[i] + u * sum_neighbour;

        Cartesian3 modified_nor = (1 - degree[i] * u) * normals[i] + u * sum_normal;

        vertices[i] = modified_pos;
    }

}

void DirectedEdgeSurface::LoopTime(int t)
{
    int time = 0;
    if(t == 0)
    {

    }
    else if(t == 1)
    {
        if(time == 0)
        {
            time++;
            LoopSubdivision();
        }
    }
    else if(t == 2)
    {
        if(time == 1)
        {
            time++;
            LoopSubdivision();
        }
    }
    else if(t == 3)
    {
        if(time == 2)
        {
            time++;
            LoopSubdivision();
        }
    }
}
// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()

