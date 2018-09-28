#ifndef LOAD_H
#define LOAD_H
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include<cmath>
#include<algorithm>
using namespace std;
class LoadModel {
public:
	int Vnum, Vnnum, Fnum;
	int VNum, VnNum, FNum;
	float sum_x, sum_y, sum_z;
	float	Max_z;
	LoadModel()
	{
	}

	LoadModel(vector<std::string>filename, vector<float> &vertices, vector<unsigned int> &indices)
	{ 
		int num = 0;
        Vnum = 0;
		Vnnum = 0;
		Fnum = 0;
		for (int i = 0; i < filename.size(); i++)
		{  
			string line;
			
			int vnum = 0;
			int vnnum = 0;
			int fnum = 0;
			vector<float> positions;
			vector<float> normals;

			//	float max_x = 0.0, max_y = 0.0, max_z = 0.0;
			sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
			ifstream inf;
			inf.open(filename[i], ios::in);
			if (!inf.is_open())
				cout << "Something went wrong when opening objfile." << endl;

			while (getline(inf, line))
			{
				vector<string>parameters;
				string ans = "";
				line = line.append(" ");
				for (int i = 0;i < line.length();i++)
				{
					if (line[i] != ' ')
						ans += line[i];
					else
					{
						parameters.push_back(ans);
						parameters.push_back(" ");
						ans = "";
					}
				}
				if (parameters[0] == "v")
				{
					float ver;
					vnum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						//cout << ver << endl;
						positions.push_back(ver);
					}
					float pos_x = atof(parameters[2].c_str());
					sum_x = sum_x + pos_x;
					float pos_y = atof(parameters[4].c_str());
					sum_y = sum_y + pos_y;
					float pos_z = atof(parameters[6].c_str());
					sum_z = sum_z + pos_z;
				}
				
				if (parameters[0] == "vn")
				{
					float ver;
					vnnum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						//cout << ver << endl;
						normals.push_back(ver);
					}
				}
				if (parameters[0] == "f")
				{
					unsigned int ver;
					fnum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						ver = ver - 1+ num;
						//cout << ver << endl;
						indices.push_back(ver);
					}
				}
				
			}
			num = num + vnum;
			Vnum = Vnum +vnum;
			Vnnum = Vnnum+ vnnum;
			Fnum= Fnum+fnum;
			//	cout << max_x << " " << max_y << " " << max_z << endl;
			for (int j = 0; j < vnum; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					vertices.push_back(positions[3 * j + i]);
				}
				for (int i = 0; i < 3; i++)
				{
					vertices.push_back(normals[3 * j + i]);
				}

			}
			inf.close();
		}
		}

	LoadModel(string filename, vector<float> &vertices, vector<unsigned int> &indices)
	{
			string line;
			VNum = 0;
			VnNum = 0;
			FNum = 0;
			vector<float> positions;
			vector<float> normals;
			vector<float> texCoords;
			float tex;
			Max_z = 0.0;
			sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
			ifstream inf;
			inf.open(filename, ios::in);
			if (!inf.is_open())
				cout << "Something went wrong when opening objfile." << endl;

			while (getline(inf, line))
			{
				vector<string>parameters;
				string ans = "";
				line = line.append(" ");
				for (int i = 0;i < line.length();i++)
				{
					if (line[i] != ' ')
						ans += line[i];
					else
					{
						parameters.push_back(ans);
						parameters.push_back(" ");
						ans = "";
					}
				}
				if (parameters[0] == "v")
				{
					float ver;
					VNum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						positions.push_back(ver);
					}			
					float pos_x = atof(parameters[2].c_str());
					sum_x = sum_x + pos_x;
					float pos_y = atof(parameters[4].c_str());
					sum_y = sum_y + pos_y;
					float pos_z = atof(parameters[6].c_str());
					sum_z = sum_z + pos_z;
					if (pos_z > Max_z)
						Max_z = pos_z;
				}

				if (parameters[0] == "vn")
				{
					float ver;
					VnNum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						normals.push_back(ver);
					}
				}
				if (parameters[0] == "f")
				{
					unsigned int ver;
					FNum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						ver = ver - 1;
						//cout << ver << endl;
						indices.push_back(ver);
					}
				}

			}
			//纹理坐标
			for (int i = 0; i < 3*VNum; i++)
			{
				if (i % 3!= 2) 
				{
					tex = positions[i] / 2.0 + 0.5;
					texCoords.push_back(tex);
				}
			
			}
			//切线和副切线
			vector<float> tangents(3*VNum, 0.0);
			vector<float> bitangents(3*VNum, 0.0);
			for (int i = 0; i < 3*FNum; i+=3)
			{
				float edge1_x = positions[3 * indices[i + 1] + 0] - positions[3 * indices[i] + 0];
				float edge1_y = positions[3 * indices[i + 1] + 1] - positions[3 * indices[i] + 1];
				float edge1_z = positions[3 * indices[i + 1] + 2] - positions[3 * indices[i] + 2];
				float uv1_x = texCoords[2 * indices[i + 1] + 0] - texCoords[2 * indices[i] + 0];
				float uv1_y = texCoords[2 * indices[i + 1] + 1] - texCoords[2 * indices[i] + 1];

				float edge2_x = positions[3 * indices[i + 2] + 0] - positions[3 * indices[i] + 0];
				float edge2_y = positions[3 * indices[i + 2] + 1] - positions[3 * indices[i] + 1];
				float edge2_z = positions[3 * indices[i + 2] + 2] - positions[3 * indices[i] + 2];
				float uv2_x = texCoords[2 * indices[i + 2] + 0] - texCoords[2 * indices[i] + 0];
				float uv2_y = texCoords[2 * indices[i + 2] + 1] - texCoords[2 * indices[i] + 1];

				float f = 1.0f / (uv1_x * uv2_y - uv2_x * uv1_y);
				tangents[3 * indices[i + 0] + 0] += f*(uv2_y*edge1_x - uv1_y*edge2_x);
				tangents[3 * indices[i + 1] + 0] += f*(uv2_y*edge1_x - uv1_y*edge2_x);
				tangents[3 * indices[i + 2] + 0] += f*(uv2_y*edge1_x - uv1_y*edge2_x);

				tangents[3 * indices[i + 0] + 1] += f*(uv2_y*edge1_y - uv1_y*edge2_y);
				tangents[3 * indices[i + 1] + 1] += f*(uv2_y*edge1_y - uv1_y*edge2_y);
				tangents[3 * indices[i + 2] + 1] += f*(uv2_y*edge1_y - uv1_y*edge2_y);

				tangents[3 * indices[i + 0] + 2] += f*(uv2_y*edge1_z - uv1_y*edge2_z);
				tangents[3 * indices[i + 1] + 2] += f*(uv2_y*edge1_z - uv1_y*edge2_z);
				tangents[3 * indices[i + 2] + 2] += f*(uv2_y*edge1_z - uv1_y*edge2_z);

				bitangents[3 * indices[i + 0] + 0] += f*(-uv2_x*edge1_x + uv1_x*edge2_x);
				bitangents[3 * indices[i + 1] + 0] += f*(-uv2_x*edge1_x + uv1_x*edge2_x);
				bitangents[3 * indices[i + 2] + 0] += f*(-uv2_x*edge1_x + uv1_x*edge2_x);

				bitangents[3 * indices[i + 0] + 1] += f*(-uv2_x*edge1_y + uv1_x*edge2_y);
				bitangents[3 * indices[i + 1] + 1] += f*(-uv2_x*edge1_y + uv1_x*edge2_y);
				bitangents[3 * indices[i + 2] + 1] += f*(-uv2_x*edge1_y + uv1_x*edge2_y);

			    bitangents[3 * indices[i + 0] + 2] += f*(-uv2_x*edge1_z + uv1_x*edge2_z);
				bitangents[3 * indices[i + 1] + 2] += f*(-uv2_x*edge1_z + uv1_x*edge2_z);
				bitangents[3 * indices[i + 2] + 2] += f*(-uv2_x*edge1_z + uv1_x*edge2_z);
			}
			for (int j = 0; j < VNum; j++)
			{
				for (int i = 0; i < 3; i++)
				{
					vertices.push_back(positions[3 * j + i]);
				}
				for (int i = 0; i < 3; i++)
				{
					vertices.push_back(normals[3 * j + i]);
				}
				//for (int i = 0; i < 2; i++)
				//{
				//vertices.push_back(texCoords[2 * j + i]);
				//}
			//	for (int i = 0; i < 3; i++)
			//	{
			//		vertices.push_back(tangents[3 * j + i]);
			//	}
			//	for (int i = 0; i < 3; i++)
			//	{
			//		vertices.push_back(bitangents[3 * j + i]);
			//	}

			}
			inf.close();
		}

	LoadModel(string filename, vector<float> &vertices)
	{
		string line;
		VNum = 0;
		VnNum = 0;
		FNum = 0;
		vector<float> positions;
		vector<float> normals;
		vector<unsigned int> indices;
		float tex;
		Max_z = 0.0;
		sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
		ifstream inf;
		inf.open(filename, ios::in);
		if (!inf.is_open())
			cout << "Something went wrong when opening objfile." << endl;

		while (getline(inf, line))
		{
			vector<string>parameters;
			string ans = "";
			line = line.append(" ");
			for (int i = 0;i < line.length();i++)
			{
				if (line[i] != ' ')
					ans += line[i];
				else
				{
					parameters.push_back(ans);
					parameters.push_back(" ");
					ans = "";
				}
			}
			if (parameters[0] == "v")
			{
				float ver;
				VNum++;
				for (int i = 2; i < 8; i = i + 2)
				{
					ver = atof(parameters[i].c_str());
					positions.push_back(ver);
				}
				float pos_x = atof(parameters[2].c_str());
				sum_x = sum_x + pos_x;
				float pos_y = atof(parameters[4].c_str());
				sum_y = sum_y + pos_y;
				float pos_z = atof(parameters[6].c_str());
				sum_z = sum_z + pos_z;
				if (pos_z > Max_z)
					Max_z = pos_z;
			}

			if (parameters[0] == "vn")
			{
				float ver;
				VnNum++;
				for (int i = 2; i < 8; i = i + 2)
				{
					ver = atof(parameters[i].c_str());
					normals.push_back(ver);
				}
			}
			if (parameters[0] == "f")
			{
				unsigned int ver;
				FNum++;
				for (int i = 2; i < 8; i = i + 2)
				{
					ver = atof(parameters[i].c_str());
					ver = ver - 1;
					//cout << ver << endl;
					indices.push_back(ver);
				}
			}			
		}
		inf.close();
		vertices = positions;
	}

	LoadModel(vector<std::string>filename, vector<float> &vertices)
	{
		int num = 0;
		Vnum = 0;
		Vnnum = 0;
		Fnum = 0;
		vector<float> positions;
		vector<float> normals;
		vector<unsigned int> indices;
		for (int i = 0; i < filename.size(); i++)
		{
			string line;

			int vnum = 0;
			int vnnum = 0;
			int fnum = 0;
			//	float max_x = 0.0, max_y = 0.0, max_z = 0.0;
			sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
			ifstream inf;
			inf.open(filename[i], ios::in);
			if (!inf.is_open())
				cout << "Something went wrong when opening objfile." << endl;

			while (getline(inf, line))
			{
				vector<string>parameters;
				string ans = "";
				line = line.append(" ");
				for (int i = 0;i < line.length();i++)
				{
					if (line[i] != ' ')
						ans += line[i];
					else
					{
						parameters.push_back(ans);
						parameters.push_back(" ");
						ans = "";
					}
				}
				if (parameters[0] == "v")
				{
					float ver;
					vnum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						//cout << ver << endl;
						positions.push_back(ver);
					}
					float pos_x = atof(parameters[2].c_str());
					sum_x = sum_x + pos_x;
					float pos_y = atof(parameters[4].c_str());
					sum_y = sum_y + pos_y;
					float pos_z = atof(parameters[6].c_str());
					sum_z = sum_z + pos_z;
				}

				if (parameters[0] == "vn")
				{
					float ver;
					vnnum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						//cout << ver << endl;
						normals.push_back(ver);
					}
				}
				if (parameters[0] == "f")
				{
					unsigned int ver;
					fnum++;
					for (int i = 2; i < 8; i = i + 2)
					{
						ver = atof(parameters[i].c_str());
						ver = ver - 1 + num;
						//cout << ver << endl;
						indices.push_back(ver);
					}
				}

			}
			
			inf.close();
		}

		vertices = positions;
	}

	
	
};
#endif
