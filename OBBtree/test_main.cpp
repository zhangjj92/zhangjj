#pragma once
#pragma execution_character_set("utf-8")
#include <iostream>
#include <vector>
#include<string>
#include<OBBtree.cpp>
#include<load.h>
#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;


int main()
{
	//obb
	vector<float> vertices;
	LoadModel load;
	OBB obb;
	load = LoadModel("E:/OpenGL/resources/model/11.obj", vertices);
	MatrixXf ver;
	ver = obb.set_Vertexs(vertices);
	cout << "ver\n" << ver << endl;
	cout << "\n\n\n" << endl;
	cout << "obb.get_Centre()" << obb.get_Centre() << endl;
	cout << "obb.get_HalfExtents()" << obb.get_HalfExtents() << endl;
	cout << "obb.get_RotationMatrix()" << obb.get_RotationMatrix() << endl;
	obb.generateOBB(ver, obb.get_Centre());
	cout << "\n\n\n" << endl;
	obb.printobb();
	cout << "obb.get_Centre()" << obb.get_Centre() << endl;
	cout << "obb.get_HalfExtents()" << obb.get_HalfExtents() << endl;
	cout << "obb.get_RotationMatrix()" << obb.get_RotationMatrix() << endl;
	cout << "\n\n\n" << endl;
	//TREE
	OBBtreeNode tree, *root;
	tree.CreateRoot(obb, ver, obb.get_Centre());
	root = tree.get_root();
	tree.CreateTree(obb, ver);
	tree.preTraverse(root, 0);
	tree.get_rightChild(root);
	OBBtreeNode *l;
	l = tree.get_rightChild(root);
	cout << "\n\n" << endl;
	//点是否在obb盒中
	Eigen::RowVector3f proPoint;
	proPoint << 0, 0, -1;
	obb.point_in_OBB(proPoint, obb);
	cout << "\n\n" << endl;
	//两盒是否碰撞
	Eigen::Vector3f centre1;
	centre1 << 0, 0, 0;
	Eigen::Vector3f centre2;
	centre2 << 1, 0, 0;
	Eigen::Matrix3f rotationMatrix;
	rotationMatrix << 1, 0, 0,
		0, 1, 0,
		0, 0, 1;
	Eigen::Vector3f halfExtents1;
	Eigen::Vector3f halfExtents2;
	halfExtents1 << 1, 1, 1;
	halfExtents2 << 1, 1, 2;
	OBB c;
	OBB box0(centre1, rotationMatrix, halfExtents1), box1(centre2, rotationMatrix, halfExtents2);
	c.intersects(box0, box1);


}










