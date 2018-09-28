#pragma once
#pragma execution_character_set("utf-8")
#ifndef OBB_H
#define OBB_H

#include<MATRIXOP.cpp>
using namespace std;
using namespace Eigen;

//一个OBB包围盒需要：
                    //一个中心点
                    //一个旋转矩阵————点集合的三个坐标轴
                    //三个二分之一边长
struct OBB_data
{
	Eigen::RowVector3f Centre;
	Eigen::Matrix3f RotationMatrix;
	Eigen::RowVector3f HalfExtents;
};

class OBB {
private:
	OBB_data obb_data;
	//顶点：Eigen::MatrixXf pVertices(numVertices, 3);
	Eigen::MatrixXf Vertices;
public:	
	//初始化包围盒
	OBB();
	OBB(Eigen::Vector3f centre, Eigen::Matrix3f rotationMatrix, Eigen::Vector3f halfExtents);
	~OBB();
	//使用PCA ——协方差矩阵构造OBB包围盒
	//产生包围盒
	OBB_data generateOBB( Eigen::MatrixXf Vertices, RowVector3f Centre);

	//变换包围盒
	OBB_data transform(Eigen::Matrix4f & transMatrix);

	// 包围盒三个方向最大的方向, 用过中心的与该方向垂直的平面//判断点与面的法向是否同向
								  //大于等于0， 在法向指向的一侧或者面上
								  //小于0，在法向指向的另一侧
	bool toNormalDirection(Vector3f Normal, RowVector3f Centre, RowVector3f Vertice);

	//设置与第一主成分方向作为法向量一致方向的顶点集合
	vector<Eigen::MatrixXf> set_PCAVertices(Eigen::MatrixXf Vertices, Vector3f Normal, RowVector3f Centre);
	
	//按照第一主成分方向和中心点，获取顶点集合
	Eigen::MatrixXf get_Vertices();
	//获取中心点
	RowVector3f get_Centre();
	//获取二分之一边长
	RowVector3f get_HalfExtents();
	//获取旋转矩阵
	Matrix3f get_RotationMatrix();
	
	
	//设置中心点
	RowVector3f set_Centre(Eigen::MatrixXf Vertices, Matrix3f RotationMatrix, Eigen::RowVector3f Centre_0);
	//设置二分之一边长
	RowVector3f set_HalfExtents(Eigen::MatrixXf Vertices, Matrix3f RotationMatrix, RowVector3f Centre_0);
	//设置旋转矩阵
	Matrix3f set_RotationMatrix(Eigen::MatrixXf Vertices);
	//第一次输入顶点//根结点的顶点集
	MatrixXf set_Vertexs(vector<float> vertex);

	//打印OBB包围盒的中心点、旋转矩阵、三个二分之一边长
	void printobb();

	//判断点是否在包围盒内
	bool point_in_OBB(RowVector3f point, OBB obb);

	//检测两个包围盒是否碰撞
	bool intersects(OBB box0,OBB box1);
	
};
#endif // !OBB_H
