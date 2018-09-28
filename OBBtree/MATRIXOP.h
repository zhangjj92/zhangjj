#pragma once
#pragma execution_character_set("utf-8")
#ifndef MATRIXOP_H
#define MATRIXOP_H
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Eigen>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include<cmath>
#include<algorithm>
#define pi 3.1415926
using namespace std;
using namespace Eigen;

class MATRIXOP {
public:
	// 计算协方差矩阵
	Eigen::Matrix3f CalculateCovariance(Eigen::MatrixXf pVertices);

	//特征值//特征向量//排序//schmidt正交
	Eigen::Matrix3f CalculateEigenvectors(Eigen::Matrix3f covariance);

	//schmidt正交的单步计算模板
	Eigen::Vector3f schmidt_template(Eigen::Vector3f vector0, Eigen::Vector3f vector1);

	//一维容器转换为矩阵
	Eigen::MatrixXf vector_trans_matrix(vector<float> vertex);

	//glm矩阵变换-缩放矩阵//size缩放尺寸,定义第四个元素为1.0
	Eigen::Matrix4f scal_matrix(Eigen::Vector4f & size);

	//glm矩阵变换-位移矩阵//size移动大小,定义第四个元素为1.0
	Eigen::Matrix4f translate_matrix(Eigen::Vector4f & size);

	//glm矩阵变换-旋转矩阵//字符ch决定旋转轴
	Eigen::Matrix4f rotate_matrix(float angle, char ch);
};
#endif
