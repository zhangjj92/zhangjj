#pragma once
#pragma execution_character_set("utf-8")
#include<MATRIXOP.h>

//协方差矩阵
Eigen::Matrix3f MATRIXOP::CalculateCovariance(Eigen::MatrixXf pVertices)
{
	Eigen::Matrix3f covariance;
	//求取列向量均值
	Eigen::MatrixXf meanVec = pVertices.colwise().mean();
	//将列向量均值从MatrixXf 转换为行向量 RowVectorXf
	Eigen::RowVectorXf meanVecRow(Eigen::RowVectorXf::Map(meanVec.data(), pVertices.cols()));
	//求取上述的零均值列向量矩阵
	Eigen::MatrixXf zeroMeanMat = pVertices;
	//计算协方差矩阵
	zeroMeanMat.rowwise() -= meanVecRow;
	if (pVertices.rows() == 1)
		covariance = (zeroMeanMat.adjoint()*zeroMeanMat) / static_cast<float>(pVertices.rows());
	else
		covariance = (zeroMeanMat.adjoint()*zeroMeanMat) / static_cast<float>(pVertices.rows() - 1);
	return covariance;
}

//特征值//特征向量//排序(按照特征值从大到小排序)//schmidt正交
Eigen::Matrix3f MATRIXOP::CalculateEigenvectors(Eigen::Matrix3f covariance)
{
	EigenSolver<Matrix3f> eig(covariance);
	//特征值
	Eigen::Matrix3f eigenvalue = eig.pseudoEigenvalueMatrix();
	//特征向量
	Eigen::Matrix3f eigenvectors = eig.pseudoEigenvectors();
	//对特征值进行从大到小排序
	Eigen::RowVector3f eigenvalueVec = eigenvalue.colwise().sum();
	unsigned int index[3] = { 0,1,2 };
	for (int i = 0; i< 2; i++)
		for (int j = 0;j< 2 - i;j++)
			if (eigenvalueVec[j] < eigenvalueVec[j + 1])
			{
				float temp = eigenvalueVec[j];
				eigenvalueVec[j] = eigenvalueVec[j + 1];
				eigenvalueVec[j + 1] = temp;
				int temp_index = index[j];
				index[j] = index[j + 1];
				index[j + 1] = temp_index;

			}
	Matrix3f Sorted_eigenvector;
	Matrix3f Sorted_eigenvalue(eigenvalueVec.asDiagonal());
	for (unsigned i = 0; i < 3; i++)
	{
		Sorted_eigenvector.col(i) = eigenvectors.col(index[i]);
	}
	//schmidt正交
	Matrix3f schmidt_eigenvector;
	for (unsigned int i = 0; i < 3; i++)
	{
		schmidt_eigenvector.col(i) = Sorted_eigenvector.col(i);
		for (unsigned int j = 1; j <= i; j++)
		{
			schmidt_eigenvector.col(i) -= schmidt_template(Sorted_eigenvector.col(i), schmidt_eigenvector.col(j - 1));
		}
		schmidt_eigenvector.col(i).normalized();
	}
	return schmidt_eigenvector;
}

//schmidt正交的单步计算模板
Eigen::Vector3f MATRIXOP::schmidt_template(Eigen::Vector3f vector0, Eigen::Vector3f vector1)
{
	Vector3f schmidt_vector = vector0.dot(vector1) / vector1.squaredNorm() * vector1;
	return schmidt_vector;
}

//一维容器转换为矩阵
Eigen::MatrixXf MATRIXOP::vector_trans_matrix(vector<float> vertex)
{
	vector<Eigen::RowVector3f> Vertices;
	Eigen::RowVector3f rowvector;
	for (unsigned int i = 0; i < vertex.size() / 3; i++)
	{
		for (unsigned int j = 0; j < 3; j++)
		{
			rowvector(j) = vertex[3 * i + j];
		}
		Vertices.push_back(rowvector);
	}
	unsigned int numVertices = Vertices.size();
	Eigen::MatrixXf pVertices(numVertices, 3);
	for (unsigned int i = 0; i < numVertices; i++)
	{
		pVertices.row(i) = Vertices[i];
	}
	return pVertices;
}

//glm对应矩阵变换-缩放矩阵//size缩放尺寸,一般定义第四个元素为1.0
Eigen::Matrix4f MATRIXOP::scal_matrix(Eigen::Vector4f & size)
{
	Eigen::Matrix4f scalMatrix;
	scalMatrix.setIdentity(4, 4);
	scalMatrix.diagonal() = size;
	return scalMatrix;
}

//glm对应矩阵变换-位移矩阵//size移动大小,定义第四个元素为1.0
Eigen::Matrix4f MATRIXOP::translate_matrix(Eigen::Vector4f & size)
{
	Eigen::Matrix4f translateMatrix;
	translateMatrix.setIdentity(4, 4);
	translateMatrix.col(3) = size;
	return  translateMatrix;

}

//glm对应矩阵变换-旋转矩阵//用字符ch确定旋转轴
//ch的值有三种，x代表以x轴旋转， y以y轴旋转，z是以z轴旋转
//fov为角度
Eigen::Matrix4f MATRIXOP::rotate_matrix(float angle, char ch)
{
	float radian = pi / 180.0 * angle;
	Eigen::Matrix4f rotateMatrix;
	rotateMatrix.setIdentity(4, 4);
	switch (ch)
	{
	//x轴为旋转轴
	case 'x':
		rotateMatrix(1, 1) = cos(radian);
		rotateMatrix(1, 2) = -sin(radian);
		rotateMatrix(2, 1) = sin(radian);
		rotateMatrix(2, 2) = cos(radian);
		break;
	//y轴为旋转轴
	case 'y':
		rotateMatrix(0, 0) = cos(radian);
		rotateMatrix(0, 2) = sin(radian);
		rotateMatrix(2, 0) = -sin(radian);
		rotateMatrix(2, 2) = cos(radian);
		break;
	//z轴为旋转轴
	case 'z':
		rotateMatrix(0, 0) = cos(radian);
		rotateMatrix(0, 1) = -sin(radian);
		rotateMatrix(1, 0) = sin(radian);
		rotateMatrix(1, 1) = cos(radian);
		break;
	default:
		cout << "MATRIXOP::ERROR INPUT! Please choose x/y/z which represents X/Y/Z - axis.\n" << endl;
		break;
	}
	return rotateMatrix;
}
