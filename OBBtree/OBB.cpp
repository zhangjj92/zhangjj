#pragma once
#pragma execution_character_set("utf-8")
#include<OBB.h>
#include <Eigen/Geometry>
//初始化包围盒
OBB::OBB()  
{
	obb_data.Centre << 0.0f, 0.0f, 0.0f;
	obb_data.HalfExtents << 1.0f, 1.0f, 1.0f;
	obb_data.RotationMatrix = Eigen::Matrix3f::Identity();
}

OBB::OBB(Eigen::Vector3f centre, Eigen::Matrix3f rotationMatrix, Eigen::Vector3f halfExtents)
{
	obb_data.Centre = centre;
	obb_data.RotationMatrix = rotationMatrix;
	obb_data.HalfExtents = halfExtents;
}

OBB::~OBB()
{
}

//产生OBB包围盒//Centre_0父节点的中心点
OBB_data OBB::generateOBB(Eigen::MatrixXf Vertices, RowVector3f Centre_0)
{
	MATRIXOP matrixOp;
	//协方差矩阵
	Eigen::Matrix3f covariance = matrixOp.CalculateCovariance(Vertices);
	//schmidt正交特征向量--旋转轴
	obb_data.RotationMatrix = matrixOp.CalculateEigenvectors(covariance);
	// 中心点计算
	float infinity = std::numeric_limits<float>::infinity();
    Eigen::RowVector3f minExtents(infinity, infinity, infinity);
	Eigen::RowVector3f maxExtents(-infinity, -infinity, -infinity);
	for (unsigned int i = 0; i < Vertices.rows(); i++)
	{
		Eigen::RowVector3f vertex = Vertices.row(i);
		Eigen::RowVector3f displacement = vertex - Centre_0;
		Eigen::RowVector3f allExtents = displacement * obb_data.RotationMatrix;
		maxExtents = maxExtents.cwiseMax(allExtents);
		minExtents = minExtents.cwiseMin(allExtents);
	}
	Eigen::RowVector3f offset = (maxExtents - minExtents) / 2.0f + minExtents;
	
	obb_data.Centre += (offset(0) * obb_data.RotationMatrix.col(0).adjoint()) +
		(offset(1) * obb_data.RotationMatrix.col(1).adjoint()) +
		(offset(2) * obb_data.RotationMatrix.col(2).adjoint());
	//二分之一边长
	obb_data.HalfExtents = (maxExtents -  minExtents) / 2.0f;
	return obb_data;
}
//变换包围盒
OBB_data OBB::transform(Eigen::Matrix4f & transMatrix)
{   
	Eigen::RowVector4f Centre4f;
	Eigen::Matrix4f RotationMatrix4f;
	Eigen::RowVector4f HalfExtents4f;
	//中心点
	Centre4f << obb_data.Centre, 1.0;
	Centre4f = (transMatrix * Centre4f.adjoint()).adjoint();
	obb_data.Centre = Centre4f.head(3);
	//二分之一边长
	HalfExtents4f << obb_data.HalfExtents, 1.0;
	HalfExtents4f = (transMatrix * HalfExtents4f.adjoint()).adjoint();
	obb_data.HalfExtents = HalfExtents4f.head(3);
	//旋转矩阵
	RotationMatrix4f.setIdentity(4, 4);
	RotationMatrix4f.block(0, 0, 3, 3) = obb_data.RotationMatrix;
	RotationMatrix4f = transMatrix * RotationMatrix4f;
	obb_data.RotationMatrix = RotationMatrix4f.block(0, 0, 3, 3);
	return obb_data;
}
//判断点与面的法向是否同向
bool OBB::toNormalDirection(Vector3f Normal, RowVector3f Centre, RowVector3f Vertice)
{
	RowVector3f direction;
	direction = Vertice - Centre;
	bool isSameDirection = direction.dot(Normal) >= 0 ? 1 : 0;
	return isSameDirection;
}
//设置与第一主成分方向作为法向量一致方向的顶点集合
vector<Eigen::MatrixXf> OBB::set_PCAVertices(Eigen::MatrixXf Vertices, Vector3f Normal, RowVector3f Centre)
{
	vector<RowVector3f> vert[2];
	for (unsigned int i = 0; i < Vertices.rows(); i++)
	{
		if (toNormalDirection(Normal, Centre, Vertices.row(i)))
		{
			vert[0].push_back(Vertices.row(i));
		}
		else
		{
			vert[1].push_back(Vertices.row(i));
		}
	}
	unsigned int numVertices0 = vert[0].size();
	unsigned int numVertices1 = vert[1].size();
	Eigen::MatrixXf p0Vertices(numVertices0, 3);
	Eigen::MatrixXf p1Vertices(numVertices1, 3);
	for (unsigned int i = 0; i < numVertices0; i++)
	{
		p0Vertices.row(i) = vert[0][i];
	}
	for (unsigned int i = 0; i < numVertices1; i++)
	{
		p1Vertices.row(i) = vert[1][i];
	}
	vector<Eigen::MatrixXf> matrix_vertix;
	matrix_vertix.push_back(p0Vertices);
	matrix_vertix.push_back(p1Vertices);
	return matrix_vertix;
}


//获取顶点集合
Eigen::MatrixXf OBB::get_Vertices()
{
	return Vertices;
}
//获取中心点
RowVector3f OBB::get_Centre()
{
	return obb_data.Centre;
}
//获取二分之一边长
RowVector3f OBB::get_HalfExtents()
{
	return obb_data.HalfExtents;
}
//获取旋转矩阵
Matrix3f OBB::get_RotationMatrix()
{
	return obb_data.RotationMatrix;
}

//设置中心点//Centre_0上一层中心点
RowVector3f OBB::set_Centre(Eigen::MatrixXf Vertices, Matrix3f RotationMatrix, Eigen::RowVector3f Centre_0)
{
	// 中心点计算
	float infinity = std::numeric_limits<float>::infinity();
	Eigen::RowVector3f minExtents(infinity, infinity, infinity);
	Eigen::RowVector3f maxExtents(-infinity, -infinity, -infinity);
	for (unsigned int i = 0; i < Vertices.rows(); i++)
	{
		Eigen::RowVector3f vertex = Vertices.row(i);
		Eigen::RowVector3f displacement = vertex - Centre_0;
		Eigen::RowVector3f allExtents = displacement * RotationMatrix;
		maxExtents = maxExtents.cwiseMax(allExtents);
		minExtents = minExtents.cwiseMin(allExtents);
	}
	Eigen::RowVector3f offset = (maxExtents - minExtents) / 2.0f + minExtents;

	obb_data.Centre += (offset(0) * RotationMatrix.col(0).adjoint()) +
		(offset(1) * RotationMatrix.col(1).adjoint()) +
		(offset(2) * RotationMatrix.col(2).adjoint());
	return obb_data.Centre;

}
//设置二分之一边长
RowVector3f OBB::set_HalfExtents(Eigen::MatrixXf Vertices, Matrix3f RotationMatrix, RowVector3f Centre_0)
{
	float infinity = std::numeric_limits<float>::infinity();
	Eigen::RowVector3f minExtents(infinity, infinity, infinity);
	Eigen::RowVector3f maxExtents(-infinity, -infinity, -infinity);
	for (unsigned int i = 0; i <  Vertices.rows(); i++)
	{
		Eigen::RowVector3f vertex = Vertices.row(i);
		Eigen::RowVector3f displacement = vertex - Centre_0;
		Eigen::RowVector3f allExtents = displacement * RotationMatrix;
		maxExtents = maxExtents.cwiseMax(allExtents);
		minExtents = minExtents.cwiseMin(allExtents);
	}
	Eigen::RowVector3f offset = (maxExtents - minExtents) / 2.0f + minExtents;
	//二分之一边长
	obb_data.HalfExtents = (maxExtents - minExtents) / 2.0f;
	return obb_data.HalfExtents;
}
//设置旋转矩阵
Matrix3f OBB::set_RotationMatrix(Eigen::MatrixXf Vertices)
{
	MATRIXOP matrixOp;
	//协方差矩阵
	Eigen::Matrix3f covariance = matrixOp.CalculateCovariance(Vertices);
	//schmidt正交特征向量--旋转轴
	obb_data.RotationMatrix = matrixOp.CalculateEigenvectors(covariance);
	return obb_data.RotationMatrix;
}

//第一次输入顶点//根结点的顶点集
MatrixXf OBB::set_Vertexs(vector<float> vertex)
{
	MATRIXOP op;
	Vertices = op.vector_trans_matrix(vertex);
	return Vertices;
}

void OBB::printobb()
{
	cout << "OBB_centre\n" << obb_data.Centre << endl;
	cout << "rotation matrix:\n" << obb_data.RotationMatrix << endl;
	cout << "1/2 edge:\n" << obb_data.HalfExtents << endl;
}

//判断点是否在包围盒内
bool OBB::point_in_OBB(RowVector3f point, OBB obb)
{
	Eigen::RowVector3f displacement = point - obb.get_Centre();
	Eigen::RowVector3f proPoint = displacement * obb.get_RotationMatrix();
	Eigen::RowVector3f reduceMax = obb.get_HalfExtents() - proPoint;
	Eigen::RowVector3f reduceMin = -obb.get_HalfExtents() - proPoint;
	if (reduceMax.minCoeff() <= 0 || reduceMin.maxCoeff() >= 0)
	{
		cout << "POINT IS NOT IN OBBBOX!" << endl;
		return false;
	}
	else
	{
		cout << "POINT IS IN OBBBOX!" << endl;
		return true;
	}	
}

//检测两个包围盒是否碰撞
//判断两个包围盒在分离轴上的投影的半径和小于包围盒中心点间距在分离轴的投影距离，那么包围盒A与包围盒B处于分离状态。
//包含3+3+3*3 = 15个分离轴，分别是两个包围盒的6个轴，分别垂直于两包围盒的轴的轴；
bool  OBB::intersects(OBB box0, OBB box1)
{
	MatrixXf Separating_Axis(3, 15);
	MatrixXf Separating_Axis_cross(3, 9);
	unsigned int k = 0;
	for (unsigned int j = 0; j < 3; j++)
	{
		for (unsigned int i = 0; i < 3; i++,k++)
		{
			Separating_Axis_cross.col(k) = box0.get_RotationMatrix().col(j).cross(box1.get_RotationMatrix().col(i));
		}
	}
	Separating_Axis << box0.get_RotationMatrix(), box1.get_RotationMatrix(), Separating_Axis_cross;
	RowVector3f centreDistance = box0.get_Centre() - box1.get_Centre();
	//两个OBB盒中心点到分离轴的投影向量
	RowVectorXf projectCentre(15);
	for (unsigned i = 0; i <15 ; i++)
	{
		projectCentre(i) = abs(centreDistance.dot(Separating_Axis.col(i)));
	}
	//两个OBB盒半边长到分离轴的投影和向量
	RowVectorXf projectHalfExtends(15);
	Vector3f extentAxis0;
	Vector3f extentAxis1;
	for (unsigned j = 0; j <15; j++)
	{
		float projectHalfExtends0 = 0, projectHalfExtends1 = 0;
		for (unsigned i = 0; i < 3; i++)
		{
			projectHalfExtends0 += abs(box0.get_RotationMatrix().col(i).dot(Separating_Axis.col(j)))*box0.get_HalfExtents()(i);
			projectHalfExtends1 += abs(box0.get_RotationMatrix().col(i).dot(Separating_Axis.col(j)))*box0.get_HalfExtents()(i);
		}
		projectHalfExtends(j) = projectHalfExtends0 + projectHalfExtends1;
	}
	//用中心点投影减半边长投影的最小值如果小于0，则发生碰撞；否则，未发生碰撞。
	RowVectorXf reduce(15);
	reduce = projectCentre - projectHalfExtends;
	if (reduce.maxCoeff() > 0)
	{
		cout << "NOT COLLISION\n" << endl;
		return true;
	}
	else
	{
		cout << "COLLISION\n" << endl;
		return false;
	}


		
}
