#pragma once
#pragma execution_character_set("utf-8")
#include<OBBtree.h>

// 无参数的构造函数 
OBBtreeNode::OBBtreeNode() 
{ 
	leftChild = rightChild = NULL; 
};  // 无参数的构造函数 

//析构函数
OBBtreeNode::~OBBtreeNode()
{
}
//初始化
OBBtreeNode::OBBtreeNode(OBBtreeNode *tree, OBBtreeNode *lChild, OBBtreeNode *rChild )
{
	tree->leftChild = lChild;
	tree->rightChild = rChild;
}

//建立根结点
OBBtreeNode* OBBtreeNode::CreateRoot(OBB &obb, Eigen::MatrixXf Vertices, Eigen::RowVector3f Centre_0)
{
	OBB_data data = obb.generateOBB(Vertices, obb.get_Centre());
	root = new OBBtreeNode;
	root->obbData.RotationMatrix = data.RotationMatrix;
	root->obbData.Centre = data.Centre;
	root->obbData.HalfExtents = data.HalfExtents;
	return root;
}
//获取根结点
OBBtreeNode* OBBtreeNode::get_root()
{
	return root;
}

//获取左孩子结点
OBBtreeNode* OBBtreeNode::get_leftChild(OBBtreeNode* father)
{
	return father->leftChild;
}
//获取右孩子结点
OBBtreeNode* OBBtreeNode::get_rightChild(OBBtreeNode* father)
{
	return father->rightChild;
}

//建立二叉树结点
OBBtreeNode* OBBtreeNode::CreateNode(OBB &obb, OBBtreeNode *BiTree, Eigen::MatrixXf Vertices, Eigen::RowVector3f Centre_0)
{	
	if (Vertices.rows() <= 10)
		return NULL;
	else
	{
		OBB_data data = obb.generateOBB(Vertices, Centre_0);
		BiTree = new OBBtreeNode;
		BiTree->obbData.RotationMatrix = data.RotationMatrix;
		BiTree->obbData.Centre = data.Centre;
	    BiTree->obbData.HalfExtents = data.HalfExtents;  
		//cout << "centre:\n" << BiTree->obbData.Centre << endl;
		//cout << "rotation matrix:\n" << BiTree->obbData.RotationMatrix << endl;
		//cout << "1/2 edge:\n" << BiTree->obbData.HalfExtents << endl;
		vector<Eigen::MatrixXf>  obb_vertexs;
		obb_vertexs = obb.set_PCAVertices(Vertices, BiTree->obbData.RotationMatrix.col(0), BiTree->obbData.Centre);
		BiTree->rightChild = CreateNode(obb, BiTree->leftChild, obb_vertexs[0], BiTree->obbData.Centre);
		BiTree->leftChild = CreateNode(obb, BiTree->rightChild, obb_vertexs[1], BiTree->obbData.Centre);
	}	
	return BiTree;
}

//按照前序顺序建立二叉树
void OBBtreeNode::CreateTree(OBB &obb, Eigen::MatrixXf Vertices)
{
	if (root == NULL)
	{
		cout << "NO OBBTREE!\n" << endl;
		return;
	}
	vector<Eigen::MatrixXf>  obb_vertexs;
	obb_vertexs = obb.set_PCAVertices(Vertices, root->obbData.RotationMatrix.col(0), root->obbData.Centre);
	//结点
	root->rightChild = CreateNode(obb, root->leftChild, obb_vertexs[0], root->obbData.Centre);
	root->leftChild = CreateNode(obb, root->rightChild, obb_vertexs[1], root->obbData.Centre);
}

//先序遍历二叉树//一般根结点为0，i初始值为0
void OBBtreeNode::preTraverse(OBBtreeNode *tree, unsigned int i)
{
	if (tree)
	{  
		cout <<"\n\n"<< i << "th layer node:\n" << endl;
		cout<<"centre:\n" << tree->obbData.Centre << endl;
		cout << "rotation matrix:\n" << tree->obbData.RotationMatrix << endl;
		cout << "1/2 edge:\n" << tree->obbData.HalfExtents << endl;
		++i;
		preTraverse(tree->leftChild,i);
		preTraverse(tree->rightChild,i);
	}
}

//打印结点
void OBBtreeNode::printNode(OBBtreeNode *node)
{
	cout << "centre:\n" << node->obbData.Centre << endl;
	cout << "rotation matrix:\n" << node->obbData.RotationMatrix << endl;
	cout << "1/2 edge:\n" << node->obbData.HalfExtents << endl;
}
