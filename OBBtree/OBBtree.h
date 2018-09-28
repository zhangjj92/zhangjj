#pragma once
#pragma execution_character_set("utf-8")
#ifndef OBBTREE_H
#define OBBTREE_H
#include<OBB.CPP>

using namespace std;
using namespace Eigen;

//建立包围盒树
struct OBBtreeNode {

public:
	// 无参数的构造函数 
	OBBtreeNode();  
	//析构函数
	~OBBtreeNode();
	//参数初始化
	OBBtreeNode(OBBtreeNode *tree, OBBtreeNode *lChild, OBBtreeNode *rChild);

	//建立根结点
	OBBtreeNode* CreateRoot(OBB &obb, Eigen::MatrixXf Vertices, Eigen::RowVector3f Centre_0);
	//获取根结点
	OBBtreeNode* get_root();
	//获取左孩子结点
	OBBtreeNode* get_leftChild(OBBtreeNode* father);
	//获取右孩子结点
	OBBtreeNode* get_rightChild(OBBtreeNode* father);

	//建立二叉树结点
	OBBtreeNode * CreateNode(OBB &obb,OBBtreeNode *BiTree, Eigen::MatrixXf Vertices, Eigen::RowVector3f Centre);
	//按照前序顺序建立二叉树
	void CreateTree(OBB &obb, Eigen::MatrixXf Vertices);
	//先序遍历二叉树
	void preTraverse(OBBtreeNode *tree, unsigned int i);
	//打印结点
	void printNode(OBBtreeNode *node);

private:
	OBB_data obbData;
	OBBtreeNode *leftChild, *rightChild;
	OBBtreeNode *root;
};

#endif // !OBBtree_h
