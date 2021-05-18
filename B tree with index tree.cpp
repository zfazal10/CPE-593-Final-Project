#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class Index_BTree;
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class Index_BTreeNode;

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
struct location_in_IndexBTree {
	Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> * node;
	int index;
};

struct student {
	int ID;	
	string firstName;
	string lastName;
	string DOB;
	string major;
	float GPA;
};

struct faculty {
	int ID;
	string firstName;
	string lastName;
	string subject;
	int salary;
	
};


template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class Index_BTreeNode
{
	Index_BTreeNode **C; // An array of child pointers	
	bool leaf; // Is true when node is leaf. Otherwise false
	int n;     // Current number of keys

public:	
	Record **keys;  // An array of keys
	Index_BTreeNode(bool _leaf);
	void traverse();
	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> search(const IndexKeyType &k);

	int findKey(const IndexKeyType &k);

	void insertNonFull(const Record &k);
	void splitChild(int i, Index_BTreeNode *y);

	void remove(const Record &k);
	void removeFromLeaf(int idx);
	void removeFromNonLeaf(int idx);
	const Record& getPred(int idx);
	const Record& getSucc(int idx);
	void fill(int idx);
	void borrowFromPrev(int idx);
	void borrowFromNext(int idx);
	void merge(int idx);

	friend Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>;
};

constexpr uint32_t blocksize = 1024; //<- should be 1024 on linux, 4096 on windows 10
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class Index_BTree
{
	Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> *root;
	
public:
	constexpr static uint32_t rec_per_block = (blocksize - 8) / (8 + 8);
	constexpr static uint32_t min_rec_per_block = (rec_per_block + 1) / 2;
	constexpr static uint32_t max_rec_per_block = min_rec_per_block * 2 - 1;

	Index_BTree() : root(new Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>(true)) {
	}

	void remove(const Record &k);
	void insert(const Record &k);
	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> search(IndexKeyType &k) {
		return root->search(k);
	}

	void traverse() {
		root->traverse();
	}

};

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::Index_BTreeNode(bool leaf) :leaf(leaf), n(0) {
	keys = new Record *[Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::max_rec_per_block];
	C = new Index_BTreeNode *[Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::max_rec_per_block + 1];

}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
int Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::findKey(const IndexKeyType &k)
{
	int idx = 0;
	while (idx < n && *keys[idx].*IndexKey < k)
		++idx;
	return idx;
}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::remove(const Record &k)
{
	Record *x = const_cast<Record*> (&k);
	int idx = findKey(*x.*IndexKey);

	if (idx < n && *keys[idx].*IndexKey == *x.*IndexKey) {
		if (leaf) {
			removeFromLeaf(idx);
		}
		else {
			removeFromNonLeaf(idx);
		}
	}
	else {
		if (leaf) {
			return;
		}
		bool flag = ((idx == n) ? true : false);

		if (C[idx]->n < Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block)
			fill(idx);
		if (flag && idx > n)
			C[idx - 1]->remove(k);
		else
			C[idx]->remove(k);
	}
	return;
}




template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::removeFromLeaf(int idx) {
	for (int i = idx + 1; i < n; ++i) {
		keys[i - 1] = keys[i];
	}
	n--;

	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::removeFromNonLeaf(int idx) {
	const Record k = (*keys[idx]);

	if (C[idx]->n >= Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block) {
		const Record& pred = getPred(idx);
		Record *x = const_cast<Record*> (&pred);
		keys[idx] = x;
		C[idx]->remove(pred);
	}

	else if (C[idx + 1]->n >= Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block) {
		const Record& succ = getSucc(idx);
		Record *x = const_cast<Record*> (&succ);
		keys[idx] = x;
		C[idx + 1]->remove(succ);
	}

	else {
		merge(idx);
		C[idx]->remove(k);
	}
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
const Record& Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::getPred(int idx) {
	Index_BTreeNode *cur = C[idx];
	while (!cur->leaf)
		cur = cur->C[cur->n];

	return *(cur->keys[cur->n - 1]);
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
const Record& Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::getSucc(int idx) {

	Index_BTreeNode *cur = C[idx + 1];
	while (!cur->leaf)
		cur = cur->C[0];

	return *(cur->keys[0]);
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::fill(int idx) {

	if (idx != 0 && C[idx - 1]->n >= Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block)
		borrowFromPrev(idx);

	else if (idx != n && C[idx + 1]->n >= Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block)
		borrowFromNext(idx);

	else {
		if (idx != n)
			merge(idx);
		else
			merge(idx - 1);
	}
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::borrowFromPrev(int idx) {
	Index_BTreeNode *child = C[idx];
	Index_BTreeNode *sibling = C[idx - 1];

	for (int i = child->n - 1; i >= 0; --i)
		child->keys[i + 1] = child->keys[i];

	if (!child->leaf) {
		for (int i = child->n; i >= 0; --i)
			child->C[i + 1] = child->C[i];
	}

	child->keys[0] = keys[idx - 1];

	if (!child->leaf)
		child->C[0] = sibling->C[sibling->n];

	keys[idx - 1] = sibling->keys[sibling->n - 1];

	child->n += 1;
	sibling->n -= 1;

	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::borrowFromNext(int idx) {
	Index_BTreeNode *child = C[idx];
	Index_BTreeNode *sibling = C[idx + 1];

	child->keys[(child->n)] = keys[idx];

	if (!(child->leaf))
		child->C[(child->n) + 1] = sibling->C[0];

	keys[idx] = sibling->keys[0];

	for (int i = 1; i < sibling->n; ++i)
		sibling->keys[i - 1] = sibling->keys[i];

	if (!sibling->leaf) {
		for (int i = 1; i <= sibling->n; ++i)
			sibling->C[i - 1] = sibling->C[i];
	}

	child->n += 1;
	sibling->n -= 1;
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::merge(int idx) {
	Index_BTreeNode *child = C[idx];
	Index_BTreeNode *sibling = C[idx + 1];

	child->keys[Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1] = keys[idx];
	for (int i = 0; i < sibling->n; ++i)
		child->keys[i + Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block] = sibling->keys[i];

	if (!child->leaf) {
		for (int i = 0; i <= sibling->n; ++i)
			child->C[i + Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block] = sibling->C[i];
	}

	for (int i = idx + 1; i < n; ++i)
		keys[i - 1] = keys[i];

	for (int i = idx + 2; i <= n; ++i)
		C[i - 1] = C[i];

	child->n += sibling->n + 1;
	n--;
	delete(sibling);
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::insert(const Record &k) {
	Record *x = const_cast<Record*> (&k);
	if (root->n == 0) {

		root->keys[0] = x;
		root->n = 1;
	}
	else {
		if (root->n == max_rec_per_block) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>*s = new Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>(false);
			s->C[0] = root;
			s->splitChild(0, root);

			int i = 0;
			if (s->keys[0]->*IndexKey < *x.*IndexKey)
				i++;
			s->C[i]->insertNonFull(k);
			root = s;
		}
		else
			root->insertNonFull(k);
	}
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::insertNonFull(const Record &k) {
	Record *x = const_cast<Record*> (&k);

	int i = n - 1;
	if (leaf == true) {
		while (i >= 0 && *keys[i].*IndexKey > *x.*IndexKey) {
			keys[i + 1] = keys[i];
			i--;
		}

		keys[i + 1] = x;
		n = n + 1;
	}
	else {
		while (i >= 0 && *keys[i].*IndexKey > *x.*IndexKey)
			i--;

		if (C[i + 1]->n == Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::max_rec_per_block) {
			splitChild(i + 1, C[i + 1]);

			if (*keys[i + 1].*IndexKey < *x.*IndexKey)
				i++;
		}
		C[i + 1]->insertNonFull(k);
	}
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::splitChild(int i, Index_BTreeNode *y) {
	Index_BTreeNode *z = new Index_BTreeNode(y->leaf);
	z->n = Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1;

	for (int j = 0; j < Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1; j++)
		z->keys[j] = y->keys[j + Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block];

	if (y->leaf == false) {
		for (int j = 0; j < Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block; j++)
			z->C[j] = y->C[j + Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block];
	}
	y->n = Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1;

	for (int j = n; j >= i + 1; j--)
		C[j + 1] = C[j];

	C[i + 1] = z;

	for (int j = n - 1; j >= i; j--)
		keys[j + 1] = keys[j];

	keys[i] = y->keys[Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1];
	n = n + 1;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::traverse() {
	int i;
	for (i = 0; i < n; i++) {
		if (leaf == false) {
			Record **that = &C[i]->keys[0];
			C[i]->traverse();
		}
		cout << " (" << (*keys[i]).*IndexKey << ") ";
	}

	if (leaf == false) {
		C[i]->traverse();
	}
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::search(const IndexKeyType &k) {
	int i = 0;
	while (i < n && k > *keys[i].*IndexKey)
		i++;
	if (i < n) {
		if (*keys[i].*IndexKey == k) {
			return { this, i };
		}
	}
	if (leaf == true) {
		return { NULL, NULL };
	}
	return C[i]->search(k);
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::remove(const Record &k) {
	if (root->n == 0) {
		cout << "The tree is empty\n";
		return;
	}
	root->remove(k);

	if (root->n == 0) {
		Index_BTreeNode <Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> *tmp = root;
		if (root->leaf)
			root->n = 0;
		else
			root = root->C[0];
		delete tmp;
	}
	return;
}

//================================================================	End of Index Classes ===================================================================


template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class BTree;
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class BTreeNode;

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
struct location_in_BTree {
	BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> * node;
	int index;
};

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class BTreeNode
{

	BTreeNode **C; // An array of child pointers
	int n;     // Current number of keys
	bool leaf; // Is true when node is leaf. Otherwise false

public:
	Record *keys;  // An array of keys
	BTreeNode(bool _leaf);

	void traverse() {
		int i;
		for (i = 0; i < n; i++) {
			if (leaf == false) {
				C[i]->traverse();
			}
			Record *ptr = &keys[i];
			cout << " (" << keys[i].*PrimaryKey << "), ";			
		}

		if (leaf == false) {
			C[i]->traverse();
		}
	}

	location_in_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> search(const PrimaryKeyType &k);
	
	int findKey(const PrimaryKeyType &k);

	void insertNonFull(const Record &k, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	void splitChild(int i, BTreeNode *y, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);

	void remove(const Record &k, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	void removeFromLeaf(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	void removeFromNonLeaf(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	const Record& getPred(int idx);
	const Record& getSucc(int idx);
	void fill(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	void borrowFromPrev(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	void borrowFromNext(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);
	void merge(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree);


	friend BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>;
};

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
class BTree
{
	BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> *root;
	Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>	i_tree;

public:
	
	constexpr static uint32_t rec_per_block = (blocksize - 8) / (sizeof(Record) + 8);
	constexpr static uint32_t min_rec_per_block = (rec_per_block + 1) / 2;
	constexpr static uint32_t max_rec_per_block = min_rec_per_block * 2 - 1;

	BTree() : root(new BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>(true)) {

	}

	void remove(const Record &k);
	void insert(const Record &k);
	location_in_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> search(PrimaryKeyType k) {
		return root->search(k);
	}
	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> search(IndexKeyType k) {
		return i_tree.search(k);
	}

	void traverse() {		
		cout << "Primary B-Tree Traversal: " << endl;
		root->traverse();
		cout << endl;
		cout << endl;
		cout << "Index B-Tree Traversal: " << endl;
		i_tree.traverse();
		cout << endl;
		cout << endl;
	}

};

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::BTreeNode(bool leaf) :leaf(leaf), n(0) {
	keys = new Record[BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::max_rec_per_block];
	C = new BTreeNode *[BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::max_rec_per_block + 1];
}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
int BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::findKey(const PrimaryKeyType &k)
{
	int idx = 0;
	while (idx < n && keys[idx].*PrimaryKey < k)
		++idx;
	return idx;
}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::remove(const Record &k, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree)
{
	int idx = findKey(k.*PrimaryKey);

	if (idx < n && keys[idx].*PrimaryKey == k.*PrimaryKey) {
		if (leaf) {
			removeFromLeaf(idx, i_tree);
		}
		else {
			removeFromNonLeaf(idx, i_tree);
		}
	}
	else {
		if (leaf) {
			return;
		}

		bool flag = ((idx == n) ? true : false);

		if (C[idx]->n < BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block) {
			fill(idx, i_tree);
		}
		if (flag && idx > n)
			C[idx - 1]->remove(k, i_tree);
		else
			C[idx]->remove(k, i_tree);
	}
	return;
}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::removeFromLeaf(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	int i = idx + 1;
	for (i = idx + 1; i < n; ++i) {

		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[i].*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &keys[i - 1];
		}
		
		keys[i - 1] = keys[i];
	}

	n--;	
	return;
}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::removeFromNonLeaf(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	const Record k = keys[idx];
	if (C[idx]->n >= BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block) {
		const Record& pred = getPred(idx);

		Record *unconst_pred = const_cast<Record*> (&pred);
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(*unconst_pred.*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &keys[idx];
		}
		


		keys[idx] = pred;		//key change!
		C[idx]->remove(pred, i_tree);
	}

	else if (C[idx + 1]->n >= BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block) {
		const Record& succ = getSucc(idx);
		
		Record *unconst_succ = const_cast<Record*> (&succ);
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(*unconst_succ.*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &keys[idx];
		}
		

		keys[idx] = succ;	//key change!
		C[idx + 1]->remove(succ, i_tree);
	}

	else {
		merge(idx, i_tree);
		C[idx]->remove(k, i_tree);
	}
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
const Record& BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::getPred(int idx) {
	BTreeNode *cur = C[idx];
	while (!cur->leaf)
		cur = cur->C[cur->n];

	return cur->keys[cur->n - 1];
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
const Record& BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::getSucc(int idx) {
	BTreeNode *cur = C[idx + 1];
	while (!cur->leaf)
		cur = cur->C[0];

	return cur->keys[0];
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::fill(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	if (idx != 0 && C[idx - 1]->n >= BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block)
		borrowFromPrev(idx, i_tree);

	else if (idx != n && C[idx + 1]->n >= BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block)
		borrowFromNext(idx, i_tree);

	else {
		if (idx != n)
			merge(idx, i_tree);
		else
			merge(idx - 1, i_tree);
	}
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::borrowFromPrev(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	BTreeNode *child = C[idx];
	BTreeNode *sibling = C[idx - 1];

	for (int i = child->n - 1; i >= 0; --i) {
	
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(child->keys[i].*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &(child->keys[i + 1]);
		}
		
		child->keys[i + 1] = child->keys[i];
	}

	if (!child->leaf) {
		for (int i = child->n; i >= 0; --i)
			child->C[i + 1] = child->C[i];
	}


	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[idx - 1].*IndexKey);
	if (i_tree_key_loc.node != nullptr) {
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &(child->keys[0]);
	}

	child->keys[0] = keys[idx - 1];

	if (!child->leaf)
		child->C[0] = sibling->C[sibling->n];

	i_tree_key_loc = i_tree.search(sibling->keys[sibling->n - 1].*IndexKey);
	if (i_tree_key_loc.node != nullptr) {
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &(keys[idx - 1]);
	}

	keys[idx - 1] = sibling->keys[sibling->n - 1];

	child->n += 1;
	sibling->n -= 1;

	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::borrowFromNext(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	BTreeNode *child = C[idx];
	BTreeNode *sibling = C[idx + 1];

	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[idx].*IndexKey);
	if (i_tree_key_loc.node != nullptr) {
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &(child->keys[(child->n)]);
	}

	child->keys[(child->n)] = keys[idx];

	if (!(child->leaf))
		child->C[(child->n) + 1] = sibling->C[0];

	i_tree_key_loc = i_tree.search(sibling->keys[0].*IndexKey);
	if (i_tree_key_loc.node != nullptr) {
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &(keys[idx]);
	}

	keys[idx] = sibling->keys[0];

	for (int i = 1; i < sibling->n; ++i) {
		i_tree_key_loc = i_tree.search(sibling->keys[i].*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &(sibling->keys[i - 1]);
		}

		sibling->keys[i - 1] = sibling->keys[i];
	}

	if (!sibling->leaf) {
		for (int i = 1; i <= sibling->n; ++i)
			sibling->C[i - 1] = sibling->C[i];
	}

	child->n += 1;
	sibling->n -= 1;
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::merge(int idx, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	BTreeNode *child = C[idx];
	BTreeNode *sibling = C[idx + 1];
	int minRec = BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block;
	
	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[idx].*IndexKey);
	if (i_tree_key_loc.node != nullptr) {
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &(child->keys[minRec - 1]);
	}

	child->keys[minRec - 1] = keys[idx];		//key change

	for (int i = 0; i < sibling->n; ++i) {
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(sibling->keys[i].*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &(child->keys[i + minRec]);
		}
		
		child->keys[i + minRec] = sibling->keys[i];		//key change
	}

	if (!child->leaf) {
		for (int i = 0; i <= sibling->n; ++i)
			child->C[i + minRec] = sibling->C[i];
	}

	for (int i = idx + 1; i < n; ++i) {
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[i].*IndexKey);
		if (i_tree_key_loc.node != nullptr) {
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &(keys[i - 1]);
		}
		
		keys[i - 1] = keys[i];		//key change
	}

	for (int i = idx + 2; i <= n; ++i)
		C[i - 1] = C[i];

	child->n += sibling->n + 1;
	n--;
	delete(sibling);
	return;
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::insert(const Record &k) {
	
	if (root->n == 0) {
		root->keys[0] = k;
		root->n = 1;
		Record* ptr = &(root->keys[0]);
		i_tree.insert(*ptr);
	}


	else {
		if (root->n == max_rec_per_block) {
			BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>*s = new BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>(false);
			s->C[0] = root;
			s->splitChild(0, root, i_tree);

			int i = 0;
			if (s->keys[0].*PrimaryKey < k.*PrimaryKey)
				i++;
			
			s->C[i]->insertNonFull(k, i_tree);
			root = s;
		}
		else {
			root->insertNonFull(k, i_tree);
		}
	}

}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::insertNonFull(const Record &k, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	int i = n - 1;
	if (leaf == true) {

		while (i >= 0 && keys[i].*PrimaryKey > k.*PrimaryKey) { //TODO: analyze complexity of functions with index tree: logn^2

			location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[i].*IndexKey);
			
			Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
			i_node.keys[i_tree_key_loc.index] = &keys[i+1];
			
			keys[i + 1] = keys[i];
			i--;
		}
		keys[i + 1] = k;
		n = n + 1;
		Record* ptr = &(keys[i + 1]);
		i_tree.insert(*ptr);
	}
	else {
		while (i >= 0 && keys[i].*PrimaryKey > k.*PrimaryKey)
			i--;

		if (C[i + 1]->n == BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::max_rec_per_block) {
			splitChild(i + 1, C[i + 1], i_tree);
			if (keys[i + 1].*PrimaryKey < k.*PrimaryKey)
				i++;
		}
		
		C[i + 1]->insertNonFull(k, i_tree);
	}
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::splitChild(int i, BTreeNode *y, Index_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> &i_tree) {
	int minRecs = BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block;
	BTreeNode *z = new BTreeNode(y->leaf);
	z->n = BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1;

	for (int j = 0; j < BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1; j++) {
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(y->keys[j + minRecs].*IndexKey);
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &(z->keys[j]);

		z->keys[j] = y->keys[j + minRecs];
	}

	if (y->leaf == false) {
		for (int j = 0; j < BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block; j++)
			z->C[j] = y->C[j + BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block];

	}
	y->n = BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::min_rec_per_block - 1;

	for (int j = n; j >= i + 1; j--)
		C[j + 1] = C[j];

	C[i + 1] = z;

	for (int j = n - 1; j >= i; j--) {
		location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(keys[j].*IndexKey);
		Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
		i_node.keys[i_tree_key_loc.index] = &keys[j + 1];
		
		keys[j + 1] = keys[j];
	}

	location_in_IndexBTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_tree_key_loc = i_tree.search(y->keys[minRecs - 1].*IndexKey);
	Index_BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> i_node = *i_tree_key_loc.node;
	i_node.keys[i_tree_key_loc.index] = &keys[i];

	keys[i] = y->keys[minRecs - 1];
	n = n + 1;
}

template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
location_in_BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> BTreeNode<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::search(const PrimaryKeyType &k) {
	
	int i = 0;
	while (i < n && k > keys[i].*PrimaryKey)
		i++;

	if (keys[i].*PrimaryKey == k)
		return { this, i };

	if (leaf == true)
		return { NULL, NULL };

	return C[i]->search(k);
}
template<typename Record, typename PrimaryKeyType, PrimaryKeyType Record::*PrimaryKey, typename IndexKeyType, IndexKeyType Record::*IndexKey>
void BTree<Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey>::remove(const Record &k) {
	if (root->n == 0) {
		cout << "The tree is empty\n";
		return;
	}	

	i_tree.remove(k);
	root->remove(k, i_tree);

	if (root->n == 0) {
		BTreeNode <Record, PrimaryKeyType, PrimaryKey, IndexKeyType, IndexKey> *tmp = root;
		if (root->leaf)
			root->n = 0;
		else
			root = root->C[0];
		delete tmp;
	}
	return;
}


//================================================================	Beginning of main function ===================================================================


int main() {

	//create and fill faculty database
	BTree<faculty, int, &faculty::ID, string, &faculty::lastName> faculties;	
	ifstream facultyCSV("Faculty Database Data Table.csv");
	string line;	
	while (getline(facultyCSV, line, '\n')) {
		if (line != "First Name,Last Name,ID,Subject,Salary") {
			vector<string> result;
			stringstream s_stream(line); //create string stream from the string
			while (s_stream.good()) {
				string substr;
				getline(s_stream, substr, ','); //get first string delimited by comma
				result.push_back(substr);
			}

			string FN = result.at(0);
			string LN = result.at(1);
			string ID_str = result.at(2);
			stringstream to_int(ID_str);
			int ID = 0;
			to_int >> ID;
			string Subject = result.at(3);
			string Sal_str = result.at(4);
			stringstream to_int2(Sal_str);
			int salary = 0;
			to_int2 >> salary;
			faculties.insert({ ID, FN, LN, Subject, salary });
		}
	}
	
	//create and fill student database
	BTree<student, int, &student::ID, string, &student::lastName> students;
	ifstream studentCSV("Student Database Data Table.csv");
	while (getline(studentCSV, line, '\n')) {
		if (line != "First Name,Last Name,ID,Birthdate,Major,GPA") {
			vector<string> result;
			stringstream s_stream(line); //create string stream from the string
			while (s_stream.good()) {
				string substr;
				getline(s_stream, substr, ','); //get first string delimited by comma
				result.push_back(substr);
			}

			string FN = result.at(0);
			string LN = result.at(1);
			string ID_str = result.at(2);
			stringstream to_int(ID_str);
			int ID = 0;
			to_int >> ID;
			string DOB = result.at(3);
			string Major = result.at(4);
			string GPA_str = result.at(5);
			float GPA = stof(GPA_str);
			students.insert({ ID, FN, LN, DOB, Major, GPA });
		}
	}


	char database = ' ';
	cout << "Welcome to the b tree database demo! Two b trees have already been created: one for students and one for faculty" << endl;
	cout << "Both of these b trees are sorted by student/faculty ID, and both have created an index b tree of pointers sorted by last name" << endl;
	cout << "Insert, Remove, Traverse, and Search functions are all available at your leisure for either database (everything is case sensitive!)" << endl;
	while (database != 's' && database != 'f') {
		cout << "If you would like to interact with the student database type 's'. If you would like to interact with the faculty database type 'f'" << endl;
		cin >> database;
	}
	
	bool stud;
	if (database == 's') stud = true;
	else   stud = false;

	string input = " ";
	while (input != "q")
	{
		cout << "type 'i' to insert, 'r' to remove, 't' to traverse, 's' to search, 'c' to change to the other database, and 'q' to quit the program" << endl;
		cin >> input;
		if (input == "i") {
			int ID;
			string firstName;
			string lastName;
			cout << "Enter first name (string): ";
			cin >> firstName;
			cout << "Enter last name (string): ";
			cin >> lastName;
			if (stud) {
				string DOB;
				string major;
				float GPA;
				bool validID = false;
				while (validID == false) {
					cout << "Enter ID (int): ";
					cin >> ID;
					if (students.search(ID).node == nullptr) {
						validID = true;
					}
					else {
						cout << "That ID is already taken!" << endl;
					}
				}
				cout << "Enter date of birth (MM/DD/YYYY)(string): ";
				cin >> DOB;
				cout << "Enter major (no spaces)(string): ";
				cin >> major;
				cout << "Enter GPA (up to three digits)(float): ";
				cin >> GPA;
				students.insert({ ID, firstName, lastName, DOB, major, GPA });
			}
			else {
				string subject;
				int salary;
				bool validID = false;
				while (validID == false) {
					cout << "Enter ID (int): ";
					cin >> ID;
					if (faculties.search(ID).node == nullptr) {
						validID = true;
					}
					else {
						cout << "That ID is already taken!" << endl;
					}
				}
				cout << "Enter subject (string): ";
				cin >> subject;
				cout << "Enter salary (int): ";
				cin >> salary;
				faculties.insert({ ID, firstName, lastName, subject, salary });
			}
		}
		else if (input == "r") {
			cout << "Instead of typing in all info for student/faculty to remove, search for them with a specific key and then opt to remove" << endl;
		}
		else if (input == "t") {
			if (stud)	students.traverse();
			else        faculties.traverse();
		}
		else if (input == "s") {
			int tree;
			cout << "Would you like to search by [0]ID or [1]last name? type 0 or 1: ";
			cin >> tree;
			if (tree == 0) {
				int search;
				cout << "type ID (int) to search for: ";
				cin >> search;
				if (stud) {
					location_in_BTree<student, int, &student::ID, string, &student::lastName> loc = students.search(search);
					if (loc.node == nullptr) {
						cout << "student not found" << endl;
					}
					else {
						cout << "Student located at node at: " << loc.node << ", at index: " << loc.index << endl;
						BTreeNode<student, int, &student::ID, string, &student::lastName> node = *loc.node;
						student found = node.keys[loc.index];
						cout << "first name: " << found.firstName << endl;
						cout << "last name: " << found.lastName << endl;
						cout << "ID: " << found.ID << endl;
						cout << "DOB: " << found.DOB << endl;
						cout << "Major: " << found.major << endl;
						cout << "GPA: " << found.GPA << endl;
						int erase;
						cout << "Would you like to remove this student? (0 for yes, 1 for no): ";
						cin >> erase;
						if (erase == 0) {
							//remove
							
							students.remove(node.keys[loc.index]);
						}
					}
				}
				else {
					location_in_BTree<faculty, int, &faculty::ID, string, &faculty::lastName> loc = faculties.search(search);
					if (loc.node == nullptr) {
						cout << "faculty not found" << endl;
					}
					else {
						cout << "Faculty located at node at: " << loc.node << ", at index: " << loc.index << endl;
						BTreeNode<faculty, int, &faculty::ID, string, &faculty::lastName> node = *loc.node;
						faculty found = node.keys[loc.index];
						cout << "first name: " << found.firstName << endl;
						cout << "last name: " << found.lastName << endl;
						cout << "ID: " << found.ID << endl;
						cout << "Subject: " << found.subject << endl;
						cout << "Salary: $" << found.salary << endl;
						int erase;
						cout << "Would you like to remove this faculty? (0 for yes, 1 for no): ";
						cin >> erase;
						if (erase == 0) {
							//remove

							faculties.remove(node.keys[loc.index]);
						}
					}
				}
			}
			else if (tree == 1) {
				string search;
				cout << "type last name (string) to search for: ";
				cin >> search;
				if (stud) {
					location_in_IndexBTree<student, int, &student::ID, string, &student::lastName> loc = students.search(search);
					if (loc.node == nullptr) {
						cout << "student not found" << endl;
					}
					else {
						cout << "Student located at node at: " << loc.node << ", at index: " << loc.index << endl;
						Index_BTreeNode<student, int, &student::ID, string, &student::lastName> node = *loc.node;
						student found = *node.keys[loc.index];
						cout << "first name: " << found.firstName << endl;
						cout << "last name: " << found.lastName << endl;
						cout << "ID: " << found.ID << endl;
						cout << "DOB: " << found.DOB << endl;
						cout << "Major: " << found.major << endl;
						cout << "GPA: " << found.GPA << endl;
						int erase;
						cout << "Would you like to remove this student? (0 for yes, 1 for no): ";
						cin >> erase;
						if (erase == 0) {
							//remove

							students.remove(*node.keys[loc.index]);
						}
					}
				}
				else {
					location_in_IndexBTree<faculty, int, &faculty::ID, string, &faculty::lastName> loc = faculties.search(search);
					if (loc.node == nullptr) {
						cout << "faculty not found" << endl;
					}
					else {
						cout << "Faculty located at node at: " << loc.node << ", at index: " << loc.index << endl;
						Index_BTreeNode<faculty, int, &faculty::ID, string, &faculty::lastName> node = *loc.node;
						faculty found = *node.keys[loc.index];
						cout << "first name: " << found.firstName << endl;
						cout << "last name: " << found.lastName << endl;
						cout << "ID: " << found.ID << endl;
						cout << "Subject: " << found.subject << endl;
						cout << "Salary: $" << found.salary << endl;
						int erase;
						cout << "Would you like to remove this faculty? (0 for yes, 1 for no): ";
						cin >> erase;
						if (erase == 0) {
							//remove

							faculties.remove(*node.keys[loc.index]);
						}
					}
				}
			}
			else {
				cout << "Invalid input" << endl;
			}
		}
		else if (input == "c") {
			stud = !stud;
		}
		else if (input == "q") {
			break;
		}
		else {
			cout << "Invalid input, please try again: " << endl;
			cin >> input;
		}
	}

	system("pause");
	return 0;
}
