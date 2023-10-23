#include <blitz/array.h>
#include <mxml.h>

#include <iostream>

void main(void)
{
	FILE *fp;
    mxml_node_t *tree;

    fp = fopen("locations.xml", "r");
    tree = mxmlLoadFile(NULL, fp, MXML_NO_CALLBACK);
    fclose(fp);

	mxml_node_t *node = mxmlFindElement(tree, tree, "zheight", NULL, NULL, MXML_DESCEND);

	
	cout << atof( node->child->value.text.string ) << endl;
}

